//! Binary PCK reader for Earth orientation (ITRF93 ↔ ICRF/J2000).
//!
//! Parses JPL binary PCK files (DAF/PCK format) to evaluate the high-precision
//! Earth rotation kernel. A PCK file contains Type 2 segments with Chebyshev
//! coefficients for three Euler angles that describe the rotation from
//! ECLIPJ2000 into the body-fixed ITRF93 frame. This module reads the binary
//! format directly — no CSPICE dependency — and combines data from multiple
//! kernel files (historical, current high-precision, long-term predict) for
//! full temporal coverage.
//!
//! # Layout summary
//!
//! * File is a DAF: 1024-byte file record (header), then summary records in a
//!   linked list starting at `FWARD`, then segment data.
//! * For a `DAF/PCK`, `ND=2` (start/stop ET per segment) and `NI=5`
//!   (body, frame, type, initial DAF addr, final DAF addr).
//! * The Earth high-precision kernels use body 3000 (ITRF93), frame 17
//!   (ECLIPJ2000), and segment type 2 (Chebyshev polynomials, fixed interval).
//! * Each Type 2 segment stores `N` equal-length records back-to-back, followed
//!   by a 4-double directory `[INIT, INTLEN, RSIZE, N]`. Each record is
//!   `[MID, RADIUS, a1_cheby[NCOEF], a2_cheby[NCOEF], a3_cheby[NCOEF]]` with
//!   `NCOEF = (RSIZE - 2) / 3`.
//!
//! The Euler angles `(a1, a2, a3) = (π/2 + RA, π/2 - DEC, W)` describe the
//! 3-1-3 passive rotation from ECLIPJ2000 into ITRF93.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

/// NAIF frame ID for ECLIPJ2000 (ecliptic of J2000).
const FRAME_ECLIPJ2000: i32 = 17;
/// NAIF body ID for ITRF93 high-precision Earth body frame.
const BODY_ITRF93: i32 = 3000;
/// Type 2 PCK segment: Chebyshev coefficients, fixed interval length.
const PCK_TYPE_CHEBYSHEV: i32 = 2;
/// DAF record size in bytes.
const DAF_RECORD_BYTES: usize = 1024;
/// MJD (days) of the J2000 epoch.
const J2000_MJD: f64 = 51544.5;
/// Seconds per day.
const SEC_PER_DAY: f64 = 86400.0;
/// J2000 obliquity of the ecliptic used by SPICE's ECLIPJ2000 frame,
/// defined as exactly `84381.448 arcsec = 84381.448 · π / 648000 rad`.
/// Carried to full f64 precision so the ECLIPJ2000 → J2000-equatorial
/// rotation matches SPICE to ~1e-13.
const OBLIQUITY_ECLIPJ2000: f64 = 0.409_092_804_222_328_97;

/// Errors from the Earth orientation reader.
#[derive(Debug)]
pub enum EarthOrientationError {
    /// File I/O failure.
    Io(std::io::Error),
    /// File is not a DAF/PCK or cannot be parsed.
    BadFormat(String),
    /// No segment covers the requested epoch.
    OutOfRange { mjd_tdb: f64 },
}

impl std::fmt::Display for EarthOrientationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::BadFormat(m) => write!(f, "bad PCK format: {m}"),
            Self::OutOfRange { mjd_tdb } => {
                write!(f, "no PCK segment covers MJD TDB {mjd_tdb}")
            }
        }
    }
}

impl std::error::Error for EarthOrientationError {}

impl From<std::io::Error> for EarthOrientationError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

type Result<T> = std::result::Result<T, EarthOrientationError>;

// ─── Type 2 segment ────────────────────────────────────────────────────────

/// An in-memory Type 2 PCK segment: all Chebyshev records decoded to `f64`.
#[derive(Debug, Clone)]
struct Type2Segment {
    /// Segment ET coverage (inclusive start, inclusive stop; seconds past J2000 TDB).
    et_start: f64,
    et_stop: f64,
    /// ET of the first interval's left edge.
    init: f64,
    /// Length of each interval in ET seconds.
    intlen: f64,
    /// Number of Chebyshev coefficients per angle (`= (rsize - 2) / 3`).
    ncoef: usize,
    /// Flat storage for `n_records` records, each laid out as
    /// `[MID, RADIUS, a1[ncoef], a2[ncoef], a3[ncoef]]` (`rsize` doubles).
    records: Vec<f64>,
    /// Record stride in doubles (`= 2 + 3 * ncoef`).
    rsize: usize,
    /// Number of records.
    n_records: usize,
}

impl Type2Segment {
    fn covers(&self, et: f64) -> bool {
        et >= self.et_start && et <= self.et_stop
    }

    /// Locate the record containing `et` and return its (mid, radius, a1, a2, a3).
    fn locate(&self, et: f64) -> (f64, f64, &[f64], &[f64], &[f64]) {
        // Interval index. Clamp to the last record in case et == et_stop exactly.
        let raw = ((et - self.init) / self.intlen).floor() as i64;
        let idx = raw.clamp(0, self.n_records as i64 - 1) as usize;

        let base = idx * self.rsize;
        let mid = self.records[base];
        let radius = self.records[base + 1];
        let off = base + 2;
        let a1 = &self.records[off..off + self.ncoef];
        let a2 = &self.records[off + self.ncoef..off + 2 * self.ncoef];
        let a3 = &self.records[off + 2 * self.ncoef..off + 3 * self.ncoef];
        (mid, radius, a1, a2, a3)
    }
}

// ─── Public interface ─────────────────────────────────────────────────────

/// Earth orientation loader for binary PCK kernel files.
///
/// Loads one or more JPL Earth orientation PCKs and evaluates the ITRF93 →
/// ICRF/J2000 equatorial rotation matrix at any epoch covered by the kernels.
///
/// # Kernel priority
///
/// Files are searched in the order provided. For epochs covered by more than
/// one file, the earlier-listed kernel wins. A typical high-precision setup is:
///
/// ```no_run
/// # use assist_rs::earth_orientation::EarthOrientation;
/// # use std::path::PathBuf;
/// # let cache: PathBuf = PathBuf::new();
/// let eo = EarthOrientation::from_paths(&[
///     cache.join("earth_latest_high_prec.bpc"),
///     cache.join("earth_620120_250826.bpc"),
///     cache.join("earth_2025_250826_2125_predict.bpc"),
/// ])?;
/// # Ok::<(), assist_rs::earth_orientation::EarthOrientationError>(())
/// ```
///
/// which serves the current file first (highest precision, ~year 2000–2026),
/// falling back to the historical file for earlier dates and the long-term
/// predict for later dates.
#[derive(Debug)]
pub struct EarthOrientation {
    /// Segments in user-supplied order (priority = earlier first).
    segments: Vec<Type2Segment>,
}

impl EarthOrientation {
    /// Load and merge one or more binary PCK kernel files.
    pub fn from_paths<P: AsRef<Path>>(paths: &[P]) -> Result<Self> {
        let mut segments = Vec::new();
        for path in paths {
            let path_ref = path.as_ref();
            let mut file_segments = load_pck(path_ref)?;
            segments.append(&mut file_segments);
        }
        if segments.is_empty() {
            return Err(EarthOrientationError::BadFormat(
                "no usable segments found in supplied PCK files".into(),
            ));
        }
        Ok(Self { segments })
    }

    /// Return the rotation matrix `R` such that `v_j2000 = R @ v_itrf93`.
    ///
    /// The output is a 3×3 row-major matrix mapping vectors from ITRF93
    /// (Earth body-fixed) into the ICRF/J2000 equatorial frame (the frame
    /// ASSIST uses internally).
    pub fn rotation_itrf_to_j2000(&self, mjd_tdb: f64) -> Result<[[f64; 3]; 3]> {
        let et = (mjd_tdb - J2000_MJD) * SEC_PER_DAY;
        self.rotation_itrf_to_j2000_et(et).map_err(|e| match e {
            EarthOrientationError::OutOfRange { .. } => {
                EarthOrientationError::OutOfRange { mjd_tdb }
            }
            other => other,
        })
    }

    /// Like [`rotation_itrf_to_j2000`] but keyed on an ephemeris time (seconds
    /// past J2000 TDB). Used internally and by callers that already have ET,
    /// bypassing the ~2.5e-7 s precision loss of an MJD round-trip.
    pub fn rotation_itrf_to_j2000_et(&self, et: f64) -> Result<[[f64; 3]; 3]> {
        let seg = self.segments.iter().find(|s| s.covers(et)).ok_or(
            EarthOrientationError::OutOfRange {
                mjd_tdb: J2000_MJD + et / SEC_PER_DAY,
            },
        )?;

        let (mid, radius, a1_c, a2_c, a3_c) = seg.locate(et);
        // Chebyshev argument normalised to [-1, 1] within the interval.
        let x = (et - mid) / radius;
        let a1 = chebyshev_eval(a1_c, x);
        let a2 = chebyshev_eval(a2_c, x);
        let a3 = chebyshev_eval(a3_c, x);
        Ok(build_rotation(a1, a2, a3))
    }

    /// Return the rotation matrix `R` such that `v_itrf93 = R @ v_j2000`.
    pub fn rotation_j2000_to_itrf(&self, mjd_tdb: f64) -> Result<[[f64; 3]; 3]> {
        let m = self.rotation_itrf_to_j2000(mjd_tdb)?;
        Ok(transpose3(&m))
    }

    /// Number of loaded segments.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Union of ET coverage across all segments, as `(et_min, et_max)`.
    pub fn coverage_et(&self) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for s in &self.segments {
            if s.et_start < lo {
                lo = s.et_start;
            }
            if s.et_stop > hi {
                hi = s.et_stop;
            }
        }
        (lo, hi)
    }
}

// ─── DAF / PCK binary parsing ──────────────────────────────────────────────

fn load_pck(path: &Path) -> Result<Vec<Type2Segment>> {
    let mut file = File::open(path).map_err(|e| {
        EarthOrientationError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;

    // ── File record (1024 bytes) ──
    let mut header = [0u8; DAF_RECORD_BYTES];
    file.read_exact(&mut header)?;

    let locidw = std::str::from_utf8(&header[0..8])
        .map_err(|_| EarthOrientationError::BadFormat("LOCIDW is not ASCII".into()))?;
    if !locidw.starts_with("DAF/PCK") {
        return Err(EarthOrientationError::BadFormat(format!(
            "expected 'DAF/PCK' identifier, got {locidw:?}"
        )));
    }

    let locfmt = std::str::from_utf8(&header[88..96])
        .map_err(|_| EarthOrientationError::BadFormat("LOCFMT is not ASCII".into()))?;
    if !locfmt.starts_with("LTL-IEEE") {
        return Err(EarthOrientationError::BadFormat(format!(
            "only LTL-IEEE byte order supported, got {locfmt:?}"
        )));
    }

    let nd = read_i32_le(&header[8..12]);
    let ni = read_i32_le(&header[12..16]);
    let fward = read_i32_le(&header[76..80]);
    if nd != 2 || ni != 5 {
        return Err(EarthOrientationError::BadFormat(format!(
            "unexpected DAF/PCK (ND,NI); got ({nd},{ni}), want (2,5)"
        )));
    }

    // Summary "size" is the number of 8-byte words per summary.
    let summary_words = nd as usize + (ni as usize + 1) / 2;
    let summary_bytes = summary_words * 8;

    // ── Walk summary records, collecting (et_start, et_stop, body, frame, type, init_addr, final_addr) ──
    let mut segment_infos: Vec<SegmentInfo> = Vec::new();
    let mut next_rec = fward as i64;
    while next_rec > 0 {
        let mut rec = [0u8; DAF_RECORD_BYTES];
        read_daf_record(&mut file, next_rec as u64, &mut rec)?;

        let next_d = f64::from_le_bytes(rec[0..8].try_into().unwrap());
        // prev_d at 8..16 (unused)
        let nsum = f64::from_le_bytes(rec[16..24].try_into().unwrap()) as usize;

        for i in 0..nsum {
            let off = 24 + i * summary_bytes;
            let s = &rec[off..off + summary_bytes];
            let et_start = f64::from_le_bytes(s[0..8].try_into().unwrap());
            let et_stop = f64::from_le_bytes(s[8..16].try_into().unwrap());
            let body = read_i32_le(&s[16..20]);
            let frame = read_i32_le(&s[20..24]);
            let dtype = read_i32_le(&s[24..28]);
            let init_addr = read_i32_le(&s[28..32]);
            let final_addr = read_i32_le(&s[32..36]);

            segment_infos.push(SegmentInfo {
                et_start,
                et_stop,
                body,
                frame,
                dtype,
                init_addr,
                final_addr,
            });
        }

        next_rec = next_d as i64;
    }

    // ── Filter to Earth / ECLIPJ2000 / Chebyshev, decode each ──
    let mut out = Vec::new();
    for info in segment_infos {
        if info.body != BODY_ITRF93
            || info.frame != FRAME_ECLIPJ2000
            || info.dtype != PCK_TYPE_CHEBYSHEV
        {
            // Silently skip segments we don't know how to evaluate.
            continue;
        }
        out.push(load_type2_segment(&mut file, &info, path)?);
    }
    Ok(out)
}

#[derive(Debug, Clone, Copy)]
struct SegmentInfo {
    et_start: f64,
    et_stop: f64,
    body: i32,
    frame: i32,
    dtype: i32,
    /// DAF 1-indexed double-word address of first segment word.
    init_addr: i32,
    /// DAF 1-indexed double-word address of last segment word.
    final_addr: i32,
}

fn load_type2_segment(file: &mut File, info: &SegmentInfo, source: &Path) -> Result<Type2Segment> {
    let _ = source;
    // DAF addresses are 1-based in 8-byte words. Byte offset = (addr - 1) * 8.
    // Directory occupies the last 4 words of the segment.
    let dir_addr = (info.final_addr - 3) as u64; // 1-based addr of INIT
    let dir_byte = (dir_addr - 1) * 8;
    file.seek(SeekFrom::Start(dir_byte))?;
    let mut dir = [0u8; 32];
    file.read_exact(&mut dir)?;
    let init = f64::from_le_bytes(dir[0..8].try_into().unwrap());
    let intlen = f64::from_le_bytes(dir[8..16].try_into().unwrap());
    let rsize_f = f64::from_le_bytes(dir[16..24].try_into().unwrap());
    let n_f = f64::from_le_bytes(dir[24..32].try_into().unwrap());

    if !(intlen > 0.0) || !(rsize_f >= 3.0) || !(n_f >= 1.0) {
        return Err(EarthOrientationError::BadFormat(format!(
            "{}: invalid Type 2 directory INIT={init} INTLEN={intlen} RSIZE={rsize_f} N={n_f}",
            source.display()
        )));
    }

    let rsize = rsize_f as usize;
    let n_records = n_f as usize;
    if (rsize - 2) % 3 != 0 {
        return Err(EarthOrientationError::BadFormat(format!(
            "{}: RSIZE - 2 = {} is not divisible by 3",
            source.display(),
            rsize - 2
        )));
    }
    let ncoef = (rsize - 2) / 3;

    // Sanity: records + 4-word directory must equal segment length.
    let seg_words = (info.final_addr - info.init_addr + 1) as usize;
    let expected = n_records * rsize + 4;
    if seg_words != expected {
        return Err(EarthOrientationError::BadFormat(format!(
            "{}: segment size mismatch: {seg_words} words, expected {expected}",
            source.display()
        )));
    }

    // Read all records in one shot.
    let records_byte = (info.init_addr - 1) as u64 * 8;
    file.seek(SeekFrom::Start(records_byte))?;
    let total_bytes = n_records * rsize * 8;
    let mut raw = vec![0u8; total_bytes];
    file.read_exact(&mut raw)?;
    let mut records = vec![0.0f64; n_records * rsize];
    for (i, chunk) in raw.chunks_exact(8).enumerate() {
        records[i] = f64::from_le_bytes(chunk.try_into().unwrap());
    }

    Ok(Type2Segment {
        et_start: info.et_start,
        et_stop: info.et_stop,
        init,
        intlen,
        ncoef,
        records,
        rsize,
        n_records,
    })
}

fn read_daf_record(file: &mut File, record_number: u64, buf: &mut [u8]) -> Result<()> {
    // DAF records are 1-indexed, 1024 bytes each.
    let byte_off = (record_number - 1) * DAF_RECORD_BYTES as u64;
    file.seek(SeekFrom::Start(byte_off))?;
    file.read_exact(buf)?;
    Ok(())
}

fn read_i32_le(bytes: &[u8]) -> i32 {
    i32::from_le_bytes(bytes[0..4].try_into().unwrap())
}

// ─── Chebyshev evaluation and rotation construction ────────────────────────

/// Clenshaw evaluation of a Chebyshev-T series `sum_k c_k T_k(x)`.
fn chebyshev_eval(coefs: &[f64], x: f64) -> f64 {
    let n = coefs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return coefs[0];
    }
    let two_x = 2.0 * x;
    let mut b_kp1 = 0.0f64;
    let mut b_kp2 = 0.0f64;
    for k in (1..n).rev() {
        let b_k = coefs[k] + two_x * b_kp1 - b_kp2;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }
    coefs[0] + x * b_kp1 - b_kp2
}

/// Build the ITRF93 → ICRF/J2000 equatorial rotation matrix from the 3 Euler
/// angles read out of a Type 2 Earth PCK segment.
///
/// The angles encode a 3-1-3 passive rotation from ECLIPJ2000 into ITRF93:
/// `a1 = π/2 + RA_pole`, `a2 = π/2 - Dec_pole`, `a3 = W`.
/// We reconstruct the equivalent *active* rotation, take its transpose to flip
/// the direction, then pre-multiply by `Rx(+ε)` to go from ECLIPJ2000 to the
/// J2000 equatorial frame (ICRF) that ASSIST/REBOUND use internally.
fn build_rotation(a1: f64, a2: f64, a3: f64) -> [[f64; 3]; 3] {
    // Active-rotation inverse of NAIF passive 3-1-3: negated angles, reversed order.
    // M_eclj_to_itrf = Rz(-a3) · Rx(-a2) · Rz(-a1)
    let m_eclj_to_itrf = matmul3(
        &matmul3(&rot_z(-a3), &rot_x(-a2)),
        &rot_z(-a1),
    );
    let m_itrf_to_eclj = transpose3(&m_eclj_to_itrf);
    // ECLIPJ2000 → J2000 equatorial is an active rotation about +X by +ε.
    let r_eclj_to_j2000 = rot_x(OBLIQUITY_ECLIPJ2000);
    matmul3(&r_eclj_to_j2000, &m_itrf_to_eclj)
}

fn rot_x(theta: f64) -> [[f64; 3]; 3] {
    let (s, c) = theta.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]
}

fn rot_z(theta: f64) -> [[f64; 3]; 3] {
    let (s, c) = theta.sin_cos();
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

fn matmul3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for k in 0..3 {
                s += a[i][k] * b[k][j];
            }
            out[i][j] = s;
        }
    }
    out
}

fn transpose3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// Reference matrices produced by `spiceypy.pxform('J2000', 'ITRF93', et)`
    /// using `earth_latest_high_prec.bpc` (queried 2026-04-14). Each tuple is
    /// `(et_seconds_past_J2000_TDB, J2000→ITRF93 rotation matrix)`.
    #[allow(clippy::type_complexity)]
    const REFERENCE_MATRICES: &[(f64, [[f64; 3]; 3])] = &[
        (
            0.0,
            [
                [1.76980593920407242e-01, -9.84214340916322761e-01, -2.25882366727270245e-05],
                [9.84214340853249436e-01, 1.76980593228519362e-01, 2.96527445451250093e-05],
                [-2.51869769025144716e-05, -2.74796268088239870e-05, 9.99999999305243081e-01],
            ],
        ),
        (
            50000.0,
            [
                [-6.30643919962785349e-01, 7.76072320204723942e-01, 5.10181939461507028e-06],
                [-7.76072319714293024e-01, -6.30643919299331168e-01, -4.02995604470435254e-05],
                [-2.80579420008189413e-05, -2.93740535854958296e-05, 9.99999999174958409e-01],
            ],
        ),
        (
            86400.0,
            [
                [1.93884145136247354e-01, -9.81024432810858049e-01, -2.22018276302415885e-05],
                [9.81024432746447350e-01, 1.93884144450126555e-01, 2.97548960264265583e-05],
                [-2.48856976423294409e-05, -2.75495379365442439e-05, 9.99999999310862586e-01],
            ],
        ),
        (
            10_000_000.0,
            [
                [5.14075082542250206e-01, -8.57745188977827921e-01, -1.71630107713705016e-05],
                [8.57745189149325959e-01, 5.14075082453315679e-01, 9.58142343673862307e-06],
                [6.04656321041702567e-07, -1.96470609645782446e-05, 9.99999999806813755e-01],
            ],
        ),
        (
            86_400_000.0,
            [
                [-9.94874022034065053e-01, -1.01121833051657423e-01, 2.34865994768539144e-04],
                [1.01121833890863072e-01, -9.94874049637020752e-01, -8.32967580344234904e-06],
                [2.34504395423373261e-04, 1.54631020408135811e-05, 9.99999972384290059e-01],
            ],
        ),
        (
            100_000_000.0,
            [
                [-6.10321136736908021e-01, 7.92154080234992897e-01, 1.52378544191200671e-04],
                [-7.92154045927301431e-01, -6.10321155354273892e-01, 2.34196564775329685e-04],
                [2.78519613505794257e-04, 2.22278333399383854e-05, 9.99999960966373336e-01],
            ],
        ),
        (
            500_000_000.0,
            [
                [-5.40465715167075533e-01, -8.41365666675846113e-01, 7.90991814210717514e-04],
                [8.41364633946004137e-01, -5.40466291909377339e-01, -1.31911130469442028e-03],
                [1.53735937505101940e-03, -4.74218964589545600e-05, 9.99998817137958129e-01],
            ],
        ),
        (
            788_896_800.0,
            [
                [-9.86134927339477052e-01, -1.65928105443120949e-01, 2.40185453777035790e-03],
                [1.65927694454056252e-01, -9.86137851794280262e-01, -3.70771794214408423e-04],
                [2.43008113556495006e-03, 3.29031695190939999e-05, 9.99997046807167278e-01],
            ],
        ),
        (
            800_000_000.0,
            [
                [-7.60910330896355958e-01, 6.48854406998714550e-01, 1.85117625969774546e-03],
                [-6.48852373067780652e-01, -7.60912581883302686e-01, 1.62502184476126521e-03],
                [2.46298589273021918e-03, 3.53558005394827291e-05, 9.99996966220627947e-01],
            ],
        ),
    ];

    fn cache_path(name: &str) -> PathBuf {
        let home = std::env::var("HOME").expect("HOME env var");
        PathBuf::from(home).join(".cache/assist-rs").join(name)
    }

    fn max_abs_diff(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> f64 {
        let mut m = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                let d = (a[i][j] - b[i][j]).abs();
                if d > m {
                    m = d;
                }
            }
        }
        m
    }

    #[test]
    fn chebyshev_clenshaw_matches_direct() {
        // f(x) = 3 T_0 + 2 T_1 + 1 T_2 = 3 + 2x + (2x^2 - 1)
        let coefs = [3.0, 2.0, 1.0];
        for &x in &[-0.9, -0.3, 0.0, 0.4, 0.95] {
            let direct = 3.0 + 2.0 * x + (2.0 * x * x - 1.0);
            let clen = chebyshev_eval(&coefs, x);
            assert!(
                (direct - clen).abs() < 1e-14,
                "x={x}: direct={direct} clen={clen}"
            );
        }
    }

    #[test]
    fn rotation_matches_spice_reference() {
        let path = cache_path("earth_latest_high_prec.bpc");
        if !path.exists() {
            eprintln!("skipping: {} not present", path.display());
            return;
        }

        let eo = EarthOrientation::from_paths(&[&path]).expect("load PCK");
        assert!(eo.segment_count() >= 1);

        for (et, expected) in REFERENCE_MATRICES {
            let itrf_to_j2000 = eo.rotation_itrf_to_j2000_et(*et).unwrap();
            // Reference is J2000 → ITRF93 (pxform). Compare against our transpose.
            let j2000_to_itrf = transpose3(&itrf_to_j2000);
            let diff = max_abs_diff(&j2000_to_itrf, expected);
            assert!(
                diff < 1e-13,
                "et={et}: diff={diff:.3e}\ngot={:#?}\nexpected={:#?}",
                j2000_to_itrf,
                expected
            );
        }
    }

    /// Reference rotation matrices from the *historical* PCK.
    #[allow(clippy::type_complexity)]
    const HISTORICAL_REFERENCE: &[(f64, [[f64; 3]; 3])] = &[
        (
            -1_000_000_000.0,
            [
                [9.94364553444000543e-01, 1.05970302576597813e-01, 3.07080216520302929e-03],
                [-1.05969692207008692e-01, 9.94369294492217670e-01, -3.61254206012628920e-04],
                [-3.09179360005633732e-03, 3.38064169663754299e-05, 9.99995219823305304e-01],
            ],
        ),
        (
            -500_000_000.0,
            [
                [7.99970199503151957e-01, -6.00038405624525262e-01, 1.26161898562854313e-03],
                [6.00037686369575995e-01, 7.99971194223579762e-01, 9.29165618249683245e-04],
                [-1.56679390272404951e-03, 1.37141322137379618e-05, 9.99998772483641107e-01],
            ],
        ),
        (
            -100_000_000.0,
            [
                [2.96054724022567550e-01, 9.55170970204011538e-01, 1.34400329954154074e-04],
                [-9.55170930067680501e-01, 2.96054751062657939e-01, -2.80583088907757072e-04],
                [-3.07794677482192702e-04, -4.53073392118197305e-05, 9.99999951604839565e-01],
            ],
        ),
    ];

    /// Reference rotation matrices from the *long-term predict* PCK.
    #[allow(clippy::type_complexity)]
    const PREDICT_REFERENCE: &[(f64, [[f64; 3]; 3])] = &[
        (
            1_000_000_000.0,
            [
                [-9.68628820561433290e-01, -2.48493871286336765e-01, 3.00065174523937717e-03],
                [2.48492558544244874e-01, -9.68633467582096008e-01, -8.08596326593691295e-04],
                [3.10746293650065878e-03, -3.75900766640402395e-05, 9.99995171118883142e-01],
            ],
        ),
        (
            2_000_000_000.0,
            [
                [-3.19030968337769172e-01, 9.47742239093496308e-01, 1.97217633886831445e-03],
                [-9.47724230103578091e-01, -3.19037064009499882e-01, 5.84255618531083487e-03],
                [6.16643462995759405e-03, -5.12294501492371879e-06, 9.99980987348114470e-01],
            ],
        ),
        (
            3_500_000_000.0,
            [
                [9.39214248749958514e-01, -3.43181457839852599e-01, -1.01529276548461694e-02],
                [3.43163334530398545e-01, 9.39269108836152888e-01, -3.53086676571068292e-03],
                [1.07480593145401586e-02, -1.67872132489699233e-04, 9.99942223850917467e-01],
            ],
        ),
    ];

    #[test]
    fn historical_kernel_matches_spice() {
        let path = cache_path("earth_620120_250826.bpc");
        if !path.exists() {
            return;
        }
        let eo = EarthOrientation::from_paths(&[&path]).unwrap();
        for (et, expected) in HISTORICAL_REFERENCE {
            let m = eo.rotation_itrf_to_j2000_et(*et).unwrap();
            let diff = max_abs_diff(&transpose3(&m), expected);
            assert!(diff < 1e-13, "historical et={et}: diff={diff:.3e}");
        }
    }

    #[test]
    fn predict_kernel_matches_spice() {
        let path = cache_path("earth_2025_250826_2125_predict.bpc");
        if !path.exists() {
            return;
        }
        let eo = EarthOrientation::from_paths(&[&path]).unwrap();
        for (et, expected) in PREDICT_REFERENCE {
            let m = eo.rotation_itrf_to_j2000_et(*et).unwrap();
            let diff = max_abs_diff(&transpose3(&m), expected);
            assert!(diff < 1e-13, "predict et={et}: diff={diff:.3e}");
        }
    }

    #[test]
    fn three_kernel_coverage_spans_historical_current_predict() {
        let current = cache_path("earth_latest_high_prec.bpc");
        let historical = cache_path("earth_620120_250826.bpc");
        let predict = cache_path("earth_2025_250826_2125_predict.bpc");
        if !current.exists() || !historical.exists() || !predict.exists() {
            return;
        }
        let eo = EarthOrientation::from_paths(&[&current, &historical, &predict]).unwrap();

        // Historical-era ET (current kernel doesn't cover, historical does).
        let hist_et = -500_000_000.0;
        let m_hist = eo.rotation_itrf_to_j2000_et(hist_et).unwrap();
        let expected_hist = HISTORICAL_REFERENCE
            .iter()
            .find(|(et, _)| *et == hist_et)
            .unwrap()
            .1;
        assert!(max_abs_diff(&transpose3(&m_hist), &expected_hist) < 1e-13);

        // Current-era ET (current kernel wins per priority).
        let cur_et = 500_000_000.0;
        let m_cur = eo.rotation_itrf_to_j2000_et(cur_et).unwrap();
        let expected_cur = REFERENCE_MATRICES.iter().find(|(et, _)| *et == cur_et).unwrap().1;
        assert!(max_abs_diff(&transpose3(&m_cur), &expected_cur) < 1e-13);

        // Predict-era ET (far future, only predict kernel covers).
        let pred_et = 2_000_000_000.0;
        let m_pred = eo.rotation_itrf_to_j2000_et(pred_et).unwrap();
        let expected_pred = PREDICT_REFERENCE.iter().find(|(et, _)| *et == pred_et).unwrap().1;
        assert!(max_abs_diff(&transpose3(&m_pred), &expected_pred) < 1e-13);
    }

    #[test]
    fn out_of_range_returns_error() {
        let path = cache_path("earth_latest_high_prec.bpc");
        if !path.exists() {
            return;
        }
        let eo = EarthOrientation::from_paths(&[&path]).unwrap();
        // Far past, before current kernel coverage.
        let err = eo.rotation_itrf_to_j2000_et(-1e10).unwrap_err();
        assert!(matches!(err, EarthOrientationError::OutOfRange { .. }));
    }

    #[test]
    fn mjd_api_matches_et_api_within_micro_second() {
        let path = cache_path("earth_latest_high_prec.bpc");
        if !path.exists() {
            return;
        }
        let eo = EarthOrientation::from_paths(&[&path]).unwrap();
        // Round-tripping MJD ↔ ET loses a few hundred ns, which limits
        // rotation agreement to ~1e-10. That's ~2 μas at 1 AU — well below
        // our 1 mas astrometric tolerance.
        let mjd = 60676.5;
        let et = (mjd - J2000_MJD) * SEC_PER_DAY;
        let m_mjd = eo.rotation_itrf_to_j2000(mjd).unwrap();
        let m_et = eo.rotation_itrf_to_j2000_et(et).unwrap();
        assert!((max_abs_diff(&m_mjd, &m_et)) < 1e-9);
    }

    #[test]
    fn rotation_is_orthogonal_and_proper() {
        let path = cache_path("earth_latest_high_prec.bpc");
        if !path.exists() {
            eprintln!("skipping: {} not present", path.display());
            return;
        }
        let eo = EarthOrientation::from_paths(&[&path]).expect("load PCK");
        // MJD around 2025-01-01 TDB.
        let m = eo.rotation_itrf_to_j2000(60676.5).unwrap();
        // R^T R = I
        let rt = transpose3(&m);
        let id = matmul3(&rt, &m);
        for i in 0..3 {
            for j in 0..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((id[i][j] - want).abs() < 1e-14, "R^T R not identity");
            }
        }
        // det(R) = 1
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        assert!((det - 1.0).abs() < 1e-14, "det(R) = {det}");
    }
}
