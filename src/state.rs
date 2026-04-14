//! `assist_get_state` — Query state of any solar system body or observatory.

use crate::coordinates::equatorial_to_ecliptic;
use crate::ffi;
use crate::observatory::ObservatoryTable;
use crate::wrappers::Ephemeris;
use crate::{Error, Result};

/// State vector result from `assist_get_state`.
#[derive(Debug, Clone)]
pub struct BodyState {
    /// Heliocentric ecliptic J2000 state [x, y, z, vx, vy, vz] (AU, AU/day).
    pub state: [f64; 6],
    /// Epoch (MJD TDB).
    pub epoch: f64,
}

/// Body identifiers that `assist_get_state` accepts.
#[derive(Debug, Clone)]
pub enum Target {
    /// Solar system body by ASSIST ID.
    Body(i32),
    /// MPC observatory code (e.g., "I11", "W84", "500").
    Observatory(String),
}

impl Target {
    /// Parse a target string: body names map to ASSIST IDs, anything else
    /// is treated as an MPC observatory code.
    pub fn parse(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "sun" => Target::Body(ffi::ASSIST_BODY_SUN),
            "mercury" => Target::Body(ffi::ASSIST_BODY_MERCURY),
            "venus" => Target::Body(ffi::ASSIST_BODY_VENUS),
            "earth" => Target::Body(ffi::ASSIST_BODY_EARTH),
            "moon" => Target::Body(ffi::ASSIST_BODY_MOON),
            "mars" => Target::Body(ffi::ASSIST_BODY_MARS),
            "jupiter" => Target::Body(ffi::ASSIST_BODY_JUPITER),
            "saturn" => Target::Body(ffi::ASSIST_BODY_SATURN),
            "uranus" => Target::Body(ffi::ASSIST_BODY_URANUS),
            "neptune" => Target::Body(ffi::ASSIST_BODY_NEPTUNE),
            "pluto" => Target::Body(ffi::ASSIST_BODY_PLUTO),
            _ => Target::Observatory(s.to_string()),
        }
    }
}

/// Get the heliocentric ecliptic J2000 state of a body or observatory.
///
/// For solar system bodies, queries the ASSIST ephemeris directly.
/// For observatories, computes the geocentric position from parallax
/// coefficients and adds it to Earth's heliocentric state.
pub fn assist_get_state(
    ephem: &Ephemeris,
    target: &str,
    epoch_mjd: f64,
    obs_table: Option<&ObservatoryTable>,
) -> Result<BodyState> {
    let jd_ref = ephem.jd_ref();
    let t = mjd_to_assist_time(epoch_mjd, jd_ref);

    match Target::parse(target) {
        Target::Body(body_id) => {
            // Get body state in barycentric equatorial ICRF
            let p = ephem.get_body_state(body_id, t)?;
            let bary_eq = [p.x, p.y, p.z, p.vx, p.vy, p.vz];

            // Convert to heliocentric: subtract Sun
            let sun = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t)?;
            let helio_eq = [
                bary_eq[0] - sun.x,
                bary_eq[1] - sun.y,
                bary_eq[2] - sun.z,
                bary_eq[3] - sun.vx,
                bary_eq[4] - sun.vy,
                bary_eq[5] - sun.vz,
            ];

            Ok(BodyState {
                state: equatorial_to_ecliptic(&helio_eq),
                epoch: epoch_mjd,
            })
        }
        Target::Observatory(code) => {
            let table = obs_table.ok_or_else(|| {
                Error::Other("Observatory table required for observatory codes".into())
            })?;
            compute_observatory_state(ephem, &code, epoch_mjd, t, table)
        }
    }
}

/// Compute heliocentric ecliptic state of a ground observatory.
fn compute_observatory_state(
    ephem: &Ephemeris,
    code: &str,
    epoch_mjd: f64,
    t: f64,
    obs_table: &ObservatoryTable,
) -> Result<BodyState> {
    // Geocentric observatory (code "500"): use Earth center directly
    let entry = obs_table.get(code).ok_or_else(|| {
        Error::InvalidObservatory(format!("Unknown MPC code: {code}"))
    })?;

    // Get Earth barycentric state
    let earth = ephem.get_body_state(ffi::ASSIST_BODY_EARTH, t)?;
    let sun = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t)?;

    if entry.is_space_based() || code == "500" {
        // Geocentric: Earth center relative to Sun
        let helio_eq = [
            earth.x - sun.x,
            earth.y - sun.y,
            earth.z - sun.z,
            earth.vx - sun.vx,
            earth.vy - sun.vy,
            earth.vz - sun.vz,
        ];
        return Ok(BodyState {
            state: equatorial_to_ecliptic(&helio_eq),
            epoch: epoch_mjd,
        });
    }

    // Compute observatory ECEF position from parallax coefficients
    let earth_radius_au = ephem.earth_radius_au();
    let lon_rad = entry.longitude_deg.to_radians();
    let gmst = compute_gmst(epoch_mjd);
    let theta = gmst + lon_rad; // local sidereal angle

    // Position in equatorial frame (AU, relative to Earth center)
    let (sin_theta, cos_theta) = theta.sin_cos();
    let obs_x = earth_radius_au * entry.cos_lat * cos_theta;
    let obs_y = earth_radius_au * entry.cos_lat * sin_theta;
    let obs_z = earth_radius_au * entry.sin_lat;

    // Velocity from Earth rotation (ω_earth ≈ 7.2921150e-5 rad/s = 6.30038809866574 rad/day)
    const OMEGA_EARTH: f64 = 6.300_388_098_665_74; // rad/day
    let obs_vx = -OMEGA_EARTH * obs_y;
    let obs_vy = OMEGA_EARTH * obs_x;
    let obs_vz = 0.0;

    // Heliocentric = (Earth_bary + obs_offset) - Sun_bary
    let helio_eq = [
        (earth.x + obs_x) - sun.x,
        (earth.y + obs_y) - sun.y,
        (earth.z + obs_z) - sun.z,
        (earth.vx + obs_vx) - sun.vx,
        (earth.vy + obs_vy) - sun.vy,
        (earth.vz + obs_vz) - sun.vz,
    ];

    Ok(BodyState {
        state: equatorial_to_ecliptic(&helio_eq),
        epoch: epoch_mjd,
    })
}

/// Greenwich Mean Sidereal Time in radians.
///
/// Uses the IAU 2006 approximation:
///   GMST = 280.46061837° + 360.98564736629° × (JD - 2451545.0)
///        + 0.000387933° × T² − T³/38710000
/// where T = Julian centuries from J2000.0.
fn compute_gmst(mjd_tdb: f64) -> f64 {
    let jd = mjd_tdb + 2_400_000.5;
    let d = jd - 2_451_545.0; // days from J2000.0
    let t = d / 36525.0; // Julian centuries

    let gmst_deg = 280.460_618_37
        + 360.985_647_366_29 * d
        + 0.000_387_933 * t * t
        - t * t * t / 38_710_000.0;

    // Normalize to [0, 360) then convert to radians
    let gmst_deg = gmst_deg.rem_euclid(360.0);
    gmst_deg.to_radians()
}

fn mjd_to_assist_time(mjd_tdb: f64, jd_ref: f64) -> f64 {
    (mjd_tdb + 2_400_000.5) - jd_ref
}
