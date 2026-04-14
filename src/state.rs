//! `assist_get_state` — Query state of any solar system body or observatory.

use crate::coordinates::equatorial_to_ecliptic;
use crate::ffi;
use crate::observatory::ObservatoryTable;
use crate::origin::Origin;
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

/// Get the heliocentric ecliptic J2000 state of a body or observatory at one or more epochs.
///
/// # Arguments
/// - `ephem`: ASSIST ephemeris data.
/// - `origin`: the body or observatory to query.
/// - `epochs`: one or more epochs (MJD TDB).
/// - `obs_table`: optional observatory table (required if origin is an `Observatory`).
///
/// # Returns
/// One `BodyState` per epoch, in the same order.
pub fn assist_get_state(
    ephem: &Ephemeris,
    origin: &Origin,
    epochs: &[f64],
    obs_table: Option<&ObservatoryTable>,
) -> Result<Vec<BodyState>> {
    let mut results = Vec::with_capacity(epochs.len());
    for &epoch_mjd in epochs {
        let state = resolve_origin_state(ephem, origin, epoch_mjd, obs_table)?;
        results.push(BodyState {
            state,
            epoch: epoch_mjd,
        });
    }
    Ok(results)
}

/// Resolve the heliocentric ecliptic J2000 state of an origin at a single epoch.
///
/// Used internally by both `assist_get_state` and `assist_generate_ephemeris`.
pub(crate) fn resolve_origin_state(
    ephem: &Ephemeris,
    origin: &Origin,
    epoch_mjd: f64,
    obs_table: Option<&ObservatoryTable>,
) -> Result<[f64; 6]> {
    let jd_ref = ephem.jd_ref();
    let t = mjd_to_assist_time(epoch_mjd, jd_ref);

    if let Some(body_id) = origin.body_id() {
        // Named body: query ephemeris
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
        Ok(equatorial_to_ecliptic(&helio_eq))
    } else if let Origin::SolarSystemBarycenter = origin {
        // SSB in heliocentric = -Sun_bary
        let sun = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t)?;
        let helio_eq = [-sun.x, -sun.y, -sun.z, -sun.vx, -sun.vy, -sun.vz];
        Ok(equatorial_to_ecliptic(&helio_eq))
    } else if let Origin::Observatory(code) = origin {
        let table = obs_table.ok_or_else(|| {
            Error::Other("Observatory table required for observatory codes".into())
        })?;
        compute_observatory_state(ephem, code, t, table)
    } else {
        unreachable!()
    }
}

/// Compute heliocentric ecliptic state of a ground observatory.
fn compute_observatory_state(
    ephem: &Ephemeris,
    code: &str,
    t: f64,
    obs_table: &ObservatoryTable,
) -> Result<[f64; 6]> {
    let entry = obs_table
        .get(code)
        .ok_or_else(|| Error::InvalidObservatory(format!("Unknown MPC code: {code}")))?;

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
        return Ok(equatorial_to_ecliptic(&helio_eq));
    }

    // Compute observatory ECEF position from parallax coefficients
    let earth_radius_au = ephem.earth_radius_au();
    let lon_rad = entry.longitude_deg.to_radians();
    let jd = t + ephem.jd_ref();
    let mjd = jd - 2_400_000.5;
    let gmst = compute_gmst(mjd);
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

    Ok(equatorial_to_ecliptic(&helio_eq))
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

    let gmst_deg =
        280.460_618_37 + 360.985_647_366_29 * d + 0.000_387_933 * t * t - t * t * t / 38_710_000.0;

    // Normalize to [0, 360) then convert to radians
    let gmst_deg = gmst_deg.rem_euclid(360.0);
    gmst_deg.to_radians()
}

fn mjd_to_assist_time(mjd_tdb: f64, jd_ref: f64) -> f64 {
    (mjd_tdb + 2_400_000.5) - jd_ref
}
