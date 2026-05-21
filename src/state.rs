//! `assist_get_state` — Query state of any solar system body or observatory.

use crate::assist_data::AssistData;
use crate::coordinates::{bary_to_helio, equatorial_to_ecliptic};
use libassist_sys::Ephemeris;
use libassist_sys::ffi;

use crate::observatory::ObservatoryTable;
use crate::origin::Origin;
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
/// - `data`: bundle holding the ASSIST ephemeris and (for observatory
///   origins) the observatory table.
/// - `origin`: the body or observatory to query.
/// - `epochs`: one or more epochs (MJD TDB).
/// - `num_threads`: follows the same `None`/`Some(1)`/`Some(n)` convention as
///   [`crate::assist_propagate`]. Each epoch is cheap (~1–5 µs for a planetary
///   body, ~2–10 µs for an ITRF-rotated ground observatory), so serial
///   (`Some(1)`) is typically the right choice; parallelism pays off when
///   `epochs.len()` is in the thousands or more.
///
/// # Returns
/// One `BodyState` per epoch, in the same order.
pub fn assist_get_state(
    data: &AssistData,
    origin: &Origin,
    epochs: &[f64],
    num_threads: Option<usize>,
) -> Result<Vec<BodyState>> {
    let ephem = &data.ephem;
    let obs_table = data.observatory.as_ref();
    let op = |&epoch_mjd: &f64| -> Result<BodyState> {
        let state = resolve_origin_state(ephem, origin, epoch_mjd, obs_table)?;
        Ok(BodyState {
            state,
            epoch: epoch_mjd,
        })
    };
    crate::propagate::map_with_threads(epochs, num_threads, op)
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
    let t = ephem.mjd_to_assist_time(epoch_mjd);

    if let Some(body_id) = origin.body_id() {
        let bary_eq = ephem.get_body_state_array(body_id, t)?;
        let sun = ephem.get_body_state_array(ffi::ASSIST_BODY_SUN, t)?;
        Ok(equatorial_to_ecliptic(&bary_to_helio(&bary_eq, &sun)))
    } else if let Origin::SolarSystemBarycenter = origin {
        // SSB in heliocentric = -sun_bary.
        let sun = ephem.get_body_state_array(ffi::ASSIST_BODY_SUN, t)?;
        let helio_eq = [-sun[0], -sun[1], -sun[2], -sun[3], -sun[4], -sun[5]];
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

    let earth = ephem.get_body_state_array(ffi::ASSIST_BODY_EARTH, t)?;
    let sun = ephem.get_body_state_array(ffi::ASSIST_BODY_SUN, t)?;

    if entry.is_space_based() || code == "500" {
        // Geocentric: Earth centre, heliocentric.
        return Ok(equatorial_to_ecliptic(&bary_to_helio(&earth, &sun)));
    }

    // Ground observatory: require Earth orientation kernel. The previous
    // GMST-based fallback was good only to ~50 mas — fine for finding-charts
    // but well above modern astrometric precision and a silent error source
    // when callers forgot to attach an EOP kernel.
    let eo = obs_table
        .earth_orientation()
        .ok_or_else(|| Error::MissingEarthOrientation(code.to_string()))?;

    // Body-fixed (ITRF93) position of the observatory, AU.
    let earth_radius_au = ephem.earth_radius_au();
    let lon_rad = entry.longitude_deg.to_radians();
    let (sin_lon, cos_lon) = lon_rad.sin_cos();
    let bf_x = earth_radius_au * entry.cos_lat * cos_lon;
    let bf_y = earth_radius_au * entry.cos_lat * sin_lon;
    let bf_z = earth_radius_au * entry.sin_lat;

    // Rotate ITRF93 → ICRF/J2000 equatorial. ET avoids the precision loss
    // of an MJD round-trip inside the rotation lookup.
    let jd = t + ephem.jd_ref();
    let et = (jd - 2_451_545.0) * 86_400.0;
    let r = eo
        .rotation_itrf_to_j2000_et(et)
        .map_err(|e| Error::Other(format!("earth orientation lookup: {e}")))?;
    // Position: v_icrf = R · v_itrf
    let obs_x = r[0][0] * bf_x + r[0][1] * bf_y + r[0][2] * bf_z;
    let obs_y = r[1][0] * bf_x + r[1][1] * bf_y + r[1][2] * bf_z;
    let obs_z = r[2][0] * bf_x + r[2][1] * bf_y + r[2][2] * bf_z;
    // Velocity in ICRF: d/dt (R · r_body) = R · (ω × r_body), where
    // ω × r_body in the body frame is (-Ω·y_bf, +Ω·x_bf, 0) with Ω =
    // Earth sidereal rate. r_body is topocentric (time independent), so
    // take the body-frame rate and rotate it.
    const OMEGA_EARTH_RAD_S: f64 = 7.292_115_0e-5; // rad / sec
    const SEC_PER_DAY: f64 = 86_400.0;
    let omega = OMEGA_EARTH_RAD_S * SEC_PER_DAY; // rad / day
    let vb_x = -omega * bf_y;
    let vb_y = omega * bf_x;
    let vb_z = 0.0;
    let obs_vx = r[0][0] * vb_x + r[0][1] * vb_y + r[0][2] * vb_z;
    let obs_vy = r[1][0] * vb_x + r[1][1] * vb_y + r[1][2] * vb_z;
    let obs_vz = r[2][0] * vb_x + r[2][1] * vb_y + r[2][2] * vb_z;

    // Heliocentric topocentric = (Earth_bary + obs_offset) - Sun_bary.
    let obs_bary = [
        earth[0] + obs_x,
        earth[1] + obs_y,
        earth[2] + obs_z,
        earth[3] + obs_vx,
        earth[4] + obs_vy,
        earth[5] + obs_vz,
    ];
    Ok(equatorial_to_ecliptic(&bary_to_helio(&obs_bary, &sun)))
}
