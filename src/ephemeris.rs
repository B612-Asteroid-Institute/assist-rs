//! `assist_generate_ephemeris` — Propagate orbit to observer epochs with
//! light-time correction and topocentric spherical output.

use crate::coordinates::{cartesian_to_spherical, ecliptic_to_equatorial};
use crate::ffi;
use crate::propagate::assist_propagate;
use crate::wrappers::Ephemeris;
use crate::Result;

/// Ephemeris result for a single observer epoch.
#[derive(Debug, Clone)]
pub struct EphemerisResult {
    /// Topocentric spherical: [rho (AU), ra (rad), dec (rad), drho, dra, ddec].
    pub spherical: [f64; 6],
    /// Aberrated (light-time corrected) heliocentric ecliptic state (AU, AU/day).
    pub aberrated_state: [f64; 6],
    /// One-way light time in days.
    pub light_time: f64,
    /// Epoch of the observation (MJD TDB).
    pub epoch: f64,
}

/// Observer specification: pre-computed heliocentric ecliptic J2000 state.
#[derive(Debug, Clone)]
pub struct Observer {
    /// Heliocentric ecliptic J2000 state [x, y, z, vx, vy, vz] (AU, AU/day).
    pub state: [f64; 6],
    /// Epoch (MJD TDB).
    pub epoch: f64,
}

/// Generate ephemeris for a test orbit as seen from a set of observers.
///
/// For each observer at epoch t_obs:
/// 1. Propagate orbit to t_obs (initial estimate)
/// 2. Iterate light-time correction to find emission time t_emit
/// 3. Compute topocentric spherical coordinates (RA, Dec, range + rates)
///
/// # Arguments
/// - `ephem`: ASSIST ephemeris data.
/// - `orbit_state`: initial heliocentric ecliptic J2000 state (AU, AU/day).
/// - `orbit_epoch`: initial epoch (MJD TDB).
/// - `observers`: pre-computed observer states and epochs.
///
/// # Returns
/// One `EphemerisResult` per observer, in the same order.
pub fn assist_generate_ephemeris(
    ephem: &Ephemeris,
    orbit_state: &[f64; 6],
    orbit_epoch: f64,
    observers: &[Observer],
) -> Result<Vec<EphemerisResult>> {
    if observers.is_empty() {
        return Ok(vec![]);
    }

    let c = ephem.c_au_per_day();

    let mut results = Vec::with_capacity(observers.len());

    for obs in observers {
        let result = compute_single_ephemeris(
            ephem,
            orbit_state,
            orbit_epoch,
            obs,
            c,
        )?;
        results.push(result);
    }

    Ok(results)
}

/// Compute ephemeris for a single observer with light-time iteration.
fn compute_single_ephemeris(
    ephem: &Ephemeris,
    orbit_state: &[f64; 6],
    orbit_epoch: f64,
    obs: &Observer,
    c: f64,
) -> Result<EphemerisResult> {
    let t_obs = obs.epoch;

    // Step 1: Propagate to observation epoch for initial distance estimate
    let states = assist_propagate(ephem, orbit_state, orbit_epoch, &[t_obs], false)?;
    let helio_ecl_obs = states[0].state;

    // Initial light-time estimate
    let dx = [
        helio_ecl_obs[0] - obs.state[0],
        helio_ecl_obs[1] - obs.state[1],
        helio_ecl_obs[2] - obs.state[2],
    ];
    let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
    let mut tau = dist / c;

    // Step 2: Light-time iteration
    const MAX_ITER: usize = 10;
    // ~10 ns: above the ~4e-14 day integrator re-integration noise floor so
    // the loop can't limit-cycle between two adjacent f64 values. Position
    // error from this tolerance is < 1 m even for the closest NEOs.
    const TOL: f64 = 1e-13;
    let mut aberrated_state = helio_ecl_obs;

    for _ in 0..MAX_ITER {
        let t_emit = t_obs - tau;
        let emit_states = assist_propagate(ephem, orbit_state, orbit_epoch, &[t_emit], false)?;
        aberrated_state = emit_states[0].state;

        let dx_new = [
            aberrated_state[0] - obs.state[0],
            aberrated_state[1] - obs.state[1],
            aberrated_state[2] - obs.state[2],
        ];
        let dist_new = (dx_new[0] * dx_new[0] + dx_new[1] * dx_new[1] + dx_new[2] * dx_new[2]).sqrt();
        let tau_new = dist_new / c;

        if (tau_new - tau).abs() < TOL {
            tau = tau_new;
            break;
        }
        tau = tau_new;
    }

    // Step 3: Topocentric vector in barycentric equatorial (inertial frame).
    //
    // The astrometric direction must be computed in an inertial frame (the
    // SSB). Using heliocentric states would introduce a spurious offset
    // equal to the Sun's barycentric displacement during the light travel
    // time, which reaches ~9 mas for TNOs.
    let jd_ref = ephem.jd_ref();
    let t_emit = t_obs - tau;
    let t_obs_a = mjd_to_assist_time(t_obs, jd_ref);
    let t_emit_a = mjd_to_assist_time(t_emit, jd_ref);

    let sun_obs = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t_obs_a)?;
    let sun_emit = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t_emit_a)?;

    let obs_helio_eq = ecliptic_to_equatorial(&obs.state);
    let obs_bary_eq = [
        obs_helio_eq[0] + sun_obs.x,
        obs_helio_eq[1] + sun_obs.y,
        obs_helio_eq[2] + sun_obs.z,
        obs_helio_eq[3] + sun_obs.vx,
        obs_helio_eq[4] + sun_obs.vy,
        obs_helio_eq[5] + sun_obs.vz,
    ];

    let ast_helio_eq = ecliptic_to_equatorial(&aberrated_state);
    let ast_bary_eq = [
        ast_helio_eq[0] + sun_emit.x,
        ast_helio_eq[1] + sun_emit.y,
        ast_helio_eq[2] + sun_emit.z,
        ast_helio_eq[3] + sun_emit.vx,
        ast_helio_eq[4] + sun_emit.vy,
        ast_helio_eq[5] + sun_emit.vz,
    ];

    let topo_eq = [
        ast_bary_eq[0] - obs_bary_eq[0],
        ast_bary_eq[1] - obs_bary_eq[1],
        ast_bary_eq[2] - obs_bary_eq[2],
        ast_bary_eq[3] - obs_bary_eq[3],
        ast_bary_eq[4] - obs_bary_eq[4],
        ast_bary_eq[5] - obs_bary_eq[5],
    ];

    let dx_eq = [topo_eq[0], topo_eq[1], topo_eq[2]];
    let dv_eq = [topo_eq[3], topo_eq[4], topo_eq[5]];

    let spherical = cartesian_to_spherical(dx_eq, dv_eq);

    Ok(EphemerisResult {
        spherical,
        aberrated_state,
        light_time: tau,
        epoch: t_obs,
    })
}

fn mjd_to_assist_time(mjd_tdb: f64, jd_ref: f64) -> f64 {
    (mjd_tdb + 2_400_000.5) - jd_ref
}
