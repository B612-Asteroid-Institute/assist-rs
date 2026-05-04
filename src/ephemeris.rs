//! `assist_generate_ephemeris_single` — Propagate orbit to observer epochs with
//! light-time correction and topocentric spherical output.

use crate::assist_data::AssistData;
use crate::coordinates::{cartesian_to_spherical, ecliptic_to_equatorial, equatorial_to_ecliptic};
use crate::ffi;
use crate::observatory::ObservatoryTable;
use crate::orbit::Orbit;
use crate::origin::Origin;
use crate::state::resolve_origin_state;
use crate::wrappers::{AssistSim, Ephemeris, IntegratorConfig, Simulation};
use crate::{Error, Result};

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

/// Observer specification: an origin and epoch.
///
/// The observer's heliocentric state is resolved internally from the origin.
#[derive(Debug, Clone)]
pub struct Observer {
    /// Where the observer is (a body or observatory).
    pub origin: Origin,
    /// Epoch (MJD TDB).
    pub epoch: f64,
}

impl Observer {
    /// Create an observer at a named body or observatory at a given epoch.
    pub fn new(origin: Origin, epoch: f64) -> Self {
        Self { origin, epoch }
    }
}

/// Generate ephemeris for a test orbit as seen from a set of observers.
///
/// For each observer at epoch t_obs:
/// 1. Propagate orbit to t_obs (initial estimate)
/// 2. Iterate light-time correction to find emission time t_emit
/// 3. Compute topocentric spherical coordinates (RA, Dec, range + rates)
///
/// Observers are independent — each gets a fresh internal simulation — so
/// this dispatches through rayon when `num_threads != 1` and the `parallel`
/// cargo feature is on.
///
/// # Arguments
/// - `data`: bundle holding the ASSIST ephemeris and (for observatory-based
///   observers) the observatory table.
/// - `orbit`: initial orbit (state, epoch, optional non-grav params).
/// - `observers`: observer origins and epochs.
/// - `num_threads`: `None` → rayon global pool (one worker per core),
///   `Some(1)` → serial (no rayon overhead), `Some(n)` for `n > 1` →
///   dedicated pool with `n` workers. `Some(0)` returns `Error::Other`.
///   Ignored without the `parallel` feature.
///
/// # Returns
/// One `EphemerisResult` per observer, in input order regardless of the
/// threading mode.
pub fn assist_generate_ephemeris_single(
    data: &AssistData,
    orbit: &Orbit,
    observers: &[Observer],
    num_threads: Option<usize>,
    integrator: &IntegratorConfig,
) -> Result<Vec<EphemerisResult>> {
    if observers.is_empty() {
        return Ok(vec![]);
    }

    let ephem = &data.ephem;
    let obs_table = data.observatory.as_ref();
    let c = ephem.c_au_per_day();
    let op = |obs: &Observer| -> Result<EphemerisResult> {
        // One simulation per observer — reused across the initial
        // propagation and every light-time iteration.
        let mut sim = EphemerisSim::new(ephem, orbit, integrator)?;
        compute_single_ephemeris(ephem, &mut sim, obs, c, obs_table)
    };
    crate::propagate::map_with_threads(observers, num_threads, op)
}

/// Generate ephemeris for many orbits viewed from a shared set of observers.
///
/// Shape: the returned `Vec<Vec<EphemerisResult>>` is indexed
/// `[orbit_index][observer_index]`.
///
/// Parallelism is across orbits. Within each orbit the observer loop runs
/// serially, reusing one `EphemerisSim` across observers (and exploiting
/// `assist_integrate_or_interpolate`'s in-step interpolation for
/// close-together observer epochs). When orbits outnumber cores — the
/// typical OD / catalog-wide case — this is the right axis to split on.
///
/// For the reverse case (one orbit, many observers that benefit from
/// splitting), use [`assist_generate_ephemeris_single`] with a non-trivial
/// `num_threads` instead.
///
/// `num_threads` follows the same `None` / `Some(1)` / `Some(n)` convention
/// as [`assist_generate_ephemeris_single`] and
/// [`crate::assist_propagate`].
pub fn assist_generate_ephemeris(
    data: &AssistData,
    orbits: &[Orbit],
    observers: &[Observer],
    num_threads: Option<usize>,
    integrator: &IntegratorConfig,
) -> Result<Vec<Vec<EphemerisResult>>> {
    if orbits.is_empty() {
        return Ok(vec![]);
    }
    // Each orbit runs the (serial-per-orbit) observer loop.
    let op = |orbit: &Orbit| -> Result<Vec<EphemerisResult>> {
        // `Some(1)` here means "no further parallelism inside this orbit"
        // — observers run serially so we don't nest rayon pools.
        assist_generate_ephemeris_single(data, orbit, observers, Some(1), integrator)
    };
    crate::propagate::map_with_threads(orbits, num_threads, op)
}

/// Compute ephemeris for a single observer with light-time iteration.
fn compute_single_ephemeris(
    ephem: &Ephemeris,
    sim: &mut EphemerisSim,
    obs: &Observer,
    c: f64,
    obs_table: Option<&ObservatoryTable>,
) -> Result<EphemerisResult> {
    let t_obs = obs.epoch;

    let obs_state = resolve_origin_state(ephem, &obs.origin, t_obs, obs_table)?;

    // Step 1: Propagate to observation epoch for initial distance estimate.
    let helio_ecl_obs = sim.integrate_to(t_obs)?;

    let dx = [
        helio_ecl_obs[0] - obs_state[0],
        helio_ecl_obs[1] - obs_state[1],
        helio_ecl_obs[2] - obs_state[2],
    ];
    let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
    let mut tau = dist / c;

    // Step 2: Light-time iteration.
    const MAX_ITER: usize = 10;
    // ~10 ns: above the ~4e-14 day integrator re-integration noise floor so
    // the loop can't limit-cycle between two adjacent f64 values. Position
    // error from this tolerance is < 1 m even for the closest NEOs.
    const TOL: f64 = 1e-13;
    let mut aberrated_state = helio_ecl_obs;
    let mut converged = false;

    for _ in 0..MAX_ITER {
        let t_emit = t_obs - tau;
        aberrated_state = sim.integrate_to(t_emit)?;

        let dx_new = [
            aberrated_state[0] - obs_state[0],
            aberrated_state[1] - obs_state[1],
            aberrated_state[2] - obs_state[2],
        ];
        let dist_new =
            (dx_new[0] * dx_new[0] + dx_new[1] * dx_new[1] + dx_new[2] * dx_new[2]).sqrt();
        let tau_new = dist_new / c;

        if (tau_new - tau).abs() < TOL {
            tau = tau_new;
            converged = true;
            break;
        }
        tau = tau_new;
    }

    if !converged {
        return Err(Error::LightTimeConvergence(MAX_ITER));
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

    let obs_helio_eq = ecliptic_to_equatorial(&obs_state);
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

/// Reusable simulation for ephemeris generation.
///
/// Owns an `AssistSim` set up for a single orbit and supports repeated
/// integration to different target epochs (forward or backward in time)
/// without re-creating the simulation. This amortizes the setup cost over
/// all light-time iterations for one observer.
struct EphemerisSim<'a> {
    asim: AssistSim,
    ephem: &'a Ephemeris,
    jd_ref: f64,
}

impl<'a> EphemerisSim<'a> {
    fn new(ephem: &'a Ephemeris, orbit: &Orbit, integrator: &IntegratorConfig) -> Result<Self> {
        let jd_ref = ephem.jd_ref();
        let t0 = mjd_to_assist_time(orbit.epoch, jd_ref);

        // Heliocentric ecliptic → barycentric equatorial ICRF
        let eq_state = ecliptic_to_equatorial(&orbit.state);
        let sun = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t0)?;
        let bary_state = [
            eq_state[0] + sun.x,
            eq_state[1] + sun.y,
            eq_state[2] + sun.z,
            eq_state[3] + sun.vx,
            eq_state[4] + sun.vy,
            eq_state[5] + sun.vz,
        ];

        let mut sim = Simulation::new()?;
        sim.set_t(t0);
        integrator.apply(&mut sim);
        let mut asim = AssistSim::new(sim, ephem)?;

        let non_grav = orbit.non_grav.as_ref();
        let mut forces = ffi::ASSIST_FORCES_DEFAULT;
        if non_grav.is_some() {
            forces |= ffi::ASSIST_FORCE_NON_GRAVITATIONAL;
        }
        asim.set_forces(forces);

        asim.sim_mut().add_test_particle(
            bary_state[0],
            bary_state[1],
            bary_state[2],
            bary_state[3],
            bary_state[4],
            bary_state[5],
        );

        if let Some(ng) = non_grav {
            if let Some(v) = ng.alpha {
                asim.set_alpha(v);
            }
            if let Some(v) = ng.nk {
                asim.set_nk(v);
            }
            if let Some(v) = ng.nm {
                asim.set_nm(v);
            }
            if let Some(v) = ng.nn {
                asim.set_nn(v);
            }
            if let Some(v) = ng.r0 {
                asim.set_r0(v);
            }
        }

        if let Some(ng) = non_grav {
            let n_total = asim.sim().n_particles();
            let mut params = vec![0.0f64; 3 * n_total];
            params[0] = ng.a1;
            params[1] = ng.a2;
            params[2] = ng.a3;
            asim.set_particle_params(params);
        }

        Ok(Self {
            asim,
            ephem,
            jd_ref,
        })
    }

    /// Integrate to the given epoch (MJD TDB) and return the heliocentric
    /// ecliptic J2000 state of the test particle.
    ///
    /// Uses `assist_integrate_or_interpolate` under the hood: the first
    /// call does a full IAS15 integration (with `exact_finish_time = 0`,
    /// potentially overshooting) and reconstructs the state at the exact
    /// target via the integrator's `br` Hermite-interpolation coefficients.
    /// Subsequent calls whose target falls inside the last completed step's
    /// interval skip integration entirely and return pure polynomial
    /// interpolation — typically ~1 µs vs ~50 µs for a short-back-step
    /// re-integration. This is the central optimization for the light-time
    /// iteration loop in [`compute_single_ephemeris`], where tau shrinks
    /// by femtoseconds-to-microseconds per iteration and successive
    /// t_emit values cluster inside the last completed step.
    fn integrate_to(&mut self, mjd_tdb: f64) -> Result<[f64; 6]> {
        let t_target = mjd_to_assist_time(mjd_tdb, self.jd_ref);
        self.asim.integrate_or_interpolate(t_target)?;

        let particles = self.asim.sim().particles();
        if particles.is_empty() {
            return Err(Error::Other("No particles after integration".into()));
        }
        let p = &particles[0];
        let bary_eq = [p.x, p.y, p.z, p.vx, p.vy, p.vz];

        let sun_t = self.ephem.get_body_state(ffi::ASSIST_BODY_SUN, t_target)?;
        let helio_eq = [
            bary_eq[0] - sun_t.x,
            bary_eq[1] - sun_t.y,
            bary_eq[2] - sun_t.z,
            bary_eq[3] - sun_t.vx,
            bary_eq[4] - sun_t.vy,
            bary_eq[5] - sun_t.vz,
        ];
        Ok(equatorial_to_ecliptic(&helio_eq))
    }
}

fn mjd_to_assist_time(mjd_tdb: f64, jd_ref: f64) -> f64 {
    (mjd_tdb + 2_400_000.5) - jd_ref
}
