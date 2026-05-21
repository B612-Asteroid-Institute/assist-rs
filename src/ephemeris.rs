//! `assist_generate_ephemeris_single` — Propagate orbit to observer epochs with
//! light-time correction and topocentric spherical output.

use crate::assist_data::AssistData;
use crate::coordinates::{
    bary_to_helio, cartesian_to_spherical, ecliptic_to_equatorial, equatorial_to_ecliptic,
    helio_to_bary,
};
use libassist_sys::ffi;
use libassist_sys::{AssistSim, Ephemeris, IntegratorConfig, Simulation};

use crate::observatory::ObservatoryTable;
use crate::orbit::Orbit;
use crate::origin::Origin;
use crate::state::resolve_origin_state;
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

    // Sort observers by epoch so each thread's chunk integrates monotonically
    // forward. `EphemerisSim::integrate_to` uses `integrate_or_interpolate`,
    // which extends the IAS15 trajectory and serves polynomial interpolation
    // for in-step queries — both depend on processing observers in time order
    // to amortize integrator state across the whole observer set.
    let mut indexed: Vec<(usize, &Observer)> = observers.iter().enumerate().collect();
    indexed.sort_by(|a, b| {
        a.1.epoch
            .partial_cmp(&b.1.epoch)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Build chunks of contiguous-in-time observers, one per worker. Each
    // chunk reuses a single `EphemerisSim` so the per-observer cost is the
    // marginal IAS15 integration to the next epoch, not a full cold-start.
    // Without this, 73 800 observers spent 90+% of wall time on
    // `Simulation::new` + force attach + first IAS15 step.
    let n_workers = match num_threads {
        Some(0) => {
            return Err(Error::Other(
                "num_threads = Some(0) is not valid; use None for the default pool".into(),
            ));
        }
        Some(n) => n,
        None => {
            #[cfg(feature = "parallel")]
            {
                rayon::current_num_threads().max(1)
            }
            #[cfg(not(feature = "parallel"))]
            {
                1
            }
        }
    };
    let n_workers = n_workers.min(indexed.len()).max(1);
    let chunk_size = indexed.len().div_ceil(n_workers);
    let chunks: Vec<&[(usize, &Observer)]> = indexed.chunks(chunk_size).collect();

    // Per-chunk worker: build one EphemerisSim, integrate observers in order,
    // emit (orig_index, EphemerisResult) so we can reorder back to caller's
    // original observer ordering.
    let process_chunk = |chunk: &&[(usize, &Observer)]| -> Result<Vec<(usize, EphemerisResult)>> {
        let mut sim = EphemerisSim::new(ephem, orbit, integrator)?;
        let mut out = Vec::with_capacity(chunk.len());
        for &(orig_idx, obs) in chunk.iter() {
            let r = compute_single_ephemeris(ephem, &mut sim, obs, c, obs_table)?;
            out.push((orig_idx, r));
        }
        Ok(out)
    };

    let chunk_results: Vec<Vec<(usize, EphemerisResult)>> =
        crate::propagate::map_with_threads(&chunks, num_threads, process_chunk)?;

    // Scatter chunk results back into caller's original ordering.
    let mut out: Vec<Option<EphemerisResult>> = (0..observers.len()).map(|_| None).collect();
    for chunk in chunk_results {
        for (idx, r) in chunk {
            out[idx] = Some(r);
        }
    }
    out.into_iter()
        .map(|opt| opt.ok_or_else(|| Error::Other("missing ephemeris result for observer".into())))
        .collect()
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
    let t_emit = t_obs - tau;
    let t_obs_a = ephem.mjd_to_assist_time(t_obs);
    let t_emit_a = ephem.mjd_to_assist_time(t_emit);

    let sun_obs = ephem.get_body_state_array(ffi::ASSIST_BODY_SUN, t_obs_a)?;
    let sun_emit = ephem.get_body_state_array(ffi::ASSIST_BODY_SUN, t_emit_a)?;

    let obs_bary_eq = helio_to_bary(&ecliptic_to_equatorial(&obs_state), &sun_obs);
    let ast_bary_eq = helio_to_bary(&ecliptic_to_equatorial(&aberrated_state), &sun_emit);
    let topo_eq = bary_to_helio(&ast_bary_eq, &obs_bary_eq);

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
}

impl<'a> EphemerisSim<'a> {
    fn new(ephem: &'a Ephemeris, orbit: &Orbit, integrator: &IntegratorConfig) -> Result<Self> {
        use crate::propagate::{
            add_particles_and_variationals, apply_nongrav_scalars, configure_forces,
            install_particle_params,
        };

        let t0 = ephem.mjd_to_assist_time(orbit.epoch);

        // Heliocentric ecliptic → barycentric equatorial ICRF.
        let sun = ephem.get_body_state_array(ffi::ASSIST_BODY_SUN, t0)?;
        let bary_state = helio_to_bary(&ecliptic_to_equatorial(&orbit.state), &sun);

        let mut sim = Simulation::new()?;
        sim.set_t(t0);
        integrator.apply(&mut sim);
        let mut asim = AssistSim::new(sim, ephem)?;

        let has_nongrav = orbit.non_grav.is_some();
        configure_forces(&mut asim, has_nongrav);
        // Ephemeris generation never uses variational equations (n_var = 0).
        add_particles_and_variationals(&mut asim, &bary_state, 0);

        if let Some(ng) = orbit.non_grav.as_ref() {
            apply_nongrav_scalars(&mut asim, ng);
            install_particle_params(&mut asim, ng, false);
        }

        Ok(Self { asim, ephem })
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
        let t_target = self.ephem.mjd_to_assist_time(mjd_tdb);
        self.asim.integrate_or_interpolate(t_target)?;

        let particles = self.asim.sim().particles();
        if particles.is_empty() {
            return Err(Error::Other("No particles after integration".into()));
        }
        let p = &particles[0];
        let bary_eq = [p.x, p.y, p.z, p.vx, p.vy, p.vz];

        let sun_t = self
            .ephem
            .get_body_state_array(ffi::ASSIST_BODY_SUN, t_target)?;
        Ok(equatorial_to_ecliptic(&bary_to_helio(&bary_eq, &sun_t)))
    }
}
