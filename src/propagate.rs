//! `assist_propagate_single` — N-body propagation of a test particle.
//!
//! Two entry points:
//!
//! * [`assist_propagate_single`] — one-shot free function. Creates a fresh
//!   `AssistSim`, propagates, and tears down on return.
//! * [`PropagatorPool`] — stateful container that reuses the underlying
//!   `AssistSim` (and its REBOUND simulation, ASSIST extras, and IAS15
//!   allocations) across many orbits with the same force-model +
//!   variational-dimension configuration.

use crate::assist_data::AssistData;
use crate::coordinates::{ecliptic_to_equatorial, equatorial_to_ecliptic, rotate_matrix_eq_to_ecl};
use crate::ffi;
use crate::orbit::{NonGravParams, Orbit};
use crate::wrappers::{AssistSim, Ephemeris, Simulation};
use crate::{Error, Result};

/// Result of propagating to a single epoch.
#[derive(Debug, Clone)]
pub struct PropagatedState {
    /// Heliocentric ecliptic J2000 state [x, y, z, vx, vy, vz] (AU, AU/day).
    pub state: [f64; 6],
    /// Epoch (MJD TDB).
    pub epoch: f64,
    /// 6×6 state transition matrix Φ(t, t₀) = ∂x(t)/∂x₀, row-major, in
    /// heliocentric ecliptic J2000. Populated when `compute_stm` is true.
    pub stm: Option<[[f64; 6]; 6]>,
    /// 6×3 partials ∂x(t)/∂(A1, A2, A3), row-major (6 state rows, 3 parameter
    /// columns), in heliocentric ecliptic J2000. Populated when `compute_stm`
    /// is true *and* the orbit carries non-gravitational parameters.
    /// Columns are ordered A1, A2, A3.
    pub nongrav_partials: Option<[[f64; 3]; 6]>,
}

impl PropagatedState {
    /// Linearly propagate a 6×6 initial-state covariance to this epoch:
    ///
    /// ```text
    ///   P(t) = Φ · P₀ · Φᵀ
    /// ```
    ///
    /// where `Φ` is `self.stm`. Returns `None` when `stm` is not populated
    /// (i.e., the orbit was propagated with `compute_stm = false`).
    ///
    /// The input and output covariances are in heliocentric ecliptic J2000,
    /// matching the frame of `state` and `stm`.
    pub fn propagate_covariance(&self, p0: &[[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
        self.stm.as_ref().map(|stm| covariance_6x6(stm, p0))
    }

    /// Linearly propagate a 9×9 initial covariance over
    /// `(x, y, z, vx, vy, vz, A1, A2, A3)` to the 6×6 state covariance at this
    /// epoch:
    ///
    /// ```text
    ///   J = [ Φ | G ]         (6×9, Φ = stm, G = nongrav_partials)
    ///   P(t) = J · P₀ · Jᵀ
    /// ```
    ///
    /// Returns `None` unless both `stm` and `nongrav_partials` are populated
    /// (i.e., the orbit has non-gravitational parameters and was propagated
    /// with `compute_stm = true`).
    pub fn propagate_covariance_with_nongrav(&self, p0: &[[f64; 9]; 9]) -> Option<[[f64; 6]; 6]> {
        let stm = self.stm.as_ref()?;
        let ng = self.nongrav_partials.as_ref()?;
        Some(covariance_9x9(stm, ng, p0))
    }
}

/// P = Φ · P₀ · Φᵀ, where Φ and P₀ are 6×6 (row-major).
fn covariance_6x6(stm: &[[f64; 6]; 6], p0: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut tmp = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut s = 0.0;
            for k in 0..6 {
                s += stm[i][k] * p0[k][j];
            }
            tmp[i][j] = s;
        }
    }
    let mut out = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut s = 0.0;
            for k in 0..6 {
                // stm[j][k] is (Φᵀ)[k][j]
                s += tmp[i][k] * stm[j][k];
            }
            out[i][j] = s;
        }
    }
    out
}

/// P = J · P₀ · Jᵀ, where J = [stm | ng] is 6×9 (stm 6×6, ng 6×3) and P₀ is 9×9.
fn covariance_9x9(stm: &[[f64; 6]; 6], ng: &[[f64; 3]; 6], p0: &[[f64; 9]; 9]) -> [[f64; 6]; 6] {
    let j = |i: usize, k: usize| -> f64 { if k < 6 { stm[i][k] } else { ng[i][k - 6] } };
    // tmp = J · P₀  (6×9)
    let mut tmp = [[0.0f64; 9]; 6];
    for i in 0..6 {
        for col in 0..9 {
            let mut s = 0.0;
            for k in 0..9 {
                s += j(i, k) * p0[k][col];
            }
            tmp[i][col] = s;
        }
    }
    // out = tmp · Jᵀ  (6×6)
    let mut out = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for col in 0..6 {
            let mut s = 0.0;
            for k in 0..9 {
                s += tmp[i][k] * j(col, k);
            }
            out[i][col] = s;
        }
    }
    out
}

/// Propagate a test particle from an initial heliocentric ecliptic J2000 orbit.
///
/// One-shot convenience: creates a fresh simulation, integrates, and drops on
/// return. For repeated propagation of many orbits that share the same force
/// model and variational-equation setup, prefer [`PropagatorPool`], which
/// amortizes the ~25 µs REBOUND/ASSIST setup across calls.
///
/// # Arguments
/// - `data`: bundle holding the ASSIST ephemeris (observatory table unused
///   here — propagation never needs it).
/// - `orbit`: initial orbit (state, epoch, optional non-grav params).
/// - `target_epochs`: sorted slice of target epochs (MJD TDB).
/// - `compute_stm`: whether to compute the state transition matrix via variational equations.
///
/// # Returns
/// One `PropagatedState` per target epoch, in the same order.
pub fn assist_propagate_single(
    data: &AssistData,
    orbit: &Orbit,
    target_epochs: &[f64],
    compute_stm: bool,
) -> Result<Vec<PropagatedState>> {
    if target_epochs.is_empty() {
        return Ok(vec![]);
    }

    let ephem = &data.ephem;
    let jd_ref = ephem.jd_ref();
    let t0 = mjd_to_assist_time(orbit.epoch, jd_ref);

    let has_nongrav = orbit.non_grav.is_some();
    let want_nongrav_partials = compute_stm && has_nongrav;

    let mut sim = Simulation::new()?;
    sim.set_t(t0);
    let mut asim = AssistSim::new(sim, ephem)?;
    configure_forces(&mut asim, has_nongrav);

    let bary_state = ecl_orbit_to_bary_eq(&orbit.state, ephem, t0)?;
    let n_var = variational_count(compute_stm, has_nongrav);
    add_particles_and_variationals(&mut asim, &bary_state, n_var);

    if let Some(ng) = orbit.non_grav.as_ref() {
        apply_nongrav_scalars(&mut asim, ng);
        install_particle_params(&mut asim, ng, want_nongrav_partials);
    }

    run_integration(
        &mut asim,
        ephem,
        target_epochs,
        compute_stm,
        want_nongrav_partials,
    )
}

/// Propagate many orbits to a shared set of target epochs.
///
/// Shape: the returned `Vec<Vec<PropagatedState>>` is indexed
/// `[orbit_index][target_epoch_index]`.
///
/// # Parallelism
///
/// `num_threads` controls how work is distributed:
///
/// * `None` — use rayon's default global pool (typically one worker per
///   CPU core). Best for single-process workloads that want maximum
///   throughput.
/// * `Some(1)` — serial loop, no rayon scheduling overhead. Equivalent
///   to `orbits.iter().map(|o| assist_propagate_single(...)).collect()`.
/// * `Some(n)` for `n > 1` — build a dedicated rayon thread pool with
///   exactly `n` workers and run the batch inside it. Useful when the
///   caller already parallelizes at a higher level and wants to bound
///   the per-call concurrency.
///
/// `Some(0)` returns `Error::Other` — use `None` for "default".
///
/// Each orbit still pays its own [`assist_propagate_single`] setup cost (~1 µs of
/// `reb_simulation_create` + `assist_attach`); the parallel win is the
/// embarrassingly-parallel map across orbits. `Ephemeris` is `Send + Sync`
/// so one ephemeris serves all worker threads.
///
/// Without the `parallel` cargo feature, `num_threads` is ignored and the
/// function always runs serially.
///
/// # Errors
/// Returns the first error encountered across any orbit's propagation.
pub fn assist_propagate(
    data: &AssistData,
    orbits: &[Orbit],
    target_epochs: &[f64],
    compute_stm: bool,
    num_threads: Option<usize>,
) -> Result<Vec<Vec<PropagatedState>>> {
    let op = |orbit: &Orbit| assist_propagate_single(data, orbit, target_epochs, compute_stm);
    map_with_threads(orbits, num_threads, op)
}

/// Dispatch a fallible per-item function across `items`:
/// `None` = rayon's global pool, `Some(1)` = serial, `Some(n)` = custom
/// pool with `n` workers. `Some(0)` is an error.
///
/// Factored out so [`assist_propagate`] and
/// [`crate::ephemeris::assist_generate_ephemeris_single`] share one dispatch
/// implementation. `Send + Sync` bounds are required on `T` and the
/// closure's captures because rayon moves work across threads.
pub(crate) fn map_with_threads<T, U, F>(
    items: &[T],
    num_threads: Option<usize>,
    f: F,
) -> Result<Vec<U>>
where
    T: Sync,
    U: Send,
    F: Fn(&T) -> Result<U> + Send + Sync,
{
    if num_threads == Some(0) {
        return Err(Error::Other(
            "num_threads = Some(0) is not valid; use None for the default pool".into(),
        ));
    }
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        match num_threads {
            None => items.par_iter().map(&f).collect(),
            Some(1) => items.iter().map(&f).collect(),
            Some(n) => {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(n)
                    .build()
                    .map_err(|e| Error::Other(format!("rayon pool build failed: {e}")))?;
                pool.install(|| items.par_iter().map(&f).collect())
            }
        }
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = num_threads;
        items.iter().map(&f).collect()
    }
}

// ─── PropagatorPool: reusable simulation for many-orbit workloads ───────────

/// Configuration locked at [`PropagatorPool`] construction. A single pool can
/// only serve orbits whose variational dimension matches this config:
/// creating a `has_nongrav = false` pool and then passing an orbit with
/// non-gravitational parameters is an error, and vice versa.
#[derive(Debug, Clone, Copy)]
pub struct PropagatorConfig {
    /// Compute the 6×6 state transition matrix via variational equations.
    pub compute_stm: bool,
    /// Whether orbits served by this pool carry non-gravitational parameters.
    /// When true, an extra 3 variational particles are added to yield the
    /// 6×3 `∂x/∂A` matrix in each result (when `compute_stm` is also true).
    pub has_nongrav: bool,
}

impl PropagatorConfig {
    /// Gravity-only orbits, no STM — the cheapest configuration.
    pub fn gravity_only() -> Self {
        Self {
            compute_stm: false,
            has_nongrav: false,
        }
    }

    /// Gravity-only orbits with STM (6×6).
    pub fn gravity_with_stm() -> Self {
        Self {
            compute_stm: true,
            has_nongrav: false,
        }
    }

    /// Non-grav orbits with STM (6×6) and parameter partials (6×3).
    pub fn nongrav_with_stm() -> Self {
        Self {
            compute_stm: true,
            has_nongrav: true,
        }
    }
}

/// Reusable propagator for many orbits sharing the same force model and
/// variational-equation configuration.
///
/// Construction attaches ASSIST to a REBOUND simulation, pre-adds the real
/// test particle + up to 9 variational particles, and (when `has_nongrav`)
/// installs a zeroed `particle_params` array of the right size. Each call
/// to [`propagate`] rewrites particle state in place, resets the IAS15
/// integrator's scratch arrays, and re-runs the integration loop —
/// avoiding the ~25 µs per-orbit setup cost that [`assist_propagate_single`] pays.
///
/// The pool is stateful and not thread-safe; clone [`Ephemeris`] (which is
/// `Send + Sync`) and build one pool per worker thread.
///
/// ```no_run
/// # use assist_rs::{AssistData, Orbit, propagate::{PropagatorPool, PropagatorConfig}};
/// # fn run(data: &AssistData, orbits: &[Orbit], targets: &[f64]) -> Result<(), Box<dyn std::error::Error>> {
/// let mut pool = PropagatorPool::new(data, PropagatorConfig::gravity_with_stm())?;
/// for orbit in orbits {
///     let result = pool.propagate(orbit, targets)?;
///     // ... use result[0].state, result[0].stm ...
/// }
/// # Ok(()) }
/// ```
///
/// [`propagate`]: PropagatorPool::propagate
pub struct PropagatorPool<'a> {
    asim: AssistSim,
    ephem: &'a Ephemeris,
    jd_ref: f64,
    config: PropagatorConfig,
    /// Cached variational count (0, 6, or 9) — derived from config.
    n_var: usize,
}

impl<'a> PropagatorPool<'a> {
    /// Construct a pool for a given force/variational configuration.
    ///
    /// Placeholder state is installed for the real test particle (zeros) and
    /// the variational particles (unit perturbations / zeros as required);
    /// these are overwritten by the first [`propagate`] call. The pool
    /// borrows `data` for the ephemeris; `data.observatory` is not used by
    /// the propagator.
    ///
    /// [`propagate`]: PropagatorPool::propagate
    pub fn new(data: &'a AssistData, config: PropagatorConfig) -> Result<Self> {
        let ephem = &data.ephem;
        let jd_ref = ephem.jd_ref();
        let mut sim = Simulation::new()?;
        sim.set_t(0.0);
        let mut asim = AssistSim::new(sim, ephem)?;
        configure_forces(&mut asim, config.has_nongrav);

        let placeholder = [0.0f64; 6];
        let n_var = variational_count(config.compute_stm, config.has_nongrav);
        add_particles_and_variationals(&mut asim, &placeholder, n_var);

        if config.has_nongrav {
            // Install a zeroed params array of the correct total length; each
            // propagate() call updates the first three slots in place.
            let placeholder_ng = NonGravParams::new(0.0, 0.0, 0.0);
            install_particle_params(
                &mut asim,
                &placeholder_ng,
                config.compute_stm && config.has_nongrav,
            );
        }

        Ok(Self {
            asim,
            ephem,
            jd_ref,
            config,
            n_var,
        })
    }

    /// The configuration this pool was built with.
    pub fn config(&self) -> PropagatorConfig {
        self.config
    }

    /// Cumulative IAS15 step count across all propagations done by this
    /// pool. Useful for diagnostics; no functional impact.
    pub fn steps_done(&self) -> u64 {
        self.asim.sim().steps_done()
    }

    /// Propagate `orbit` to `target_epochs`, reusing the internal simulation.
    ///
    /// The orbit's non-grav flag must match the pool's `config.has_nongrav`
    /// (otherwise `Error::Other` is returned — swapping variational
    /// dimensions mid-pool would require re-adding particles and is out of
    /// scope).
    pub fn propagate(
        &mut self,
        orbit: &Orbit,
        target_epochs: &[f64],
    ) -> Result<Vec<PropagatedState>> {
        if target_epochs.is_empty() {
            return Ok(vec![]);
        }

        if orbit.non_grav.is_some() != self.config.has_nongrav {
            return Err(Error::Other(format!(
                "orbit non-grav flag ({}) does not match pool config \
                 has_nongrav ({}) — rebuild the pool with a matching config",
                orbit.non_grav.is_some(),
                self.config.has_nongrav,
            )));
        }

        let want_nongrav_partials = self.config.compute_stm && self.config.has_nongrav;
        let t0 = mjd_to_assist_time(orbit.epoch, self.jd_ref);

        // 1) Clear IAS15 scratch arrays so no compensated-sum or predictor
        //    state leaks from the previous orbit. In-place memset, no
        //    allocator churn.
        self.asim.reset_integrator();

        // 2) Re-seed epoch and timestep. `dt = 0` would lock IAS15 into a
        //    zero-step infinite loop (copysign preserves zero); use REBOUND's
        //    default of 0.001 days so IAS15 has a fresh, adaptive starting
        //    point independent of the previous orbit's final dt.
        self.asim.sim_mut().set_t(t0);
        self.asim.sim_mut().set_dt(0.001);

        // 3) Overwrite real particle + re-init variationals to their unit-
        //    perturbation ICs (STM identity basis, param variationals zero).
        let bary_state = ecl_orbit_to_bary_eq(&orbit.state, self.ephem, t0)?;
        overwrite_particles(&mut self.asim, &bary_state, self.n_var);

        // 4) Non-grav scalar params + in-place A1/A2/A3 update (no realloc).
        if let Some(ng) = orbit.non_grav.as_ref() {
            apply_nongrav_scalars(&mut self.asim, ng);
            self.asim.update_nongrav_coeffs(ng.a1, ng.a2, ng.a3);
        }

        run_integration(
            &mut self.asim,
            self.ephem,
            target_epochs,
            self.config.compute_stm,
            want_nongrav_partials,
        )
    }
}

// ─── Private helpers shared by assist_propagate_single and PropagatorPool ──────────

/// Heliocentric ecliptic J2000 → barycentric equatorial ICRF at ASSIST time `t`.
fn ecl_orbit_to_bary_eq(ecl_state: &[f64; 6], ephem: &Ephemeris, t: f64) -> Result<[f64; 6]> {
    let eq_state = ecliptic_to_equatorial(ecl_state);
    let sun = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t)?;
    Ok([
        eq_state[0] + sun.x,
        eq_state[1] + sun.y,
        eq_state[2] + sun.z,
        eq_state[3] + sun.vx,
        eq_state[4] + sun.vy,
        eq_state[5] + sun.vz,
    ])
}

/// Set ASSIST force flags (default + optional non-grav).
fn configure_forces(asim: &mut AssistSim, has_nongrav: bool) {
    let mut forces = ffi::ASSIST_FORCES_DEFAULT;
    if has_nongrav {
        forces |= ffi::ASSIST_FORCE_NON_GRAVITATIONAL;
    }
    asim.set_forces(forces);
}

/// Variational particle count: 0, 6, or 9.
fn variational_count(compute_stm: bool, has_nongrav: bool) -> usize {
    match (compute_stm, has_nongrav) {
        (false, _) => 0,
        (true, false) => 6,
        (true, true) => 9,
    }
}

/// Apply the non-grav g(r) scalar parameters (α, nk, nm, nn, r0).
fn apply_nongrav_scalars(asim: &mut AssistSim, ng: &NonGravParams) {
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

/// Add the real test particle + `n_var` variational particles (0, 6, or 9)
/// to an empty simulation. Variational particles 0..5 get unit state
/// perturbations; 6..8 (when present) get zero state (the correct IC for
/// parameter sensitivities).
fn add_particles_and_variationals(asim: &mut AssistSim, bary_state: &[f64; 6], n_var: usize) {
    asim.sim_mut().add_test_particle(
        bary_state[0],
        bary_state[1],
        bary_state[2],
        bary_state[3],
        bary_state[4],
        bary_state[5],
    );
    for _ in 0..n_var {
        asim.sim_mut().add_variation_1st_order(0);
    }
    init_variational_state_perturbations(asim, n_var);
}

/// Overwrite particle[0] with `bary_state` and re-initialize variational
/// particles to their unit-perturbation / zero ICs. Used by
/// `PropagatorPool::propagate` between orbits.
fn overwrite_particles(asim: &mut AssistSim, bary_state: &[f64; 6], n_var: usize) {
    unsafe {
        let ptr = ffi::assist_rs_sim_get_particles(asim.sim().ptr);
        let p = &mut *ptr;
        p.x = bary_state[0];
        p.y = bary_state[1];
        p.z = bary_state[2];
        p.vx = bary_state[3];
        p.vy = bary_state[4];
        p.vz = bary_state[5];
        // ax/ay/az are set by the force evaluator on the first step.
    }
    init_variational_state_perturbations(asim, n_var);
}

/// Write the six identity-column state perturbations to particles[1..=6]
/// and zero the parameter variational particles (7..=9 if present).
///
/// `n_var` must be 0, 6, or 9; the caller is responsible for having
/// added exactly that many variational particles beforehand. When
/// `n_var == 0` this is a no-op (nothing to initialize).
fn init_variational_state_perturbations(asim: &mut AssistSim, n_var: usize) {
    debug_assert!(matches!(n_var, 0 | 6 | 9));
    if n_var == 0 {
        return;
    }
    unsafe {
        let ptr = ffi::assist_rs_sim_get_particles(asim.sim().ptr);
        for i in 0..n_var {
            *ptr.add(1 + i) = ffi::reb_particle::default();
        }
        // Six state-sensitivity columns: ∂/∂x₀, ∂/∂y₀, ... ∂/∂vz₀.
        for d in 0..6 {
            let p = &mut *ptr.add(1 + d);
            match d {
                0 => p.x = 1.0,
                1 => p.y = 1.0,
                2 => p.z = 1.0,
                3 => p.vx = 1.0,
                4 => p.vy = 1.0,
                5 => p.vz = 1.0,
                _ => unreachable!(),
            }
        }
        // Parameter-sensitivity columns (if present) remain zero — the
        // correct IC for ∂x/∂A evaluated at the orbit's epoch.
    }
}

/// Build and install the `particle_params` array of length `3 * n_total`.
///
/// Layout: `[A1, A2, A3 | 0…0 for 6 state variationals | 100, 010, 001 for 3 param variationals]`
/// where the last block only exists when `want_nongrav_partials` is true.
fn install_particle_params(asim: &mut AssistSim, ng: &NonGravParams, want_nongrav_partials: bool) {
    let n_total = asim.sim().n_particles();
    let mut params = vec![0.0f64; 3 * n_total];
    params[0] = ng.a1;
    params[1] = ng.a2;
    params[2] = ng.a3;
    if want_nongrav_partials {
        let n_real = 1usize;
        for k in 0..3 {
            params[3 * (n_real + 6 + k) + k] = 1.0;
        }
    }
    asim.set_particle_params(params);
}

/// Run the integration loop: step to each target epoch in order, then
/// extract helio-ecliptic state + optional STM / non-grav partials.
fn run_integration(
    asim: &mut AssistSim,
    ephem: &Ephemeris,
    target_epochs: &[f64],
    compute_stm: bool,
    want_nongrav_partials: bool,
) -> Result<Vec<PropagatedState>> {
    let jd_ref = ephem.jd_ref();
    let mut results = Vec::with_capacity(target_epochs.len());

    for &target_mjd in target_epochs {
        let t_target = mjd_to_assist_time(target_mjd, jd_ref);
        asim.integrate(t_target)?;

        let particles = asim.sim().particles();
        if particles.is_empty() {
            return Err(Error::Other("No particles after integration".into()));
        }
        let p = &particles[0];
        let bary_eq = [p.x, p.y, p.z, p.vx, p.vy, p.vz];

        // Convert barycentric equatorial → heliocentric ecliptic
        let sun_t = ephem.get_body_state(ffi::ASSIST_BODY_SUN, t_target)?;
        let helio_eq = [
            bary_eq[0] - sun_t.x,
            bary_eq[1] - sun_t.y,
            bary_eq[2] - sun_t.z,
            bary_eq[3] - sun_t.vx,
            bary_eq[4] - sun_t.vy,
            bary_eq[5] - sun_t.vz,
        ];
        let helio_ecl = equatorial_to_ecliptic(&helio_eq);

        let (stm, nongrav_partials) = if compute_stm {
            extract_stm_and_partials(particles, want_nongrav_partials)
        } else {
            (None, None)
        };

        results.push(PropagatedState {
            state: helio_ecl,
            epoch: target_mjd,
            stm,
            nongrav_partials,
        });
    }

    Ok(results)
}

/// Read variational particles out of `particles[1..]` and rotate the
/// columns from barycentric equatorial into heliocentric ecliptic.
fn extract_stm_and_partials(
    particles: &[ffi::reb_particle],
    want_nongrav_partials: bool,
) -> (Option<[[f64; 6]; 6]>, Option<[[f64; 3]; 6]>) {
    let n_real = 1usize;
    let mut stm_eq = [[0.0f64; 6]; 6];
    for (d, vp) in particles[n_real..n_real + 6].iter().enumerate() {
        stm_eq[0][d] = vp.x;
        stm_eq[1][d] = vp.y;
        stm_eq[2][d] = vp.z;
        stm_eq[3][d] = vp.vx;
        stm_eq[4][d] = vp.vy;
        stm_eq[5][d] = vp.vz;
    }
    let stm = Some(rotate_matrix_eq_to_ecl(&stm_eq));

    let nongrav = if want_nongrav_partials {
        let mut ng_eq = [[0.0f64; 3]; 6];
        for (k, vp) in particles[n_real + 6..n_real + 9].iter().enumerate() {
            ng_eq[0][k] = vp.x;
            ng_eq[1][k] = vp.y;
            ng_eq[2][k] = vp.z;
            ng_eq[3][k] = vp.vx;
            ng_eq[4][k] = vp.vy;
            ng_eq[5][k] = vp.vz;
        }
        let mut ng_ecl = [[0.0f64; 3]; 6];
        for k in 0..3 {
            let col = [
                ng_eq[0][k],
                ng_eq[1][k],
                ng_eq[2][k],
                ng_eq[3][k],
                ng_eq[4][k],
                ng_eq[5][k],
            ];
            let rotated = equatorial_to_ecliptic(&col);
            for r in 0..6 {
                ng_ecl[r][k] = rotated[r];
            }
        }
        Some(ng_ecl)
    } else {
        None
    };
    (stm, nongrav)
}

/// Convert MJD TDB to ASSIST time (days from jd_ref).
///
/// ASSIST time = JD - jd_ref, where JD = MJD + 2400000.5.
fn mjd_to_assist_time(mjd_tdb: f64, jd_ref: f64) -> f64 {
    (mjd_tdb + 2_400_000.5) - jd_ref
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity6() -> [[f64; 6]; 6] {
        let mut m = [[0.0; 6]; 6];
        for i in 0..6 {
            m[i][i] = 1.0;
        }
        m
    }

    fn identity9() -> [[f64; 9]; 9] {
        let mut m = [[0.0; 9]; 9];
        for i in 0..9 {
            m[i][i] = 1.0;
        }
        m
    }

    /// A non-trivial "STM-like" 6×6: entries chosen so both rotation and
    /// scaling happen, and it's not symmetric.
    fn sample_stm() -> [[f64; 6]; 6] {
        [
            [1.0, 0.0, 0.0, 30.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 30.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 30.0],
            [0.0001, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0001, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0001, 0.0, 0.0, 1.0],
        ]
    }

    fn max_abs_diff_6x6(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> f64 {
        let mut m: f64 = 0.0;
        for i in 0..6 {
            for j in 0..6 {
                m = m.max((a[i][j] - b[i][j]).abs());
            }
        }
        m
    }

    #[test]
    fn covariance_6x6_identity_p0() {
        // With P₀ = I: P(t) = Φ · Φᵀ (manually checkable).
        let stm = sample_stm();
        let got = covariance_6x6(&stm, &identity6());
        // Build the reference Φ · Φᵀ by hand.
        let mut want = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let mut s = 0.0;
                for k in 0..6 {
                    s += stm[i][k] * stm[j][k];
                }
                want[i][j] = s;
            }
        }
        assert!(max_abs_diff_6x6(&got, &want) < 1e-14);
        // Result is symmetric.
        for i in 0..6 {
            for j in 0..6 {
                assert!((got[i][j] - got[j][i]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn covariance_6x6_zero_stm_is_zero() {
        let stm = [[0.0f64; 6]; 6];
        let p0 = identity6();
        let got = covariance_6x6(&stm, &p0);
        for i in 0..6 {
            for j in 0..6 {
                assert_eq!(got[i][j], 0.0);
            }
        }
    }

    #[test]
    fn propagate_covariance_method_wraps_helper() {
        let state = PropagatedState {
            state: [0.0; 6],
            epoch: 0.0,
            stm: Some(sample_stm()),
            nongrav_partials: None,
        };
        let p0 = identity6();
        let via_method = state.propagate_covariance(&p0).unwrap();
        let via_helper = covariance_6x6(state.stm.as_ref().unwrap(), &p0);
        assert!(max_abs_diff_6x6(&via_method, &via_helper) < 1e-14);
    }

    #[test]
    fn propagate_covariance_none_when_no_stm() {
        let state = PropagatedState {
            state: [0.0; 6],
            epoch: 0.0,
            stm: None,
            nongrav_partials: None,
        };
        assert!(state.propagate_covariance(&identity6()).is_none());
    }

    #[test]
    fn covariance_9x9_reduces_to_6x6_when_nongrav_block_is_zero() {
        // If the 9×9 P₀ has only the upper-left 6×6 block populated and the
        // A-parameter block is zero, the result must equal the 6×6
        // propagation.
        let stm = sample_stm();
        let ng: [[f64; 3]; 6] = [
            [1e-3, 2e-3, -1e-3],
            [0.0, 1e-3, 0.0],
            [0.0, 0.0, 1e-3],
            [1e-5, 0.0, 0.0],
            [0.0, 1e-5, 0.0],
            [0.0, 0.0, 1e-5],
        ];

        // Build a P₀ with P_xx = diag(0.1, 0.1, 0.1, 1e-4, 1e-4, 1e-4) and
        // zero everywhere else.
        let p_xx = {
            let mut m = [[0.0; 6]; 6];
            for (i, v) in [0.1, 0.1, 0.1, 1e-4, 1e-4, 1e-4].iter().enumerate() {
                m[i][i] = *v;
            }
            m
        };
        let mut p0_9 = [[0.0; 9]; 9];
        for i in 0..6 {
            for j in 0..6 {
                p0_9[i][j] = p_xx[i][j];
            }
        }

        let got = covariance_9x9(&stm, &ng, &p0_9);
        let want = covariance_6x6(&stm, &p_xx);
        assert!(
            max_abs_diff_6x6(&got, &want) < 1e-14,
            "9×9 with zero A-block should match 6×6 path; diff={:.3e}",
            max_abs_diff_6x6(&got, &want)
        );
    }

    #[test]
    fn covariance_9x9_picks_up_pure_nongrav_covariance() {
        // If only the lower-right 3×3 A-block of P₀ is non-zero, the result
        // must equal G · P_AA · Gᵀ.
        let stm = sample_stm();
        let ng: [[f64; 3]; 6] = [
            [10.0, 5.0, 2.0],
            [3.0, 8.0, 1.0],
            [1.0, 2.0, 6.0],
            [0.1, 0.05, 0.02],
            [0.03, 0.08, 0.01],
            [0.01, 0.02, 0.06],
        ];
        let p_aa = [[1.0, 0.2, 0.1], [0.2, 1.0, 0.3], [0.1, 0.3, 1.0]];
        let mut p0_9 = [[0.0; 9]; 9];
        for i in 0..3 {
            for j in 0..3 {
                p0_9[6 + i][6 + j] = p_aa[i][j];
            }
        }

        let got = covariance_9x9(&stm, &ng, &p0_9);

        // Expected: G · P_AA · Gᵀ, computed the long way.
        let mut tmp = [[0.0; 3]; 6]; // G · P_AA
        for i in 0..6 {
            for j in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += ng[i][k] * p_aa[k][j];
                }
                tmp[i][j] = s;
            }
        }
        let mut want = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += tmp[i][k] * ng[j][k];
                }
                want[i][j] = s;
            }
        }
        assert!(max_abs_diff_6x6(&got, &want) < 1e-12);
    }

    #[test]
    fn covariance_9x9_identity_p0_includes_both_blocks() {
        // With P₀ = I_9 the result must equal the 6×6 "identity state P₀"
        // result PLUS the pure-nongrav contribution. Verifies cross-term
        // handling doesn't double- or drop-count.
        let stm = sample_stm();
        let ng: [[f64; 3]; 6] = [
            [1.0, 0.5, 0.2],
            [0.3, 0.8, 0.1],
            [0.1, 0.2, 0.6],
            [0.01, 0.005, 0.002],
            [0.003, 0.008, 0.001],
            [0.001, 0.002, 0.006],
        ];
        let got = covariance_9x9(&stm, &ng, &identity9());
        let state_part = covariance_6x6(&stm, &identity6());
        // Pure-nongrav part: G · G^T (since P_AA = I_3)
        let mut ng_part = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += ng[i][k] * ng[j][k];
                }
                ng_part[i][j] = s;
            }
        }
        let mut want = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                want[i][j] = state_part[i][j] + ng_part[i][j];
            }
        }
        assert!(max_abs_diff_6x6(&got, &want) < 1e-12);
    }
}
