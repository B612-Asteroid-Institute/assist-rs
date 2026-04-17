//! Safe RAII wrappers around REBOUND and ASSIST C objects.

use std::ffi::CString;
use std::path::Path;

use crate::ffi;
use crate::{Error, Result};

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

/// Owned REBOUND simulation. Freed on drop.
///
/// Not `Send`/`Sync` — each thread must create its own simulation.
pub struct Simulation {
    pub(crate) ptr: *mut ffi::reb_simulation,
}

impl Simulation {
    /// Create a new, empty REBOUND simulation.
    pub fn new() -> Result<Self> {
        let ptr = unsafe { ffi::reb_simulation_create() };
        if ptr.is_null() {
            return Err(Error::Other("reb_simulation_create returned null".into()));
        }
        Ok(Self { ptr })
    }

    pub fn t(&self) -> f64 {
        unsafe { ffi::assist_rs_sim_get_t(self.ptr) }
    }
    pub fn set_t(&mut self, t: f64) {
        unsafe { ffi::assist_rs_sim_set_t(self.ptr, t) }
    }

    pub fn dt(&self) -> f64 {
        unsafe { ffi::assist_rs_sim_get_dt(self.ptr) }
    }
    pub fn set_dt(&mut self, dt: f64) {
        unsafe { ffi::assist_rs_sim_set_dt(self.ptr, dt) }
    }

    pub fn n_particles(&self) -> usize {
        unsafe { ffi::assist_rs_sim_get_N(self.ptr) as usize }
    }

    /// Total IAS15 steps (accepted + rejected) since this simulation was
    /// created. Useful for diagnosing adaptive-timestep behavior.
    pub fn steps_done(&self) -> u64 {
        unsafe { ffi::assist_rs_sim_get_steps_done(self.ptr) }
    }

    /// Raw mutable pointer to the underlying REBOUND simulation. Intended for
    /// benchmark probes and low-level FFI workarounds; prefer the safe API
    /// whenever possible.
    #[doc(hidden)]
    pub fn raw_ptr_mut(&mut self) -> *mut ffi::reb_simulation {
        self.ptr
    }

    pub fn n_var(&self) -> i32 {
        unsafe { ffi::assist_rs_sim_get_N_var(self.ptr) }
    }

    pub fn status(&self) -> i32 {
        unsafe { ffi::assist_rs_sim_get_status(self.ptr) }
    }

    pub fn set_exact_finish_time(&mut self, v: bool) {
        unsafe { ffi::assist_rs_sim_set_exact_finish_time(self.ptr, v as i32) }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: ffi::reb_particle) {
        unsafe { ffi::reb_simulation_add(self.ptr, p) }
    }

    /// Add a test particle with given position and velocity (mass = 0).
    pub fn add_test_particle(&mut self, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) {
        let p = ffi::reb_particle {
            x,
            y,
            z,
            vx,
            vy,
            vz,
            m: 0.0,
            ..Default::default()
        };
        self.add_particle(p);
    }

    /// Read-only access to the particle array.
    pub fn particles(&self) -> &[ffi::reb_particle] {
        let n = self.n_particles();
        if n == 0 {
            return &[];
        }
        let ptr = unsafe { ffi::assist_rs_sim_get_particles(self.ptr) };
        if ptr.is_null() {
            return &[];
        }
        unsafe { std::slice::from_raw_parts(ptr, n) }
    }

    /// Integrate to target time. Returns the status code.
    pub fn integrate(&mut self, tmax: f64) -> Result<()> {
        let status = unsafe { ffi::reb_simulation_integrate(self.ptr, tmax) };
        match status {
            ffi::REB_STATUS_SUCCESS | ffi::REB_STATUS_RUNNING => Ok(()),
            ffi::REB_STATUS_NO_PARTICLES => Err(Error::NoParticles),
            ffi::REB_STATUS_ENCOUNTER => Err(Error::CloseEncounter),
            ffi::REB_STATUS_ESCAPE => Err(Error::Escape),
            ffi::REB_STATUS_COLLISION => Err(Error::Collision),
            other => Err(Error::IntegrationFailed(other)),
        }
    }

    /// Add first-order variational particles for a test particle.
    /// Returns the index of the first variational particle.
    pub fn add_variation_1st_order(&mut self, testparticle: i32) -> i32 {
        unsafe { ffi::reb_simulation_add_variation_1st_order(self.ptr, testparticle) }
    }
}

impl Drop for Simulation {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::reb_simulation_free(self.ptr) }
        }
    }
}

// ---------------------------------------------------------------------------
// Ephemeris
// ---------------------------------------------------------------------------

/// Owned ASSIST ephemeris data. Freed on drop.
///
/// Read-only after creation — safe to share across threads.
pub struct Ephemeris {
    pub(crate) ptr: *mut ffi::assist_ephem,
}

impl Ephemeris {
    /// Load ephemeris from SPK files.
    pub fn from_paths(planets: &Path, asteroids: &Path) -> Result<Self> {
        let planets_cstr = path_to_cstring(planets)?;
        let asteroids_cstr = path_to_cstring(asteroids)?;
        let ptr =
            unsafe { ffi::assist_ephem_create(planets_cstr.as_ptr(), asteroids_cstr.as_ptr()) };
        if ptr.is_null() {
            return Err(Error::EphemerisError(
                "assist_ephem_create returned null — check file paths".into(),
            ));
        }
        Ok(Self { ptr })
    }

    /// Raw pointer to the underlying `assist_ephem`. Useful for direct FFI calls.
    ///
    /// Returns a `*const` pointer because `Ephemeris` implements `Sync` on the
    /// premise that the underlying data is read-only after construction. Use
    /// [`pointer::cast_mut`] at the call site if the target FFI signature
    /// requires `*mut`; that cast is the caller's assertion of unique access.
    pub fn as_ptr(&self) -> *const ffi::assist_ephem {
        self.ptr
    }

    /// Reference Julian date for the ephemeris (typically 2451545.0 = J2000.0 TDB).
    pub fn jd_ref(&self) -> f64 {
        unsafe { ffi::assist_rs_ephem_get_jd_ref(self.ptr) }
    }

    /// Override the reference Julian date.
    ///
    /// Requires `&mut self`, which prevents concurrent mutation when the
    /// `Ephemeris` is shared across threads via `Arc`. Must be called before
    /// any `AssistSim` is attached.
    pub fn set_jd_ref(&mut self, jd: f64) {
        unsafe { ffi::assist_rs_ephem_set_jd_ref(self.ptr, jd) }
    }

    /// Speed of light in AU/day.
    pub fn c_au_per_day(&self) -> f64 {
        unsafe { ffi::assist_rs_ephem_get_c_au_per_day(self.ptr) }
    }

    /// Earth equatorial radius in AU.
    pub fn earth_radius_au(&self) -> f64 {
        unsafe { ffi::assist_rs_ephem_get_re_eq(self.ptr) }
    }

    /// Earth/Moon mass ratio.
    pub fn emrat(&self) -> f64 {
        unsafe { ffi::assist_rs_ephem_get_emrat(self.ptr) }
    }

    /// Get a solar system body's state at time `t` (days from `jd_ref`).
    pub fn get_body_state(&self, body_id: i32, t: f64) -> Result<ffi::reb_particle> {
        let mut error: i32 = 0;
        let p = unsafe { ffi::assist_get_particle_with_error(self.ptr, body_id, t, &mut error) };
        if error != 0 {
            return Err(Error::EphemerisError(format!(
                "assist_get_particle failed for body {body_id} at t={t}: error code {error}"
            )));
        }
        Ok(p)
    }
}

impl Drop for Ephemeris {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ffi::assist_ephem_free(self.ptr) }
        }
    }
}

// SAFETY: Ephemeris data is read-only after construction. The underlying
// assist_ephem struct is only mutated during init (assist_ephem_create),
// and all subsequent access (assist_get_particle*) is const-correct in C.
// `set_jd_ref` takes &mut self so Rust's aliasing rules prevent races through
// Arc<Ephemeris>, and `as_ptr` returns *const to forbid back-door mutation.
//
// Note on process-wide concurrency with AssistSim: REBOUND's IAS15 integrator
// and ASSIST force routines use only `static const` tables (no mutable global
// state). The only shared process-global is the SIGINT handler / flag
// (reb_sigint), which all concurrent simulations legitimately share.
// Concurrent AssistSim instances on separate threads are therefore safe.
unsafe impl Send for Ephemeris {}
unsafe impl Sync for Ephemeris {}

// ---------------------------------------------------------------------------
// AssistSim — Simulation with ASSIST forces attached
// ---------------------------------------------------------------------------

/// A REBOUND simulation with ASSIST ephemeris forces attached.
///
/// Owns the simulation. Borrows the ephemeris (caller must keep it alive).
/// ASSIST extras are freed on drop, then the simulation is freed.
pub struct AssistSim {
    pub(crate) sim: Simulation,
    ax: *mut ffi::assist_extras,
    /// Backing storage for ASSIST's `particle_params` pointer. Kept alive
    /// here so its heap buffer lives as long as the simulation.
    particle_params: Option<Vec<f64>>,
}

impl AssistSim {
    /// Create a new ASSIST-powered simulation.
    ///
    /// The `ephem` must outlive this `AssistSim`. ASSIST stores a raw pointer
    /// to the ephemeris data internally.
    pub fn new(mut sim: Simulation, ephem: &Ephemeris) -> Result<Self> {
        let ax = unsafe { ffi::assist_attach(sim.ptr, ephem.ptr) };
        if ax.is_null() {
            return Err(Error::Other("assist_attach returned null".into()));
        }
        // ASSIST sets integrator=IAS15, gravity=NONE, registers force callbacks.
        // Ensure exact finish time is on.
        sim.set_exact_finish_time(true);
        Ok(Self {
            sim,
            ax,
            particle_params: None,
        })
    }

    /// Set the ASSIST force model flags.
    pub fn set_forces(&mut self, flags: i32) {
        unsafe { ffi::assist_rs_extras_set_forces(self.ax, flags) }
    }

    /// Get current force model flags.
    pub fn forces(&self) -> i32 {
        unsafe { ffi::assist_rs_extras_get_forces(self.ax) }
    }

    /// Access the underlying simulation.
    pub fn sim(&self) -> &Simulation {
        &self.sim
    }

    /// Mutable access to the underlying simulation.
    pub fn sim_mut(&mut self) -> &mut Simulation {
        &mut self.sim
    }

    // --- Non-gravitational force model parameters ---

    /// Set the g(r) model exponent α. Default: 1.0.
    pub fn set_alpha(&mut self, v: f64) {
        unsafe { ffi::assist_rs_extras_set_alpha(self.ax, v) }
    }
    pub fn alpha(&self) -> f64 {
        unsafe { ffi::assist_rs_extras_get_alpha(self.ax) }
    }

    /// Set the g(r) model exponent k. Default: 0.0 (pure inverse-power law).
    pub fn set_nk(&mut self, v: f64) {
        unsafe { ffi::assist_rs_extras_set_nk(self.ax, v) }
    }
    pub fn nk(&self) -> f64 {
        unsafe { ffi::assist_rs_extras_get_nk(self.ax) }
    }

    /// Set the g(r) model exponent m. Default: 2.0 (inverse-square).
    pub fn set_nm(&mut self, v: f64) {
        unsafe { ffi::assist_rs_extras_set_nm(self.ax, v) }
    }
    pub fn nm(&self) -> f64 {
        unsafe { ffi::assist_rs_extras_get_nm(self.ax) }
    }

    /// Set the g(r) model exponent n. Default: 5.093 (Marsden-Sekanina water ice).
    pub fn set_nn(&mut self, v: f64) {
        unsafe { ffi::assist_rs_extras_set_nn(self.ax, v) }
    }
    pub fn nn(&self) -> f64 {
        unsafe { ffi::assist_rs_extras_get_nn(self.ax) }
    }

    /// Set the g(r) model scale distance r₀ in AU. Default: 1.0.
    pub fn set_r0(&mut self, v: f64) {
        unsafe { ffi::assist_rs_extras_set_r0(self.ax, v) }
    }
    pub fn r0(&self) -> f64 {
        unsafe { ffi::assist_rs_extras_get_r0(self.ax) }
    }

    /// Install ASSIST's `particle_params` array (3 doubles per particle:
    /// `[A1, A2, A3]`, in `[real | variational]` order).
    ///
    /// Takes ownership of the `Vec`; its heap buffer lives for as long as the
    /// `AssistSim`, matching the lifetime ASSIST requires for the pointer it
    /// stashes internally. Must be called *after* all particles (real +
    /// variational) have been added; `params.len()` must equal
    /// `3 * n_particles`.
    ///
    /// Replacing a previously installed array drops the old storage; the
    /// previous pointer ASSIST held is already overwritten at that point.
    pub(crate) fn set_particle_params(&mut self, mut params: Vec<f64>) {
        let n = self.sim.n_particles();
        assert_eq!(
            params.len(),
            3 * n,
            "particle_params length must equal 3 * n_particles (got {}, expected {})",
            params.len(),
            3 * n
        );
        let ptr = params.as_mut_ptr();
        unsafe { ffi::assist_rs_extras_set_particle_params(self.ax, ptr) }
        self.particle_params = Some(params);
    }

    /// Integrate to target time.
    pub fn integrate(&mut self, tmax: f64) -> Result<()> {
        self.sim.integrate(tmax)
    }

    /// Raw mutable pointer to the underlying REBOUND simulation.
    /// See [`Simulation::raw_ptr_mut`].
    #[doc(hidden)]
    pub fn raw_sim_ptr_mut(&mut self) -> *mut ffi::reb_simulation {
        self.sim.raw_ptr_mut()
    }

    /// Zero IAS15's compensated-summation and predictor state (csx, csv,
    /// csa0, b, e, br, er, g) *in place*, leaving the allocations intact.
    /// Also invalidates the ASSIST ephemeris-lookup cache (sets every slot's
    /// `t` to a sentinel so matched-t comparisons miss on the first
    /// post-reset call). Required between two unrelated orbits integrated on
    /// the same simulation — otherwise stale b/e predictor state seeds the
    /// new orbit's first step and causes extra corrector iterations, and a
    /// populated ephem cache causes per-lookup LRU work that adds up across
    /// ~7000 lookups per 30-day integrate (≈190 µs regression).
    ///
    /// Cheaper than [`reb_integrator_ias15_reset`] (no free/malloc), and
    /// faster in practice: a pool-style benchmark with this helper matches
    /// or beats the unpooled free-function path.
    pub(crate) fn reset_integrator(&mut self) {
        unsafe {
            ffi::assist_rs_ias15_zero_state(self.sim.ptr);
            ffi::assist_rs_ephem_cache_reset(self.ax);
        }
    }

    /// Rewrite the first three slots of the installed `particle_params`
    /// array (the real test particle's A1, A2, A3) without reallocating.
    /// Does nothing if no params array has been installed yet.
    ///
    /// The variational-particle parameter columns (indices 3 onward) are
    /// orbit-invariant IC perturbations (identity for parameter
    /// variationals, zero for state variationals) and are left untouched.
    pub(crate) fn update_nongrav_coeffs(&mut self, a1: f64, a2: f64, a3: f64) {
        if let Some(params) = self.particle_params.as_mut() {
            params[0] = a1;
            params[1] = a2;
            params[2] = a3;
        }
    }
}

impl Drop for AssistSim {
    fn drop(&mut self) {
        if !self.ax.is_null() {
            // Detach ASSIST first, then assist_free, then sim drops automatically.
            unsafe {
                ffi::assist_detach(self.sim.ptr, self.ax);
                ffi::assist_free(self.ax);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn path_to_cstring(path: &Path) -> Result<CString> {
    let s = path.to_str().ok_or_else(|| {
        Error::Other(format!(
            "path contains non-UTF8 characters: {}",
            path.display()
        ))
    })?;
    CString::new(s).map_err(|e| Error::Other(format!("path contains null byte: {e}")))
}
