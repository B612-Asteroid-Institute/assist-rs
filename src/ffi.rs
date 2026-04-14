//! Raw FFI bindings to REBOUND and ASSIST C libraries.
//!
//! `reb_simulation`, `assist_ephem`, and `assist_extras` are treated as opaque
//! types. Field access goes through thin C helpers compiled in `helpers.c`.

use libc::{c_char, c_double, c_int, c_uint, c_void};

// ---------------------------------------------------------------------------
// reb_particle — the only REBOUND struct we reproduce in full.
// Passed by value in reb_simulation_add and returned by assist_get_particle.
// ---------------------------------------------------------------------------

/// REBOUND particle: position, velocity, acceleration, mass, radius, hash.
///
/// Layout matches `struct reb_particle` in `rebound.h` on LP64 platforms.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct reb_particle {
    pub x: c_double,
    pub y: c_double,
    pub z: c_double,
    pub vx: c_double,
    pub vy: c_double,
    pub vz: c_double,
    pub ax: c_double,
    pub ay: c_double,
    pub az: c_double,
    pub m: c_double,
    pub r: c_double,
    pub last_collision: c_double,
    pub c: *mut c_void, // reb_treecell*
    pub hash: u32,
    pub ap: *mut c_void,
    pub sim: *mut c_void, // reb_simulation*
}

impl Default for reb_particle {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0,
            m: 0.0,
            r: 0.0,
            last_collision: 0.0,
            c: std::ptr::null_mut(),
            hash: 0,
            ap: std::ptr::null_mut(),
            sim: std::ptr::null_mut(),
        }
    }
}

// ---------------------------------------------------------------------------
// Opaque types — we never access their fields from Rust directly.
// ---------------------------------------------------------------------------

/// Opaque REBOUND simulation.
#[repr(C)]
pub struct reb_simulation {
    _opaque: [u8; 0],
}

/// Opaque ASSIST ephemeris data.
#[repr(C)]
pub struct assist_ephem {
    _opaque: [u8; 0],
}

/// Opaque ASSIST extras (attaches ASSIST forces to a REBOUND simulation).
#[repr(C)]
pub struct assist_extras {
    _opaque: [u8; 0],
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

// REBOUND integrator IDs (from enum in reb_simulation)
pub const REB_INTEGRATOR_IAS15: c_int = 0;
pub const REB_INTEGRATOR_WHFAST: c_int = 1;

// REBOUND gravity modes
pub const REB_GRAVITY_NONE: c_int = 0;
pub const REB_GRAVITY_BASIC: c_int = 1;

// REBOUND status codes
pub const REB_STATUS_SUCCESS: c_int = 0;
pub const REB_STATUS_RUNNING: c_int = -1;
pub const REB_STATUS_GENERIC_ERROR: c_int = 1;
pub const REB_STATUS_NO_PARTICLES: c_int = 2;
pub const REB_STATUS_ENCOUNTER: c_int = 3;
pub const REB_STATUS_ESCAPE: c_int = 4;

// ASSIST body IDs
pub const ASSIST_BODY_SUN: c_int = 0;
pub const ASSIST_BODY_MERCURY: c_int = 1;
pub const ASSIST_BODY_VENUS: c_int = 2;
pub const ASSIST_BODY_EARTH: c_int = 3;
pub const ASSIST_BODY_MOON: c_int = 4;
pub const ASSIST_BODY_MARS: c_int = 5;
pub const ASSIST_BODY_JUPITER: c_int = 6;
pub const ASSIST_BODY_SATURN: c_int = 7;
pub const ASSIST_BODY_URANUS: c_int = 8;
pub const ASSIST_BODY_NEPTUNE: c_int = 9;
pub const ASSIST_BODY_PLUTO: c_int = 10;
pub const ASSIST_BODY_NPLANETS: c_int = 11;

// ASSIST force flags
pub const ASSIST_FORCE_SUN: c_int = 0x01;
pub const ASSIST_FORCE_PLANETS: c_int = 0x02;
pub const ASSIST_FORCE_ASTEROIDS: c_int = 0x04;
pub const ASSIST_FORCE_NON_GRAVITATIONAL: c_int = 0x08;
pub const ASSIST_FORCE_EARTH_HARMONICS: c_int = 0x10;
pub const ASSIST_FORCE_SUN_HARMONICS: c_int = 0x20;
pub const ASSIST_FORCE_GR_EIH: c_int = 0x40;
pub const ASSIST_FORCE_GR_SIMPLE: c_int = 0x80;
pub const ASSIST_FORCE_GR_POTENTIAL: c_int = 0x100;

/// Default force flags: Sun + planets + asteroids + Earth J2/J3/J4 + Sun J2 + GR (EIH).
pub const ASSIST_FORCES_DEFAULT: c_int = ASSIST_FORCE_SUN
    | ASSIST_FORCE_PLANETS
    | ASSIST_FORCE_ASTEROIDS
    | ASSIST_FORCE_EARTH_HARMONICS
    | ASSIST_FORCE_SUN_HARMONICS
    | ASSIST_FORCE_GR_EIH;

// ASSIST status codes
pub const ASSIST_SUCCESS: c_int = 0;
pub const ASSIST_ERROR_EPHEM_FILE: c_int = 1;
pub const ASSIST_ERROR_AST_FILE: c_int = 2;

// ---------------------------------------------------------------------------
// REBOUND functions
// ---------------------------------------------------------------------------

unsafe extern "C" {
    // Simulation lifecycle
    pub fn reb_simulation_create() -> *mut reb_simulation;
    pub fn reb_simulation_free(r: *mut reb_simulation);

    // Particle management
    pub fn reb_simulation_add(r: *mut reb_simulation, pt: reb_particle);

    // Integration
    pub fn reb_simulation_integrate(r: *mut reb_simulation, tmax: c_double) -> c_int;
    pub fn reb_simulation_step(r: *mut reb_simulation);
    pub fn reb_simulation_synchronize(r: *mut reb_simulation);

    // Variational equations
    pub fn reb_simulation_add_variation_1st_order(
        r: *mut reb_simulation,
        testparticle: c_int,
    ) -> c_int;
}

// ---------------------------------------------------------------------------
// ASSIST functions
// ---------------------------------------------------------------------------

unsafe extern "C" {
    pub fn assist_ephem_create(
        planets_path: *const c_char,
        asteroids_path: *const c_char,
    ) -> *mut assist_ephem;
    pub fn assist_ephem_free(ephem: *mut assist_ephem);

    pub fn assist_attach(sim: *mut reb_simulation, ephem: *mut assist_ephem) -> *mut assist_extras;
    pub fn assist_free(ax: *mut assist_extras);
    pub fn assist_detach(sim: *mut reb_simulation, ax: *mut assist_extras);

    pub fn assist_get_particle(
        ephem: *const assist_ephem,
        particle_id: c_int,
        t: c_double,
    ) -> reb_particle;
    pub fn assist_get_particle_with_error(
        ephem: *const assist_ephem,
        particle_id: c_int,
        t: c_double,
        error: *mut c_int,
    ) -> reb_particle;

    pub fn assist_integrate_or_interpolate(ax: *mut assist_extras, t: c_double);
}

// ---------------------------------------------------------------------------
// C helper functions (from src/helpers.c)
// ---------------------------------------------------------------------------

unsafe extern "C" {
    // reb_simulation field accessors
    pub fn assist_rs_sim_get_t(r: *const reb_simulation) -> c_double;
    pub fn assist_rs_sim_set_t(r: *mut reb_simulation, t: c_double);
    pub fn assist_rs_sim_get_dt(r: *const reb_simulation) -> c_double;
    pub fn assist_rs_sim_set_dt(r: *mut reb_simulation, dt: c_double);
    pub fn assist_rs_sim_get_N(r: *const reb_simulation) -> c_uint;
    pub fn assist_rs_sim_get_N_var(r: *const reb_simulation) -> c_int;
    pub fn assist_rs_sim_get_N_active(r: *const reb_simulation) -> c_int;
    pub fn assist_rs_sim_set_N_active(r: *mut reb_simulation, n: c_int);
    pub fn assist_rs_sim_get_particles(r: *const reb_simulation) -> *mut reb_particle;
    pub fn assist_rs_sim_get_exact_finish_time(r: *const reb_simulation) -> c_int;
    pub fn assist_rs_sim_set_exact_finish_time(r: *mut reb_simulation, v: c_int);
    pub fn assist_rs_sim_get_force_is_velocity_dependent(r: *const reb_simulation) -> c_uint;
    pub fn assist_rs_sim_get_status(r: *const reb_simulation) -> c_int;
    pub fn assist_rs_sim_get_extras(r: *const reb_simulation) -> *mut c_void;
    pub fn assist_rs_sim_get_integrator(r: *const reb_simulation) -> c_int;
    pub fn assist_rs_sim_set_integrator(r: *mut reb_simulation, i: c_int);
    pub fn assist_rs_sim_get_gravity(r: *const reb_simulation) -> c_int;
    pub fn assist_rs_sim_set_gravity(r: *mut reb_simulation, g: c_int);
    pub fn assist_rs_sim_get_ias15_epsilon(r: *const reb_simulation) -> c_double;
    pub fn assist_rs_sim_set_ias15_epsilon(r: *mut reb_simulation, eps: c_double);

    // assist_extras field accessors
    pub fn assist_rs_extras_get_forces(ax: *const assist_extras) -> c_int;
    pub fn assist_rs_extras_set_forces(ax: *mut assist_extras, f: c_int);
    pub fn assist_rs_extras_get_geocentric(ax: *const assist_extras) -> c_int;
    pub fn assist_rs_extras_set_geocentric(ax: *mut assist_extras, g: c_int);
    pub fn assist_rs_extras_get_particle_params(ax: *const assist_extras) -> *mut c_double;
    pub fn assist_rs_extras_set_particle_params(ax: *mut assist_extras, p: *mut c_double);

    // Non-gravitational force model parameters
    pub fn assist_rs_extras_get_alpha(ax: *const assist_extras) -> c_double;
    pub fn assist_rs_extras_set_alpha(ax: *mut assist_extras, v: c_double);
    pub fn assist_rs_extras_get_nk(ax: *const assist_extras) -> c_double;
    pub fn assist_rs_extras_set_nk(ax: *mut assist_extras, v: c_double);
    pub fn assist_rs_extras_get_nm(ax: *const assist_extras) -> c_double;
    pub fn assist_rs_extras_set_nm(ax: *mut assist_extras, v: c_double);
    pub fn assist_rs_extras_get_nn(ax: *const assist_extras) -> c_double;
    pub fn assist_rs_extras_set_nn(ax: *mut assist_extras, v: c_double);
    pub fn assist_rs_extras_get_r0(ax: *const assist_extras) -> c_double;
    pub fn assist_rs_extras_set_r0(ax: *mut assist_extras, v: c_double);

    // assist_ephem field accessors
    pub fn assist_rs_ephem_get_jd_ref(ephem: *const assist_ephem) -> c_double;
    pub fn assist_rs_ephem_set_jd_ref(ephem: *mut assist_ephem, jd: c_double);
    pub fn assist_rs_ephem_get_au(ephem: *const assist_ephem) -> c_double;
    pub fn assist_rs_ephem_get_clight(ephem: *const assist_ephem) -> c_double;
    pub fn assist_rs_ephem_get_c_au_per_day(ephem: *const assist_ephem) -> c_double;
    pub fn assist_rs_ephem_get_re(ephem: *const assist_ephem) -> c_double;
    pub fn assist_rs_ephem_get_re_eq(ephem: *const assist_ephem) -> c_double;
    pub fn assist_rs_ephem_get_emrat(ephem: *const assist_ephem) -> c_double;
}
