// helpers.c — Thin C accessors for reb_simulation fields.
//
// reb_simulation is enormous (~5KB with nested integrator structs).
// Rather than reproducing it in Rust, we treat it as opaque and
// expose only the handful of fields the FFI layer needs.

#include <string.h>
#include "rebound.h"
#include "assist.h"
#include "spk.h"

// --- reb_simulation field accessors ---

double assist_rs_sim_get_t(const struct reb_simulation* r) { return r->t; }
void   assist_rs_sim_set_t(struct reb_simulation* r, double t) { r->t = t; }

double assist_rs_sim_get_dt(const struct reb_simulation* r) { return r->dt; }
void   assist_rs_sim_set_dt(struct reb_simulation* r, double dt) { r->dt = dt; }

unsigned long long assist_rs_sim_get_steps_done(const struct reb_simulation* r) { return r->steps_done; }

// Reset the ephemeris cache to "all slots invalid (-1e306)". Emulates the
// assist_init cache-initialization behavior between propagations.
//
// WHY THIS MATTERS FOR PERFORMANCE: `assist_all_ephem` does a 7-slot LRU on
// the cache; when cache is fresh (all slots = -1e306), the "find oldest" loop
// always picks slot 0 (single branch), but when slots are populated with
// stale-but-realistic t values from a previous propagation, the loop does real
// comparisons every iteration. Over ~7000 ephem lookups per 30-day Ceres-type
// integrate, this compounds to ~190 µs of overhead — the full PropagatorPool
// vs assist_propagate regression. Invalidating the cache between pool
// propagations restores fresh-sim performance.
void assist_rs_ephem_cache_reset(struct assist_extras* ax) {
    if (!ax || !ax->ephem_cache || !ax->ephem_cache->t) return;
    int N_total = ASSIST_BODY_NPLANETS;
    if (ax->ephem && ax->ephem->spk_asteroids) {
        N_total += ax->ephem->spk_asteroids->num;
    }
    for (int i = 0; i < 7 * N_total; i++) {
        ax->ephem_cache->t[i] = -1e306;
    }
}

// Zero IAS15's compensated-summation accumulators (csx, csv) and predictor
// arrays (b, e, br, er, g) in place, leaving their allocations intact. This
// is equivalent to `reb_integrator_ias15_reset` for correctness — it clears
// all state that carries across `reb_simulation_integrate` calls — but
// avoids the 13 free/malloc pairs that reset incurs. Required between
// propagations of two unrelated orbits that share the same simulation.
void assist_rs_ias15_zero_state(struct reb_simulation* r) {
    const int N_allocated = r->ri_ias15.N_allocated;
    if (N_allocated == 0) { return; }   // never stepped; nothing to zero

    // Compensated-summation accumulators for positions and velocities.
    memset(r->ri_ias15.csx,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.csv,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.csa0, 0, sizeof(double) * N_allocated);

    // b/e/br/er/g: 7-component tables of length N3 each. Zero all p0..p6.
    memset(r->ri_ias15.b.p0,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.b.p1,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.b.p2,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.b.p3,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.b.p4,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.b.p5,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.b.p6,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.e.p0,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.e.p1,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.e.p2,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.e.p3,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.e.p4,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.e.p5,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.e.p6,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.br.p0, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.br.p1, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.br.p2, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.br.p3, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.br.p4, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.br.p5, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.br.p6, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.er.p0, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.er.p1, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.er.p2, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.er.p3, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.er.p4, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.er.p5, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.er.p6, 0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.g.p0,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.g.p1,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.g.p2,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.g.p3,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.g.p4,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.g.p5,  0, sizeof(double) * N_allocated);
    memset(r->ri_ias15.g.p6,  0, sizeof(double) * N_allocated);
}

unsigned int assist_rs_sim_get_N(const struct reb_simulation* r) { return r->N; }
int assist_rs_sim_get_N_var(const struct reb_simulation* r) { return r->N_var; }
int assist_rs_sim_get_N_active(const struct reb_simulation* r) { return r->N_active; }
void assist_rs_sim_set_N_active(struct reb_simulation* r, int n) { r->N_active = n; }

struct reb_particle* assist_rs_sim_get_particles(const struct reb_simulation* r) {
    return r->particles;
}

int assist_rs_sim_get_exact_finish_time(const struct reb_simulation* r) {
    return r->exact_finish_time;
}
void assist_rs_sim_set_exact_finish_time(struct reb_simulation* r, int v) {
    r->exact_finish_time = v;
}

unsigned int assist_rs_sim_get_force_is_velocity_dependent(const struct reb_simulation* r) {
    return r->force_is_velocity_dependent;
}

int assist_rs_sim_get_status(const struct reb_simulation* r) { return r->status; }

void* assist_rs_sim_get_extras(const struct reb_simulation* r) { return r->extras; }

// Integrator enum (inline in reb_simulation, not a separate field)
int assist_rs_sim_get_integrator(const struct reb_simulation* r) { return r->integrator; }
void assist_rs_sim_set_integrator(struct reb_simulation* r, int i) { r->integrator = i; }

// Gravity enum
int assist_rs_sim_get_gravity(const struct reb_simulation* r) { return r->gravity; }
void assist_rs_sim_set_gravity(struct reb_simulation* r, int g) { r->gravity = g; }

// IAS15 epsilon
double assist_rs_sim_get_ias15_epsilon(const struct reb_simulation* r) {
    return r->ri_ias15.epsilon;
}
void assist_rs_sim_set_ias15_epsilon(struct reb_simulation* r, double eps) {
    r->ri_ias15.epsilon = eps;
}

// IAS15 minimum timestep floor. When the adaptive step would shrink below
// this value, the integrator clamps it instead of grinding toward zero.
// Default 0 = no floor.
double assist_rs_sim_get_ias15_min_dt(const struct reb_simulation* r) {
    return r->ri_ias15.min_dt;
}
void assist_rs_sim_set_ias15_min_dt(struct reb_simulation* r, double min_dt) {
    r->ri_ias15.min_dt = min_dt;
}

// IAS15 adaptive-mode selector: 0=Individual, 1=Global, 2=PRS23 (REBOUND
// default since 2024-01), 3=Aarseth85.
int assist_rs_sim_get_ias15_adaptive_mode(const struct reb_simulation* r) {
    return (int)r->ri_ias15.adaptive_mode;
}
void assist_rs_sim_set_ias15_adaptive_mode(struct reb_simulation* r, int mode) {
    r->ri_ias15.adaptive_mode = mode;
}

// Diagnostic counter: number of IAS15 steps where the predictor-corrector
// loop hit the iteration cap without converging. Monotone-increasing across
// the simulation lifetime.
unsigned long long assist_rs_sim_get_ias15_iterations_max_exceeded(const struct reb_simulation* r) {
    return r->ri_ias15.iterations_max_exceeded;
}

// --- assist_extras field accessors ---

int assist_rs_extras_get_forces(const struct assist_extras* ax) { return ax->forces; }
void assist_rs_extras_set_forces(struct assist_extras* ax, int f) { ax->forces = f; }

int assist_rs_extras_get_geocentric(const struct assist_extras* ax) { return ax->geocentric; }
void assist_rs_extras_set_geocentric(struct assist_extras* ax, int g) { ax->geocentric = g; }

double* assist_rs_extras_get_particle_params(const struct assist_extras* ax) {
    return ax->particle_params;
}
void assist_rs_extras_set_particle_params(struct assist_extras* ax, double* p) {
    ax->particle_params = p;
}

// --- assist_ephem field accessors ---

double assist_rs_ephem_get_jd_ref(const struct assist_ephem* ephem) { return ephem->jd_ref; }
void   assist_rs_ephem_set_jd_ref(struct assist_ephem* ephem, double jd) { ephem->jd_ref = jd; }

double assist_rs_ephem_get_au(const struct assist_ephem* ephem) { return ephem->AU; }
double assist_rs_ephem_get_clight(const struct assist_ephem* ephem) { return ephem->CLIGHT; }
double assist_rs_ephem_get_c_au_per_day(const struct assist_ephem* ephem) { return ephem->c_AU_per_day; }
double assist_rs_ephem_get_re(const struct assist_ephem* ephem) { return ephem->RE; }
double assist_rs_ephem_get_re_eq(const struct assist_ephem* ephem) { return ephem->Re_eq; }
double assist_rs_ephem_get_emrat(const struct assist_ephem* ephem) { return ephem->EMRAT; }

// --- assist_extras non-gravitational force parameters ---

double assist_rs_extras_get_alpha(const struct assist_extras* ax) { return ax->alpha; }
void   assist_rs_extras_set_alpha(struct assist_extras* ax, double v) { ax->alpha = v; }

double assist_rs_extras_get_nk(const struct assist_extras* ax) { return ax->nk; }
void   assist_rs_extras_set_nk(struct assist_extras* ax, double v) { ax->nk = v; }

double assist_rs_extras_get_nm(const struct assist_extras* ax) { return ax->nm; }
void   assist_rs_extras_set_nm(struct assist_extras* ax, double v) { ax->nm = v; }

double assist_rs_extras_get_nn(const struct assist_extras* ax) { return ax->nn; }
void   assist_rs_extras_set_nn(struct assist_extras* ax, double v) { ax->nn = v; }

double assist_rs_extras_get_r0(const struct assist_extras* ax) { return ax->r0; }
void   assist_rs_extras_set_r0(struct assist_extras* ax, double v) { ax->r0 = v; }
