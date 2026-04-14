// helpers.c — Thin C accessors for reb_simulation fields.
//
// reb_simulation is enormous (~5KB with nested integrator structs).
// Rather than reproducing it in Rust, we treat it as opaque and
// expose only the handful of fields the FFI layer needs.

#include "rebound.h"
#include "assist.h"

// --- reb_simulation field accessors ---

double assist_rs_sim_get_t(const struct reb_simulation* r) { return r->t; }
void   assist_rs_sim_set_t(struct reb_simulation* r, double t) { r->t = t; }

double assist_rs_sim_get_dt(const struct reb_simulation* r) { return r->dt; }
void   assist_rs_sim_set_dt(struct reb_simulation* r, double dt) { r->dt = dt; }

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
