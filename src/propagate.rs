//! `assist_propagate` — N-body propagation of a test particle.

use crate::coordinates::{ecliptic_to_equatorial, equatorial_to_ecliptic, rotate_matrix_eq_to_ecl};
use crate::ffi;
use crate::orbit::Orbit;
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

/// Propagate a test particle from an initial heliocentric ecliptic J2000 orbit.
///
/// # Arguments
/// - `ephem`: ASSIST ephemeris data.
/// - `orbit`: initial orbit (state, epoch, optional non-grav params).
/// - `target_epochs`: sorted slice of target epochs (MJD TDB).
/// - `compute_stm`: whether to compute the state transition matrix via variational equations.
///
/// # Returns
/// One `PropagatedState` per target epoch, in the same order.
pub fn assist_propagate(
    ephem: &Ephemeris,
    orbit: &Orbit,
    target_epochs: &[f64],
    compute_stm: bool,
) -> Result<Vec<PropagatedState>> {
    if target_epochs.is_empty() {
        return Ok(vec![]);
    }

    let jd_ref = ephem.jd_ref();
    let t0 = mjd_to_assist_time(orbit.epoch, jd_ref);

    // Convert heliocentric ecliptic → barycentric equatorial ICRF
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

    // Create simulation
    let mut sim = Simulation::new()?;
    sim.set_t(t0);
    let mut asim = AssistSim::new(sim, ephem)?;

    // Set force model: default + non-gravitational if requested
    let non_grav = orbit.non_grav.as_ref();
    let mut forces = ffi::ASSIST_FORCES_DEFAULT;
    if non_grav.is_some() {
        forces |= ffi::ASSIST_FORCE_NON_GRAVITATIONAL;
    }
    asim.set_forces(forces);

    // Add test particle (index 0)
    asim.sim_mut().add_test_particle(
        bary_state[0],
        bary_state[1],
        bary_state[2],
        bary_state[3],
        bary_state[4],
        bary_state[5],
    );

    // Set non-gravitational model parameters if provided.
    // particle_params is allocated after variational particles are added
    // (since the array must cover N_real + N_var particles).
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

    // Add variational particles if STM requested:
    //   d=0..5: unit state perturbations (column d of ∂x/∂x₀)
    //   d=6..8: unit perturbations of A1, A2, A3 (columns of ∂x/∂A)
    // Parameter variational particles only exist when the orbit carries
    // non-gravitational parameters; otherwise there's nothing to differentiate.
    let want_nongrav_partials = compute_stm && non_grav.is_some();
    let n_var = if want_nongrav_partials {
        9
    } else if compute_stm {
        6
    } else {
        0
    };

    if compute_stm {
        for _ in 0..n_var {
            asim.sim_mut().add_variation_1st_order(0);
        }

        // State variational particles: unit perturbation in one dimension.
        unsafe {
            let ptr = ffi::assist_rs_sim_get_particles(asim.sim().ptr);
            for d in 0..6 {
                let mut p = ffi::reb_particle::default();
                match d {
                    0 => p.x = 1.0,
                    1 => p.y = 1.0,
                    2 => p.z = 1.0,
                    3 => p.vx = 1.0,
                    4 => p.vy = 1.0,
                    5 => p.vz = 1.0,
                    _ => unreachable!(),
                }
                *ptr.add(1 + d) = p;
            }
            // Parameter variational particles already start at zero state
            // (reb_particle::default()) — that's the correct IC for ∂x/∂A.
        }
    }

    // Install particle_params for non-gravitational forces.
    // Must be done after variational particles are added since the array
    // covers all particles: 3 doubles per particle (real + variational).
    // Layout: [A1, A2, A3 | dA1_var0, dA2_var0, dA3_var0 | ...]
    //
    // For state variational particles (d=0..5), the parameter perturbation
    // is zero — ∂A_j/∂x₀_i = 0. For the three parameter variational particles
    // (d=6..8) we set a unit perturbation in one of A1, A2, A3 so ASSIST
    // integrates the corresponding sensitivity.
    if let Some(ng) = non_grav {
        let n_total = asim.sim().n_particles(); // real + variational
        let mut params = vec![0.0f64; 3 * n_total];
        params[0] = ng.a1;
        params[1] = ng.a2;
        params[2] = ng.a3;
        if want_nongrav_partials {
            // Variational particles are at indices 1..=9; param variationals
            // are the last three.
            let n_real = 1usize;
            for k in 0..3 {
                params[3 * (n_real + 6 + k) + k] = 1.0;
            }
        }
        asim.set_particle_params(params);
    }

    // Integrate to each target epoch sequentially (they should be sorted)
    let mut results = Vec::with_capacity(target_epochs.len());

    for &target_mjd in target_epochs {
        let t_target = mjd_to_assist_time(target_mjd, jd_ref);
        asim.integrate(t_target)?;

        // Read test particle state (index 0)
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

        // Extract STM and (optional) non-grav partials if requested.
        let (stm, nongrav_partials) = if compute_stm {
            let n_real = 1usize;
            let mut stm_eq = [[0.0f64; 6]; 6];
            for (d, vp) in particles[n_real..n_real + 6].iter().enumerate() {
                // Column d of the STM (in barycentric equatorial)
                stm_eq[0][d] = vp.x;
                stm_eq[1][d] = vp.y;
                stm_eq[2][d] = vp.z;
                stm_eq[3][d] = vp.vx;
                stm_eq[4][d] = vp.vy;
                stm_eq[5][d] = vp.vz;
            }
            // Rotate 6×6 STM from equatorial to ecliptic: R × STM_eq × R^T.
            let stm = Some(rotate_matrix_eq_to_ecl(&stm_eq));

            let nongrav = if want_nongrav_partials {
                let mut ng_eq = [[0.0f64; 3]; 6]; // columns: A1, A2, A3
                for (k, vp) in particles[n_real + 6..n_real + 9].iter().enumerate() {
                    ng_eq[0][k] = vp.x;
                    ng_eq[1][k] = vp.y;
                    ng_eq[2][k] = vp.z;
                    ng_eq[3][k] = vp.vx;
                    ng_eq[4][k] = vp.vy;
                    ng_eq[5][k] = vp.vz;
                }
                // The helio-ecl state is a *linear* function of bary-eq state
                // (subtract Sun(t), rotate by R); Sun(t) doesn't depend on A,
                // so the parameter partial transforms like a plain 6-vector.
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

/// Convert MJD TDB to ASSIST time (days from jd_ref).
///
/// ASSIST time = JD - jd_ref, where JD = MJD + 2400000.5.
fn mjd_to_assist_time(mjd_tdb: f64, jd_ref: f64) -> f64 {
    (mjd_tdb + 2_400_000.5) - jd_ref
}
