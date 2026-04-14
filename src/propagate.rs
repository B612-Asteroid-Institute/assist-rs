//! `assist_propagate` — N-body propagation of a test particle.

use crate::coordinates::{ecliptic_to_equatorial, equatorial_to_ecliptic, rotate_matrix_eq_to_ecl};
use crate::ffi;
use crate::wrappers::{AssistSim, Ephemeris, Simulation};
use crate::{Error, Result};

/// Result of propagating to a single epoch.
#[derive(Debug, Clone)]
pub struct PropagatedState {
    /// Heliocentric ecliptic J2000 state [x, y, z, vx, vy, vz] (AU, AU/day).
    pub state: [f64; 6],
    /// Epoch (MJD TDB).
    pub epoch: f64,
    /// 6x6 state transition matrix Phi(t, t0), row-major.
    /// Maps initial state perturbations to propagated state perturbations.
    /// Only populated when `compute_stm` is true.
    pub stm: Option<[[f64; 6]; 6]>,
}

/// Propagate a test particle from an initial heliocentric ecliptic J2000 state.
///
/// # Arguments
/// - `ephem`: ASSIST ephemeris data.
/// - `state`: initial [x, y, z, vx, vy, vz] in heliocentric ecliptic J2000 (AU, AU/day).
/// - `epoch`: initial epoch (MJD TDB).
/// - `target_epochs`: sorted slice of target epochs (MJD TDB).
/// - `compute_stm`: whether to compute the state transition matrix via variational equations.
///
/// # Returns
/// One `PropagatedState` per target epoch, in the same order.
pub fn assist_propagate(
    ephem: &Ephemeris,
    state: &[f64; 6],
    epoch: f64,
    target_epochs: &[f64],
    compute_stm: bool,
) -> Result<Vec<PropagatedState>> {
    if target_epochs.is_empty() {
        return Ok(vec![]);
    }

    let jd_ref = ephem.jd_ref();
    let t0 = mjd_to_assist_time(epoch, jd_ref);

    // Convert heliocentric ecliptic -> barycentric equatorial ICRF
    let eq_state = ecliptic_to_equatorial(state);
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
    asim.set_forces(ffi::ASSIST_FORCES_DEFAULT);

    // Add test particle (index 0)
    asim.sim_mut().add_test_particle(
        bary_state[0], bary_state[1], bary_state[2],
        bary_state[3], bary_state[4], bary_state[5],
    );

    // Add variational particles if STM requested
    let var_start_idx = if compute_stm {
        // Add 6 first-order variational equation sets, one per initial state dimension.
        // Each adds one variational particle (since we have 1 test particle).
        let idx = asim.sim_mut().add_variation_1st_order(0);
        // Initialize: each variational particle gets a unit perturbation
        // in one dimension. After integration, particle[idx+d] gives
        // the d-th column of the STM.
        let n = asim.sim().n_particles();
        for d in 0..6 {
            if d > 0 {
                asim.sim_mut().add_variation_1st_order(0);
            }
            // Set the initial perturbation: the variational particle
            // corresponding to dimension `d` gets 1.0 in that component.
            let vi = (idx as usize) + d;
            if vi < n {
                // Variational particles were just added, but we need to
                // set their initial conditions via the particle array.
                // REBOUND initializes variational particles to zero; we
                // need to set the unit perturbation.
            }
        }
        // Actually, REBOUND's add_variation_1st_order adds a set of
        // variational particles equal to N_real. For 1 test particle,
        // each call adds 1 variational particle. We need 6 calls total.
        // The first call was already made above; we made 5 more in the loop.
        // Now set initial conditions.
        //
        // REBOUND variational particles start at index N_real (after the test particle).
        // For 1 test particle, N_real = 1. After adding 6 variation sets,
        // we have particles at indices 1..6, each representing dx/dx0_d.
        //
        // We need to initialize each variational particle with a unit
        // perturbation in one phase-space dimension.
        let particles = asim.sim().particles();
        let n_real = 1usize; // our one test particle
        for d in 0..6 {
            let vi = n_real + d;
            if vi < particles.len() {
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
                // Write directly to the particle array
                unsafe {
                    let ptr = ffi::assist_rs_sim_get_particles(asim.sim().ptr);
                    *ptr.add(vi) = p;
                }
            }
        }
        idx
    } else {
        -1
    };

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

        // Convert barycentric equatorial -> heliocentric ecliptic
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

        // Extract STM if requested
        let stm = if compute_stm && var_start_idx >= 0 {
            let mut stm = [[0.0f64; 6]; 6];
            let n_real = 1usize;
            for d in 0..6 {
                let vi = n_real + d;
                if vi < particles.len() {
                    let vp = &particles[vi];
                    // Column d of the STM (in barycentric equatorial)
                    stm[0][d] = vp.x;
                    stm[1][d] = vp.y;
                    stm[2][d] = vp.z;
                    stm[3][d] = vp.vx;
                    stm[4][d] = vp.vy;
                    stm[5][d] = vp.vz;
                }
            }
            // Rotate STM from equatorial to ecliptic: R x STM_eq x R^T
            Some(rotate_matrix_eq_to_ecl(&stm))
        } else {
            None
        };

        results.push(PropagatedState {
            state: helio_ecl,
            epoch: target_mjd,
            stm,
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
