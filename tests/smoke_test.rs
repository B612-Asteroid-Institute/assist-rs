//! Hermetic + env-gated smoke tests for the lean binding surface.
//!
//! High-level propagation/ephemeris behavior is validated downstream in
//! `adam-assist` (`adam_assist_rs`), which owns that orchestration.

use std::path::PathBuf;

fn ephem_paths() -> Option<(PathBuf, PathBuf)> {
    let planets = std::env::var("ASSIST_PLANETS_PATH").ok()?;
    let asteroids = std::env::var("ASSIST_ASTEROIDS_PATH").ok()?;
    Some((PathBuf::from(planets), PathBuf::from(asteroids)))
}

fn load_ephem() -> Option<assist_rs::Ephemeris> {
    let (planets, asteroids) = ephem_paths()?;
    assist_rs::Ephemeris::from_paths(&planets, &asteroids).ok()
}

#[test]
fn create_simulation_and_free_cleanly() {
    let sim = assist_rs::Simulation::new().unwrap();
    assert_eq!(sim.n_particles(), 0);
}

#[test]
fn simulation_time_and_dt_roundtrip() {
    let mut sim = assist_rs::Simulation::new().unwrap();
    sim.set_t(123.5);
    sim.set_dt(0.25);
    assert_eq!(sim.t(), 123.5);
    assert_eq!(sim.dt(), 0.25);
}

#[test]
fn particles_are_added_and_readable() {
    let mut sim = assist_rs::Simulation::new().unwrap();
    sim.add_test_particle(1.0, 2.0, 3.0, 0.1, 0.2, 0.3);
    sim.add_test_particle(-1.0, 0.0, 0.5, 0.0, -0.1, 0.0);
    assert_eq!(sim.n_particles(), 2);
    let particles = sim.particles();
    assert_eq!(particles[0].x, 1.0);
    assert_eq!(particles[1].vy, -0.1);
}

#[test]
fn ias15_knobs_roundtrip() {
    let mut sim = assist_rs::Simulation::new().unwrap();
    sim.set_ias15_epsilon(1e-8);
    sim.set_ias15_min_dt(1e-9);
    sim.set_ias15_adaptive_mode(assist_rs::Ias15AdaptiveMode::Global);
    assert_eq!(sim.ias15_epsilon(), 1e-8);
    assert_eq!(sim.ias15_min_dt(), 1e-9);
    assert_eq!(
        sim.ias15_adaptive_mode(),
        assist_rs::Ias15AdaptiveMode::Global
    );
}

#[test]
fn integrator_config_applies_every_knob() {
    let mut sim = assist_rs::Simulation::new().unwrap();
    let config = assist_rs::IntegratorConfig {
        initial_dt: Some(1e-6),
        min_dt: Some(1e-9),
        epsilon: Some(1e-6),
        adaptive_mode: Some(assist_rs::Ias15AdaptiveMode::Global),
    };
    config.apply(&mut sim);
    assert_eq!(sim.dt(), 1e-6);
    assert_eq!(sim.ias15_min_dt(), 1e-9);
    assert_eq!(sim.ias15_epsilon(), 1e-6);
    assert_eq!(
        sim.ias15_adaptive_mode(),
        assist_rs::Ias15AdaptiveMode::Global
    );
}

#[test]
fn variational_particles_attach_to_test_particle() {
    let mut sim = assist_rs::Simulation::new().unwrap();
    sim.add_test_particle(1.0, 0.0, 0.0, 0.0, 0.017, 0.0);
    let first_var = sim.add_variation_1st_order(0);
    assert_eq!(first_var, 1);
    assert_eq!(sim.n_particles(), 2);
    assert_eq!(sim.n_var(), 1);
}

// --- env-gated: require real DE440/SB441 kernels ---------------------------

#[test]
fn ephemeris_loads_and_serves_body_states_with_env_kernels() {
    let Some(ephem) = load_ephem() else {
        eprintln!("skipping: set ASSIST_PLANETS_PATH and ASSIST_ASTEROIDS_PATH");
        return;
    };
    assert!(ephem.jd_ref() > 2_400_000.0);
    assert!(ephem.c_au_per_day() > 173.0 && ephem.c_au_per_day() < 174.0);
    let sun = ephem
        .get_body_state(assist_rs::ffi::ASSIST_BODY_SUN, 0.0)
        .unwrap();
    assert!(sun.x.abs() < 0.02, "sun stays near the barycenter");
}

#[test]
fn assist_sim_attaches_and_integrates_with_env_kernels() {
    let Some(ephem) = load_ephem() else {
        eprintln!("skipping: set ASSIST_PLANETS_PATH and ASSIST_ASTEROIDS_PATH");
        return;
    };
    let mut sim = assist_rs::Simulation::new().unwrap();
    sim.set_t(0.0);
    let mut asim = assist_rs::AssistSim::new(sim, &ephem).unwrap();
    asim.set_forces(assist_rs::ffi::ASSIST_FORCES_DEFAULT);
    // Roughly 1 AU circular-ish test particle; integrate 10 days.
    asim.sim_mut()
        .add_test_particle(1.0, 0.0, 0.0, 0.0, 0.0172, 0.0);
    asim.integrate(10.0).unwrap();
    let particle = asim.sim().particles()[0];
    let radius =
        (particle.x * particle.x + particle.y * particle.y + particle.z * particle.z).sqrt();
    assert!(
        (0.5..1.5).contains(&radius),
        "radius stayed bounded: {radius}"
    );
    assert!(asim.sim().steps_done() > 0);
}
