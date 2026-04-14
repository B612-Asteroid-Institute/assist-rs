//! Smoke tests for assist-rs.
//!
//! These tests require ephemeris data files. Set environment variables:
//!   ASSIST_PLANETS_PATH  — path to de440.bsp
//!   ASSIST_ASTEROIDS_PATH — path to sb441-n16.bsp
//!
//! Or the tests will be skipped.

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
fn test_create_simulation() {
    // No ephemeris needed — just test that REBOUND creates/frees cleanly.
    let sim = assist_rs::Simulation::new().unwrap();
    assert_eq!(sim.n_particles(), 0);
    assert!((sim.t() - 0.0).abs() < 1e-14);
    drop(sim);
}

#[test]
fn test_load_ephemeris() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ASSIST_PLANETS_PATH / ASSIST_ASTEROIDS_PATH not set");
        return;
    };
    // jd_ref should be J2000.0 = 2451545.0
    let jd_ref = ephem.jd_ref();
    assert!(
        (jd_ref - 2_451_545.0).abs() < 1.0,
        "Unexpected jd_ref: {jd_ref}"
    );
    let c = ephem.c_au_per_day();
    assert!(c > 170.0 && c < 175.0, "c_au_per_day out of range: {c}");
}

#[test]
fn test_get_earth_state() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };
    // Get Earth state at J2000.0 (t = 0)
    let earth = ephem.get_body_state(assist_rs::ffi::ASSIST_BODY_EARTH, 0.0).unwrap();
    // Earth should be ~1 AU from SSB
    let dist = (earth.x * earth.x + earth.y * earth.y + earth.z * earth.z).sqrt();
    assert!(
        dist > 0.9 && dist < 1.1,
        "Earth distance from SSB at J2000: {dist} AU"
    );
}

#[test]
fn test_create_assist_sim() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };
    let sim = assist_rs::Simulation::new().unwrap();
    let mut asim = assist_rs::AssistSim::new(sim, &ephem).unwrap();

    // Add a test particle near Earth (offset by 0.01 AU to avoid singularity)
    let earth = ephem.get_body_state(assist_rs::ffi::ASSIST_BODY_EARTH, 0.0).unwrap();
    asim.sim_mut().add_test_particle(
        earth.x + 0.01, earth.y, earth.z,
        earth.vx, earth.vy, earth.vz,
    );
    assert_eq!(asim.sim().n_particles(), 1);

    // Integrate 1 day
    asim.integrate(1.0).unwrap();
    let p = &asim.sim().particles()[0];
    let dist = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    assert!(
        dist > 0.5 && dist < 1.5,
        "Particle distance after 1 day: {dist} AU"
    );
}

#[test]
fn test_propagate_ceres() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    // Ceres heliocentric ecliptic J2000 state at MJD 60000.0 TDB
    // (approximate — from JPL Horizons query for epoch 2023-Feb-25)
    // x, y, z in AU; vx, vy, vz in AU/day
    let ceres_state = [
        -1.938_169_72,   // x
         2.289_213_79,   // y
         1.094_048_30,   // z
        -0.008_744_54,   // vx
        -0.005_523_16,   // vy
         0.001_174_22,   // vz
    ];
    let epoch = 60000.0; // MJD TDB

    // Propagate 30 days forward
    let target = [epoch + 30.0];
    let results = assist_rs::assist_propagate(
        &ephem,
        &ceres_state,
        epoch,
        &target,
        false,
    ).unwrap();

    assert_eq!(results.len(), 1);
    let r = &results[0];

    // The state should still be ~2-3 AU from the Sun
    let dist = (r.state[0] * r.state[0] + r.state[1] * r.state[1] + r.state[2] * r.state[2]).sqrt();
    assert!(
        dist > 1.5 && dist < 4.0,
        "Ceres heliocentric distance after 30 days: {dist} AU"
    );

    // Position should have changed by ~0.2-0.3 AU over 30 days
    let dx = r.state[0] - ceres_state[0];
    let dy = r.state[1] - ceres_state[1];
    let dz = r.state[2] - ceres_state[2];
    let displacement = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(
        displacement > 0.05 && displacement < 1.0,
        "Ceres displacement over 30 days: {displacement} AU"
    );

    eprintln!("Ceres after 30 days: state={:?}, dist={dist:.4} AU", r.state);
}

#[test]
fn test_propagate_with_stm() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    // Simple test: propagate a particle and get STM
    let state = [
        -1.938_169_72,
         2.289_213_79,
         1.094_048_30,
        -0.008_744_54,
        -0.005_523_16,
         0.001_174_22,
    ];
    let epoch = 60000.0;
    let target = [epoch + 10.0];

    let results = assist_rs::assist_propagate(
        &ephem,
        &state,
        epoch,
        &target,
        true,
    ).unwrap();

    assert_eq!(results.len(), 1);
    let stm = results[0].stm.expect("STM should be populated");

    // STM should be close to identity for short propagation
    // Diagonal elements should be approximately 1.0
    for i in 0..6 {
        assert!(
            stm[i][i].abs() > 0.5 && stm[i][i].abs() < 2.0,
            "STM diagonal [{i}][{i}] = {} (expected ~1.0)",
            stm[i][i]
        );
    }

    eprintln!("STM diagonal: [{:.4}, {:.4}, {:.4}, {:.4}, {:.4}, {:.4}]",
        stm[0][0], stm[1][1], stm[2][2], stm[3][3], stm[4][4], stm[5][5]);
}

#[test]
fn test_get_state() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    // Get Sun's heliocentric state (should be ~0)
    let sun = assist_rs::assist_get_state(&ephem, "sun", 60000.0, None).unwrap();
    let sun_dist = (sun.state[0] * sun.state[0] + sun.state[1] * sun.state[1] + sun.state[2] * sun.state[2]).sqrt();
    assert!(
        sun_dist < 0.01,
        "Sun heliocentric distance should be ~0, got {sun_dist} AU"
    );

    // Get Earth's heliocentric state
    let earth = assist_rs::assist_get_state(&ephem, "earth", 60000.0, None).unwrap();
    let earth_dist = (earth.state[0] * earth.state[0] + earth.state[1] * earth.state[1] + earth.state[2] * earth.state[2]).sqrt();
    assert!(
        earth_dist > 0.98 && earth_dist < 1.02,
        "Earth heliocentric distance: {earth_dist} AU"
    );

    eprintln!("Earth heliocentric ecliptic at MJD 60000: {:?}", earth.state);
}

#[test]
fn test_generate_ephemeris() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    // Use Ceres state
    let orbit_state = [
        -1.938_169_72, 2.289_213_79, 1.094_048_30,
        -0.008_744_54, -0.005_523_16, 0.001_174_22,
    ];
    let orbit_epoch = 60000.0;

    // Observer at geocenter
    let earth = assist_rs::assist_get_state(&ephem, "earth", 60010.0, None).unwrap();
    let observer = assist_rs::ephemeris::Observer {
        state: earth.state,
        epoch: 60010.0,
    };

    let results = assist_rs::assist_generate_ephemeris(
        &ephem,
        &orbit_state,
        orbit_epoch,
        &[observer],
    ).unwrap();

    assert_eq!(results.len(), 1);
    let eph = &results[0];

    // Range should be reasonable (Ceres is 1.5-4 AU from Earth)
    let rho = eph.spherical[0];
    assert!(
        rho > 1.0 && rho < 5.0,
        "Geocentric range to Ceres: {rho} AU"
    );

    // RA should be in [0, 2pi)
    let ra = eph.spherical[1];
    assert!(
        ra >= 0.0 && ra < std::f64::consts::TAU,
        "RA out of range: {ra} rad"
    );

    // Dec should be in [-pi/2, pi/2]
    let dec = eph.spherical[2];
    assert!(
        dec >= -std::f64::consts::FRAC_PI_2 && dec <= std::f64::consts::FRAC_PI_2,
        "Dec out of range: {dec} rad"
    );

    // Light time should be positive and reasonable
    assert!(
        eph.light_time > 0.001 && eph.light_time < 0.05,
        "Light time: {} days",
        eph.light_time
    );

    eprintln!("Ceres from geocenter at MJD 60010:");
    eprintln!("  rho={rho:.4} AU, RA={:.4} deg, Dec={:.4} deg",
        ra.to_degrees(), dec.to_degrees());
    eprintln!("  light_time={:.6} days ({:.1} min)",
        eph.light_time, eph.light_time * 24.0 * 60.0);
}
