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

/// Compile-time assertion that `Ephemeris` is `Send + Sync` — the property
/// relied on by callers that share a single loaded ephemeris across threads
/// via `Arc<Ephemeris>`.
const _EPHEMERIS_IS_SEND_SYNC: fn() = || {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<assist_rs::Ephemeris>();
};

#[test]
fn test_concurrent_assistsim_matches_serial() {
    // The audit in assist-rs-753 concluded that concurrent `AssistSim`
    // instances on separate threads are safe (REBOUND's hot paths use only
    // `static const` tables; the only process-shared state is the SIGINT
    // handler). This test actually exercises that: propagate the same orbit
    // many times in parallel and serially, and assert bit-for-bit equal
    // results. A hidden race on global state would show up as nondeterministic
    // output across runs or between the serial and parallel paths.
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    let orbit = assist_rs::Orbit::new(
        [
            -1.938_169_72,
            2.289_213_79,
            1.094_048_30,
            -0.008_744_54,
            -0.005_523_16,
            0.001_174_22,
        ],
        60000.0,
    );
    let targets = [60030.0];
    const N: usize = 32;

    // Serial baseline.
    let serial: Vec<[f64; 6]> = (0..N)
        .map(|_| assist_rs::assist_propagate(&ephem, &orbit, &targets, false).unwrap()[0].state)
        .collect();
    // All serial runs must be bitwise identical (same IC, same integrator).
    for (i, s) in serial.iter().enumerate().skip(1) {
        assert_eq!(
            serial[0], *s,
            "serial run {i} disagrees with serial run 0 — integrator is non-deterministic?"
        );
    }

    // Parallel runs of the same orbit must also agree bitwise with serial.
    use rayon::prelude::*;
    let parallel: Vec<[f64; 6]> = (0..N)
        .into_par_iter()
        .map(|_| assist_rs::assist_propagate(&ephem, &orbit, &targets, false).unwrap()[0].state)
        .collect();
    for (i, p) in parallel.iter().enumerate() {
        assert_eq!(
            serial[0], *p,
            "parallel run {i} disagrees with serial — REBOUND global state race?"
        );
    }
}

#[test]
fn test_integrate_empty_sim_maps_to_no_particles() {
    // REBOUND returns REB_STATUS_NO_PARTICLES when integrating a simulation
    // with zero particles; verify our wrapper maps that to the named
    // `Error::NoParticles` variant (not the opaque IntegrationFailed fallback).
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };
    let sim = assist_rs::Simulation::new().unwrap();
    let mut asim = assist_rs::AssistSim::new(sim, &ephem).unwrap();
    // No particles added — integrate should fail with NoParticles.
    let err = asim.integrate(1.0).unwrap_err();
    assert!(
        matches!(err, assist_rs::Error::NoParticles),
        "expected Error::NoParticles, got {err:?}"
    );
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
    let earth = ephem
        .get_body_state(assist_rs::ffi::ASSIST_BODY_EARTH, 0.0)
        .unwrap();
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
    let earth = ephem
        .get_body_state(assist_rs::ffi::ASSIST_BODY_EARTH, 0.0)
        .unwrap();
    asim.sim_mut().add_test_particle(
        earth.x + 0.01,
        earth.y,
        earth.z,
        earth.vx,
        earth.vy,
        earth.vz,
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
    let orbit = assist_rs::Orbit::new(
        [
            -1.938_169_72, // x
            2.289_213_79,  // y
            1.094_048_30,  // z
            -0.008_744_54, // vx
            -0.005_523_16, // vy
            0.001_174_22,  // vz
        ],
        60000.0,
    );

    // Propagate 30 days forward
    let target = [orbit.epoch + 30.0];
    let results = assist_rs::assist_propagate(&ephem, &orbit, &target, false).unwrap();

    assert_eq!(results.len(), 1);
    let r = &results[0];

    // The state should still be ~2-3 AU from the Sun
    let dist = (r.state[0] * r.state[0] + r.state[1] * r.state[1] + r.state[2] * r.state[2]).sqrt();
    assert!(
        dist > 1.5 && dist < 4.0,
        "Ceres heliocentric distance after 30 days: {dist} AU"
    );

    // Position should have changed by ~0.2-0.3 AU over 30 days
    let dx = r.state[0] - orbit.state[0];
    let dy = r.state[1] - orbit.state[1];
    let dz = r.state[2] - orbit.state[2];
    let displacement = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(
        displacement > 0.05 && displacement < 1.0,
        "Ceres displacement over 30 days: {displacement} AU"
    );

    eprintln!(
        "Ceres after 30 days: state={:?}, dist={dist:.4} AU",
        r.state
    );
}

#[test]
fn test_propagate_with_stm() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    let orbit = assist_rs::Orbit::new(
        [
            -1.938_169_72,
            2.289_213_79,
            1.094_048_30,
            -0.008_744_54,
            -0.005_523_16,
            0.001_174_22,
        ],
        60000.0,
    );
    let target = [orbit.epoch + 10.0];

    let results = assist_rs::assist_propagate(&ephem, &orbit, &target, true).unwrap();

    assert_eq!(results.len(), 1);
    let stm = results[0].stm.expect("STM should be populated");

    // STM should be close to identity for short propagation
    for (i, row) in stm.iter().enumerate() {
        assert!(
            row[i].abs() > 0.5 && row[i].abs() < 2.0,
            "STM diagonal [{i}][{i}] = {} (expected ~1.0)",
            row[i]
        );
    }

    eprintln!(
        "STM diagonal: [{:.4}, {:.4}, {:.4}, {:.4}, {:.4}, {:.4}]",
        stm[0][0], stm[1][1], stm[2][2], stm[3][3], stm[4][4], stm[5][5]
    );
}

#[test]
fn test_get_state() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    // Get Sun's heliocentric state (should be ~0)
    let sun =
        assist_rs::assist_get_state(&ephem, &assist_rs::Origin::Sun, &[60000.0], None).unwrap();
    let sun_dist = (sun[0].state[0] * sun[0].state[0]
        + sun[0].state[1] * sun[0].state[1]
        + sun[0].state[2] * sun[0].state[2])
        .sqrt();
    assert!(
        sun_dist < 0.01,
        "Sun heliocentric distance should be ~0, got {sun_dist} AU"
    );

    // Get Earth's heliocentric state
    let earth =
        assist_rs::assist_get_state(&ephem, &assist_rs::Origin::Earth, &[60000.0], None).unwrap();
    let earth_dist = (earth[0].state[0] * earth[0].state[0]
        + earth[0].state[1] * earth[0].state[1]
        + earth[0].state[2] * earth[0].state[2])
        .sqrt();
    assert!(
        earth_dist > 0.98 && earth_dist < 1.02,
        "Earth heliocentric distance: {earth_dist} AU"
    );

    eprintln!(
        "Earth heliocentric ecliptic at MJD 60000: {:?}",
        earth[0].state
    );
}

#[test]
fn test_generate_ephemeris() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    // Use Ceres orbit
    let orbit = assist_rs::Orbit::new(
        [
            -1.938_169_72,
            2.289_213_79,
            1.094_048_30,
            -0.008_744_54,
            -0.005_523_16,
            0.001_174_22,
        ],
        60000.0,
    );

    // Observer at geocenter
    let observer = assist_rs::Observer::new(assist_rs::Origin::Earth, 60010.0);

    let results = assist_rs::assist_generate_ephemeris(&ephem, &orbit, &[observer], None).unwrap();

    assert_eq!(results.len(), 1);
    let eph = &results[0];

    // Range should be reasonable (Ceres is 1.5-4 AU from Earth)
    let rho = eph.spherical[0];
    assert!(
        rho > 1.0 && rho < 5.0,
        "Geocentric range to Ceres: {rho} AU"
    );

    // RA should be in [0, 2π)
    let ra = eph.spherical[1];
    assert!(
        (0.0..std::f64::consts::TAU).contains(&ra),
        "RA out of range: {ra} rad"
    );

    // Dec should be in [-π/2, π/2]
    let dec = eph.spherical[2];
    assert!(
        (-std::f64::consts::FRAC_PI_2..=std::f64::consts::FRAC_PI_2).contains(&dec),
        "Dec out of range: {dec} rad"
    );

    // Light time should be positive and reasonable
    assert!(
        eph.light_time > 0.001 && eph.light_time < 0.05,
        "Light time: {} days",
        eph.light_time
    );

    eprintln!("Ceres from geocenter at MJD 60010:");
    eprintln!(
        "  rho={rho:.4} AU, RA={:.4}°, Dec={:.4}°",
        ra.to_degrees(),
        dec.to_degrees()
    );
    eprintln!(
        "  light_time={:.6} days ({:.1} min)",
        eph.light_time,
        eph.light_time * 24.0 * 60.0
    );
}

#[test]
fn test_propagate_with_non_grav() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    let ceres_state = [
        -1.938_169_72,
        2.289_213_79,
        1.094_048_30,
        -0.008_744_54,
        -0.005_523_16,
        0.001_174_22,
    ];
    let epoch = 60000.0;
    let target = [epoch + 30.0];

    // Propagate without non-grav forces (baseline)
    let orbit_grav = assist_rs::Orbit::new(ceres_state, epoch);
    let baseline = assist_rs::assist_propagate(&ephem, &orbit_grav, &target, false).unwrap();

    // Propagate with a small transverse non-grav acceleration (A2)
    let ng = assist_rs::NonGravParams::new(0.0, 1e-10, 0.0);
    let orbit_ng = assist_rs::Orbit::with_non_grav(ceres_state, epoch, ng);
    let with_ng = assist_rs::assist_propagate(&ephem, &orbit_ng, &target, false).unwrap();

    // The states should differ
    let dx: f64 = (0..6)
        .map(|i| (baseline[0].state[i] - with_ng[0].state[i]).powi(2))
        .sum::<f64>()
        .sqrt();

    assert!(
        dx > 1e-15,
        "Non-grav force had no effect: state difference = {dx}"
    );

    // But the difference should be small for this weak force over 30 days
    let pos_diff = ((baseline[0].state[0] - with_ng[0].state[0]).powi(2)
        + (baseline[0].state[1] - with_ng[0].state[1]).powi(2)
        + (baseline[0].state[2] - with_ng[0].state[2]).powi(2))
    .sqrt();
    assert!(
        pos_diff < 0.01,
        "Non-grav position difference too large: {pos_diff} AU"
    );

    eprintln!("Non-grav effect over 30 days: pos_diff={pos_diff:.2e} AU, state_diff={dx:.2e}");
}

#[test]
fn test_nongrav_partials_match_finite_differences() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    let ceres_state = [
        -1.938_169_72,
        2.289_213_79,
        1.094_048_30,
        -0.008_744_54,
        -0.005_523_16,
        0.001_174_22,
    ];
    let epoch = 60000.0;
    let target = [epoch + 30.0];

    // Baseline run: variational STM + nongrav partials.
    let a = [2e-10, 1e-10, -5e-11];
    let ng = assist_rs::NonGravParams::new(a[0], a[1], a[2]);
    let orbit = assist_rs::Orbit::with_non_grav(ceres_state, epoch, ng);
    let with_partials = assist_rs::assist_propagate(&ephem, &orbit, &target, true).unwrap();
    assert!(with_partials[0].stm.is_some(), "STM not populated");
    let partials = with_partials[0]
        .nongrav_partials
        .expect("nongrav_partials not populated for non-grav orbit");

    // Finite-difference each A_k in turn. Step chosen well above the
    // adaptive integrator's per-step noise but small enough to stay in
    // the linear regime of A perturbations.
    let h = 1e-12;
    for k in 0..3 {
        let mut a_plus = a;
        let mut a_minus = a;
        a_plus[k] += h;
        a_minus[k] -= h;
        let ng_plus = assist_rs::NonGravParams::new(a_plus[0], a_plus[1], a_plus[2]);
        let ng_minus = assist_rs::NonGravParams::new(a_minus[0], a_minus[1], a_minus[2]);
        let o_plus = assist_rs::Orbit::with_non_grav(ceres_state, epoch, ng_plus);
        let o_minus = assist_rs::Orbit::with_non_grav(ceres_state, epoch, ng_minus);
        let s_plus = assist_rs::assist_propagate(&ephem, &o_plus, &target, false).unwrap();
        let s_minus = assist_rs::assist_propagate(&ephem, &o_minus, &target, false).unwrap();

        for row in 0..6 {
            let fd = (s_plus[0].state[row] - s_minus[0].state[row]) / (2.0 * h);
            let analytic = partials[row][k];
            // Within 0.5% of the FD estimate — loose enough to absorb the
            // adaptive-step integrator's per-run noise, tight enough to
            // detect a missing factor of 2 or a sign flip.
            let tol = 5e-3 * analytic.abs().max(1e-6);
            let err = (fd - analytic).abs();
            assert!(
                err < tol,
                "∂x[{row}]/∂A{} mismatch: analytic={analytic:.6e} FD={fd:.6e} (tol={tol:.2e})",
                k + 1,
            );
        }
    }
}

#[test]
fn test_nongrav_partials_absent_without_nongrav() {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping: ephemeris not available");
        return;
    };

    let ceres_state = [
        -1.938_169_72,
        2.289_213_79,
        1.094_048_30,
        -0.008_744_54,
        -0.005_523_16,
        0.001_174_22,
    ];
    let orbit = assist_rs::Orbit::new(ceres_state, 60000.0);
    let result = assist_rs::assist_propagate(&ephem, &orbit, &[60030.0], true).unwrap();
    assert!(result[0].stm.is_some());
    assert!(
        result[0].nongrav_partials.is_none(),
        "nongrav_partials should be None for gravity-only orbit"
    );
}
