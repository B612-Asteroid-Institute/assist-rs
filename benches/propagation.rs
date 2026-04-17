//! Benchmarks for assist-rs propagation.
//!
//! Requires ephemeris data files. Set environment variables:
//!   ASSIST_PLANETS_PATH  — path to de440.bsp
//!   ASSIST_ASTEROIDS_PATH — path to sb441-n16.bsp
//!
//! Run with:
//!   cargo bench
//!
//! Compares:
//! - Rust high-level API (`assist_propagate`) which includes coordinate
//!   transforms, RAII setup/teardown, and error handling
//! - Raw C FFI calls (same integration, no coordinate transforms or wrappers)
//!   to measure the overhead of the Rust API layer

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::path::PathBuf;

// ---------------------------------------------------------------------------
// Ephemeris loading (shared across benchmarks)
// ---------------------------------------------------------------------------

fn ephem_paths() -> Option<(PathBuf, PathBuf)> {
    let planets = std::env::var("ASSIST_PLANETS_PATH").ok()?;
    let asteroids = std::env::var("ASSIST_ASTEROIDS_PATH").ok()?;
    Some((PathBuf::from(planets), PathBuf::from(asteroids)))
}

fn load_ephem() -> Option<assist_rs::Ephemeris> {
    let (planets, asteroids) = ephem_paths()?;
    assist_rs::Ephemeris::from_paths(&planets, &asteroids).ok()
}

// Ceres heliocentric ecliptic J2000 at MJD 60000.0 TDB
const CERES_STATE: [f64; 6] = [
    -1.938_169_72,
    2.289_213_79,
    1.094_048_30,
    -0.008_744_54,
    -0.005_523_16,
    0.001_174_22,
];
const EPOCH: f64 = 60000.0;

fn ceres_orbit() -> assist_rs::Orbit {
    assist_rs::Orbit::new(CERES_STATE, EPOCH)
}

// ---------------------------------------------------------------------------
// Rust high-level API benchmarks
// ---------------------------------------------------------------------------

fn bench_propagate_single(c: &mut Criterion) {
    let Some(ephem) = load_ephem() else {
        eprintln!("Skipping benchmarks: ASSIST_PLANETS_PATH / ASSIST_ASTEROIDS_PATH not set");
        return;
    };

    let orbit = ceres_orbit();
    let mut group = c.benchmark_group("propagate_single");

    // Propagate to N epochs over 30 days
    for n_epochs in [1, 10, 100] {
        let targets: Vec<f64> = (1..=n_epochs)
            .map(|i| EPOCH + 30.0 * (i as f64) / (n_epochs as f64))
            .collect();

        group.bench_with_input(
            BenchmarkId::new("rust_api", n_epochs),
            &targets,
            |b, targets| {
                b.iter(|| assist_rs::assist_propagate(&ephem, &orbit, targets, false).unwrap());
            },
        );
    }

    group.finish();
}

fn bench_propagate_with_stm(c: &mut Criterion) {
    let Some(ephem) = load_ephem() else { return };

    let orbit = ceres_orbit();
    let targets = vec![EPOCH + 30.0];

    c.benchmark_group("propagate_stm")
        .bench_function("without_stm", |b| {
            b.iter(|| assist_rs::assist_propagate(&ephem, &orbit, &targets, false).unwrap());
        })
        .bench_function("with_stm", |b| {
            b.iter(|| assist_rs::assist_propagate(&ephem, &orbit, &targets, true).unwrap());
        });
}

fn bench_propagate_with_nongrav(c: &mut Criterion) {
    let Some(ephem) = load_ephem() else { return };

    let orbit_grav = ceres_orbit();
    let ng = assist_rs::NonGravParams::new(0.0, 1e-10, 0.0);
    let orbit_ng = assist_rs::Orbit::with_non_grav(CERES_STATE, EPOCH, ng);
    let targets = vec![EPOCH + 30.0];

    c.benchmark_group("propagate_nongrav")
        .bench_function("gravity_only", |b| {
            b.iter(|| assist_rs::assist_propagate(&ephem, &orbit_grav, &targets, false).unwrap());
        })
        .bench_function("with_a2", |b| {
            b.iter(|| assist_rs::assist_propagate(&ephem, &orbit_ng, &targets, false).unwrap());
        });
}

// ---------------------------------------------------------------------------
// Raw C FFI benchmark — same integration, no Rust wrapper overhead
// ---------------------------------------------------------------------------

/// Propagate using raw FFI calls: no coordinate transforms, no RAII overhead
/// beyond what's needed for memory safety. This isolates the C integration time.
fn raw_c_propagate(
    ephem: &assist_rs::Ephemeris,
    bary_eq_state: &[f64; 6],
    t0: f64,
    t_target: f64,
) -> [f64; 6] {
    use assist_rs::ffi;
    unsafe {
        let sim = ffi::reb_simulation_create();
        ffi::assist_rs_sim_set_t(sim, t0);
        ffi::assist_rs_sim_set_exact_finish_time(sim, 1);

        let ax = ffi::assist_attach(sim, ephem.as_ptr().cast_mut());

        let p = ffi::reb_particle {
            x: bary_eq_state[0],
            y: bary_eq_state[1],
            z: bary_eq_state[2],
            vx: bary_eq_state[3],
            vy: bary_eq_state[4],
            vz: bary_eq_state[5],
            ..Default::default()
        };
        ffi::reb_simulation_add(sim, p);

        ffi::reb_simulation_integrate(sim, t_target);

        let particles = ffi::assist_rs_sim_get_particles(sim);
        let pp = &*particles;
        let result = [pp.x, pp.y, pp.z, pp.vx, pp.vy, pp.vz];

        ffi::assist_detach(sim, ax);
        ffi::assist_free(ax);
        ffi::reb_simulation_free(sim);

        result
    }
}

fn bench_rust_vs_raw_c(c: &mut Criterion) {
    let Some(ephem) = load_ephem() else { return };

    // Pre-compute the barycentric equatorial state for the raw C path
    let jd_ref = ephem.jd_ref();
    let t0 = (EPOCH + 2_400_000.5) - jd_ref;
    let eq_state = assist_rs::coordinates::ecliptic_to_equatorial(&CERES_STATE);
    let sun = ephem
        .get_body_state(assist_rs::ffi::ASSIST_BODY_SUN, t0)
        .unwrap();
    let bary_eq = [
        eq_state[0] + sun.x,
        eq_state[1] + sun.y,
        eq_state[2] + sun.z,
        eq_state[3] + sun.vx,
        eq_state[4] + sun.vy,
        eq_state[5] + sun.vz,
    ];
    let t_target = ((EPOCH + 30.0) + 2_400_000.5) - jd_ref;

    let orbit = ceres_orbit();
    let targets = vec![EPOCH + 30.0];

    let mut group = c.benchmark_group("rust_vs_raw_c");

    group.bench_function("rust_api", |b| {
        b.iter(|| assist_rs::assist_propagate(&ephem, &orbit, &targets, false).unwrap());
    });

    group.bench_function("raw_c_ffi", |b| {
        b.iter(|| raw_c_propagate(&ephem, &bary_eq, t0, t_target));
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Parallel propagation benchmark (rayon)
// ---------------------------------------------------------------------------

fn bench_parallel_propagation(c: &mut Criterion) {
    let Some(ephem) = load_ephem() else { return };

    // 28 slightly different orbits (perturb velocity slightly)
    let orbits: Vec<assist_rs::Orbit> = (0..28)
        .map(|i| {
            let mut s = CERES_STATE;
            s[3] += (i as f64) * 1e-6;
            assist_rs::Orbit::new(s, EPOCH)
        })
        .collect();
    let targets = vec![EPOCH + 30.0];

    let mut group = c.benchmark_group("parallel");

    group.bench_function("serial_28_orbits", |b| {
        b.iter(|| {
            orbits
                .iter()
                .map(|orbit| assist_rs::assist_propagate(&ephem, orbit, &targets, false).unwrap())
                .collect::<Vec<_>>()
        });
    });

    group.bench_function("rayon_28_orbits", |b| {
        use rayon::prelude::*;
        b.iter(|| {
            orbits
                .par_iter()
                .map(|orbit| assist_rs::assist_propagate(&ephem, orbit, &targets, false).unwrap())
                .collect::<Vec<_>>()
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Propagation duration scaling
// ---------------------------------------------------------------------------

fn bench_duration_scaling(c: &mut Criterion) {
    let Some(ephem) = load_ephem() else { return };

    let orbit = ceres_orbit();
    let mut group = c.benchmark_group("duration_scaling");

    for days in [1, 10, 30, 100, 365] {
        let targets = vec![EPOCH + days as f64];
        group.bench_with_input(BenchmarkId::new("days", days), &targets, |b, targets| {
            b.iter(|| assist_rs::assist_propagate(&ephem, &orbit, targets, false).unwrap());
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// PropagatorPool vs assist_propagate: the pool amortizes REBOUND/ASSIST setup
// (~25 µs / orbit) across many orbits with matching force-model config.
// Measured by looping N orbits through each path. The per-iteration work of
// each bench is N full propagations; criterion divides out to per-loop time.
// ---------------------------------------------------------------------------

/// Generate N distinct orbits by perturbing Ceres's initial state. The
/// perturbation is small (~1e-6 AU on position, ~1e-8 AU/day on velocity) so
/// every orbit stays in roughly the main belt — IAS15 step counts stay
/// comparable and we're timing the setup overhead, not dynamical variation.
fn varied_orbits(n: usize) -> Vec<assist_rs::Orbit> {
    (0..n)
        .map(|i| {
            let f = i as f64;
            let mut state = CERES_STATE;
            state[0] += f * 1e-6;
            state[1] += f * 1.3e-6;
            state[2] += f * -0.7e-6;
            state[3] += f * 1e-8;
            state[4] += f * -1.1e-8;
            state[5] += f * 0.5e-8;
            assist_rs::Orbit::new(state, EPOCH)
        })
        .collect()
}

fn bench_pool_vs_unpooled(c: &mut Criterion) {
    let Some(ephem) = load_ephem() else { return };

    // 128 is small enough to keep the criterion inner-loop fast (<100 ms per
    // sample at dt=365) but large enough to wash out pool-construction cost.
    // The bench reports per-batch-of-128 time; divide for per-orbit.
    const N: usize = 128;
    let orbits = varied_orbits(N);

    for &days in &[30.0f64, 365.0] {
        let targets = [EPOCH + days];

        let mut group = c.benchmark_group(format!("pool_vs_unpooled_{}d", days as u32));

        group.bench_function("unpooled", |b| {
            b.iter(|| {
                for o in &orbits {
                    criterion::black_box(
                        assist_rs::assist_propagate(&ephem, o, &targets, false).unwrap(),
                    );
                }
            });
        });

        group.bench_function("pooled", |b| {
            b.iter_with_setup(
                || {
                    assist_rs::PropagatorPool::new(
                        &ephem,
                        assist_rs::PropagatorConfig::gravity_only(),
                    )
                    .unwrap()
                },
                |mut pool| {
                    for o in &orbits {
                        criterion::black_box(pool.propagate(o, &targets).unwrap());
                    }
                },
            );
        });

        // STM variant: same orbits, both paths with compute_stm=true. Pool
        // win should be a bit smaller here since the STM integration itself
        // costs much more than the setup it amortizes.
        group.bench_function("unpooled_with_stm", |b| {
            b.iter(|| {
                for o in &orbits {
                    criterion::black_box(
                        assist_rs::assist_propagate(&ephem, o, &targets, true).unwrap(),
                    );
                }
            });
        });

        group.bench_function("pooled_with_stm", |b| {
            b.iter_with_setup(
                || {
                    assist_rs::PropagatorPool::new(
                        &ephem,
                        assist_rs::PropagatorConfig::gravity_with_stm(),
                    )
                    .unwrap()
                },
                |mut pool| {
                    for o in &orbits {
                        criterion::black_box(pool.propagate(o, &targets).unwrap());
                    }
                },
            );
        });

        group.finish();
    }
}

// ---------------------------------------------------------------------------
// Criterion harness
// ---------------------------------------------------------------------------

criterion_group!(
    benches,
    bench_propagate_single,
    bench_propagate_with_stm,
    bench_propagate_with_nongrav,
    bench_rust_vs_raw_c,
    bench_parallel_propagation,
    bench_duration_scaling,
    bench_pool_vs_unpooled,
);
criterion_main!(benches);
