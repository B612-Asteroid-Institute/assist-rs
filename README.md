# assist-rs: Domain layer for ASSIST + REBOUND solar-system propagation
#### A Rust crate by the Asteroid Institute, a program of the B612 Foundation

[![REBOUND 4.6.0](https://img.shields.io/badge/REBOUND-4.6.0-orange?style=flat-square)](https://github.com/hannorein/rebound)
[![ASSIST 1.2.0](https://img.shields.io/badge/ASSIST-1.2.0-orange?style=flat-square)](https://github.com/matthewholman/assist)<br/>
[![CI](https://github.com/B612-Asteroid-Institute/assist-rs/actions/workflows/rust.yml/badge.svg?style=flat-square)](https://github.com/B612-Asteroid-Institute/assist-rs/actions/workflows/rust.yml)
[![crates.io](https://img.shields.io/crates/v/assist-rs.svg?style=flat-square)](https://crates.io/crates/assist-rs)
[![docs.rs](https://img.shields.io/docsrs/assist-rs?style=flat-square)](https://docs.rs/assist-rs)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg?style=flat-square)](LICENSE)<br/>
[![GitHub](https://img.shields.io/badge/GitHub-B612--Asteroid--Institute-181717?style=flat-square&logo=github&logoColor=white)](https://github.com/B612-Asteroid-Institute)
[![Website](https://img.shields.io/badge/Website-asteroid.institute-1f6feb?style=flat-square&logo=googlechrome&logoColor=white)](https://asteroid.institute/)

High-level Rust API for ephemeris-quality integration of test particles in the solar system. Wraps [ASSIST](https://github.com/matthewholman/assist) + [REBOUND](https://github.com/hannorein/rebound) and adds orbits, observatories, light-time iteration, STM propagation, and data-file management.

ASSIST uses REBOUND's IAS15 integrator (15th-order, adaptive step-size) to propagate test particle trajectories under the gravitational influence of the Sun, planets, Moon, and 16 massive asteroids, with positions sourced from JPL DE440/DE441 ephemerides.

## Crate hierarchy

```text
assist-rs            ← domain types: Orbit, Origin, Observer, DataManager (this crate)
  └── libassist-sys  ← raw ASSIST FFI + AssistSim/Ephemeris RAII
        └── librebound-sys  ← raw REBOUND FFI + Simulation RAII
```

Most users want **`assist-rs`**. The `-sys` crates are useful only if you need raw FFI access or want to call REBOUND directly without ASSIST forces.

## API

### Core types

#### `Orbit`

Bundles a state vector, epoch, and optional non-gravitational parameters:

```rust
// Gravity-only orbit
let orbit = Orbit::new(
    [x, y, z, vx, vy, vz],  // heliocentric ecliptic J2000 (AU, AU/day)
    epoch_mjd,                // MJD TDB
);

// Orbit with non-gravitational forces
let ng = NonGravParams::new(0.0, 1e-10, 0.0);  // A1, A2, A3
let orbit = Orbit::with_non_grav([x, y, z, vx, vy, vz], epoch_mjd, ng);
```

#### `Origin`

Identifies a solar system body or ground observatory:

```rust
Origin::Sun
Origin::Earth
Origin::Moon
Origin::JupiterBarycenter    // giant planets are system barycenters
Origin::SaturnBarycenter
Origin::Observatory("I11".into())  // MPC observatory code
Origin::SolarSystemBarycenter      // SSB

// Parse from string. Body names are case-insensitive; MPC codes must be
// exactly 3 alphanumeric characters. Returns Err on typos.
Origin::parse("earth")?      // → Origin::Earth
Origin::parse("jupiter")?    // → Origin::JupiterBarycenter
Origin::parse("W84")?        // → Origin::Observatory("W84")
Origin::parse("earty").is_err()   // typo → error, not silent observatory
```

### `AssistData`

Bundle of the data resources every high-level entry point needs: the JPL SPK `Ephemeris`, and optionally an `ObservatoryTable` (which itself may carry an `EarthOrientation` for sub-mas observatory rotation). Load once, pass by reference to every `assist_*` call:

```rust
let ephem = Ephemeris::from_paths(planets, asteroids)?;

// Minimal: propagation / named-body state queries only
let data = AssistData::new(ephem);

// With observatory support (required for Origin::Observatory queries and
// ground-based observers)
let data = AssistData::new(ephem).with_observatory(obs_table);
```

### `assist_propagate` / `assist_propagate_single`

N-body propagation of a test particle with optional state transition matrix (STM) via variational equations. `assist_propagate` takes many orbits and parallelizes across them; `_single` takes one orbit.

```rust
let results = assist_propagate_single(
    &data,
    &orbit,              // Orbit (state + epoch + optional non-grav)
    &[t1, t2, t3],      // target epochs (MJD TDB, sorted)
    true,                // compute STM
)?;
// results[i].state             -> [f64; 6]
// results[i].stm               -> Option<[[f64; 6]; 6]>            // ∂x(t)/∂x₀
// results[i].nongrav_partials  -> Option<[[f64; 3]; 6]>            // ∂x(t)/∂(A1,A2,A3)
```

`nongrav_partials` is populated only when `compute_stm = true` *and* the orbit carries non-gravitational parameters.

Each result can linearly propagate an initial covariance to its epoch:

```rust
// 6×6 state covariance → propagated 6×6
let p_t = results[i].propagate_covariance(&p0_6x6);

// 9×9 over (state, A1, A2, A3) → propagated 6×6 state covariance
let p_t = results[i].propagate_covariance_with_nongrav(&p0_9x9);
```

Both return `None` when the required partials are absent.

### `assist_get_state`

Query the state of any solar system body or ground observatory at one or more epochs. The fourth argument is `num_threads` (`None` = rayon global pool, `Some(1)` = serial).

```rust
let earth = assist_get_state(&data, &Origin::Earth, &[60000.0, 60001.0], Some(1))?;
// earth[0].state -> [f64; 6] heliocentric ecliptic J2000

let obs = assist_get_state(&data, &Origin::Observatory("I11".into()), &[60000.0], Some(1))?;
// `data` must have been built with `.with_observatory(..)` for this
```

### `assist_generate_ephemeris` / `assist_generate_ephemeris_single`

Propagate an orbit to observer epochs with light-time correction, returning topocentric spherical coordinates (range, RA, Dec + rates). `assist_generate_ephemeris` takes many orbits and parallelizes across them; `_single` takes one.

```rust
let results = assist_generate_ephemeris_single(
    &data,               // must include .with_observatory(..) for ground observatories
    &orbit,
    &[
        Observer::new(Origin::Earth, 60010.0),
        Observer::new(Origin::Observatory("I11".into()), 60011.0),
    ],
    Some(1),             // num_threads: None = global pool, Some(1) = serial
)?;
// results[i].spherical      -> [rho, ra, dec, drho, dra, ddec]
// results[i].aberrated_state -> [f64; 6] light-time-corrected heliocentric
// results[i].light_time      -> f64 (days)
```

### `NonGravParams`

Marsden-Sekanina non-gravitational force model with configurable g(r):

```text
a_ng = g(r) * (A1*r_hat + A2*t_hat + A3*n_hat)
g(r) = alpha * (r/r0)^(-m) * (1 + (r/r0)^n)^(-k)
```

```rust
// Default model: g(r) = r^-2 (inverse-square law)
let ng = NonGravParams::new(0.0, 1e-10, 0.0);  // A1, A2, A3

// Custom model (e.g. Marsden-Sekanina water ice sublimation)
let ng = NonGravParams {
    a1: 0.0, a2: 1e-10, a3: 0.0,
    alpha: Some(0.111262),
    nm: Some(2.15),
    nn: Some(5.093),
    nk: Some(4.6142),
    r0: Some(2.808),
};

let orbit = Orbit::with_non_grav(state, epoch, ng);
let results = assist_propagate_single(&data, &orbit, &targets, false)?;
```

## Setup

### Ephemeris data

Requires two JPL SPK files, available from the [B612 Asteroid Institute data packages](https://b612.ai/opensource/data_packages/) or directly from JPL:

| File | Size | Source |
|------|------|--------|
| `de440.bsp` | 114 MB | Planetary ephemeris (JPL DE440) |
| `sb441-n16.bsp` | 616 MB | 16 massive asteroid perturbers |

```rust
let ephem = Ephemeris::from_paths(
    Path::new("/path/to/de440.bsp"),
    Path::new("/path/to/sb441-n16.bsp"),
)?;
```

`Ephemeris` is `Send + Sync` -- load once, share across threads.

### Observatory table (optional)

For ground observatory lookups, load the MPC observatory codes from the [`mpc-obscodes`](https://pypi.org/project/mpc-obscodes/) package:

```rust
let obs_table = ObservatoryTable::from_json(
    Path::new("/path/to/obscodes_extended.json"),
)?;
```

### Earth orientation

Ground-based observatory positions are rotated from ITRF93 into ICRF via a binary PCK kernel (NAIF `earth_*.bpc` files). Looking up a ground observatory without an attached kernel returns `Error::MissingEarthOrientation` — no silent fallback to a low-precision approximation.

```rust
use assist_rs::earth_orientation::EarthOrientation;

let eo = EarthOrientation::from_paths(&[
    "/path/to/earth_latest_high_prec.bpc",
    "/path/to/earth_620120_250826.bpc",
])?;
let obs_table = obs_table.with_earth_orientation(eo);
```

### Automatic data download (optional)

The `data` feature (on by default) provides `DataManager`, which downloads and caches everything — ephemerides, observatory codes, and Earth orientation PCKs — directly from NAIF / JPL / MPC:

```rust
use assist_rs::data::DataManager;
use assist_rs::earth_orientation::EarthOrientation;

let dm = DataManager::new();           // default: ~/.cache/assist-rs
let paths = dm.ensure_ready()?;        // downloads missing files; re-fetches stale

let ephem = Ephemeris::from_paths(&paths.planets, &paths.asteroids)?;
let eo = EarthOrientation::from_paths(&paths.eop_kernels())?;
let obs_table = ObservatoryTable::from_json(&paths.obscodes)?
    .with_earth_orientation(eo);
let data = AssistData::new(ephem).with_observatory(obs_table);
```

`paths.eop_kernels()` returns the three EOP kernel paths in SPICE-idiomatic order (predict, historical, current) so the high-precision kernel wins at epochs it covers.

## Building

The crate vendors REBOUND and ASSIST as git submodules and compiles them from C source via the `cc` crate. No system-level dependencies beyond a C compiler and `libm`.

```bash
git clone --recursive <repo-url>
cd assist-rs
cargo build --release
```

For **maximum performance**, build with `CC=clang`:

```bash
CC=clang cargo build --release
```

`build.rs` sets `-flto=thin` (LLVM ThinLTO), which is only picked up by clang. On clang this enables cross-translation-unit inlining of ASSIST's hot ephemeris/force-evaluation routines and yields **~6–7 % faster propagation** than the default GCC build. GCC silently ignores the flag; GCC builds are still correct and work fine.

## Testing

Tests require ephemeris data files. Set environment variables and run single-threaded to keep test output readable:

```bash
export ASSIST_PLANETS_PATH=/path/to/de440.bsp
export ASSIST_ASTEROIDS_PATH=/path/to/sb441-n16.bsp
cargo test -- --test-threads=1
```

Concurrent `AssistSim` instances on separate threads are actually safe — see the thread-safety notes on `Ephemeris` in `src/wrappers.rs`. The `--test-threads=1` flag is preserved because stdout/stderr from multiple tests interleaves badly, not because REBOUND has races we're hiding.

The Horizons validation suite (`tests/horizons_v2_test.rs`) additionally needs the observatory table and (optionally) Earth orientation kernels. Set `ASSIST_DATA_DIR` to the directory containing `obscodes_extended.json` and `earth_*.bpc` files, or use the default `~/.cache/assist-rs`.

## Versioning

`assist-rs` follows **standard Rust semver** — version reflects the Rust API, not the upstream C library. Pre-1.0 (`0.x.y`) means breaking changes can happen on any minor bump; the crate will move to 1.0.0 once the domain-layer API is settled.

This differs from the wrapped `-sys` crates ([`librebound-sys`](https://crates.io/crates/librebound-sys), [`libassist-sys`](https://crates.io/crates/libassist-sys)), whose `major.minor` mirrors the upstream C library version. The two schemes are independent: an ASSIST 1.3 release bumps `libassist-sys` to 1.3.0, but `assist-rs` only bumps if its Rust API also changes.

| assist-rs | ASSIST | REBOUND | libassist-sys | librebound-sys |
|---|---|---|---|---|
| 0.1.x | [1.2.0](https://github.com/matthewholman/assist/releases/tag/v1.2.0) | [4.6.0](https://github.com/hannorein/rebound/releases/tag/4.6.0) | 1.2.x | 4.6.x |

Pinned upstream versions are also recorded in `[package.metadata.vendored]` in each crate's Cargo.toml.

## License

GPL-3.0 — required by the vendored ASSIST and REBOUND sources (pulled in transitively via [`libassist-sys`](https://github.com/B612-Asteroid-Institute/libassist-sys) and [`librebound-sys`](https://github.com/B612-Asteroid-Institute/librebound-sys), each of which preserves its upstream LICENSE). See [LICENSE](LICENSE).

## References

- Holman et al. 2023, "ASSIST: An Ephemeris-Quality Test Particle Integrator", [PSJ 4 69](https://doi.org/10.3847/PSJ/acc9a9) ([arXiv:2303.16246](https://arxiv.org/abs/2303.16246))
- Rein & Liu 2012, "REBOUND: An open-source multi-purpose N-body code for collisional dynamics", [A&A 537 A128](https://doi.org/10.1051/0004-6361/201118085)
- Rein & Spiegel 2015, "IAS15: a fast, adaptive, high-order integrator for gravitational dynamics", [MNRAS 446 1424](https://doi.org/10.1093/mnras/stu2164)

## Acknowledgments

This crate is a Rust wrapper plus a domain-layer API. All credit for the underlying physics and integrator implementations belongs to the upstream projects:

- **ASSIST** — Matthew Holman and contributors. Source: <https://github.com/matthewholman/assist>. Vendored via the companion [`libassist-sys`](https://github.com/B612-Asteroid-Institute/libassist-sys) crate ([crates.io](https://crates.io/crates/libassist-sys)), which preserves the upstream LICENSE.
- **REBOUND** — Hanno Rein and contributors. Source: <https://github.com/hannorein/rebound>. Vendored via the companion [`librebound-sys`](https://github.com/B612-Asteroid-Institute/librebound-sys) crate ([crates.io](https://crates.io/crates/librebound-sys)), which preserves the upstream LICENSE.

If you use this crate in published work, please cite the ASSIST and REBOUND papers listed in [References](#references) — not this crate.

The Rust wrapper and domain-layer API are developed by the [Asteroid Institute](https://asteroid.institute/), a program of the [B612 Foundation](https://b612foundation.org/).
