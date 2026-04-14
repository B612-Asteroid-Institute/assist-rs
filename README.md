# assist-rs

Rust FFI bindings and safe wrappers for [ASSIST](https://github.com/matthewholman/assist) + [REBOUND](https://github.com/hannorein/rebound), providing ephemeris-quality integration of test particles in the solar system.

ASSIST uses REBOUND's IAS15 integrator (15th-order, adaptive step-size) to propagate test particle trajectories under the gravitational influence of the Sun, planets, Moon, and 16 massive asteroids, with positions sourced from JPL DE440/DE441 ephemerides.

## API

Three high-level functions, all operating in **heliocentric ecliptic J2000** coordinates (AU, AU/day) with **MJD TDB** epochs:

### `assist_propagate`

N-body propagation of a test particle with optional state transition matrix (STM) via variational equations, and optional non-gravitational forces.

```rust
let results = assist_propagate(
    &ephem,
    &[x, y, z, vx, vy, vz],  // initial state
    epoch_mjd,                 // initial epoch
    &[t1, t2, t3],            // target epochs (sorted)
    true,                      // compute STM
    None,                      // non-gravitational params (or Some(&ng))
)?;
// results[i].state  -> [f64; 6]
// results[i].stm    -> Option<[[f64; 6]; 6]>
```

### `assist_get_state`

Query the state of any solar system body (by name) or ground observatory (by MPC code).

```rust
let earth = assist_get_state(&ephem, "earth", 60000.0, None)?;
let obs   = assist_get_state(&ephem, "I11", 60000.0, Some(&obs_table))?;
// .state -> [f64; 6] heliocentric ecliptic J2000
```

### `assist_generate_ephemeris`

Propagate an orbit to observer epochs with light-time correction, returning topocentric spherical coordinates (range, RA, Dec + rates).

```rust
let results = assist_generate_ephemeris(
    &ephem,
    &orbit_state,
    orbit_epoch,
    &[Observer { state: obs_state, epoch: obs_epoch }],
    None,  // non-gravitational params
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

let results = assist_propagate(&ephem, &state, epoch, &targets, false, Some(&ng))?;
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

## Building

The crate vendors REBOUND and ASSIST as git submodules and compiles them from C source via the `cc` crate. No system-level dependencies beyond a C compiler and `libm`.

```bash
git clone --recursive <repo-url>
cd assist-rs
cargo build
```

## Testing

Tests require ephemeris data files. Set environment variables and run single-threaded (multiple `Ephemeris` instances in parallel can conflict over file descriptors):

```bash
export ASSIST_PLANETS_PATH=/path/to/de440.bsp
export ASSIST_ASTEROIDS_PATH=/path/to/sb441-n16.bsp
cargo test -- --test-threads=1
```

## Versioning

The crate version tracks the vendored ASSIST version: **assist-rs 1.2.x wraps ASSIST 1.2.0**. The patch version (the `x`) is reserved for Rust-side fixes that don't change the underlying C library. A new ASSIST release (e.g., 1.3.0) would become assist-rs 1.3.0.

The corresponding REBOUND version is pinned by the git submodule and recorded in `[package.metadata.vendored]` in Cargo.toml.

| assist-rs | ASSIST | REBOUND | License |
|-----------|--------|---------|---------|
| 1.2.x | [1.2.0](https://github.com/matthewholman/assist/releases/tag/v1.2.0) | [4.6.0](https://github.com/hannorein/rebound/releases/tag/4.6.0) | GPL-3.0 |

## License

GPL-3.0 -- required by the vendored REBOUND and ASSIST libraries.

## References

- Holman et al. 2023, "ASSIST: An Ephemeris-Quality Test Particle Integrator", [PSJ 4 69](https://doi.org/10.3847/PSJ/acc9a9) ([arXiv:2303.16246](https://arxiv.org/abs/2303.16246))
- Rein & Liu 2012, "REBOUND: An open-source multi-purpose N-body code for collisional dynamics", [A&A 537 A128](https://doi.org/10.1051/0004-6361/201118085)
- Rein & Spiegel 2015, "IAS15: a fast, adaptive, high-order integrator for gravitational dynamics", [MNRAS 446 1424](https://doi.org/10.1093/mnras/stu2164)

## Acknowledgments

Developed by the [Asteroid Institute](https://asteroid.institute/), a program of the [B612 Foundation](https://b612foundation.org/).
