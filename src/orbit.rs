//! Core orbit type used across the API.

/// Marsden-Sekanina non-gravitational force parameters.
///
/// The non-gravitational acceleration uses the RTN (radial/transverse/normal)
/// decomposition with a g(r) model:
///
/// ```text
/// a_ng = g(r) * (A1·r̂ + A2·t̂ + A3·n̂)
/// g(r) = α * (r/r₀)^(-m) * (1 + (r/r₀)^n)^(-k)
/// ```
///
/// Default model parameters (α=1, m=2, n=5.093, k=0, r₀=1 AU) give
/// g(r) = r⁻² — a pure inverse-square law. The Marsden-Sekanina water ice
/// sublimation model uses α=0.111262, m=2.15, n=5.093, k=4.6142, r₀=2.808 AU.
#[derive(Debug, Clone)]
pub struct NonGravParams {
    /// Radial non-grav coefficient A1 (AU/day²).
    pub a1: f64,
    /// Transverse non-grav coefficient A2 (AU/day²).
    pub a2: f64,
    /// Normal non-grav coefficient A3 (AU/day²).
    pub a3: f64,
    /// g(r) normalization. Default: 1.0.
    pub alpha: Option<f64>,
    /// g(r) exponent k. Default: 0.0.
    pub nk: Option<f64>,
    /// g(r) exponent m. Default: 2.0.
    pub nm: Option<f64>,
    /// g(r) exponent n. Default: 5.093.
    pub nn: Option<f64>,
    /// g(r) scale distance in AU. Default: 1.0.
    pub r0: Option<f64>,
}

impl NonGravParams {
    /// Create non-grav params with the default g(r) = r⁻² model.
    pub fn new(a1: f64, a2: f64, a3: f64) -> Self {
        Self {
            a1,
            a2,
            a3,
            alpha: None,
            nk: None,
            nm: None,
            nn: None,
            r0: None,
        }
    }
}

/// A test particle orbit: state vector + epoch + optional non-gravitational parameters.
///
/// All coordinates are **heliocentric ecliptic J2000** (AU, AU/day) with **MJD TDB** epochs.
#[derive(Debug, Clone)]
pub struct Orbit {
    /// Heliocentric ecliptic J2000 state [x, y, z, vx, vy, vz] (AU, AU/day).
    pub state: [f64; 6],
    /// Epoch (MJD TDB).
    pub epoch: f64,
    /// Optional non-gravitational force parameters.
    pub non_grav: Option<NonGravParams>,
}

impl Orbit {
    /// Create an orbit with no non-gravitational forces.
    pub fn new(state: [f64; 6], epoch: f64) -> Self {
        Self {
            state,
            epoch,
            non_grav: None,
        }
    }

    /// Create an orbit with non-gravitational force parameters.
    pub fn with_non_grav(state: [f64; 6], epoch: f64, non_grav: NonGravParams) -> Self {
        Self {
            state,
            epoch,
            non_grav: Some(non_grav),
        }
    }
}
