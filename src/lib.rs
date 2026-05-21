//! High-level domain layer for ASSIST + REBOUND solar-system propagation.
//!
//! The raw FFI bindings and safe RAII wrappers live in the companion
//! `libassist-sys` crate (which itself depends on `librebound-sys` for the
//! REBOUND C ABI). This crate adds:
//!
//! - [`Orbit`] + [`NonGravParams`] + [`Origin`] + [`AssistData`]: domain
//!   types for orbital states, non-gravitational coefficients, coordinate
//!   origins, and the bundled (ephemeris + observatory table) data handle.
//! - [`ObservatoryTable`] + [`earth_orientation`]: MPC observatory lookups
//!   and Earth-orientation kernel handling for topocentric observations.
//! - [`propagate`]: single-orbit + batched + rayon-parallel propagators,
//!   [`PropagatorPool`] for reusing simulations across many orbits, and STM
//!   + covariance propagation.
//! - [`ephemeris`]: topocentric ephemeris generation with light-time
//!   iteration.
//! - [`data`] (feature `data`): a [`data::DataManager`] that fetches +
//!   caches ephemeris and Earth-orientation kernels.
//!
//! All of `libassist-sys`'s public surface is re-exported here so downstream
//! consumers can keep a single `assist_rs` import.

mod assist_data;
pub mod coordinates;
#[cfg(feature = "data")]
pub mod data;
pub mod earth_orientation;
pub mod ephemeris;
mod observatory;
mod orbit;
mod origin;
pub mod propagate;
mod state;

pub use assist_data::AssistData;
pub use coordinates::{ecliptic_to_equatorial, equatorial_to_ecliptic};
pub use ephemeris::{
    EphemerisResult, Observer, assist_generate_ephemeris, assist_generate_ephemeris_single,
};
pub use observatory::ObservatoryTable;
pub use orbit::{NonGravParams, Orbit};
pub use origin::Origin;
pub use propagate::{
    PropagatedState, PropagatorConfig, PropagatorPool, assist_propagate, assist_propagate_single,
};
pub use state::{BodyState, assist_get_state};

// Re-export the FFI + RAII layer so downstream code keeps using `assist_rs::*`
// for both domain types and low-level sim/ephemeris handles.
pub use libassist_sys::{
    AssistSim, Ephemeris, Ias15AdaptiveMode, IntegratorConfig, Simulation, ffi,
};
// Re-export the sys crates themselves so downstream code can pattern-match on
// nested error variants (Error::Sys(libassist_sys::Error::Reb(_))) and access
// raw FFI symbols by their canonical paths if needed.
pub use {libassist_sys, librebound_sys};

/// Error type for assist-rs operations.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Wrapped FFI-layer error from libassist-sys (which itself wraps
    /// `librebound_sys::Error` for REBOUND integration-exit conditions).
    #[error(transparent)]
    Sys(#[from] libassist_sys::Error),

    #[error("light-time iteration did not converge after {0} iterations")]
    LightTimeConvergence(usize),

    #[error("invalid body identifier: {0}")]
    InvalidBody(String),

    #[error("invalid observatory code: {0}")]
    InvalidObservatory(String),

    #[error(
        "observatory {0} requires Earth orientation kernel; \
         attach via ObservatoryTable::with_earth_orientation"
    )]
    MissingEarthOrientation(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("{0}")]
    Other(String),
}

impl From<librebound_sys::Error> for Error {
    fn from(e: librebound_sys::Error) -> Self {
        Error::Sys(libassist_sys::Error::Reb(e))
    }
}

pub type Result<T> = std::result::Result<T, Error>;
