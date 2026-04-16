//! Rust FFI bindings and safe wrappers for ASSIST + REBOUND.
//!
//! ASSIST is a C library for ephemeris-quality integration of test particles
//! in the solar system, built on top of the REBOUND N-body code. This crate
//! provides:
//!
//! - [`ffi`]: Raw `extern "C"` bindings to REBOUND and ASSIST functions.
//! - [`Simulation`], [`Ephemeris`], [`AssistSim`]: Safe RAII wrappers.
//! - Three high-level functions matching the THOR propagator interface:
//!   [`assist_propagate`], [`assist_get_state`], [`assist_generate_ephemeris`].

pub mod ffi;
mod wrappers;

pub use wrappers::{AssistSim, Ephemeris, Simulation};

pub mod coordinates;
#[cfg(feature = "data")]
pub mod data;
pub mod earth_orientation;
pub mod ephemeris;
mod observatory;
mod orbit;
mod origin;
mod propagate;
mod state;

pub use coordinates::{ecliptic_to_equatorial, equatorial_to_ecliptic};
pub use ephemeris::{EphemerisResult, Observer, assist_generate_ephemeris};
pub use observatory::ObservatoryTable;
pub use orbit::{NonGravParams, Orbit};
pub use origin::Origin;
pub use propagate::{PropagatedState, assist_propagate};
pub use state::{BodyState, assist_get_state};

/// Error type for assist-rs operations.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Integration ended early because no particles remain (`REB_STATUS_NO_PARTICLES`).
    #[error("integration ended: no particles remain in the simulation")]
    NoParticles,

    /// Integration ended early because two particles had a close encounter
    /// (`REB_STATUS_ENCOUNTER`; triggered by `exit_min_distance`).
    #[error("integration ended: close encounter")]
    CloseEncounter,

    /// Integration ended early because a particle escaped
    /// (`REB_STATUS_ESCAPE`; triggered by `exit_max_distance`).
    #[error("integration ended: particle escape")]
    Escape,

    /// Integration ended early because two particles collided
    /// (`REB_STATUS_COLLISION`).
    #[error("integration ended: collision")]
    Collision,

    /// REBOUND returned a generic/unknown error status.
    ///
    /// Holds the raw `REB_STATUS` code for diagnostics; use the named variants
    /// above to match on the common integration-exit conditions.
    #[error("REBOUND integration failed with status {0}")]
    IntegrationFailed(i32),

    #[error("ASSIST ephemeris error: {0}")]
    EphemerisError(String),

    #[error("light-time iteration did not converge after {0} iterations")]
    LightTimeConvergence(usize),

    #[error("invalid body identifier: {0}")]
    InvalidBody(String),

    #[error("invalid observatory code: {0}")]
    InvalidObservatory(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("{0}")]
    Other(String),
}

pub type Result<T> = std::result::Result<T, Error>;
