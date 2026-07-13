//! Lean Rust FFI bindings and safe wrappers for ASSIST + REBOUND.
//!
//! ASSIST is a C library for ephemeris-quality integration of test particles
//! in the solar system, built on top of the REBOUND N-body code. This crate
//! deliberately stays at the binding layer:
//!
//! - [`ffi`]: Raw `extern "C"` bindings to REBOUND and ASSIST functions.
//! - [`Simulation`]: Owned REBOUND simulation (particles, IAS15 knobs,
//!   integrate/step, variational particles).
//! - [`Ephemeris`]: Owned ASSIST ephemeris (`Send + Sync`; load once, share
//!   across threads).
//! - [`AssistSim`]: A simulation with ASSIST forces attached (force flags,
//!   non-gravitational model scalars, `particle_params`, integrator reset).
//! - [`IntegratorConfig`] / [`Ias15AdaptiveMode`]: IAS15 configuration.
//!
//! Everything above this layer — frame conventions, orbit containers, batch
//! or pooled propagation, STM/covariance handling, observatory tables, Earth
//! orientation, ephemeris generation, and kernel-file discovery — is
//! intentionally out of scope and owned by consumers (see `adam-assist`'s
//! `adam_assist_rs`, which hosts the previous high-level orchestration).

pub mod ffi;
mod wrappers;

pub use wrappers::{AssistSim, Ephemeris, Ias15AdaptiveMode, IntegratorConfig, Simulation};

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

    #[error("{0}")]
    Other(String),
}

pub type Result<T> = std::result::Result<T, Error>;
