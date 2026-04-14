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

mod propagate;
mod state;
pub mod ephemeris;
mod observatory;
pub mod coordinates;

pub use propagate::{assist_propagate, NonGravParams, PropagatedState};
pub use state::{assist_get_state, BodyState, Target};
pub use ephemeris::{assist_generate_ephemeris, EphemerisResult, Observer};
pub use observatory::ObservatoryTable;
pub use coordinates::{equatorial_to_ecliptic, ecliptic_to_equatorial};

/// Error type for assist-rs operations.
#[derive(Debug, thiserror::Error)]
pub enum Error {
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
