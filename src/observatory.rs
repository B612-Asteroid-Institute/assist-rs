//! MPC observatory code table.
//!
//! Parses `obscodes_extended.json` from the `mpc-obscodes` package into a
//! lookup table mapping 3-character MPC codes to geodetic parallax
//! coefficients.
//!
//! An optional [`EarthOrientation`] can be attached via
//! [`ObservatoryTable::with_earth_orientation`]. When present, the table
//! drives the ITRF93 → ICRF rotation from a binary PCK kernel (matching
//! JPL Horizons to ~μas); when absent, `compute_observatory_state` falls
//! back to a simplified IAU GMST formula good to ~50 mas.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use crate::earth_orientation::EarthOrientation;
use crate::{Error, Result};

/// A single observatory's parallax coefficients.
#[derive(Debug, Clone)]
pub struct ObservatoryEntry {
    /// Longitude in degrees (NaN for space-based observatories).
    pub longitude_deg: f64,
    /// cos(geocentric latitude) × (ρ/R_eq).
    pub cos_lat: f64,
    /// sin(geocentric latitude) × (ρ/R_eq).
    pub sin_lat: f64,
    /// Observatory name.
    pub name: String,
}

impl ObservatoryEntry {
    /// True if this is a space-based or geocentric observatory (no surface coordinates).
    pub fn is_space_based(&self) -> bool {
        self.longitude_deg.is_nan() || (self.cos_lat == 0.0 && self.sin_lat == 0.0)
    }
}

/// Lookup table from MPC observatory code to parallax coefficients.
///
/// Optionally carries an [`EarthOrientation`] so that ground-based
/// observatory states are rotated from ITRF93 into ICRF using a binary
/// PCK kernel rather than a simplified GMST approximation.
#[derive(Debug, Clone)]
pub struct ObservatoryTable {
    entries: HashMap<String, ObservatoryEntry>,
    earth_orientation: Option<Arc<EarthOrientation>>,
}

impl ObservatoryTable {
    /// Load from `obscodes_extended.json`.
    ///
    /// Expected format: JSON object where keys are 3-char MPC codes and values
    /// are objects with fields: `Longitude`, `cos`, `sin`, `Name`.
    pub fn from_json(path: &Path) -> Result<Self> {
        let data = std::fs::read_to_string(path).map_err(|e| {
            Error::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {e}", path.display()),
            ))
        })?;
        let raw: HashMap<String, serde_json::Value> = serde_json::from_str(&data)
            .map_err(|e| Error::Other(format!("Failed to parse observatory JSON: {e}")))?;

        let mut entries = HashMap::with_capacity(raw.len());
        for (code, val) in &raw {
            let longitude_deg = val
                .get("Longitude")
                .and_then(|v| v.as_f64())
                .unwrap_or(f64::NAN);
            let cos_lat = val.get("cos").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let sin_lat = val.get("sin").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let name = val
                .get("Name")
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string();

            entries.insert(
                code.clone(),
                ObservatoryEntry {
                    longitude_deg,
                    cos_lat,
                    sin_lat,
                    name,
                },
            );
        }

        Ok(Self {
            entries,
            earth_orientation: None,
        })
    }

    /// Attach an [`EarthOrientation`] so ground-observer states use the
    /// binary PCK ITRF93 → ICRF rotation instead of the GMST fallback.
    pub fn with_earth_orientation(mut self, eo: EarthOrientation) -> Self {
        self.earth_orientation = Some(Arc::new(eo));
        self
    }

    /// The attached Earth orientation kernel, if any.
    pub(crate) fn earth_orientation(&self) -> Option<&EarthOrientation> {
        self.earth_orientation.as_deref()
    }

    /// Look up an observatory by its MPC code.
    pub fn get(&self, code: &str) -> Option<&ObservatoryEntry> {
        self.entries.get(code)
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
