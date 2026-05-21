//! MPC observatory code table.
//!
//! Parses `obscodes_extended.json` from the `mpc-obscodes` package into a
//! lookup table mapping 3-character MPC codes to geodetic parallax
//! coefficients.
//!
//! An [`EarthOrientation`] kernel must be attached via
//! [`ObservatoryTable::with_earth_orientation`] before ground-based
//! observatory states can be computed. The kernel drives the ITRF93 → ICRF
//! rotation from a binary PCK file (matching JPL Horizons to ~μas).
//! Looking up a ground-based observatory without an attached kernel
//! returns `Error::MissingEarthOrientation`.

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
/// Optionally carries an [`EarthOrientation`] kernel that rotates ground-
/// based observatory positions from ITRF93 into ICRF. Without one, queries
/// for ground observatories return `Error::MissingEarthOrientation` rather
/// than falling back to an approximation.
#[derive(Debug, Clone)]
pub struct ObservatoryTable {
    entries: HashMap<String, ObservatoryEntry>,
    earth_orientation: Option<Arc<EarthOrientation>>,
}

impl ObservatoryTable {
    /// Load from `obscodes_extended.json`.
    ///
    /// Expected format: JSON object whose keys are MPC observatory codes.
    /// Each value is one of:
    ///
    /// * Ground/geocentric — has all four fields: `Longitude` (deg), `cos`
    ///   (cos(geocentric_lat) × ρ/R_eq), `sin` (same with sin), `Name`.
    /// * Space-based — has only `Name`; no surface coordinates.
    ///
    /// Partial entries (e.g. `Longitude` without `cos`/`sin`) are rejected
    /// rather than silently filled with NaN/zero defaults.
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
            entries.insert(code.clone(), parse_entry(code, val)?);
        }

        Ok(Self {
            entries,
            earth_orientation: None,
        })
    }

    /// Attach an [`EarthOrientation`] so ground-observer states can be
    /// rotated from ITRF93 into ICRF.
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

fn parse_entry(code: &str, val: &serde_json::Value) -> Result<ObservatoryEntry> {
    let obj = val
        .as_object()
        .ok_or_else(|| Error::Other(format!("observatory entry {code:?}: expected JSON object")))?;

    let name = obj
        .get("Name")
        .ok_or_else(|| Error::Other(format!("observatory entry {code:?}: missing Name")))?
        .as_str()
        .ok_or_else(|| Error::Other(format!("observatory entry {code:?}: Name is not a string")))?
        .to_string();

    let has_lon = obj.contains_key("Longitude");
    let has_cos = obj.contains_key("cos");
    let has_sin = obj.contains_key("sin");

    match (has_lon, has_cos, has_sin) {
        (false, false, false) => Ok(ObservatoryEntry {
            longitude_deg: f64::NAN,
            cos_lat: 0.0,
            sin_lat: 0.0,
            name,
        }),
        (true, true, true) => {
            let longitude_deg = parse_f64(obj, "Longitude", code)?;
            let cos_lat = parse_f64(obj, "cos", code)?;
            let sin_lat = parse_f64(obj, "sin", code)?;
            Ok(ObservatoryEntry {
                longitude_deg,
                cos_lat,
                sin_lat,
                name,
            })
        }
        _ => Err(Error::Other(format!(
            "observatory entry {code:?}: partial surface coordinates \
             (Longitude={has_lon}, cos={has_cos}, sin={has_sin}); \
             ground entries need all three, space-based entries need none"
        ))),
    }
}

fn parse_f64(
    obj: &serde_json::Map<String, serde_json::Value>,
    key: &str,
    code: &str,
) -> Result<f64> {
    obj.get(key)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| Error::Other(format!("observatory entry {code:?}: {key} is not a number")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_json(body: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(body.as_bytes()).unwrap();
        f
    }

    #[test]
    fn parses_ground_and_space_based_entries() {
        let f = write_json(
            r#"{
              "500": {"Longitude": 0.0, "cos": 0.0, "sin": 0.0, "Name": "Geocentric"},
              "I11": {"Longitude": 289.26345, "cos": 0.86502, "sin": -0.500901, "Name": "Gemini South"},
              "250": {"Name": "HST"}
            }"#,
        );
        let t = ObservatoryTable::from_json(f.path()).unwrap();
        assert_eq!(t.len(), 3);
        assert!(!t.get("I11").unwrap().is_space_based());
        assert!(t.get("250").unwrap().is_space_based());
        assert!(t.get("500").unwrap().is_space_based()); // geocentric flagged via zero cos/sin
    }

    #[test]
    fn rejects_partial_surface_coordinates() {
        let f = write_json(r#"{"X": {"Longitude": 1.0, "cos": 0.5, "Name": "missing sin"}}"#);
        let err = ObservatoryTable::from_json(f.path()).unwrap_err();
        assert!(err.to_string().contains("partial surface coordinates"));
    }

    #[test]
    fn rejects_missing_name() {
        let f = write_json(r#"{"X": {"Longitude": 1.0, "cos": 0.5, "sin": 0.5}}"#);
        let err = ObservatoryTable::from_json(f.path()).unwrap_err();
        assert!(err.to_string().contains("missing Name"));
    }

    #[test]
    fn rejects_non_numeric_coords() {
        let f =
            write_json(r#"{"X": {"Longitude": "deg", "cos": 0.5, "sin": 0.5, "Name": "broken"}}"#);
        let err = ObservatoryTable::from_json(f.path()).unwrap_err();
        assert!(err.to_string().contains("Longitude is not a number"));
    }
}
