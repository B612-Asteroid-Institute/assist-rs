//! Origin types for state queries and observer positions.

use libassist_sys::ffi;

use crate::{Error, Result};

/// An origin point in the solar system: a named body, barycenter, or MPC observatory.
///
/// Bodies sourced from JPL DE440 ephemeris. The giant planets and Pluto are
/// system barycenters (not planet centers) in the underlying ephemeris data.
/// Earth and Moon are resolved individually by ASSIST.
#[derive(Debug, Clone)]
pub enum Origin {
    /// Solar system barycenter (SSB).
    SolarSystemBarycenter,
    Sun,
    /// Mercury system barycenter (≈ Mercury center; no significant satellites).
    MercuryBarycenter,
    /// Venus system barycenter (≈ Venus center; no significant satellites).
    VenusBarycenter,
    /// Earth center (resolved from Earth-Moon barycenter by ASSIST).
    Earth,
    /// Moon center (resolved from Earth-Moon barycenter by ASSIST).
    Moon,
    /// Mars system barycenter (≈ Mars center; Phobos/Deimos are negligible).
    MarsBarycenter,
    /// Jupiter system barycenter (not Jupiter center).
    JupiterBarycenter,
    /// Saturn system barycenter (not Saturn center).
    SaturnBarycenter,
    /// Uranus system barycenter (not Uranus center).
    UranusBarycenter,
    /// Neptune system barycenter (not Neptune center).
    NeptuneBarycenter,
    /// Pluto system barycenter (not Pluto center; Charon is massive).
    PlutoBarycenter,
    /// Ground-based or space-based observatory by MPC code (e.g., "I11", "W84", "500").
    Observatory(String),
}

impl Origin {
    /// Parse a string into an `Origin`.
    ///
    /// Recognises named bodies (case-insensitive) and 3-character alphanumeric
    /// MPC observatory codes (case-preserving). Anything else is rejected —
    /// previously a typo'd body name silently became `Origin::Observatory(s)`
    /// and failed later at lookup time with a misleading "unknown MPC code"
    /// error.
    ///
    /// Body names accept short and explicit forms ("earth", "jupiter",
    /// "jupiter_barycenter", "ssb").
    pub fn parse(s: &str) -> Result<Self> {
        if let Some(o) = parse_body_name(s) {
            return Ok(o);
        }
        if is_mpc_code(s) {
            return Ok(Origin::Observatory(s.to_string()));
        }
        Err(Error::InvalidBody(format!(
            "{s:?} is neither a known body name nor a valid 3-character MPC code"
        )))
    }

    /// Return the ASSIST body ID if this is a named body (not SSB or observatory).
    pub(crate) fn body_id(&self) -> Option<i32> {
        match self {
            Origin::Sun => Some(ffi::ASSIST_BODY_SUN),
            Origin::MercuryBarycenter => Some(ffi::ASSIST_BODY_MERCURY),
            Origin::VenusBarycenter => Some(ffi::ASSIST_BODY_VENUS),
            Origin::Earth => Some(ffi::ASSIST_BODY_EARTH),
            Origin::Moon => Some(ffi::ASSIST_BODY_MOON),
            Origin::MarsBarycenter => Some(ffi::ASSIST_BODY_MARS),
            Origin::JupiterBarycenter => Some(ffi::ASSIST_BODY_JUPITER),
            Origin::SaturnBarycenter => Some(ffi::ASSIST_BODY_SATURN),
            Origin::UranusBarycenter => Some(ffi::ASSIST_BODY_URANUS),
            Origin::NeptuneBarycenter => Some(ffi::ASSIST_BODY_NEPTUNE),
            Origin::PlutoBarycenter => Some(ffi::ASSIST_BODY_PLUTO),
            Origin::SolarSystemBarycenter | Origin::Observatory(_) => None,
        }
    }

    /// Return the MPC observatory code if this is an observatory.
    pub fn observatory_code(&self) -> Option<&str> {
        match self {
            Origin::Observatory(code) => Some(code),
            _ => None,
        }
    }
}

fn parse_body_name(s: &str) -> Option<Origin> {
    Some(match s.to_lowercase().as_str() {
        "ssb" | "solar_system_barycenter" => Origin::SolarSystemBarycenter,
        "sun" => Origin::Sun,
        "mercury" | "mercury_barycenter" => Origin::MercuryBarycenter,
        "venus" | "venus_barycenter" => Origin::VenusBarycenter,
        "earth" => Origin::Earth,
        "moon" => Origin::Moon,
        "mars" | "mars_barycenter" => Origin::MarsBarycenter,
        "jupiter" | "jupiter_barycenter" => Origin::JupiterBarycenter,
        "saturn" | "saturn_barycenter" => Origin::SaturnBarycenter,
        "uranus" | "uranus_barycenter" => Origin::UranusBarycenter,
        "neptune" | "neptune_barycenter" => Origin::NeptuneBarycenter,
        "pluto" | "pluto_barycenter" => Origin::PlutoBarycenter,
        _ => return None,
    })
}

/// MPC observatory codes are exactly three ASCII alphanumeric characters
/// (every entry in `obscodes_extended.json` follows this shape).
fn is_mpc_code(s: &str) -> bool {
    s.len() == 3 && s.bytes().all(|b| b.is_ascii_alphanumeric())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_known_body_names() {
        assert!(matches!(Origin::parse("earth").unwrap(), Origin::Earth));
        assert!(matches!(
            Origin::parse("JUPITER_BARYCENTER").unwrap(),
            Origin::JupiterBarycenter
        ));
        assert!(matches!(
            Origin::parse("ssb").unwrap(),
            Origin::SolarSystemBarycenter
        ));
    }

    #[test]
    fn parses_three_char_mpc_codes() {
        match Origin::parse("I11").unwrap() {
            Origin::Observatory(c) => assert_eq!(c, "I11"),
            other => panic!("expected observatory, got {other:?}"),
        }
        // Numeric and mixed are fine.
        assert!(matches!(Origin::parse("500"), Ok(Origin::Observatory(_))));
        assert!(matches!(Origin::parse("W84"), Ok(Origin::Observatory(_))));
    }

    #[test]
    fn rejects_typoed_body_names() {
        let err = Origin::parse("earty").unwrap_err();
        assert!(err.to_string().contains("neither a known body name"));
    }

    #[test]
    fn rejects_non_three_char_strings() {
        assert!(Origin::parse("ABCD").is_err());
        assert!(Origin::parse("AB").is_err());
        assert!(Origin::parse("").is_err());
        // Non-alphanumeric.
        assert!(Origin::parse("I-1").is_err());
    }
}
