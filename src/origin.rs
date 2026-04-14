//! Origin types for state queries and observer positions.

use crate::ffi;

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
    /// Parse a string into an `Origin`. Named bodies are matched case-insensitively;
    /// anything else is treated as an MPC observatory code.
    ///
    /// Accepts both short names ("earth", "jupiter") and explicit barycenter names
    /// ("jupiter_barycenter", "ssb").
    pub fn parse(s: &str) -> Self {
        match s.to_lowercase().as_str() {
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
            _ => Origin::Observatory(s.to_string()),
        }
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
