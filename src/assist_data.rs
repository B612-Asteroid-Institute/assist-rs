//! [`AssistData`] — bundle of the data resources every high-level entry
//! point needs: the JPL SPK [`Ephemeris`], and optionally an MPC
//! [`ObservatoryTable`] (which itself may carry an
//! [`EarthOrientation`](crate::earth_orientation::EarthOrientation) for
//! sub-mas observatory rotation).
//!
//! Replacing the old `(ephem, obs_table)` argument pair across the public
//! API with a single `&AssistData` argument keeps call sites short and
//! matches how real workflows load these resources: once at startup, all
//! together, held for the lifetime of the process.

use crate::observatory::ObservatoryTable;
use crate::wrappers::Ephemeris;

/// Bundle of the data resources the high-level `assist_*` functions need.
///
/// Minimal usage — propagation only:
///
/// ```no_run
/// # use assist_rs::{AssistData, Ephemeris};
/// # use std::path::Path;
/// # fn f() -> assist_rs::Result<()> {
/// let ephem = Ephemeris::from_paths(Path::new("de440.bsp"), Path::new("sb441-n16.bsp"))?;
/// let data = AssistData::new(ephem);
/// # Ok(())
/// # }
/// ```
///
/// With observatory support (needed for `Origin::Observatory(..)` queries
/// and for ground-based `Observer`s in `assist_generate_ephemeris`):
///
/// ```no_run
/// # use assist_rs::{AssistData, Ephemeris, ObservatoryTable};
/// # use std::path::Path;
/// # fn f() -> assist_rs::Result<()> {
/// let ephem = Ephemeris::from_paths(Path::new("de440.bsp"), Path::new("sb441-n16.bsp"))?;
/// let obs = ObservatoryTable::from_json(Path::new("obscodes_extended.json"))?;
/// let data = AssistData::new(ephem).with_observatory(obs);
/// # Ok(())
/// # }
/// ```
///
/// With Earth-orientation precision for ground observatories, attach the
/// EOP kernels to the observatory table before bundling:
///
/// ```no_run
/// # use assist_rs::{AssistData, Ephemeris, ObservatoryTable};
/// # use assist_rs::earth_orientation::EarthOrientation;
/// # use std::path::Path;
/// # fn f() -> Result<(), Box<dyn std::error::Error>> {
/// # let ephem = Ephemeris::from_paths(Path::new("de440.bsp"), Path::new("sb441-n16.bsp"))?;
/// let eo = EarthOrientation::from_paths(&[Path::new("earth_latest_high_prec.bpc")])?;
/// let obs = ObservatoryTable::from_json(Path::new("obscodes_extended.json"))?
///     .with_earth_orientation(eo);
/// let data = AssistData::new(ephem).with_observatory(obs);
/// # Ok(())
/// # }
/// ```
pub struct AssistData {
    /// JPL SPK ephemeris (loaded once from `de440.bsp` + `sb441-n16.bsp`).
    pub ephem: Ephemeris,
    /// MPC observatory codes + (optional) Earth orientation.
    ///
    /// `None` when the caller only needs propagation or state queries for
    /// named bodies / the SSB. Required for
    /// [`crate::Origin::Observatory`] state queries and for
    /// observatory-based [`crate::Observer`]s.
    pub observatory: Option<ObservatoryTable>,
}

impl AssistData {
    /// Bundle with no observatory support. Sufficient for
    /// [`crate::assist_propagate`] / [`crate::assist_propagate_single`],
    /// [`crate::PropagatorPool`], and any [`crate::assist_get_state`] /
    /// [`crate::assist_generate_ephemeris`] call whose origins / observers
    /// are all named bodies or the SSB.
    pub fn new(ephem: Ephemeris) -> Self {
        Self {
            ephem,
            observatory: None,
        }
    }

    /// Attach an [`ObservatoryTable`]. Enables `Origin::Observatory(..)`
    /// state queries and ground-based observers in ephemeris generation.
    pub fn with_observatory(mut self, observatory: ObservatoryTable) -> Self {
        self.observatory = Some(observatory);
        self
    }
}
