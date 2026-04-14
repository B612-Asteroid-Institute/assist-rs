//! Data file management for ASSIST.
//!
//! ASSIST needs three files to run end-to-end:
//!
//! 1. `de440.bsp` — NAIF DE440 planetary ephemeris SPK
//! 2. `sb441-n16.bsp` — JPL SB441 N=16 small-body perturber ephemeris
//! 3. `obscodes_extended.json` — MPC observatory codes (for observer lookups)
//!
//! [`DataManager`] downloads these to a local cache directory on demand and
//! returns resolved paths for [`crate::Ephemeris::from_paths`] and
//! [`crate::ObservatoryTable::from_json`].
//!
//! # Default data directory
//!
//! `$ASSIST_DATA_DIR` if set, otherwise `$XDG_CACHE_HOME/assist-rs/` or
//! `~/.cache/assist-rs/`.
//!
//! # Example
//!
//! ```no_run
//! use assist_rs::data::DataManager;
//!
//! let dm = DataManager::new();
//! let paths = dm.ensure_ready()?;
//! # Ok::<(), assist_rs::data::DataError>(())
//! ```
//!
//! Each downloaded file gets a `<filename>.meta.json` sidecar with MD5 hash,
//! Content-Length, and Last-Modified from the HTTP response. On subsequent
//! runs, non-static files are checked via HEAD and re-downloaded if changed.

use std::fmt;
use std::fs::{self, File};
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

// ─── Kernel catalog ─────────────────────────────────────────────────────────

struct KernelEntry {
    filename: &'static str,
    url: &'static str,
    gzipped: bool,
    is_static: bool,
}

const DEFAULT_KERNELS: &[KernelEntry] = &[
    KernelEntry {
        filename: "de440.bsp",
        url: "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
        gzipped: false,
        is_static: true,
    },
    KernelEntry {
        filename: "sb441-n16.bsp",
        url: "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp",
        gzipped: false,
        is_static: true,
    },
    KernelEntry {
        filename: "obscodes_extended.json",
        url: "https://minorplanetcenter.net/Extended_Files/obscodes_extended.json.gz",
        gzipped: true,
        is_static: false,
    },
];

// ─── Sidecar metadata ──────────────────────────────────────────────────────

#[derive(Debug, Serialize, Deserialize)]
struct FileMeta {
    url: String,
    downloaded_at: u64,
    content_length: Option<u64>,
    last_modified: Option<String>,
    md5: String,
}

// ─── AssistDataPaths ────────────────────────────────────────────────────────

/// Resolved paths to the three files ASSIST needs.
#[derive(Debug, Clone)]
pub struct AssistDataPaths {
    /// DE440 planetary ephemeris SPK.
    pub planets: PathBuf,
    /// SB441 N=16 small-body perturber SPK.
    pub asteroids: PathBuf,
    /// MPC observatory codes JSON (decompressed).
    pub obscodes: PathBuf,
}

// ─── DataError ──────────────────────────────────────────────────────────────

/// Errors from the data manager.
#[derive(Debug)]
pub enum DataError {
    /// Required files are missing (offline mode).
    MissingFiles(Vec<String>),
    /// HTTP request failed.
    Http(String),
    /// I/O error.
    Io(io::Error),
}

impl fmt::Display for DataError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingFiles(files) => write!(f, "missing data files: {}", files.join(", ")),
            Self::Http(msg) => write!(f, "HTTP error: {msg}"),
            Self::Io(e) => write!(f, "I/O error: {e}"),
        }
    }
}

impl std::error::Error for DataError {}

// ─── DataManager ────────────────────────────────────────────────────────────

/// Manages downloading and caching of the files ASSIST needs.
///
/// Files are stored in a local data directory. Default location is
/// `$ASSIST_DATA_DIR` if set, otherwise `$XDG_CACHE_HOME/assist-rs/` or
/// `~/.cache/assist-rs/`.
pub struct DataManager {
    data_dir: PathBuf,
}

impl Default for DataManager {
    fn default() -> Self {
        Self::new()
    }
}

impl DataManager {
    /// Create a `DataManager` with the default data directory.
    pub fn new() -> Self {
        let data_dir = std::env::var("ASSIST_DATA_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|_| {
                if let Ok(xdg) = std::env::var("XDG_CACHE_HOME") {
                    PathBuf::from(xdg).join("assist-rs")
                } else {
                    let home = std::env::var("HOME").expect("HOME environment variable not set");
                    PathBuf::from(home).join(".cache").join("assist-rs")
                }
            });
        Self { data_dir }
    }

    /// Create a `DataManager` with a custom data directory.
    pub fn with_dir(dir: impl Into<PathBuf>) -> Self {
        Self {
            data_dir: dir.into(),
        }
    }

    /// The data directory path.
    pub fn data_dir(&self) -> &Path {
        &self.data_dir
    }

    fn paths(&self) -> AssistDataPaths {
        AssistDataPaths {
            planets: self.data_dir.join("de440.bsp"),
            asteroids: self.data_dir.join("sb441-n16.bsp"),
            obscodes: self.data_dir.join("obscodes_extended.json"),
        }
    }

    /// Return paths if all three files exist. No network access.
    pub fn offline(&self) -> Result<AssistDataPaths, DataError> {
        let missing: Vec<String> = DEFAULT_KERNELS
            .iter()
            .filter(|e| !self.data_dir.join(e.filename).exists())
            .map(|e| e.filename.to_string())
            .collect();
        if !missing.is_empty() {
            return Err(DataError::MissingFiles(missing));
        }
        Ok(self.paths())
    }

    /// Ensure all three files exist, downloading any that are missing or stale.
    ///
    /// - Missing → download.
    /// - Exists and static (`de440.bsp`, `sb441-n16.bsp`) → skip.
    /// - Exists and non-static with sidecar (`obscodes_extended.json`) → HEAD
    ///   to check Content-Length / Last-Modified; re-download if either differs.
    /// - Exists without sidecar → skip (assume valid).
    pub fn ensure_ready(&self) -> Result<AssistDataPaths, DataError> {
        fs::create_dir_all(&self.data_dir).map_err(DataError::Io)?;

        for entry in DEFAULT_KERNELS {
            let path = self.data_dir.join(entry.filename);
            let meta_path = self.data_dir.join(format!("{}.meta.json", entry.filename));

            if path.exists() {
                if entry.is_static {
                    continue;
                }
                if meta_path.exists() {
                    if let Ok(meta) = read_meta(&meta_path) {
                        match is_stale(entry.url, &meta) {
                            Ok(true) => {
                                eprintln!("Updating {} (remote changed)...", entry.filename);
                                download(entry, &path, &meta_path)?;
                            }
                            Ok(false) => {}
                            Err(e) => {
                                eprintln!(
                                    "Warning: staleness check failed for {}: {e}",
                                    entry.filename
                                );
                            }
                        }
                    }
                }
            } else {
                eprintln!("Downloading {}...", entry.filename);
                download(entry, &path, &meta_path)?;
            }
        }

        Ok(self.paths())
    }

    /// Remove the data directory and all its contents.
    pub fn clean(&self) -> Result<(), DataError> {
        if self.data_dir.exists() {
            fs::remove_dir_all(&self.data_dir).map_err(DataError::Io)?;
        }
        Ok(())
    }
}

// ─── Private helpers ────────────────────────────────────────────────────────

fn is_stale(url: &str, meta: &FileMeta) -> Result<bool, DataError> {
    let response = ureq::head(url)
        .call()
        .map_err(|e| DataError::Http(format!("HEAD {url}: {e}")))?;

    let remote_length: Option<u64> = response
        .headers()
        .get("content-length")
        .and_then(|v| v.to_str().ok())
        .and_then(|v| v.parse().ok());

    let remote_modified: Option<&str> = response
        .headers()
        .get("last-modified")
        .and_then(|v| v.to_str().ok());

    if let (Some(remote), Some(local)) = (remote_length, meta.content_length) {
        if remote != local {
            return Ok(true);
        }
    }

    if let (Some(remote), Some(local)) = (remote_modified, meta.last_modified.as_deref()) {
        if remote != local {
            return Ok(true);
        }
    }

    Ok(false)
}

fn download(entry: &KernelEntry, path: &Path, meta_path: &Path) -> Result<(), DataError> {
    let response = ureq::get(entry.url)
        .call()
        .map_err(|e| DataError::Http(format!("GET {}: {e}", entry.url)))?;

    let content_length: Option<u64> = response
        .headers()
        .get("content-length")
        .and_then(|v| v.to_str().ok())
        .and_then(|v| v.parse().ok());

    let last_modified: Option<String> = response
        .headers()
        .get("last-modified")
        .and_then(|v| v.to_str().ok())
        .map(|v| v.to_string());

    if let Some(size) = content_length {
        eprintln!("  {} ({:.1} MB)", entry.filename, size as f64 / 1_048_576.0);
    }

    let tmp_path = path.with_extension("tmp");
    {
        let mut body = response.into_body();
        let file = File::create(&tmp_path).map_err(DataError::Io)?;
        let mut writer = BufWriter::new(file);
        if entry.gzipped {
            let mut decoder = flate2::read::GzDecoder::new(body.as_reader());
            io::copy(&mut decoder, &mut writer).map_err(DataError::Io)?;
        } else {
            io::copy(&mut body.as_reader(), &mut writer).map_err(DataError::Io)?;
        }
        writer.flush().map_err(DataError::Io)?;
    }

    let md5_hex = compute_md5(&tmp_path)?;

    fs::rename(&tmp_path, path).map_err(DataError::Io)?;

    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs();
    let meta = FileMeta {
        url: entry.url.to_string(),
        downloaded_at: now,
        content_length,
        last_modified,
        md5: md5_hex,
    };
    let json =
        serde_json::to_string_pretty(&meta).map_err(|e| DataError::Io(io::Error::other(e)))?;
    fs::write(meta_path, json).map_err(DataError::Io)?;

    Ok(())
}

fn read_meta(path: &Path) -> Result<FileMeta, DataError> {
    let content = fs::read_to_string(path).map_err(DataError::Io)?;
    serde_json::from_str(&content)
        .map_err(|e| DataError::Io(io::Error::new(io::ErrorKind::InvalidData, e)))
}

fn compute_md5(path: &Path) -> Result<String, DataError> {
    let mut file = File::open(path).map_err(DataError::Io)?;
    let mut context = md5::Context::new();
    let mut buffer = [0u8; 65536];
    loop {
        let n = file.read(&mut buffer).map_err(DataError::Io)?;
        if n == 0 {
            break;
        }
        context.consume(&buffer[..n]);
    }
    Ok(format!("{:x}", context.compute()))
}
