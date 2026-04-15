#!/usr/bin/env python3
"""Fetch Horizons reference data anchored at each object's SBDB fit epoch.

For each of 28 sample orbits:
  1. Query SBDB for the osculating-element epoch (JD TDB) and non-grav params.
  2. Query Horizons VECTORS from the solar-system barycenter at epoch + {0, 30, 90, 365} days,
     in ICRF equatorial (REF_PLANE='FRAME'), TIME_TYPE='TDB', VEC_CORR='NONE'.
  3. Query Horizons OBSERVER at the same epochs from I41, X05, 500 —
     QUANTITIES='1,3,20' → RA, Dec, dRA*cos(Dec), dDEC, delta, deldot.

Output: validation/horizons_reference_v2.json
"""

import csv
import datetime as dt
import json
import re
import sys
import time
from pathlib import Path

import requests

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"

DT_DAYS = [0, 30, 90, 365]
OBSERVERS = ["X05", "500", "I41"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def horizons_designation(targetname: str) -> str:
    """MPC/Horizons-compatible designation with trailing ';' for small-body."""
    m = re.match(r"^(\d+)\s", targetname)
    if m:
        return m.group(1) + ";"
    m = re.match(r"^(\d+[A-Z])", targetname)
    if m:
        return m.group(1) + ";"
    m = re.search(r"\(([^)]+)\)", targetname)
    if m:
        return m.group(1) + ";"
    return targetname + ";"


def http_get_json(url: str, params: dict, retries: int = 3) -> dict:
    last_err = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=90)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_err = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"GET {url} failed after {retries}: {last_err}")


def _parse_csv_block(result: str) -> list[list[str]]:
    """Extract CSV rows between $$SOE and $$EOE."""
    lines = result.splitlines()
    try:
        start = lines.index("$$SOE")
        end = lines.index("$$EOE")
    except ValueError as e:
        raise RuntimeError(
            f"Missing $$SOE/$$EOE:\n{result[:500]}"
        ) from e
    rows = []
    for line in lines[start + 1 : end]:
        line = line.strip()
        if not line:
            continue
        rows.append([c.strip() for c in line.split(",")])
    return rows


# ---------------------------------------------------------------------------
# SBDB
# ---------------------------------------------------------------------------

def fetch_sbdb(des: str) -> dict:
    """Return {epoch_jd_tdb, non_grav} from SBDB."""
    sstr = des.rstrip(";")
    data = http_get_json(SBDB_URL, {"sstr": sstr, "full-prec": "true"})
    orbit = data.get("orbit") or {}

    epoch_str = orbit.get("epoch")
    if epoch_str is None:
        raise RuntimeError(f"SBDB has no epoch for {sstr}")
    epoch_jd_tdb = float(epoch_str)

    pars = orbit.get("model_pars") or []
    by_name = {p["name"]: p for p in pars}

    non_grav = None
    if any(k in by_name for k in ("A1", "A2", "A3")):
        def fv(name):
            p = by_name.get(name)
            if p is None or p.get("value") in (None, ""):
                return None
            return float(p["value"])

        non_grav = {
            "a1": fv("A1") or 0.0,
            "a2": fv("A2") or 0.0,
            "a3": fv("A3") or 0.0,
            "aln": fv("ALN"),
            "nk": fv("NK"),
            "nm": fv("NM"),
            "nn": fv("NN"),
            "r0": fv("R0"),
        }

    return {"epoch_jd_tdb": epoch_jd_tdb, "non_grav": non_grav}


# ---------------------------------------------------------------------------
# Horizons VECTORS
# ---------------------------------------------------------------------------

def fetch_vectors(des: str, jd_tdbs: list[float]) -> list[dict]:
    """Barycentric ICRF equatorial state at each epoch. Returns list of dicts."""
    tlist = " ".join(f"{jd:.9f}" for jd in jd_tdbs)
    params = {
        "format": "json",
        "COMMAND": f"'{des}'",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": "'@0'",
        "REF_PLANE": "FRAME",
        "REF_SYSTEM": "ICRF",
        "OUT_UNITS": "AU-D",
        "VEC_TABLE": "2",
        "VEC_LABELS": "NO",
        "VEC_CORR": "NONE",
        "CSV_FORMAT": "YES",
        "TLIST": f"'{tlist}'",
        "TIME_TYPE": "TDB",
        "OBJ_DATA": "NO",
    }
    resp = http_get_json(HORIZONS_URL, params)
    rows = _parse_csv_block(resp["result"])
    if len(rows) != len(jd_tdbs):
        raise RuntimeError(
            f"VECTORS {des}: expected {len(jd_tdbs)} rows, got {len(rows)}"
        )
    out = []
    for jd, row in zip(jd_tdbs, rows):
        nums = [c for c in row if c != ""]
        floats = []
        for c in nums:
            try:
                floats.append(float(c))
            except ValueError:
                pass
        # Expect: JDTDB, Cal, X, Y, Z, VX, VY, VZ  → 8 numeric-ish but only 7 parseable
        # Actually with CSV_FORMAT and VEC_LABELS=NO: JDTDB, Calendar(str), X, Y, Z, VX, VY, VZ
        if len(floats) < 7:
            raise RuntimeError(f"VECTORS parse {des}: floats={floats} row={row}")
        # floats[0] = JDTDB, floats[1..7] = state
        state = floats[1:7]
        out.append({
            "jd_tdb": jd,
            "x": state[0], "y": state[1], "z": state[2],
            "vx": state[3], "vy": state[4], "vz": state[5],
        })
    return out


# ---------------------------------------------------------------------------
# Horizons OBSERVER
# ---------------------------------------------------------------------------

def fetch_observer(des: str, center: str, jd_tdbs: list[float]) -> list[dict]:
    """Astrometric RA/Dec + rates + range from Horizons OBSERVER mode."""
    tlist = " ".join(f"{jd:.9f}" for jd in jd_tdbs)
    params = {
        "format": "json",
        "COMMAND": f"'{des}'",
        "EPHEM_TYPE": "OBSERVER",
        "CENTER": f"'{center}'",
        "QUANTITIES": "'1,3,20'",
        "ANG_FORMAT": "DEG",
        "EXTRA_PREC": "YES",
        "CSV_FORMAT": "YES",
        "TLIST": f"'{tlist}'",
        "TIME_TYPE": "TT",
        "OBJ_DATA": "NO",
    }
    resp = http_get_json(HORIZONS_URL, params)
    rows = _parse_csv_block(resp["result"])
    if len(rows) != len(jd_tdbs):
        raise RuntimeError(
            f"OBSERVER {des}@{center}: expected {len(jd_tdbs)} rows, got {len(rows)}"
        )
    out = []
    for jd, row in zip(jd_tdbs, rows):
        # Parse numeric values from the row
        floats = []
        for c in row:
            try:
                floats.append(float(c))
            except ValueError:
                pass
        # Expected floats: RA, DEC, dRA*cos(Dec), d(DEC)/dt, delta, deldot
        # (may have extra leading JDTDB — check count)
        if len(floats) < 6:
            raise RuntimeError(
                f"OBSERVER parse {des}@{center}: floats={floats} row={row}"
            )
        # Last 6 floats are the ones we want
        ra_deg, dec_deg, dra_cosdec, ddec, delta_au, deldot_kms = floats[-6:]
        out.append({
            "jd_tdb": jd,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "dra_cosdec_arcsec_hr": dra_cosdec,
            "ddec_arcsec_hr": ddec,
            "delta_au": delta_au,
            "deldot_km_s": deldot_kms,
        })
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_targets() -> list[dict]:
    csv_path = Path(__file__).parent / "data" / "elements_sun_ec.csv"
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "object_id": row["targetname"],
                "horizons_des": horizons_designation(row["targetname"]),
            })
    return rows


def main():
    targets = load_targets()
    print(f"Fetching Horizons data for {len(targets)} orbits")
    print(f"DT offsets: {DT_DAYS} days")
    print(f"Observers: {OBSERVERS}")

    orbits = []
    for i, target in enumerate(targets, 1):
        name = target["object_id"]
        des = target["horizons_des"]
        print(f"\n[{i}/{len(targets)}] {name}  (des={des})")

        # --- SBDB ---
        try:
            sbdb = fetch_sbdb(des)
        except Exception as e:
            print(f"  SBDB FAILED: {e}")
            continue
        epoch = sbdb["epoch_jd_tdb"]
        non_grav = sbdb["non_grav"]
        print(f"  SBDB epoch: JD {epoch} TDB  (MJD {epoch - 2400000.5:.1f})")
        if non_grav:
            print(f"  non-grav: A1={non_grav['a1']:.3e} A2={non_grav['a2']:.3e} A3={non_grav['a3']:.3e}")
        time.sleep(0.3)

        # Target JDs
        target_jds = [epoch + dt for dt in DT_DAYS]

        # --- VECTORS (barycentric ICRF equatorial) ---
        try:
            vectors = fetch_vectors(des, target_jds)
            ic = vectors[0]
            print(f"  IC bary eq: [{ic['x']:+.6f}, {ic['y']:+.6f}, {ic['z']:+.6f}]")
        except Exception as e:
            print(f"  VECTORS FAILED: {e}")
            continue
        time.sleep(0.3)

        # --- OBSERVER ephemeris ---
        obs_data = {}
        ok = True
        for obs in OBSERVERS:
            try:
                eph = fetch_observer(des, obs, target_jds)
                obs_data[obs] = eph
                print(f"  {obs}: {len(eph)} epochs  (T0 RA/Dec = {eph[0]['ra_deg']:.6f}, {eph[0]['dec_deg']:.6f})")
            except Exception as e:
                print(f"  OBSERVER {obs} FAILED: {e}")
                ok = False
                break
            time.sleep(0.3)
        if not ok:
            continue

        orbits.append({
            "object_id": name,
            "horizons_des": des,
            "sbdb_epoch_jd_tdb": epoch,
            "non_grav": non_grav,
            "vectors_bary_eq": vectors,
            "observer_ephemeris": obs_data,
        })

    queried_at = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    output = {
        "metadata": {
            "queried_at": queried_at,
            "horizons_source": HORIZONS_URL,
            "sbdb_source": SBDB_URL,
            "dt_days": DT_DAYS,
            "observers": OBSERVERS,
        },
        "orbits": orbits,
    }

    out_path = Path(__file__).parent / "horizons_reference_v2.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nWrote {len(orbits)}/{len(targets)} orbits to {out_path}")
    print(f"Queried at: {queried_at}")


if __name__ == "__main__":
    main()
