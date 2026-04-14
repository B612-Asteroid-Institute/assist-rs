#!/usr/bin/env python3
"""Generate reference propagation data using Python REBOUND + ASSIST.

This script propagates the 27 sample orbits from adam_core helpers using
the same ASSIST/REBOUND setup that assist-rs uses, then writes the
results to JSON for comparison in the Rust validation test.

Usage:
    python generate_reference.py

Requires:
    - rebound, assist Python packages
    - adam_core (for sample orbits and coordinate transforms)
    - ASSIST_PLANETS_PATH and ASSIST_ASTEROIDS_PATH environment variables
"""

import json
import math
import os
import sys

import numpy as np
import rebound
import assist

# J2000 obliquity (IAU 2006) — must match assist-rs constants
OBLIQUITY_J2000 = 0.40909280422  # radians
COS_EPS = math.cos(OBLIQUITY_J2000)
SIN_EPS = math.sin(OBLIQUITY_J2000)

# ASSIST default force flags (must match assist-rs)
# SUN=0x01 | PLANETS=0x02 | ASTEROIDS=0x04 | EARTH_HARMONICS=0x10 | SUN_HARMONICS=0x20 | GR_EIH=0x40
ASSIST_FORCES_DEFAULT = 0x01 | 0x02 | 0x04 | 0x10 | 0x20 | 0x40


def ecliptic_to_equatorial(state):
    """Rotate a 6-element state from ecliptic to equatorial (matching assist-rs)."""
    x, y, z, vx, vy, vz = state
    return [
        x,
        COS_EPS * y - SIN_EPS * z,
        SIN_EPS * y + COS_EPS * z,
        vx,
        COS_EPS * vy - SIN_EPS * vz,
        SIN_EPS * vy + COS_EPS * vz,
    ]


def equatorial_to_ecliptic(state):
    """Rotate a 6-element state from equatorial to ecliptic (matching assist-rs)."""
    x, y, z, vx, vy, vz = state
    return [
        x,
        COS_EPS * y + SIN_EPS * z,
        -SIN_EPS * y + COS_EPS * z,
        vx,
        COS_EPS * vy + SIN_EPS * vz,
        -SIN_EPS * vy + COS_EPS * vz,
    ]


def mjd_to_assist_time(mjd_tdb, jd_ref):
    """Convert MJD TDB to ASSIST time (days from jd_ref)."""
    return (mjd_tdb + 2400000.5) - jd_ref


def propagate_orbit(ephem, state_ecl, epoch_mjd, target_mjds, non_grav=None):
    """Propagate a single orbit using REBOUND + ASSIST.

    Matches the exact setup in assist-rs::assist_propagate:
    1. Convert heliocentric ecliptic J2000 → barycentric equatorial ICRF
    2. Create simulation with IAS15 + ASSIST forces
    3. Integrate to each target epoch
    4. Convert back to heliocentric ecliptic J2000

    Parameters
    ----------
    ephem : assist.Ephem
        ASSIST ephemeris data.
    state_ecl : list[float]
        Initial state [x, y, z, vx, vy, vz] in heliocentric ecliptic J2000 (AU, AU/day).
    epoch_mjd : float
        Initial epoch in MJD TDB.
    target_mjds : list[float]
        Target epochs in MJD TDB (must be sorted).
    non_grav : tuple[float, float, float] | None
        Optional (A1, A2, A3) non-gravitational acceleration components
        (Marsden-Sekanina model). When provided, the NON_GRAVITATIONAL
        force flag is enabled and particle_params is set.

    Returns
    -------
    list[dict]
        One dict per target epoch with 'epoch', 'state' keys.
    """
    jd_ref = ephem.jd_ref
    t0 = mjd_to_assist_time(epoch_mjd, jd_ref)

    # Convert heliocentric ecliptic → barycentric equatorial ICRF
    eq_state = ecliptic_to_equatorial(state_ecl)

    # Get Sun's barycentric position at t0
    sun = ephem.get_particle("sun", t0)
    bary_state = [
        eq_state[0] + sun.x,
        eq_state[1] + sun.y,
        eq_state[2] + sun.z,
        eq_state[3] + sun.vx,
        eq_state[4] + sun.vy,
        eq_state[5] + sun.vz,
    ]

    # Create simulation — match assist-rs setup exactly
    sim = rebound.Simulation()
    sim.t = t0
    extras = assist.Extras(sim, ephem)

    if non_grav is not None:
        extras.forces = extras.forces + ["NON_GRAVITATIONAL"]
        extras.particle_params = np.array(list(non_grav), dtype=np.float64)

    # Add test particle
    sim.add(
        x=bary_state[0],
        y=bary_state[1],
        z=bary_state[2],
        vx=bary_state[3],
        vy=bary_state[4],
        vz=bary_state[5],
        m=0.0,
    )

    sim.N_active = 0
    sim.exact_finish_time = 1

    results = []
    for target_mjd in target_mjds:
        t_target = mjd_to_assist_time(target_mjd, jd_ref)
        sim.integrate(t_target)

        p = sim.particles[0]
        bary_eq = [p.x, p.y, p.z, p.vx, p.vy, p.vz]

        # Convert barycentric equatorial → heliocentric ecliptic
        sun_t = ephem.get_particle("sun", t_target)
        helio_eq = [
            bary_eq[0] - sun_t.x,
            bary_eq[1] - sun_t.y,
            bary_eq[2] - sun_t.z,
            bary_eq[3] - sun_t.vx,
            bary_eq[4] - sun_t.vy,
            bary_eq[5] - sun_t.vz,
        ]
        helio_ecl = equatorial_to_ecliptic(helio_eq)

        results.append({
            "epoch": target_mjd,
            "state": helio_ecl,
        })

    return results


def load_sample_orbits():
    """Load the 27 sample orbits from adam_core helpers data.

    Returns list of dicts with 'object_id', 'epoch_mjd', 'state' keys.
    State is heliocentric ecliptic J2000 [x, y, z, vx, vy, vz] (AU, AU/day).
    """
    import csv

    # Read directly from the CSV (elements_sun_ec.csv) to avoid
    # adam_core import complexity. This has the Cartesian states.
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    csv_path = os.path.join(data_dir, "elements_sun_ec.csv")

    if not os.path.exists(csv_path):
        # Try adam_core installed package location
        from importlib.resources import files as pkg_files
        data_pkg = pkg_files("adam_core.utils.helpers.data")
        csv_path = str(data_pkg.joinpath("elements_sun_ec.csv"))

    orbits = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            orbits.append({
                "object_id": row["targetname"],
                "epoch_mjd": float(row["mjd_tdb"]),
                "state": [
                    float(row["x"]),
                    float(row["y"]),
                    float(row["z"]),
                    float(row["vx"]),
                    float(row["vy"]),
                    float(row["vz"]),
                ],
            })
    return orbits


def main():
    planets_path = os.environ.get("ASSIST_PLANETS_PATH")
    asteroids_path = os.environ.get("ASSIST_ASTEROIDS_PATH")

    if not planets_path or not asteroids_path:
        print("Error: Set ASSIST_PLANETS_PATH and ASSIST_ASTEROIDS_PATH", file=sys.stderr)
        sys.exit(1)

    ephem = assist.Ephem(planets_path, asteroids_path)
    print(f"Loaded ephemeris: jd_ref={ephem.jd_ref}")
    print(f"REBOUND version: {rebound.__version__}")
    print(f"ASSIST version: {assist.__version__}")

    orbits = load_sample_orbits()
    print(f"Loaded {len(orbits)} sample orbits")

    # Propagation settings: 30 days forward, 10 epochs
    dt_days = 30.0
    n_epochs = 10

    # Non-gravitational acceleration parameters for the non-grav test suite.
    # Comet-like transverse drag (A2 ≠ 0) — arbitrary but numerically exercises
    # the NON_GRAVITATIONAL force path in ASSIST.
    nongrav_params = (0.0, 1e-10, 0.0)

    reference_data = {
        "metadata": {
            "rebound_version": rebound.__version__,
            "assist_version": assist.__version__,
            "jd_ref": ephem.jd_ref,
            "cos_eps": COS_EPS,
            "sin_eps": SIN_EPS,
            "propagation_days": dt_days,
            "n_epochs": n_epochs,
            "nongrav_params": list(nongrav_params),
        },
        "orbits": [],
        "nongrav_orbits": [],
    }

    # ── Gravity-only pass ────────────────────────────────────────────────
    print("Gravity-only propagation:")
    for i, orbit in enumerate(orbits):
        epoch = orbit["epoch_mjd"]

        # Check that epoch is within ephemeris coverage
        t0 = mjd_to_assist_time(epoch, ephem.jd_ref)
        t_end = mjd_to_assist_time(epoch + dt_days, ephem.jd_ref)

        # Skip orbits whose epochs are outside the ephemeris range
        # DE440 covers ~1550 to ~2650 CE, which is MJD ~-100840 to ~301558
        # But the asteroid ephemeris may have different bounds
        try:
            ephem.get_particle("sun", t0)
            ephem.get_particle("sun", t_end)
        except Exception as e:
            print(f"  Skipping {orbit['object_id']}: epoch outside ephemeris range ({e})")
            continue

        target_epochs = [epoch + dt_days * (j + 1) / n_epochs for j in range(n_epochs)]

        try:
            results = propagate_orbit(ephem, orbit["state"], epoch, target_epochs)
        except Exception as e:
            print(f"  Skipping {orbit['object_id']}: propagation failed ({e})")
            continue

        entry = {
            "object_id": orbit["object_id"],
            "epoch_mjd": epoch,
            "initial_state": orbit["state"],
            "propagated": results,
        }
        reference_data["orbits"].append(entry)
        print(f"  [{i+1}/{len(orbits)}] {orbit['object_id']}: {len(results)} epochs OK")

    # ── Non-gravitational pass ───────────────────────────────────────────
    print(f"\nNon-gravitational propagation (A1, A2, A3) = {nongrav_params}:")
    for i, orbit in enumerate(orbits):
        epoch = orbit["epoch_mjd"]
        t0 = mjd_to_assist_time(epoch, ephem.jd_ref)
        t_end = mjd_to_assist_time(epoch + dt_days, ephem.jd_ref)

        try:
            ephem.get_particle("sun", t0)
            ephem.get_particle("sun", t_end)
        except Exception as e:
            print(f"  Skipping {orbit['object_id']}: epoch outside ephemeris range ({e})")
            continue

        target_epochs = [epoch + dt_days * (j + 1) / n_epochs for j in range(n_epochs)]

        try:
            results = propagate_orbit(
                ephem, orbit["state"], epoch, target_epochs, non_grav=nongrav_params
            )
        except Exception as e:
            print(f"  Skipping {orbit['object_id']}: propagation failed ({e})")
            continue

        entry = {
            "object_id": orbit["object_id"],
            "epoch_mjd": epoch,
            "initial_state": orbit["state"],
            "propagated": results,
        }
        reference_data["nongrav_orbits"].append(entry)
        print(f"  [{i+1}/{len(orbits)}] {orbit['object_id']}: {len(results)} epochs OK")

    output_path = os.path.join(os.path.dirname(__file__), "reference_data.json")
    with open(output_path, "w") as f:
        json.dump(reference_data, f, indent=2)

    print(
        f"\nWrote {len(reference_data['orbits'])} gravity-only orbits and "
        f"{len(reference_data['nongrav_orbits'])} non-grav orbits to {output_path}"
    )


if __name__ == "__main__":
    main()
