"""Build the LS4 field grid and read/write its scheduling database."""

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from ls4.config import (EXTRAGALACTIC, FOV_g, GRID_CSV, GRID_DB, SOLE, min_declination, sole_bulge,
                        sole_plane)


def program_ids(l_deg, b_deg):
    """SOLE for fields centered in the Galactic plane or bulge, extragalactic otherwise."""
    l = (np.asarray(l_deg) + 180) % 360 - 180
    b = np.abs(np.asarray(b_deg))
    plane = (sole_plane["l"][0] <= l) & (l <= sole_plane["l"][1]) & (b <= sole_plane["b"])
    bulge = (np.abs(l) <= sole_bulge["l"]) & (b <= sole_bulge["b"])
    return np.where(plane | bulge, SOLE, EXTRAGALACTIC)


def build_grid():
    """Tessellate the sky and keep every other field in RA and Dec."""
    from m4opt import skygrid

    vertices = skygrid.sinusoidal(FOV_g)
    df = pd.DataFrame({"ra_deg": vertices.ra.to_value(u.deg),
                       "dec_deg": vertices.dec.to_value(u.deg)})

    # Name fields by declination ring and RA position within the ring.
    dec_idx = df["dec_deg"].rank(method="dense").astype(int) - 1
    ra_idx = df.groupby("dec_deg")["ra_deg"].rank(method="first").astype(int) - 1

    galactic = SkyCoord(df["ra_deg"].to_numpy() * u.deg, df["dec_deg"].to_numpy() * u.deg).galactic
    l, b = galactic.l.to_value(u.deg), galactic.b.to_value(u.deg)
    # TODO still need to assign high cadence extragalactic fields, and assign program IDs to them.

    df["Field Name"] = ra_idx.astype(str) + "_" + dec_idx.astype(str)
    df["l"] = l
    df["b"] = b
    df["program_id"] = program_ids(l, b)
    df["last_scheduled_mjd"] = 0.0
    keep = (dec_idx % 2 == 0) & (ra_idx % 2 == 0)
    return df[keep].reset_index(drop=True)


def write_db(grid, db_path=GRID_DB, overwrite=False):
    """Store the grid in the `grid` table used by the scheduler."""
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as db:
        grid.to_sql("grid", db, if_exists="replace" if overwrite else "fail", index=False)


def fresh_db(db_path, grid_csv=GRID_CSV):
    """Create an isolated database with no scheduling history (for simulations)."""
    grid = pd.read_csv(grid_csv)
    grid["last_scheduled_mjd"] = 0.0
    write_db(grid, db_path, overwrite=True)
    return Path(db_path)


def load_schedulable(db_path=GRID_DB):
    """Fields above the declination limit, with their scheduling history."""
    with sqlite3.connect(db_path) as db:
        grid = pd.read_sql_query("SELECT * FROM grid WHERE dec_deg > ?", db,
                                 params=(min_declination.to_value(u.deg),))
    if grid.empty:
        raise ValueError(f"No schedulable fields in {db_path}")
    # Derive programs from the footprint so databases built before it changed stay valid.
    grid["program_id"] = program_ids(grid["l"], grid["b"])
    return grid


def load_history(db_path):
    with sqlite3.connect(db_path) as db:
        return pd.read_sql_query('SELECT "Field Name" AS base_field, last_scheduled_mjd FROM grid', db)


def record_pairs(last_mjd_by_field, db_path=GRID_DB):
    """Advance `last_scheduled_mjd` for the given {field name: mjd} mapping."""
    with sqlite3.connect(db_path) as db:
        db.executemany('UPDATE grid SET last_scheduled_mjd = ? WHERE "Field Name" = ?',
                       [(float(mjd), name) for name, mjd in last_mjd_by_field.items()])


def load_grid_csv(path=GRID_CSV):
    """Read the complete field universe, including unscheduled fields."""
    grid = pd.read_csv(path).rename(columns={"Field Name": "base_field", "ra_deg": "ra", "dec_deg": "dec"})
    if grid["base_field"].isna().any() or grid["base_field"].duplicated().any():
        raise ValueError(f"Grid field names must be unique and nonempty: {path}")
    return grid
