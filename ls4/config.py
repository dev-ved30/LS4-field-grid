"""Telescope, field-of-view, and scheduling parameters.

WARNING: these values drive the whole workflow (grid, scheduler, analysis).
Regenerate dependent outputs after changing them.
"""

from pathlib import Path

from astroplan import Observer
from astroplan.constraints import AirmassConstraint, AtNightConstraint, MoonSeparationConstraint
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.utils import iers
from regions import RectangleSkyRegion

# Use the bundled IERS tables so runs are reproducible and work offline.
iers.conf.auto_download = False
iers.conf.auto_max_age = None

# Same coordinates as Astropy's La Silla registry entry, fixed here so a clean
# production runner does not need a site-registry download just to import LS4.
LS4 = Observer(name="La Silla Observatory", location=EarthLocation.from_geodetic(
    lon=-70.73 * u.deg, lat=-29.256666666666668 * u.deg, height=2347 * u.m))

# Telescope specs
slew_rate = 0.8 * u.deg / u.second  # Assumed telescope slew rate; replace with a measured value when available.
read_out = 40 * u.second  # Readout time for LS4. This is often synced with the slews.
exp = 60 * u.second  # Exposure time per pointing.

FOV_length = 17826.13 * u.arcsec
FOV_width = 16924.67 * u.arcsec
FOV = FOV_width**2
FOV_g = FOV / 4

# Observing constraints
min_declination = -60 * u.deg  # Minimum declination limit for fields
max_airmass = 1.6
default_block_duration = 30 * u.minute  # Also the minimum initial-to-dither revisit time.
max_cluster_radius = 16 * u.deg  # Keep a block's field centers near its first pointing.

# Survey programs (grid `program_id`).
SOLE = 0  # Stellar Oscillations, Lensing, and Eruptions: the Galactic plane and bulge.
EXTRAGALACTIC = 1
sole_plane = {"l": (-72, 36), "b": 4}  # degrees; l is wrapped to [-180, 180).
sole_bulge = {"l": 20, "b": 8}  # |l| <= 20, |b| <= 8
sole_fraction = 0.25  # Target share of each night's pairs.
sole_cadence = 1 * u.day  # Nominal cadence for SOLE fields.

global_constraints = [AirmassConstraint(max=max_airmass),
                      AtNightConstraint.twilight_civil(),
                      MoonSeparationConstraint(min=15 * u.deg)]

# Focal-plane bands in 2x2 quadrants: i on one diagonal, z top right, g bottom left.
# Keys are (RA side, Dec side), +1 toward increasing RA / Dec. The figure's left/right
# orientation on the sky is assumed; the alternation below works either way, but
# per-band maps (ls4.bands) mirror if the camera is flipped in RA.
BANDS = ("g", "i", "z")
QUADRANT_BANDS = {(-1, 1): "i", (1, 1): "z", (-1, -1): "g", (1, -1): "i"}
# An initial pointing plus a dither half a grid step in RA gives every point i and one of
# g or z. Shifting the pair half a grid step in Dec on alternate nights (dither phase 1)
# swaps g and z, so every part of a field gets i nightly and g and z every other night.

region = RectangleSkyRegion(SkyCoord(0 * u.deg, 0 * u.deg), FOV_length, FOV_width, angle=90 * u.deg)

# Paths are relative to the repository root so commands work from any directory.
ROOT = Path(__file__).resolve().parent.parent
GRID_CSV = ROOT / "assets" / "LS4_field_grid.csv"
GRID_DB = ROOT / "assets" / "LS4_field_grid.db"
PLANS_DIR = ROOT / "plans"
OUTPUTS_DIR = ROOT / "outputs"

# Plan rows with these targets are schedule gaps, not science pointings.
NON_SCIENCE_TARGETS = ("Unused Time", "TransitionBlock")
