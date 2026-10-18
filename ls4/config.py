"""Telescope, field-of-view, and scheduling parameters.

WARNING: these values drive the whole workflow (grid, scheduler, analysis).
Regenerate dependent outputs after changing them.
"""

from pathlib import Path

from astroplan import Observer
from astroplan.constraints import AirmassConstraint, AtNightConstraint, MoonSeparationConstraint
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils import iers
from regions import RectangleSkyRegion

# Use the bundled IERS tables so runs are reproducible and work offline.
iers.conf.auto_download = False
iers.conf.auto_max_age = None

LS4 = Observer.at_site("La Silla Observatory")

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

global_constraints = [AirmassConstraint(max=max_airmass),
                      AtNightConstraint.twilight_civil(),
                      MoonSeparationConstraint(min=15 * u.deg)]

half_field_offset_ra = FOV_length / 2
half_field_offset_dec = FOV_width / 2

region = RectangleSkyRegion(SkyCoord(0 * u.deg, 0 * u.deg), FOV_length, FOV_width, angle=90 * u.deg)

# Paths are relative to the repository root so commands work from any directory.
ROOT = Path(__file__).resolve().parent.parent
GRID_CSV = ROOT / "assets" / "LS4_field_grid.csv"
GRID_DB = ROOT / "assets" / "LS4_field_grid.db"
PLANS_DIR = ROOT / "plans"
OUTPUTS_DIR = ROOT / "outputs"

# Plan rows with these targets are schedule gaps, not science pointings.
NON_SCIENCE_TARGETS = ("Unused Time", "TransitionBlock")
