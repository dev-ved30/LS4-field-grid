from astroplan import Observer
from astropy import units as u
from astropy.coordinates import SkyCoord
from regions import RectangleSkyRegion
from astroplan.constraints import AirmassConstraint, TimeConstraint, AtNightConstraint

# WARNING: This file contains constants and parameters for the LS4 field grid and observing plan. Be cautious when modifying these values, as they will affect the entire workflow.

LS4 = Observer.at_site("La Silla Observatory")

# Telescope specs
slew_rate = 0.8*u.deg/u.second # Slew rate in degrees per second
read_out = 40 * u.second # Readout time for LS4. This is often synced with the slews.
exp = 60*u.second # Exposure time per pointing.

FOV_length = 17826.13*u.arcsec
FOV_width = 16924.67*u.arcsec
FOV = FOV_width**2
FOV_g = FOV/4

# Observing constraints
min_declination = -60*u.deg  # Minimum declination limit for fields
max_airmass = 1.6
block_duration = 30*u.minute

tolerance = 4 # adding some time to account for slew overheads and other inefficiencies. This is a simple heuristic and can be optimized further.
fields_per_block = int(block_duration/(exp + read_out)) - tolerance # Number of fields that can be observed in a 30 minute block, accounting for exposure time and readout time. This is a simple calculation and can be optimized further by considering slews and other factors.

# Universal constraints for all targets
global_constraints = [AirmassConstraint(max=max_airmass),
                      AtNightConstraint.twilight_civil()]

half_field_offset_ra = FOV_length/2
half_field_offset_dec = FOV_width/2


region = RectangleSkyRegion(SkyCoord(0 * u.deg, 0 * u.deg), FOV_length, FOV_width, angle=90*u.deg)

LS4_field_grid_path = "assets/LS4_field_grid.csv"