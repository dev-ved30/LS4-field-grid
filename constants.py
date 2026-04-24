import astroplan
import pandas as pd

from astroplan import Observer
from astroplan.constraints import AirmassConstraint, TimeConstraint, AtNightConstraint
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan.scheduling import Transitioner, SequentialScheduler, Schedule
from astroplan import ObservingBlock
from regions import RectangleSkyRegion


LS4 = Observer.at_site("La Silla Observatory")

# Telescope specs
slew_rate = .8*u.deg/u.second # Slew rate in degrees per second
read_out = 40 * u.second
exp = 60*u.second

FOV_length = 17826.13*u.arcsec
FOV_width = 16924.67*u.arcsec
FOV = FOV_width**2
FOV_g = FOV/4

region = RectangleSkyRegion(SkyCoord(0 * u.deg, 0 * u.deg), FOV_length, FOV_width, angle=90*u.deg)


# Observing constraints
min_declination = -70*u.deg  # Minimum declination limit for fields