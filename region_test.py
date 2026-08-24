from astropy.coordinates import SkyCoord
from astropy import units as u
import ligo.skymap.plot
from matplotlib import pyplot as plt
from m4opt.fov import footprint
import numpy as np
from regions import Regions

# Read shapes from a DS9 region file
roman_wfi = Regions.read('LS4.reg.reg', format='ds9')

# Create a grid of target coordinates
ra, dec = np.meshgrid(
    np.linspace(-1.5, 0.5, 4) * u.deg,
    np.linspace(-1, 1, 3) * u.deg
)
target_coords = SkyCoord(ra, dec).ravel()

# Transform shapes to target coordinates
footprints = footprint(roman_wfi, target_coords, -30 * u.deg)

# Plot footprints
ax = plt.axes(
    projection='astro zoom',
    center=SkyCoord(0 * u.deg, 0 * u.deg),
    radius=2 * u.deg)
for regions in footprints:
    for region in regions:
        ax.add_patch(region.to_pixel(ax.wcs).as_artist())
ax.grid()
plt.show()