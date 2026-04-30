import numpy as np
import healpy as hp

import astropy.units as u
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord

from constants import *

def compute_union_area(pointings_ra,
                       pointings_dec,
                       nside=1024):
    """
    Compute the union area (deg^2) of rectangular FoVs on the sky.

    Parameters
    ----------
    pointings_ra : array-like
        RA of pointings (degrees)
    pointings_dec : array-like
        Dec of pointings (degrees)
    fov_width_deg : float
        Width of FoV (degrees, along RA direction in tangent plane)
    fov_height_deg : float
        Height of FoV (degrees, along Dec direction in tangent plane)
    nside : int
        HEALPix resolution (higher = more accurate, slower)

    Returns
    -------
    area_deg2 : float
        Total union area in square degrees
    """

    npix = hp.nside2npix(nside)
    observed = np.zeros(npix, dtype=int)

    half_w =  FOV_width.to(u.deg).value/ 2
    half_h = FOV_length.to(u.deg).value/ 2

    for ra, dec in zip(pointings_ra, pointings_dec):
        center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

        # Define rectangle corners in tangent plane
        offsets = [
            (-half_w, -half_h),
            ( half_w, -half_h),
            ( half_w,  half_h),
            (-half_w,  half_h),
        ]

        corners = []
        for dx, dy in offsets:
            corner = center.spherical_offsets_by(dx*u.deg, dy*u.deg)
            corners.append(corner)


        for i, c in enumerate(corners):
            theta = 90 - c.dec.deg
            if not (0 <= theta <= 180):
                print("BAD CORNER", i, c.ra.deg, c.dec.deg)

        # Convert to HEALPix vectors
        vecs = np.array([
            hp.ang2vec(np.radians(90 - c.dec.deg),
                       np.radians(c.ra.deg))
            for c in corners
        ])

        # Fill polygon
        pixels = hp.query_polygon(nside, vecs)
        observed[pixels] = True

    # Compute area
    pixel_area_sr = hp.nside2pixarea(nside)
    total_area_sr = observed.sum() * pixel_area_sr
    total_area_deg2 = total_area_sr * (180/np.pi)**2

    plot_coverage_map(observed)

    return total_area_deg2




def plot_coverage_map(observed_mask, title="Observed Sky Area"):
    """
    Plot a HEALPix coverage map.

    Parameters
    ----------
    observed_mask : array (bool)
        HEALPix boolean mask (True = observed)
    """

    # Convert boolean → numeric for plotting
    m = observed_mask.astype(float)

    hp.mollview(
        m,
        title=title,
        unit="Observed",
        cmap="viridis"
    )

    plt.show()

if __name__ == "__main__":

    import pandas as pd

    # Example usage
    df = pd.read_csv('plans/20260429.csv')
    df = df[(df['target'] != 'Unused Time') & (df['target'] != 'TransitionBlock')]

    area = compute_union_area(df['ra'].to_numpy(), df['dec'].to_numpy())
    print(f"Total observed area: {area:.2f} deg^2")