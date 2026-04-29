import numpy as np
import pandas as pd
import numpy as np

from astropy import units as u
from astropy import units as u
from m4opt.fov import footprint
from m4opt import skygrid

from visualizations import plot_field_grid
from constants import *


def main():

    # tessellate the sky with a grid of points separated by FOV_g, and assign field names based on their position in the grid
    vertices = skygrid.sinusoidal(FOV_g)

    ra_deg = vertices.ra.to(u.deg).value
    dec_deg = vertices.dec.to(u.deg).value
    df = pd.DataFrame({"ra_deg": ra_deg, "dec_deg": dec_deg})

    # Sort by declination first, then by right ascension within each declination band. This ensures that the field names are assigned in a consistent grid pattern.
    sorted_dec = np.sort(np.unique(dec_deg))
    for row in df.itertuples():

        row_ra, row_dec = row.ra_deg, row.dec_deg
        dec_idx = np.where(sorted_dec == row_dec)[0][0]

        sorted_ras_at_dec = np.sort(df[df["dec_deg"] == row_dec]["ra_deg"])
        ra_idx = np.where(sorted_ras_at_dec == row_ra)[0][0]

        df.loc[row.Index, "Field Name"] = f"{ra_idx}_{dec_idx}"

    # Save the field grid to a CSV file for later use in generating the observing plan
    df.to_csv("LS4_field_grid.csv", index=False)

    # the script below is just for testing the footprints and visualizations. It can be removed later.
    print(df)
    plot_field_grid(df)

if __name__ == "__main__":
    main()