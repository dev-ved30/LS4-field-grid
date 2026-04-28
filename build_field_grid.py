import numpy as np
import pandas as pd
import numpy as np

from astropy import units as u
from regions import RectangleSkyRegion
from astropy.coordinates import SkyCoord
from astropy import units as u
from m4opt.fov import footprint
from m4opt import skygrid

from visualizations import plot_field_grid
from constants import *


vertices = skygrid.sinusoidal(FOV_g)
target_coord = vertices

region = RectangleSkyRegion(SkyCoord(0 * u.deg, 0 * u.deg), FOV_length, FOV_width, angle=90*u.deg)

def main():
    footprints = footprint(region, target_coord)

    ra_deg = vertices.ra.to(u.deg).value
    dec_deg = vertices.dec.to(u.deg).value
    df = pd.DataFrame({"ra_deg": ra_deg, "dec_deg": dec_deg})

    sorted_dec = np.sort(np.unique(dec_deg))

    for row in df.itertuples():

        row_ra, row_dec = row.ra_deg, row.dec_deg
        dec_idx = np.where(sorted_dec == row_dec)[0][0]

        sorted_ras_at_dec = np.sort(df[df["dec_deg"] == row_dec]["ra_deg"])
        ra_idx = np.where(sorted_ras_at_dec == row_ra)[0][0]

        df.loc[row.Index, "Field Name"] = f"{ra_idx}_{dec_idx}"

    df.to_csv("LS4_field_grid.csv", index=False)
    print(df)


    plot_field_grid(vertices, [footprints[1000] ,footprints[1001], footprints[1004]])

if __name__ == "__main__":
    main()