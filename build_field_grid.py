import sqlite3
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

        coord = SkyCoord(row_ra * u.deg, row_dec * u.deg).galactic
        l, b = coord.l.to(u.deg).value, coord.b.to(u.deg).value

        if (-72 < l < 36) and np.abs(b) < 4:
            program_id = 0 # Galactic Plane
        elif np.abs(l) < 20 and np.abs(b) < 8:
            program_id = 0 # Galactic Bulge
        else:
            program_id = 1 # Extragalactic

        # TODO still need to assign high cadence extragalactic fields, and assign program IDs to them. 

        if dec_idx % 2 == 0 and ra_idx % 2 == 0:
            keep = True
        else:
            keep = False

        df.loc[row.Index, "Field Name"] = f"{ra_idx}_{dec_idx}"
        df.loc[row.Index, "l"] = l
        df.loc[row.Index, "b"] = b
        df.loc[row.Index, "program_id"] = program_id
        df.loc[row.Index, "keep"] = keep
    
    df['program_id'] = df['program_id'].astype(int)

    df = df[df["keep"] == True].reset_index().drop(columns=['keep', 'index'])

    df['last_scheduled_mjd'] = 0.0

    # Save the field grid to a CSV file for later use in generating the observing plan
    df.to_csv(LS4_field_grid_path, index=False)

    # 1. Create your database connection (creates the file if it doesn't exist)
    conn = sqlite3.connect(LS4_field_grid_db_path)

    # 2. Write the dataframe to a table named 'grid' in the database. 
    df.to_sql("grid", conn, if_exists="fail", index=False)

    # 3. Always close the connection when finished
    conn.close()

    # the script below is just for testing the footprints and visualizations. It can be removed later.
    print(df)
    plot_field_grid(df)

if __name__ == "__main__":
    main()