"""Convert a plan CSV to the telescope's .obsplan text format (after Kenneth Lin's build_plan.py).

Each line is: RA (hours) Dec (deg) Y exposure time_between num_exposures priority # target
"""

from pathlib import Path

import pandas as pd

from ls4.config import exp

exposure_time = exp.to_value("s")  # Column 4
time_between = 100.0  # this is now ignored by David's scheduler Column 5
num_exposures = 1  # this is now ignored by David's scheduler Column 6
priority_code = 0  # this is now ignored by David's scheduler Column 7


def clean_angle(angle):
    if isinstance(angle, str):
        angle = angle.strip()
        if angle == "" or angle.lower() == "nan":
            return None
        angle = angle.replace("deg", "").strip()
    try:
        return float(angle)
    except (TypeError, ValueError):
        return None


def write_obsplan(plan_csv, output_file):
    df = pd.read_csv(plan_csv)
    df["ra_deg"] = df["ra"].apply(clean_angle)
    df["dec_deg"] = df["dec"].apply(clean_angle)
    # Rows without coordinates are gaps (Unused Time), not pointings.
    df = df.dropna(subset=["ra_deg", "dec_deg"])

    lines = [
        f"{row.ra_deg / 15.0:.4f} {row.dec_deg:.4f} "
        f"Y  {exposure_time:.1f} {time_between:.1f} "
        f"{num_exposures} {priority_code} # {row.target}"
        for row in df.itertuples()
    ]
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text("".join(line + "\n" for line in lines))
    print("Wrote:", output_file)
    return output_file
