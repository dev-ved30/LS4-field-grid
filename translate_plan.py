#!/usr/bin/env python3

##############################
# -*- coding: utf-8 -*-
#
# @Author: Kenneth Lin
# @Date: 2026-04-30
# @Filename: build_plan.py
#
#
###############################

import sys
import pandas as pd
from datetime import datetime, timezone, timedelta
import argparse

# ------------------------------
# Define parameters here
# ------------------------------
exposure_time = 60.0       # Column 4
time_between = 100.0       # this is now ignored by David's scheduler Column 5
num_exposures = 1          # this is now ignored by David's scheduler Column 6
priority_code = 0          # this is now ignored by David's scheduler Column 7
output_file = "20260425.obsplan"
# ------------------------------

if len(sys.argv) < 2:
    print("Usage: python rawplan_csv.py <input_file.csv>")
    sys.exit(1)

input_csv = sys.argv[1]

def clean_angle(angle):
    if isinstance(angle, str):
        angle = angle.strip()
        if angle == "" or angle.lower() == "nan":
            return None
        if "deg" in angle:
            angle = angle.replace("deg", "").strip()
    try:
        return float(angle)
    except (TypeError, ValueError):
        return None

def build(input_csv, output_file):
    # Read CSV
    df = pd.read_csv(input_csv)

    # Clean RA/Dec values
    df["ra_deg"] = df["ra"].apply(clean_angle)
    df["dec_deg"] = df["dec"].apply(clean_angle)

    # Drop rows missing either RA or Dec
    df = df.dropna(subset=["ra_deg", "dec_deg"])

    # Convert RA decimal degrees → hours
    df["ra_hours"] = df["ra_deg"] / 15.0

    # Write output file
    with open(output_file, "w") as f:
        for idx, row in df.iterrows():
            ra  = row["ra_hours"]     # now in hours
            dec = row["dec_deg"]      # still degrees
            comment = f"# {row['target']}"

            line = (
                f"{ra:.4f} {dec:.4f} "
                f"Y  {exposure_time:.1f} {time_between:.1f} "
                f"{num_exposures} {priority_code} {comment}"
            )

            f.write(line + "\n")

    print("Wrote:", output_file)

def generate_parser():
    """
    Argument parser to run as CLI script
    """
    parser = argparse.ArgumentParser(description="Build the obsplan from CSV file from Ved")
    parser.add_argument('input_csv', metavar='CSV_FILE', type=str, help="Path to input CSV file")
    parser.add_argument('-d', '--date', metavar='DATE', type=int, nargs='+', help="Insert custom date")
    return parser

def main():
    parser = generate_parser()
    args = parser.parse_args()
    
    if args.date:
        current_utc = str(args.date[0])
    else:
        current_utc = (datetime.now(timezone.utc) + timedelta(days=1)).strftime("%Y%m%d")
    
    output_file = current_utc + '.obsplan'
    build(args.input_csv, output_file)
    
if __name__ == '__main__':
    main()

