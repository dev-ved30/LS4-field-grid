import argparse

import pandas as pd

from visualizations import plot_coverage_map
from pathlib import Path
from generate_obs_plan import compute_union_area

def argument_parser():

    parser = argparse.ArgumentParser(description="Generate an observing plan for LS4 based on the field grid and current visibility.")
    parser.add_argument('obs_plan', type=Path, help="Path to the observing plan CSV file to visualize.")
    return parser.parse_args()

def main():

    args = argument_parser()
    df = pd.read_csv(args.obs_plan)
    mask = compute_union_area(df)
    plot_coverage_map(mask)
    
if __name__ == "__main__":
    main()