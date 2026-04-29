import argparse

import pandas as pd

from visualizations import plot_obs_plan
from pathlib import Path

def argument_parser():

    parser = argparse.ArgumentParser(description="Generate an observing plan for LS4 based on the field grid and current visibility.")
    parser.add_argument('obs_plan', type=Path, help="Path to the observing plan CSV file to visualize.")
    return parser.parse_args()

def main():

    args = argument_parser()
    df = pd.read_csv(args.obs_plan)
    plot_obs_plan(df)
    
if __name__ == "__main__":
    main()