import time
import astroplan
import argparse

from networkx import efficiency
import pandas as pd

from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan.scheduling import Transitioner, SequentialScheduler, Schedule
from astroplan import ObservingBlock

from visualizations import plot_obs_plan
from constants import *
from tqdm import tqdm

# Get the current time
current_time = Time.now() # Fix this to a specific time for testing, e.g. Time("2024-06-01 00:00:00")

def argument_parser():

    parser = argparse.ArgumentParser(description="Generate an observing plan for LS4 based on the field grid and current visibility.")
    parser.add_argument('--mjd', required=False, type=float, default=current_time, help="MJD for which to generate the observing plan. If not provided, the current time will be used.")
    parser.add_argument("--output", required=False, type=str, default=None, help="Path to save the generated observing plan CSV file. Default will just save to plans/yyyymmdd.csv")
    return parser.parse_args()

def skip_alternate_fields_in_dec(fields):

    unskipped_fields = []
    for f in fields:
        field_name = f.name
        ra_idx, dec_idx = field_name.split("_")
        ra_idx = int(ra_idx)
        dec_idx = int(dec_idx)
        if dec_idx % 2 == 0 and ra_idx % 2 == 0:
            unskipped_fields.append(f)

    return unskipped_fields

def skip_polar_fields(fields):

    unskipped_fields = []
    for f in fields:
        
        if f.coord.dec > min_declination:
            unskipped_fields.append(f)
    
    return unskipped_fields

def compute_theoretical_max_images_per_night(night_start, night_end):

    night_duration = (night_end - night_start).to_value(u.second) * u.second
    total_time_per_field = exp + read_out
    theoretical_max_fields = (night_duration / total_time_per_field).to_value()

    return theoretical_max_fields

def compute_time_efficiency(obs_plan, night_start, night_end):

    total_night_duration = (night_end - night_start).to_value(u.minute) * u.minute

    # unused time
    unused_time = obs_plan[obs_plan['target'] == 'Unused Time']
    total_unused_time = unused_time['duration (minutes)'].sum() * u.minute

    efficiency = (total_night_duration - total_unused_time) / total_night_duration

    print(f"Total night duration: {total_night_duration.to_value(u.hour):.2f} hours")
    print(f"Total unused time: {total_unused_time.to_value(u.hour):.2f} hours")
    print(f"Time efficiency: {efficiency:.2%}")

    return efficiency

def get_visible_fields(night_start, night_end):
    """Get the list of fields that are currently visible from La Silla Observatory."""

    # Read in the field grid
    field_grid = pd.read_csv(LS4_field_grid_path)

    time_range = [night_start, night_end]

    # Get a list of visible fields for tonight
    all_fields = []
    for _, row in field_grid.iterrows():

        field_coords = SkyCoord(ra=row["ra_deg"]*u.deg, dec=row["dec_deg"]*u.deg)
        field_target = FixedTarget(name=row["Field Name"], coord=field_coords)
        all_fields.append(field_target)

    # only get primary pointing on the grid
    alternate_fields = skip_alternate_fields_in_dec(all_fields) 

    # skip the fields too close to the pole. These will probably be skipped due to high air mass anyway.
    non_polar_fields = skip_polar_fields(alternate_fields)

    # Check which fields are visible tonight
    is_visible = astroplan.is_observable(global_constraints, 
                                         LS4, 
                                         non_polar_fields,
                                         time_range=time_range)

     # Print visible fields
    visible_fields = [non_polar_fields[i] for i in range(len(non_polar_fields)) if is_visible[i]]

    return visible_fields

def get_obs_blocks(night_start, night_end):

    visible_fields = get_visible_fields(night_start, night_end)
    print("==============================================================")
    print(len(visible_fields), "fields are visible tonight.")

    num_exposures = 1
    blocks = []
    times = []

    # Chop up the night into 30 minute black and assign 18 field to each block, alternating between the two configs to maximize the number of fields observed while minimizing slews. This is a simple heuristic that can be improved with more sophisticated scheduling algorithms.
    current_time = night_start

    block_number = 0
    while current_time < night_end:

        # TODO: fix this since it currently wastes the last bit if the night
        if night_end - current_time < 2*block_duration:
            break

        block = []
        block_dither = []

        # Check the observability of the fields for the duration of the block and the revisit with dither.
        block_durations = [current_time, current_time + block_duration]
        revisit_duration = [current_time + block_duration, current_time + 2*block_duration]

        # Only schedule fields that are observable for the entire block duration and the revisit duration.
        fields_visible_for_revisit = []
        for f in visible_fields:
            if astroplan.is_observable(global_constraints, LS4, [f], time_range=block_durations)[0] and astroplan.is_observable(global_constraints, LS4, [f], time_range=revisit_duration)[0]:
                fields_visible_for_revisit.append(f)    
        
        # sort by angular distance from the first field to minimize slews
        if len(blocks) > 0:
            starting_point = blocks[-1][0].target # for subsequent blocks, sort by distance from the last scheduled field to minimize slews
        else:
            fields_visible_for_revisit.sort(key=lambda f: f.coord.ra)
            starting_point = fields_visible_for_revisit[0] # for the first block, just sort by RA
        
        fields_visible_for_revisit.sort(key=lambda f: f.coord.separation(starting_point.coord))

        for target in fields_visible_for_revisit[:fields_per_block]:
            b = ObservingBlock.from_exposures(target, 
                                            priority=1, 
                                            time_per_exposure=exp, 
                                            number_exposures=num_exposures, 
                                            readout_time=read_out,
                                            constraints=global_constraints,
                                            configuration={"dither": False}) # add a configuration parameter to indicate whether this block is a dither or not. This can be used later for visualization and analysis.
            
            # add target to this block
            block.append(b)

            # add a dither to the block by offsetting the target coordinates by half a field in RA and Dec. This is a simple dither pattern that can be improved with more sophisticated patterns.
            dithered_target = FixedTarget(name=target.name + "_dither",
                                        coord=SkyCoord(ra=target.coord.ra + half_field_offset_ra,
                                                        dec=target.coord.dec + half_field_offset_dec))
            b_dither = ObservingBlock.from_exposures(dithered_target,
                                            priority=1,
                                            time_per_exposure=exp,
                                            number_exposures=num_exposures,
                                            readout_time=read_out,
                                            constraints=global_constraints,
                                            configuration={"dither": True})
            block_dither.append(b_dither)

            # remove target from the list of visible fields to avoid scheduling it again
            visible_fields.remove(target)

        print(f"{len(block) + len(block_dither)} images slated for block {block_number} and {block_number + 1}. {len(visible_fields)} fields remaining to schedule.")

        current_time += 2 * block_duration
        block_number += 2

        # Add blocks for both the original pointing and the revisit with the dither.
        blocks.append(block)
        blocks.append(block_dither)

        times.append(block_durations)
        times.append(revisit_duration)

        if len(visible_fields) == 0:
            break
    
    return blocks, times
    
    
def get_obs_plan(night_start, night_end, output_path):

    print("Generating observing plan for the night of", night_start.to_datetime().strftime("%Y-%m-%d"), "to", night_end.to_datetime().strftime("%Y-%m-%d"))
    blocks, times = get_obs_blocks(night_start, night_end)
    transitioner = Transitioner(slew_rate)
    combined_obs_plan = []

    for i, (b, t) in tqdm(enumerate(zip(blocks, times)), total=len(blocks), desc="Scheduling blocks"):

        sequential_schedule = Schedule(t[0], t[1])
        seq_scheduler = SequentialScheduler(constraints = global_constraints,
                                            observer = LS4,
                                            transitioner = transitioner)
        seq_scheduler(b, sequential_schedule)   
        combined_obs_plan.append(sequential_schedule.to_table(show_unused=True).to_pandas())
        combined_obs_plan[-1]["block_number"] = i


    combined_obs_plan = pd.concat(combined_obs_plan, ignore_index=True)
    #print(combined_obs_plan)

    # remove the deg from ra and dec columns while preserving blank rows
    combined_obs_plan['ra'] = pd.to_numeric(
        combined_obs_plan['ra'].astype(str).str.replace(' deg', '', regex=False),
        errors='coerce',
    )
    combined_obs_plan['dec'] = pd.to_numeric(
        combined_obs_plan['dec'].astype(str).str.replace(' deg', '', regex=False),
        errors='coerce',
    )

    combined_obs_plan['ra_hr'] = combined_obs_plan['ra'] * u.deg.to(u.hourangle)

    print("Saving observing plan to", output_path)
    combined_obs_plan.to_csv(output_path, index=False)

    print("Done!\n==============================================================")
    images_scheduled = sum(len(b) for b in blocks)
    max_images_possible = compute_theoretical_max_images_per_night(night_start, night_end)
    print(f"{images_scheduled} images scheduled out of {max_images_possible:.0f} possible images for the night based on exposure time and readout time.")
    print("Imaging efficiency: {:.2f}%".format(100 * images_scheduled / max_images_possible))

    compute_time_efficiency(combined_obs_plan, night_start, night_end)


if __name__ == "__main__":

    start_time = time.time()

    args = argument_parser()
    mjd = args.mjd
    night_start = LS4.twilight_evening_astronomical(Time(mjd, format='mjd'), which='next')
    night_end = LS4.twilight_morning_astronomical(Time(mjd, format='mjd'), which='next')
    
    if args.output is not None:
        output_path = args.output
    else:
        output_path = f"plans/{Time(mjd, format='mjd').to_datetime().strftime('%Y%m%d')}.csv"

    get_obs_plan(night_start, night_end, output_path)

    end_time = time.time()
    print(f"Scheduling took {end_time - start_time:.2f} seconds.")