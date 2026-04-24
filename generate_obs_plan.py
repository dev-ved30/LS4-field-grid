import astroplan
import astropy
import pandas as pd

from astroplan import Observer
from astroplan.constraints import AirmassConstraint, TimeConstraint, AtNightConstraint
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan.scheduling import Transitioner, SequentialScheduler, Schedule
from astroplan import ObservingBlock


from visualizations import plot_visible_fields, plot_obs_plan
from constants import *
from build_field_grid import region


# Universal constraints for all targets
global_constraints = [AirmassConstraint(max=2),
                      AtNightConstraint.twilight_civil()]

# Get the current time
current_time = Time.now() # Fix this to a specific time for testing, e.g. Time("2024-06-01 00:00:00")
night_start = LS4.twilight_evening_astronomical(current_time, which='next')
night_end = LS4.twilight_morning_astronomical(current_time, which='next')

def skip_alternate_fields_in_dec(fields, config):

    unskipped_fields = []
    for f in fields:
        field_name = f.name
        ra_idx, dec_idx = field_name.split("_")
        ra_idx = int(ra_idx)
        dec_idx = int(dec_idx)
        if config == 1:
            if dec_idx % 2 == 0 and ra_idx % 2 == 0:
                unskipped_fields.append(f)
        elif config == 2:
            if dec_idx % 2 == 1 and ra_idx % 2 == 1:
                unskipped_fields.append(f)
    return unskipped_fields

def skip_polar_fields(fields):

    unskipped_fields = []
    for f in fields:
        
        if f.coord.dec > min_declination:
            unskipped_fields.append(f)
    
    return unskipped_fields


def get_visible_fields():
    """Get the list of fields that are currently visible from La Silla Observatory."""

    # Read in the field grid
    field_grid = pd.read_csv("LS4_field_grid.csv")


    time_range = [night_start, night_end]
    print(time_range)

    # Get a list of visible fields for tonight
    all_fields = []
    for _, row in field_grid.iterrows():

        field_coords = SkyCoord(ra=row["ra_deg"]*u.deg, dec=row["dec_deg"]*u.deg)
        field_target = FixedTarget(name=row["Field Name"], coord=field_coords)
        all_fields.append(field_target)

    alternate_fields_g1 = skip_alternate_fields_in_dec(all_fields, config=1) 
    alternate_fields_g2 = skip_alternate_fields_in_dec(all_fields, config=2)

    non_polar_fields_g1 = skip_polar_fields(alternate_fields_g1)
    non_polar_fields_g2 = skip_polar_fields(alternate_fields_g2)


    # Check which fields are visible tonight
    is_visible_g1 = astroplan.is_observable(global_constraints, 
                                         LS4, 
                                         non_polar_fields_g1,
                                         time_range=time_range)
    
    is_visible_g2 = astroplan.is_observable(global_constraints, 
                                         LS4, 
                                         non_polar_fields_g2,
                                         time_range=time_range)

     # Print visible fields
    visible_fields_g1 = [non_polar_fields_g1[i] for i in range(len(non_polar_fields_g1)) if is_visible_g1[i]]
    visible_fields_g2 = [non_polar_fields_g2[i] for i in range(len(non_polar_fields_g2)) if is_visible_g2[i]]

    return visible_fields_g1, visible_fields_g2

def get_obs_blocks():

    visible_fields_g1, visible_fields_g2 = get_visible_fields()
    print(len(visible_fields_g1), "fields are visible tonight in config 1.")
    print(len(visible_fields_g2), "fields are visible tonight in config 2.")

    num_exposures = 1
    blocks = []
    times = []

    # Chop up the night into 30 minute black and assign 18 field to each block, alternating between the two configs to maximize the number of fields observed while minimizing slews. This is a simple heuristic that can be improved with more sophisticated scheduling algorithms.
    block_duration = 30*u.minute
    fields_per_block = block_duration / (exp + read_out)
    current_time = night_start

    block_number = 0
    while current_time < night_end:

        block = []
        block_durations = [current_time, current_time + block_duration]

        print(f"Scheduling block {block_number} from {block_durations[0]} to {block_durations[1]} with config {1 if block_number % 2 == 0 else 2}...")

        # find which fields are visible during this block
        time_constraint = TimeConstraint(block_durations[0], block_durations[1])
        block_constraints = global_constraints + [time_constraint]\
        
        if block_number % 2 == 0:
            visible_fields_in_block = [f for f in visible_fields_g1 if astroplan.is_observable(block_constraints, LS4, [f], time_range=block_durations)[0]]
        else:
            visible_fields_in_block = [f for f in visible_fields_g2 if astroplan.is_observable(block_constraints, LS4, [f], time_range=block_durations)[0]]

        for target in visible_fields_in_block[:int(fields_per_block)]:
            b = ObservingBlock.from_exposures(target, 
                                            priority=1, 
                                            time_per_exposure=exp, 
                                            number_exposures=num_exposures, 
                                            readout_time=read_out,
                                            constraints=block_constraints)
            
            # add target to this block
            block.append(b)

            # remove target from the list of visible fields to avoid scheduling it again
            if block_number % 2 == 0:
                visible_fields_g1.remove(target)
            else:
                visible_fields_g2.remove(target)    

        current_time += block_duration
        block_number += 1
        blocks.append(block)
        times.append(block_durations)

        if len(visible_fields_g1) == 0 and len(visible_fields_g2) == 0:
            break

    print(f'{len(visible_fields_g1) + len(visible_fields_g2)} fields were not scheduled due to time constraints.')
    
    return blocks, times
    
    



    # # Create ObservingBlocks for each filter and target with our time
    # # constraint, and durations determined by the exposures needed
    # for target in visible_fields:
    #     # We want each filter to have separate priority (so that target
    #     # and reference are both scheduled)
    #     b = ObservingBlock.from_exposures(target, 
    #                                       priority=1, 
    #                                       time_per_exposure=exp, 
    #                                       number_exposures=num_exposures, 
    #                                       readout_time=read_out,
    #                                       constraints=global_constraints)

    #     if len(blocks) > 100:
    #         break

    #     blocks.append(b)

    # return blocks

def get_obs_plan():

    blocks, times = get_obs_blocks()
    transitioner = Transitioner(slew_rate)
    combined_obs_plan = []

    for b, t in zip(blocks, times):

        sequential_schedule = Schedule(t[0], t[1])
        seq_scheduler = SequentialScheduler(constraints = global_constraints,
                                            observer = LS4,
                                            transitioner = transitioner)
        seq_scheduler(b, sequential_schedule)   
        combined_obs_plan.append(sequential_schedule.to_table(show_unused=False).to_pandas())

    
    # concatenate the astropy tables
    combined_obs_plan = pd.concat(combined_obs_plan, ignore_index=True)
    print(combined_obs_plan)
    combined_obs_plan.to_csv("obs_plan.csv", index=False)


    # # Get the current time
    # current_time = Time.now()
    # night_start = LS4.twilight_evening_astronomical(current_time, which='next')
    # night_end = LS4.twilight_morning_astronomical(current_time, which='next')
    
    # # Initialize a Schedule object, to contain the new schedule
    # sequential_schedule = Schedule(night_start, night_end)

    # transitioner = Transitioner(slew_rate)

    # # Initialize the sequential scheduler with the constraints and transitioner
    # seq_scheduler = SequentialScheduler(constraints = global_constraints,
    #                                     observer = LS4,
    #                                     transitioner = transitioner)
    


    # # Call the schedule with the observing blocks and schedule to schedule the blocks
    # seq_scheduler(blocks, sequential_schedule)
    # print(sequential_schedule.to_table(show_unused=True))

    # # plot the schedule over the field grid
    plot_obs_plan(combined_obs_plan, region)
    



if __name__ == "__main__":
    import time
    start_time = time.time()
    get_obs_plan()
    end_time = time.time()
    print(f"Scheduling took {end_time - start_time:.2f} seconds.")