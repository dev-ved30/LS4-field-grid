import astroplan
import pandas as pd

from astroplan import Observer
from astroplan.constraints import AirmassConstraint, TimeConstraint, AtNightConstraint
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan.scheduling import Transitioner, SequentialScheduler, Schedule
from astroplan import ObservingBlock

from visualizations import plot_visible_fields

# Define the preliminaries
LS4 = Observer.at_site("La Silla Observatory")
slew_rate = .8*u.deg/u.second

# Universal constraints for all targets
global_constraints = [AirmassConstraint(max=2),
                      AtNightConstraint.twilight_civil()]


def skip_alternate_fields_in_dec(fields):

    unskipped_fields = []
    for f in fields:
        field_name = f.name
        ra_idx, dec_idx = field_name.split("_")
        ra_idx = int(ra_idx)
        dec_idx = int(dec_idx)
        if dec_idx % 2 == 0:
            unskipped_fields.append(f)
    return unskipped_fields

def skip_polar_fields(fields):

    unskipped_fields = []
    for f in fields:
        
        if f.coord.dec > -70*u.deg:
            unskipped_fields.append(f)
    
    return unskipped_fields

def get_visible_fields():
    """Get the list of fields that are currently visible from La Silla Observatory."""

    # Read in the field grid
    field_grid = pd.read_csv("LS4_field_grid.csv")

    # Get the current time
    current_time = Time.now()
    night_start = LS4.twilight_evening_astronomical(current_time, which='next')
    night_end = LS4.twilight_morning_astronomical(current_time, which='next')
    time_range = [night_start, night_end]
    print(time_range)

    # Get a list of visible fields for tonight
    all_fields = []
    for _, row in field_grid.iterrows():

        field_coords = SkyCoord(ra=row["ra_deg"]*u.deg, dec=row["dec_deg"]*u.deg)
        field_target = FixedTarget(name=row["Field Name"], coord=field_coords)
        all_fields.append(field_target)

    alternate_fields = skip_alternate_fields_in_dec(all_fields) 
    non_polar_fields = skip_polar_fields(alternate_fields)


    # Check which fields are visible tonight
    is_visible = astroplan.is_observable(global_constraints, 
                                         LS4, 
                                         non_polar_fields,
                                         time_range=time_range)
    
     # Print visible fields
    visible_fields = [non_polar_fields[i] for i in range(len(non_polar_fields)) if is_visible[i]]
    
    plot_visible_fields(non_polar_fields, is_visible)
    return visible_fields

def get_obs_blocks():

    visible_fields = get_visible_fields()
    print(len(visible_fields), "fields are visible tonight.")

    read_out = 20 * u.second
    exp = 15*u.second
    n = 1
    blocks = []


    # Create ObservingBlocks for each filter and target with our time
    # constraint, and durations determined by the exposures needed
    for target in visible_fields:
        # We want each filter to have separate priority (so that target
        # and reference are both scheduled)
        b = ObservingBlock.from_exposures(target, 
                                          priority=1, 
                                          time_per_exposure=exp, 
                                          number_exposures=1, 
                                          readout_time=read_out,
                                          constraints = global_constraints)
        if len(blocks) >= 200:  # Limit to 10 blocks for testing
            break
        blocks.append(b)

    return blocks

def get_obs_plan():

    blocks = get_obs_blocks()

    # Get the current time
    current_time = Time.now()
    night_start = LS4.twilight_evening_astronomical(current_time, which='previous')
    night_end = LS4.twilight_morning_astronomical(current_time, which='next')
    
    # Initialize a Schedule object, to contain the new schedule
    sequential_schedule = Schedule(night_start, night_end)

    transitioner = Transitioner(slew_rate)

    # Initialize the sequential scheduler with the constraints and transitioner
    seq_scheduler = SequentialScheduler(constraints = global_constraints,
                                        observer = LS4,
                                        transitioner = transitioner)
    


    # Call the schedule with the observing blocks and schedule to schedule the blocks
    seq_scheduler(blocks, sequential_schedule)
    print(sequential_schedule.to_table(show_unused=True))


if __name__ == "__main__":
    import time
    start_time = time.time()
    get_obs_plan()
    end_time = time.time()
    print(f"Scheduling took {end_time - start_time:.2f} seconds.")