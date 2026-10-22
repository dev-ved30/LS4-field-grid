"""End-to-end smoke test on a short window and a throwaway database.

Run with `pytest` or `python tests/test_smoke.py`.
"""

import sys
import tempfile
from pathlib import Path

import pandas as pd
from astropy import units as u

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ls4.cli import main  # noqa: E402
from ls4.config import EXTRAGALACTIC, SOLE  # noqa: E402
from ls4.grid import fresh_db, load_grid_csv, load_history, program_ids  # noqa: E402
from ls4.metrics import night_metrics  # noqa: E402
from ls4.scheduler import dither_phase, night_window, pair_coords, plan_night  # noqa: E402


def test_short_night_schedules_valid_pairs():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        db = fresh_db(tmp / "grid.db")
        start, _ = night_window(pd.Timestamp("2026-07-17").date())
        plan = plan_night(start, start + 2 * u.hour, tmp / "plan.csv", db_path=db)

        result = night_metrics(plan)
        assert result["images"] > 0
        assert result["unpaired_images"] == 0
        assert result["short_pairs"] == 0
        assert result["images"] == 2 * result["complete_pairs"]

        # Only fields with a complete pair get their history advanced.
        history = load_history(db)
        assert (history["last_scheduled_mjd"] > 0).sum() == result["complete_pairs"]

        main(["export", str(tmp / "plan.csv"), "--output", str(tmp / "plan.obsplan")])
        assert len((tmp / "plan.obsplan").read_text().splitlines()) == result["images"]


def test_night_window_is_one_night():
    start, end = night_window(pd.Timestamp("2026-07-17").date())
    assert start.iso.startswith("2026-07-17")
    assert 8 < (end - start).to_value(u.hour) < 12


def test_sole_footprint():
    # l is stored in [0, 360); the plane spans -72 <= l <= 36 across l = 0.
    assert list(program_ids([300, 350, 30, 330, 40, 280], [0, -7, 3.5, 5, 0, 0])) == \
        [SOLE, SOLE, SOLE, EXTRAGALACTIC, EXTRAGALACTIC, EXTRAGALACTIC]
    grid = load_grid_csv()
    assert (grid["program_id"] == program_ids(grid["l"], grid["b"])).all()
    assert (grid["program_id"] == SOLE).sum() == 55  # ~1280 square degrees


def test_dither_alternates_by_night():
    nights = [night_window(pd.Timestamp(f"2026-07-{day}").date())[0] for day in (17, 18, 19)]
    assert [dither_phase(start) for start in nights] in ([0, 1, 0], [1, 0, 1])

    grid = load_grid_csv()
    field = grid[grid["base_field"] == "56_14"]
    ring = grid.loc[grid["dec"] == field["dec"].iloc[0], "ra"].sort_values()
    ra_step = ring.diff().median()
    for phase in (0, 1):
        initial, dither = pair_coords(field["ra"], field["dec"], phase)
        assert abs(initial.ra.deg[0] - field["ra"].iloc[0]) < 1e-9
        assert abs(initial.dec.deg[0] - field["dec"].iloc[0] - phase * 180 / 77) < 1e-6  # half the 360/77-degree ring spacing
        assert abs((initial.ra.deg[0] - dither.ra.deg[0]) - ra_step / 2) < 1e-9
        assert dither.dec.deg[0] == initial.dec.deg[0]


def test_band_quadrants():
    from ls4.bands import band_visits

    exposure = pd.DataFrame({"night": ["n"], "start_time_mjd": [61000.0], "ra": [30.0], "dec": [-20.0]})
    counts = band_visits(exposure, nside=256)["band"].value_counts()
    # i covers two quadrants, g and z one each.
    assert abs(counts["i"] / counts["g"] - 2) < 0.15 and abs(counts["z"] / counts["g"] - 1) < 0.1


if __name__ == "__main__":
    test_band_quadrants()
    test_dither_alternates_by_night()
    test_sole_footprint()
    test_night_window_is_one_night()
    test_short_night_schedules_valid_pairs()
    print("ok")
