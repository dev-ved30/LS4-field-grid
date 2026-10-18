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
from ls4.grid import fresh_db, load_history  # noqa: E402
from ls4.metrics import night_metrics  # noqa: E402
from ls4.scheduler import night_window, plan_night  # noqa: E402


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


if __name__ == "__main__":
    test_night_window_is_one_night()
    test_short_night_schedules_valid_pairs()
    print("ok")
