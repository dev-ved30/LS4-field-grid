"""Observing-plan website: `ls4 site` builds a consistent static site from saved plans."""

import json
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ls4 import automation  # noqa: E402
from ls4.config import PLANS_DIR  # noqa: E402
from ls4.website import build_site  # noqa: E402


@pytest.fixture
def plans_dir(tmp_path):
    plans = tmp_path / "plans"
    forecast = plans / "forecast_20260718"
    forecast.mkdir(parents=True)
    for name in ("20260717.csv", "20260718.csv"):
        shutil.copyfile(PLANS_DIR / name, plans / name)
    shutil.copyfile(PLANS_DIR / "20260719.csv", forecast / "20260719.csv")
    (forecast / "README.md").write_text("# Forecast\n\nOne `test` night.\n")
    (forecast / "summary.csv").write_text("not a plan\n")
    return plans


def test_site_is_consistent(tmp_path, plans_dir):
    out = build_site(tmp_path / "site", plans_dir=plans_dir, nside=32)
    index = json.loads((out / "data" / "index.json").read_text())
    assert [c["key"] for c in index["collections"]] == ["production", "forecast_20260718"]
    assert index["collections"][1]["description"] == "One test night."
    assert (out / "index.html").exists() and (out / "app.js").exists()

    production = index["collections"][0]
    assert [n["night"] for n in production["nights"]] == ["20260717", "20260718"]
    survey = json.loads((out / "data" / "production" / "survey.json").read_text())
    assert [n["label"] for n in survey["nights"]] == ["20260717", "20260718"]
    assert len(survey["cumulative"]) == 2 and survey["pixels"]

    night = json.loads((out / "data" / "production" / "20260717.json").read_text())
    assert len(night["pointings"]) == night["metrics"]["images"] == production["nights"][0]["images"]
    assert sum(p[4] for p in night["pointings"]) == night["metrics"]["complete_pairs"]
    assert night["window"][0] - 1e-5 <= night["pointings"][0][0] < night["pointings"][-1][1] <= night["window"][1] + 1e-5
    obsplan = out / night["files"]["obsplan"]
    assert obsplan.name == "20260718.obsplan"
    assert len(obsplan.read_text().splitlines()) == night["metrics"]["images"]
    assert (out / night["files"]["csv"]).read_bytes() == (plans_dir / "20260717.csv").read_bytes()


def test_password_protected_site_has_no_plaintext(tmp_path, plans_dir):
    import base64
    import gzip

    from cryptography.exceptions import InvalidTag

    from ls4.website import Sealer

    out = build_site(tmp_path / "site", plans_dir=plans_dir, nside=32, password="correct horse")
    assert not [p for p in (out / "data").rglob("*") if p.is_file() and p.suffix != ".enc"]
    assert not [p for p in (out / "files").rglob("*") if p.is_file() and p.suffix != ".enc"]
    auth = json.loads((out / "auth.json").read_text())

    def unseal(raw, password):
        aead = Sealer(password, "dev-ved30/LS4-field-grid").aead
        assert base64.b64decode(auth["salt"]) == Sealer(password, "dev-ved30/LS4-field-grid").salt
        return gzip.decompress(aead.decrypt(raw[:12], raw[12:], None))

    assert unseal(base64.b64decode(auth["check"]), "correct horse") == b"LS4"
    index = json.loads(unseal((out / "data" / "index.json.enc").read_bytes(), "correct horse"))
    assert [c["key"] for c in index["collections"]] == ["production", "forecast_20260718"]
    with pytest.raises(InvalidTag):
        unseal(base64.b64decode(auth["check"]), "wrong")
    # A rebuild reuses the key (fixed salt), so remembered sign-ins survive nightly updates.
    again = build_site(tmp_path / "site", plans_dir=plans_dir, nside=32, password="correct horse")
    assert json.loads((again / "auth.json").read_text())["salt"] == auth["salt"]


def test_site_refuses_to_replace_other_directories(tmp_path, plans_dir):
    other = tmp_path / "notes"
    other.mkdir()
    (other / "keep.txt").write_text("important")
    with pytest.raises(SystemExit):
        build_site(other, plans_dir=plans_dir)
    assert (other / "keep.txt").exists()


def test_nightly_releases_keep_the_newest_revision_of_each_night():
    releases = [
        {"tag_name": "ls4-20261021-100-1", "draft": False, "published_at": "2026-10-21T20:00:00Z"},
        {"tag_name": "ls4-20261021-90-2", "draft": False, "published_at": "2026-10-21T22:00:00Z"},
        {"tag_name": "ls4-20261022-200-1", "draft": True, "published_at": None},
        {"tag_name": "ls4-20261023-300-1", "draft": False, "published_at": "2026-10-23T20:00:00Z"},
    ]
    client = automation.GitHub("owner/repo", None)
    assert "Authorization" not in client.headers
    client.api = lambda path: releases
    newest = client.nightly_releases()
    assert {stamp: r["tag_name"] for stamp, r in newest.items()} == {
        "20261021": "ls4-20261021-90-2", "20261023": "ls4-20261023-300-1"}
