"""Production invariants: history isolation, package integrity, and failure delivery."""

from datetime import date, datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from ls4 import automation, nightly, weather
from ls4.config import GRID_CSV, GRID_DB
from ls4.grid import fresh_db, load_history


def forecast_payload(start, end):
    first = start.replace(minute=0, second=0, microsecond=0)
    last = end.replace(minute=0, second=0, microsecond=0) + timedelta(hours=1)
    stamps = [first + timedelta(hours=i) for i in range(int((last-first).total_seconds()/3600)+1)]
    return {"utc_offset_seconds": 0, "latitude": -29.25, "longitude": -70.73, "elevation": 2400,
            "hourly_units": weather.VARIABLES,
            "hourly": {"time": [t.strftime("%Y-%m-%dT%H:%M") for t in stamps],
                       **{key: [10.0] * len(stamps) for key in weather.VARIABLES}}}


@pytest.fixture
def fake_weather(monkeypatch):
    original = weather.fetch_forecast
    def get_forecast(start, end, location):
        monkeypatch.setattr(weather, "request_json", lambda url: forecast_payload(start, end))
        return original(start, end, location)

    return get_forecast


@pytest.fixture
def package_setup(monkeypatch, fake_weather):
    from astropy import units as u
    from ls4 import animate, scheduler
    from PIL import Image

    original = scheduler.night_window
    monkeypatch.setattr(scheduler, "night_window", lambda day: (original(day)[0], original(day)[0] + 2*u.hour))
    monkeypatch.setattr(weather, "fetch_forecast", fake_weather)
    # Keep unit runs small; a full-size real GIF is also exercised during manual QA.
    def gif(plan, path, **kwargs):
        frames = [Image.new("RGB", (20, 20), c) for c in ("black", "white")]
        frames[0].save(path, save_all=True, append_images=frames[1:], duration=200, loop=0)
        return path

    monkeypatch.setattr(animate, "make_gif", gif)


def test_forecast_window_units_and_site(monkeypatch):
    start = datetime(2026, 10, 21, 23, 37, tzinfo=timezone.utc)
    end = start + timedelta(hours=8, minutes=13)
    payload = forecast_payload(start, end)
    urls = []
    monkeypatch.setattr(weather, "request_json", lambda url: urls.append(url) or payload)
    result = weather.fetch_forecast(start, end, {"latitude": -29.25, "longitude": -70.73, "elevation_m": 2400})
    assert result["hourly"][0]["utc"] == "2026-10-21T23:00:00+00:00"
    assert result["hourly"][-1]["utc"] == "2026-10-22T08:00:00+00:00"
    assert "elevation=2400" in urls[0] and "wind_speed_unit=ms" in urls[0]
    assert result["model_issue_time_utc"] is None
    payload["hourly"]["wind_gusts_10m"][0] = None
    with pytest.raises(ValueError, match="Missing forecast"):
        weather.fetch_forecast(start, end, {"latitude": -29, "longitude": -70, "elevation_m": 2400})


@pytest.mark.parametrize("fault", ["missing-hour", "wrong-units", "nonfinite"])
def test_forecast_rejects_incomplete_or_misleading_data(monkeypatch, fault):
    start = datetime(2026, 10, 21, 23, 30, tzinfo=timezone.utc)
    end = start + timedelta(hours=8)
    payload = forecast_payload(start, end)
    payload["hourly_units"] = dict(payload["hourly_units"])
    if fault == "missing-hour":
        payload["hourly"]["time"].pop(2)
    elif fault == "wrong-units":
        payload["hourly_units"]["wind_speed_10m"] = "km/h"
    else:
        payload["hourly"]["temperature_2m"][0] = float("nan")
    monkeypatch.setattr(weather, "request_json", lambda url: payload)
    with pytest.raises(ValueError):
        weather.fetch_forecast(start, end, {"latitude": -29, "longitude": -70, "elevation_m": 2400})


def test_safe_rerun_and_local_database_isolation(tmp_path, package_setup):
    original_db = GRID_DB.read_bytes() if GRID_DB.exists() else None
    first = nightly.generate(date(2026, 10, 21), tmp_path / "first", refresh_iers=False)
    second = nightly.generate(date(2026, 10, 21), tmp_path / "second", refresh_iers=False)
    assert first["generation"] == "success"
    assert first["images"] > 0 and first["images"] == first["complete_pairs"] * 2
    assert first["history_source"] == "generated-plans" and first["gif"] == "night.gif"
    assert first["weather"]["status"] == "available" and not first["warnings"]
    assert first["obsplan"] == "20261022.obsplan"
    assert (tmp_path / "first" / first["plan_csv"]).read_bytes() == (tmp_path / "second" / second["plan_csv"]).read_bytes()
    for name, digest in first["files"].items():
        assert automation.sha256(tmp_path / "first" / name) == digest
    assert (GRID_DB.read_bytes() if GRID_DB.exists() else None) == original_db
    with pytest.raises(FileExistsError):
        nightly.generate(date(2026, 10, 21), tmp_path / "first", refresh_iers=False)


def test_weather_and_gif_failures_preserve_valid_plan(tmp_path, package_setup, monkeypatch):
    from ls4 import animate

    def fail(*args, **kwargs):
        raise RuntimeError("service unavailable")

    monkeypatch.setattr(weather, "fetch_forecast", fail)
    monkeypatch.setattr(animate, "make_gif", fail)
    result = nightly.generate(date(2026, 10, 21), tmp_path / "package", refresh_iers=False)
    assert result["generation"] == "success"
    assert result["weather"]["status"] == "unavailable"
    assert len(result["warnings"]) == 2 and "gif" not in result
    assert "Weather unavailable" in (tmp_path / "package" / "REPORT.md").read_text()


def test_empty_plan_cannot_be_published(tmp_path, package_setup, monkeypatch):
    from ls4 import scheduler

    monkeypatch.setattr(scheduler, "plan_night", lambda *args, **kwargs: pd.DataFrame())
    with pytest.raises(Exception):
        nightly.generate(date(2026, 10, 21), tmp_path / "invalid", refresh_iers=False)
    result = json.loads((tmp_path / "invalid" / "manifest.json").read_text())
    assert result["generation"] == "failed" and result["stage"] == "validation"
    with pytest.raises(ValueError, match="invalid plan"):
        automation.GitHub("owner/repo", "unused").publish(tmp_path / "invalid", "unused", "a"*40)


def test_bootstrap_ignores_same_night_and_future_plans(tmp_path, monkeypatch):
    db = fresh_db(tmp_path / "grid.db")
    field = load_history(db)["base_field"].iloc[0]
    (tmp_path / "plans").mkdir()
    for stamp, mjd in [("20261020", 61333.0), ("20261021", 61334.0), ("20261022", 61335.0)]:
        pd.DataFrame({"target": [field, field + "_dither"],
                      "start_time_mjd": [mjd, mjd + .03]}).to_csv(tmp_path / "plans" / f"{stamp}.csv", index=False)
    monkeypatch.setattr(nightly, "ROOT", tmp_path)
    nightly.restore_history(date(2026, 10, 21), db, GRID_CSV, use_github=False)
    values = load_history(db).set_index("base_field")["last_scheduled_mjd"]
    assert values[field] == 61333.03


def test_release_history_excludes_drafts_and_uses_evening_then_revision():
    releases = [
        {"tag_name": "ls4-20261022-999-1", "draft": False},
        {"tag_name": "ls4-20261021-100-1", "draft": False},
        {"tag_name": "ls4-20261020-150-1", "draft": False},
        {"tag_name": "ls4-20261020-150-2", "draft": False},
        {"tag_name": "ls4-20261020-200-1", "draft": True},
        {"tag_name": "v2.0", "draft": False},
    ]
    client = automation.GitHub("owner/repo", "unused")
    client.api = lambda path: releases
    assert client.latest_nightly(date(2026, 10, 21))["tag_name"] == "ls4-20261020-150-2"
    assert client.latest_nightly(date(2026, 10, 21), before=False)["tag_name"] == "ls4-20261021-100-1"


def test_late_rerun_becomes_newest_revision_even_with_lower_run_id():
    client = automation.GitHub("owner/repo", "unused")
    client.api = lambda path: [
        {"tag_name": "ls4-20261020-200-1", "draft": False, "published_at": "2026-10-20T17:00:00Z"},
        {"tag_name": "ls4-20261020-100-2", "draft": False, "published_at": "2026-10-20T18:00:00Z"},
    ]
    assert client.latest_nightly(date(2026, 10, 21))["tag_name"] == "ls4-20261020-100-2"


@pytest.mark.parametrize("fault", ["checksum", "grid", "future", "field-universe"])
def test_bad_release_history_never_resets_silently(tmp_path, monkeypatch, fault):
    db = fresh_db(tmp_path / "grid.db")
    fields = dict(zip(load_history(db)["base_field"], load_history(db)["last_scheduled_mjd"]))
    if fault == "future":
        fields[next(iter(fields))] = 999999.0
    if fault == "field-universe":
        fields.pop(next(iter(fields)))
    snapshot = {"source": "generated-plans", "through_night": "2026-10-20",
                "grid_sha256": "wrong" if fault == "grid" else automation.sha256(GRID_CSV), "fields": fields}
    raw = json.dumps(snapshot).encode()
    manifest = {"generation": "success", "evening_date": "2026-10-20",
                "files": {"history.json": "wrong" if fault == "checksum" else hashlib.sha256(raw).hexdigest()}}
    client = automation.GitHub("owner/repo", "unused")
    client.latest_nightly = lambda day: {"html_url": "https://example.com/night"}
    client.read_asset = lambda release, name: raw if name == "history.json" else json.dumps(manifest).encode()
    monkeypatch.setattr(nightly, "github_client", lambda: client)
    with pytest.raises(ValueError):
        nightly.restore_history(date(2026, 10, 21), db, GRID_CSV, use_github=True)
    assert (load_history(db)["last_scheduled_mjd"] == 0).all()


def test_publication_is_draft_until_every_asset_uploads(tmp_path, monkeypatch):
    (tmp_path / "REPORT.md").write_text("Nightly report")
    manifest = {"generation": "success", "evening_date": "2026-10-21",
                "provenance": {"commit": "a"*40, "working_tree_dirty": False},
                "files": {"REPORT.md": automation.sha256(tmp_path / "REPORT.md")}}
    automation.write_json(tmp_path / "manifest.json", manifest)
    calls = []
    client = automation.GitHub("owner/repo", "unused")
    def api(path, **kwargs):
        calls.append((path, kwargs))
        return {"id": 1, "upload_url": "https://uploads.github.com/example{?name}", "html_url": "https://example.com/night"}

    client.api = api
    def upload(url, **kwargs):
        assert len(calls) == 1 and calls[0][1]["payload"]["draft"] is True
        return json.dumps({"state": "uploaded", "size": len(kwargs["data"]),
                           "browser_download_url": "https://example.com/file"}).encode()

    monkeypatch.setattr(automation, "request_bytes", upload)
    result = client.publish(tmp_path, "ls4-20261021-1-1", "a"*40)
    assert result["url"] == "https://example.com/night"
    assert calls[-1][1]["payload"]["draft"] is False
    (tmp_path / "REPORT.md").write_text("tampered")
    with pytest.raises(ValueError, match="checksum"):
        client.publish(tmp_path, "ls4-20261021-2-1", "a"*40)


def test_failed_asset_upload_never_exposes_release(tmp_path, monkeypatch):
    (tmp_path / "REPORT.md").write_text("Nightly report")
    automation.write_json(tmp_path / "manifest.json", {"generation": "success", "evening_date": "2026-10-21",
        "provenance": {"commit": "a"*40, "working_tree_dirty": False},
        "files": {"REPORT.md": automation.sha256(tmp_path / "REPORT.md")}})
    calls = []
    client = automation.GitHub("owner/repo", "unused")
    client.api = lambda path, **kwargs: calls.append(kwargs) or {"id": 1, "upload_url": "https://uploads.github.com/x{?name}"}
    monkeypatch.setattr(automation, "request_bytes", lambda *a, **kw: b'{"state":"failed"}')
    with pytest.raises(RuntimeError, match="did not complete"):
        client.publish(tmp_path, "ls4-20261021-1-1", "a"*40)
    assert len(calls) == 1 and calls[0]["payload"]["draft"] is True


def test_slack_private_gif_upload_and_failure_fallback(tmp_path, monkeypatch):
    gif = tmp_path / "night.gif"
    gif.write_bytes(b"GIF89a")
    client = automation.Slack("unused", "C0CB1J9LZPF")
    calls = []
    client.api = lambda method, payload: calls.append((method, payload)) or {
        "upload_url": "https://files.slack.com/upload", "file_id": "F123"}
    monkeypatch.setattr(automation, "request_bytes", lambda *args, **kwargs: b"ok")
    client.send("plan ready", gif)
    assert calls[-1][0] == "files.completeUploadExternal"
    assert calls[-1][1]["channel_id"] == "C0CB1J9LZPF"
    assert calls[-1][1]["initial_comment"] == "plan ready"
    def fail(*args, **kwargs):
        raise RuntimeError("upload unavailable")

    monkeypatch.setattr(automation, "request_bytes", fail)
    with pytest.raises(RuntimeError):
        client.send("plan ready", gif)
    assert calls[-1][0] == "chat.postMessage" and "GIF delivery" in calls[-1][1]["text"]


def test_slack_blocks_layout():
    weather = {"status": "available", "summary": {
        "temperature_min_c": 8.0, "temperature_max_c": 12.5, "cloud_cover_max_percent": 45,
        "precipitation_probability_max_percent": 3, "humidity_max_percent": 60,
        "wind_max_ms": 6.2, "gust_max_ms": 9.8}}
    manifest = {"generation": "success", "evening_date": "2026-10-21", "images": 200,
                "complete_pairs": 100, "warnings": ["IERS <stale> & old"], "weather": weather}
    publication = {"url": "https://github.com/x/releases/1",
                   "assets": {"20261022.obsplan": "https://github.com/x/20261022.obsplan"}}
    blocks = automation.slack_blocks(manifest, publication, run_url="https://github.com/x/runs/1")
    rendered = json.dumps(blocks)
    assert blocks[0]["type"] == "header" and "2026-10-21" in blocks[0]["text"]["text"]
    assert ":warning: *Partial Success*" in rendered and ":partly_sunny: *Max cloud*" in rendered
    assert "IERS &lt;stale&gt; &amp; old" in rendered
    assert "<https://github.com/x/20261022.obsplan|Download 20261022.obsplan>" in rendered
    footer = blocks[-1]["elements"]
    assert footer[0]["image_url"] == automation.GITHUB_LOGO and "Run logs" in footer[1]["text"]

    failed = json.dumps(automation.slack_blocks({"stage": "history", "error": "bad `tick`"}))
    assert ":x: *Failed* at `history`" in failed and "bad 'tick'" in failed
    assert "image_url" not in failed  # no links, so no GitHub footer


def test_slack_post_falls_back_to_plain_text_on_invalid_blocks():
    client = automation.Slack("unused", "C0CB1J9LZPF")
    calls = []

    def api(method, payload):
        calls.append(payload)
        if "blocks" in payload:
            raise RuntimeError("Slack chat.postMessage: invalid_blocks")
        return {"ok": True}

    client.api = api
    client.send("plain", blocks=[automation.section("rich")])
    assert "blocks" in calls[0] and "blocks" not in calls[1] and calls[1]["text"] == "plain"


def test_slack_file_upload_methods_are_form_encoded(monkeypatch):
    from urllib.parse import parse_qs

    sent = []
    monkeypatch.setattr(automation, "request_bytes", lambda url, **kwargs: sent.append((url, kwargs)) or b'{"ok": true}')
    client = automation.Slack("xoxb-test", "C0CB1J9LZPF")
    client.api("files.getUploadURLExternal", {"filename": "night.gif", "length": 6})
    client.api("files.completeUploadExternal", {"files": [{"id": "F123"}], "channel_id": "C0CB1J9LZPF"})
    client.api("chat.postMessage", {"channel": "C0CB1J9LZPF", "text": "hi"})
    (_, upload), (_, complete), (_, post) = sent
    assert upload["headers"]["Content-Type"] == "application/x-www-form-urlencoded"
    assert parse_qs(upload["data"].decode()) == {"filename": ["night.gif"], "length": ["6"]}
    assert json.loads(parse_qs(complete["data"].decode())["files"][0]) == [{"id": "F123"}]
    assert post["headers"]["Content-Type"] == "application/json"


def test_notification_distinguishes_generation_publication_and_degradation():
    manifest = {"generation": "success", "evening_date": "2026-10-21", "images": 200,
                "complete_pairs": 100, "warnings": []}
    publication = {"url": "https://example.com/night", "assets": {"night.obsplan": "https://example.com/plan"}}
    assert "— SUCCESS" in automation.slack_message(manifest, publication)
    assert "PARTIAL SUCCESS" in automation.slack_message(manifest, publication, job_status="failure")
    manifest["warnings"].append("Weather unavailable")
    assert "PARTIAL SUCCESS" in automation.slack_message(manifest, publication)
    assert "Publication failed" in automation.slack_message(manifest)
    assert "Generation failed at history" in automation.slack_message({"stage": "history"})
    text = automation.slack_message({"stage": "validation", "error": "ValueError: <Quantity 3 m> & more"},
                                    run_url="https://example.com/run")
    assert "&lt;Quantity 3 m&gt; &amp; more" in text and "<https://example.com/run|Run logs>" in text


def test_asset_redirect_drops_credentials_cross_host():
    from urllib.request import Request

    request = Request("https://api.github.com/asset", headers={"Authorization": "Bearer private"})
    redirected = automation.SafeRedirect().redirect_request(
        request, None, 302, "Found", {}, "https://release-assets.githubusercontent.com/asset")
    assert not redirected.has_header("Authorization")


def test_http_errors_do_not_expose_upload_signatures_or_retry_writes(monkeypatch):
    from urllib.error import HTTPError

    calls = []
    class BrokenOpener:
        def open(self, request, timeout):
            calls.append(request)
            raise HTTPError(request.full_url, 503, "unavailable", {}, None)

    monkeypatch.setattr(automation, "build_opener", lambda *args: BrokenOpener())
    monkeypatch.setattr(automation.time, "sleep", lambda seconds: None)
    with pytest.raises(RuntimeError) as error:
        automation.request_bytes("https://files.slack.com/upload?secret=private", method="POST", data=b"gif")
    assert "private" not in str(error.value) and len(calls) == 1
    with pytest.raises(RuntimeError):
        automation.request_bytes("https://files.slack.com/upload?secret=private")
    assert len(calls) == 4


def test_published_snapshot_restores_real_history_without_adding_today(tmp_path, monkeypatch):
    db = fresh_db(tmp_path / "grid.db")
    values = dict(zip(load_history(db)["base_field"], load_history(db)["last_scheduled_mjd"]))
    field = next(iter(values))
    values[field] = 61333.25
    snapshot = {"source": "generated-plans", "through_night": "2026-10-20",
                "grid_sha256": automation.sha256(GRID_CSV), "fields": values}
    raw = json.dumps(snapshot).encode()
    manifest = {"generation": "success", "evening_date": "2026-10-20",
                "files": {"history.json": hashlib.sha256(raw).hexdigest()}}
    client = automation.GitHub("owner/repo", "unused")
    client.latest_nightly = lambda day: {"html_url": "https://example.com/history"}
    client.read_asset = lambda release, name: raw if name == "history.json" else json.dumps(manifest).encode()
    monkeypatch.setattr(nightly, "github_client", lambda: client)
    for _ in range(2):
        nightly.restore_history(date(2026, 10, 21), db, GRID_CSV, use_github=True)
    assert load_history(db).set_index("base_field").loc[field, "last_scheduled_mjd"] == 61333.25


def test_github_history_error_fails_generation_instead_of_bootstrapping(tmp_path, package_setup, monkeypatch):
    client = automation.GitHub("owner/repo", "unused")
    def fail(day):
        raise RuntimeError("GitHub unavailable")

    client.latest_nightly = fail
    monkeypatch.setattr(nightly, "github_client", lambda: client)
    with pytest.raises(RuntimeError, match="GitHub unavailable"):
        nightly.generate(date(2026, 10, 21), tmp_path / "failed", use_github=True, refresh_iers=False)
    manifest = json.loads((tmp_path / "failed" / "manifest.json").read_text())
    assert manifest["stage"] == "history" and manifest["generation"] == "failed"
