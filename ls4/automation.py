"""Standard-library HTTP clients for durable releases and Slack notifications.

Keep this module independent of science dependencies so failure messages can be
sent even when the dependency-install or scheduling step fails.
"""

import hashlib
import json
import mimetypes
from pathlib import Path
import re
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

RELEASE_TAG = re.compile(r"^ls4-(\d{8})-(\d+)-(\d+)$")


class SafeRedirect(HTTPRedirectHandler):
    """Do not forward a GitHub credential to an asset-storage redirect host."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected and urlsplit(req.full_url).netloc != urlsplit(newurl).netloc:
            redirected.remove_header("Authorization")
        return redirected


def request_bytes(url, *, method="GET", data=None, headers=None):
    """Bound GET retries; never blindly retry a write that may have succeeded."""
    for attempt in range(3):
        try:
            req = Request(url, data=data, method=method, headers=headers or {})
            with build_opener(SafeRedirect()).open(req, timeout=30) as response:
                return response.read()
        except HTTPError as exc:
            retry = exc.code in (429, 500, 502, 503, 504)
            status = f"HTTP {exc.code}"
            exc.close()
        except (URLError, TimeoutError, OSError):
            retry, status = True, "network error"
        if method != "GET" or not retry or attempt == 2:
            # URLs can contain private upload signatures. Keep them out of logs.
            raise RuntimeError(f"{urlsplit(url).hostname}: {status}") from None
        time.sleep(2 ** attempt)


def request_json(url, *, method="GET", payload=None, headers=None):
    headers = dict(headers or {})
    data = None
    if payload is not None:
        data = json.dumps(payload, allow_nan=False).encode()
        headers["Content-Type"] = "application/json"
    return json.loads(request_bytes(url, method=method, data=data, headers=headers))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


class GitHub:
    def __init__(self, repository, token):
        if not re.fullmatch(r"[\w.-]+/[\w.-]+", repository):
            raise ValueError("Expected an owner/repository")
        self.base = f"https://api.github.com/repos/{repository}"
        self.headers = {"Authorization": f"Bearer {token}",
                        "Accept": "application/vnd.github+json",
                        "X-GitHub-Api-Version": "2022-11-28"}

    def api(self, path, *, method="GET", payload=None):
        return request_json(self.base + path, method=method, payload=payload, headers=self.headers)

    def latest_nightly(self, evening, *, before=True):
        """Find the newest revision of the most recent applicable evening.

        Scan every page: creation order differs from night order after backfills.
        Drafts are incomplete packages and must never become cadence history.
        """
        cutoff = evening.strftime("%Y%m%d")
        candidates = []
        for page in range(1, 101):
            releases = self.api(f"/releases?per_page=100&page={page}")
            for release in releases:
                match = RELEASE_TAG.fullmatch(release["tag_name"])
                if not match or release["draft"]:
                    continue
                stamp, run, attempt = match.groups()
                if ((stamp < cutoff) if before else (stamp == cutoff)):
                    # A rerun of an older run ID can be published after a newer
                    # run. Publication time determines the current revision.
                    published = release.get("published_at") or release.get("created_at") or ""
                    candidates.append(((stamp, published, int(run), int(attempt)), release))
            if len(releases) < 100:
                return max(candidates, key=lambda item: item[0])[1] if candidates else None
        raise RuntimeError("Release listing exceeded 10,000 entries; archive or revise history lookup")

    def read_asset(self, release, name):
        asset = next((a for a in release["assets"] if a["name"] == name), None)
        if asset is None:
            raise ValueError(f"Published package lacks {name}")
        headers = {**self.headers, "Accept": "application/octet-stream"}
        return request_bytes(self.base + f"/releases/assets/{asset['id']}", headers=headers)

    def publish(self, directory, tag, commit):
        """Upload to a draft; expose it only after every verified asset is present."""
        directory = Path(directory)
        manifest = json.loads((directory / "manifest.json").read_text())
        if manifest["generation"] != "success":
            raise ValueError("Cannot publish an invalid plan")
        if not RELEASE_TAG.fullmatch(tag) or not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise ValueError("Publication requires a nightly tag and exact commit SHA")
        provenance = manifest["provenance"]
        if provenance["commit"] != commit or provenance["working_tree_dirty"]:
            raise ValueError("Publication requires generation from the exact, clean release commit")
        files = manifest["files"]
        for name, digest in files.items():
            if Path(name).name != name or sha256(directory / name) != digest:
                raise ValueError(f"Package checksum mismatch: {name}")
        report = (directory / "REPORT.md").read_text()
        release = self.api("/releases", method="POST", payload={
            "tag_name": tag, "target_commitish": commit,
            "name": f"LS4 night {manifest['evening_date']} ({tag})",
            "body": report, "draft": True, "prerelease": False, "make_latest": "false",
        })
        uploads = release["upload_url"].split("{")[0]
        assets = {}
        for name in [*files, "manifest.json"]:
            path = directory / name
            result = json.loads(request_bytes(
                uploads + "?" + urlencode({"name": name}), method="POST", data=path.read_bytes(),
                headers={**self.headers, "Content-Type": mimetypes.guess_type(name)[0] or "application/octet-stream"}))
            if result.get("state") != "uploaded" or result.get("size") != path.stat().st_size:
                raise RuntimeError(f"Release upload did not complete: {name}")
            assets[name] = result["browser_download_url"]
        links = "\n\n## Downloads\n\n" + "\n".join(f"- [{name}]({url})" for name, url in assets.items())
        if "night.gif" in assets:
            links += f"\n\n![Nightly observing plan]({assets['night.gif']})\n"
        published = self.api(f"/releases/{release['id']}", method="PATCH",
                             payload={"draft": False, "make_latest": "false", "body": report + links})
        return {"url": published["html_url"], "tag": tag, "assets": assets}


class Slack:
    def __init__(self, token, channel):
        if not token or not re.fullmatch(r"[CG][A-Z0-9]+", channel):
            raise ValueError("SLACK_BOT_TOKEN and a valid Slack channel ID are required")
        self.channel = channel
        self.headers = {"Authorization": f"Bearer {token}"}

    # Despite the docs, these reject JSON bodies (invalid_arguments); send them
    # form-encoded as slack_sdk does, with structured values as JSON strings.
    FORM_METHODS = {"files.getUploadURLExternal", "files.completeUploadExternal"}

    def api(self, method, payload):
        url = f"https://slack.com/api/{method}"
        if method in self.FORM_METHODS:
            form = {key: value if isinstance(value, str) else json.dumps(value)
                    for key, value in payload.items()}
            result = json.loads(request_bytes(url, method="POST", data=urlencode(form).encode(), headers={
                **self.headers, "Content-Type": "application/x-www-form-urlencoded"}))
        else:
            result = request_json(url, method="POST", payload=payload, headers=self.headers)
        if not result.get("ok"):
            # Slack's error code contains no token or request body.
            raise RuntimeError(f"Slack {method}: {result.get('error', 'unknown_error')}")
        return result

    def post(self, text, blocks=None):
        """Post rich blocks, but never let a formatting rejection silence an alert."""
        message = {"channel": self.channel, "text": text, "unfurl_links": False, "unfurl_media": False}
        if blocks:
            try:
                return self.api("chat.postMessage", {**message, "blocks": blocks})
            except RuntimeError as exc:
                if "invalid_blocks" not in str(exc):
                    raise
                print(f"{exc}; resending as plain text", file=sys.stderr)
        return self.api("chat.postMessage", message)

    def send(self, text, gif=None, blocks=None):
        """Share a real GIF upload, including in private channels/repositories.

        `text` is the notification preview and plain-text fallback for `blocks`.
        """
        if gif is None:
            return self.post(text, blocks)
        gif = Path(gif)
        try:
            upload = self.api("files.getUploadURLExternal", {
                "filename": gif.name, "length": gif.stat().st_size,
                "alt_txt": "LS4 scheduled pointings and dithered revisits through the night",
            })
            request_bytes(upload["upload_url"], method="POST", data=gif.read_bytes(),
                          headers={"Content-Type": "application/octet-stream"})
            # Slack ignores initial_comment when blocks are supplied.
            message = {"blocks": blocks} if blocks else {"initial_comment": text}
            return self.api("files.completeUploadExternal", {
                "files": [{"id": upload["file_id"], "title": "LS4 nightly observing plan"}],
                "channel_id": self.channel, **message,
            })
        except Exception:
            note = ":warning: GIF delivery to Slack failed; use the package link."
            self.post(text + "\n" + note, blocks and [*blocks, context(note)])
            raise


GITHUB_LOGO = "https://github.githubassets.com/favicons/favicon.png"
STATUS_EMOJI = {"SUCCESS": ":white_check_mark:", "PARTIAL SUCCESS": ":warning:", "FAILED": ":x:"}


def slack_escape(text):
    """Slack mrkdwn treats &, <, > as control characters (e.g. astropy's <Quantity ...>)."""
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def section(text=None, fields=None):
    block = {"type": "section"}
    if text:
        block["text"] = {"type": "mrkdwn", "text": text[:3000]}
    if fields:
        block["fields"] = [{"type": "mrkdwn", "text": field[:2000]} for field in fields[:10]]
    return block


def context(*texts):
    return {"type": "context", "elements": [{"type": "mrkdwn", "text": t[:3000]} for t in texts]}


def header(text):
    return {"type": "header", "text": {"type": "plain_text", "text": text[:150], "emoji": True}}


def github_links(links):
    """Context footer: the GitHub mark, then `(url, label)` links; empty URLs are skipped."""
    links = [f"<{url}|{label}>" for url, label in links if url]
    if not links:
        return []
    return [{"type": "context", "elements": [
        {"type": "image", "image_url": GITHUB_LOGO, "alt_text": "GitHub"},
        {"type": "mrkdwn", "text": "   ·   ".join(links)}]}]


def notification_status(manifest, publication, job_status):
    if manifest.get("generation") == "success" and (publication or {}).get("url"):
        return "PARTIAL SUCCESS" if manifest.get("warnings") or job_status != "success" else "SUCCESS"
    return "FAILED"


def obsplan_link(publication):
    return next(((name, url) for name, url in (publication or {}).get("assets", {}).items()
                 if name.endswith(".obsplan")), None)


def slack_message(manifest, publication=None, *, job_status="success", run_url=""):
    """Plain-text version: the notification preview and fallback for `slack_blocks`."""
    publication = publication or {}
    valid = manifest.get("generation") == "success"
    published = bool(publication.get("url"))
    degraded = manifest.get("warnings", [])
    status = notification_status(manifest, publication, job_status)
    lines = [f"{STATUS_EMOJI[status]} *LS4 night {manifest.get('evening_date', 'unknown')} — {status}*"]
    if valid:
        lines.append(f"Plan generated and validated: {manifest['images']} exposures.")
    else:
        lines.append(f"Generation failed at {manifest.get('stage', 'setup')}. "
                     "No new operational plan was published.")
        if manifest.get("error"):
            lines.append(slack_escape(str(manifest["error"])[:500]))
    if valid and published and job_status != "success":
        lines.append("The plan is published, but another workflow step failed; inspect the run logs.")
    if valid and not published:
        lines.append("Publication failed or was disabled; no new operational package link is available.")
    weather = manifest.get("weather")
    if weather and weather.get("status") == "available":
        w = weather["summary"]
        lines.append(f"Temperature {w['temperature_min_c']:.1f}–{w['temperature_max_c']:.1f} °C; "
                     f"max cloud {w['cloud_cover_max_percent']:.0f}%; "
                     f"max rain chance {w['precipitation_probability_max_percent']:.0f}%; "
                     f"max humidity {w['humidity_max_percent']:.0f}%.")
    elif weather:
        lines.append("Weather unavailable: " + slack_escape(weather.get("error", "no forecast retrieved")))
    for warning in degraded:
        lines.append("Attention: " + slack_escape(warning))
    if published:
        lines.append(f"<{publication['url']}|Nightly package and weather report>")
        if obsplan_link(publication):
            lines.append(f"<{obsplan_link(publication)[1]}|Download telescope .obsplan>")
    if run_url:
        lines.append(f"<{run_url}|Run logs>")
    return "\n".join(lines)


def weather_block(weather):
    if weather.get("status") != "available":
        return section(":grey_question: *Weather unavailable* — "
                       + slack_escape(weather.get("error", "no forecast retrieved"))[:500])
    s = weather["summary"]
    cloud = s["cloud_cover_max_percent"]
    sky = ":sunny:" if cloud < 20 else ":partly_sunny:" if cloud < 60 else ":cloud:"
    return section(fields=[
        f":thermometer: *Temperature*\n{s['temperature_min_c']:.1f} to {s['temperature_max_c']:.1f} °C",
        f"{sky} *Max cloud*\n{cloud:.0f}%",
        f":umbrella_with_rain_drops: *Max rain chance*\n{s['precipitation_probability_max_percent']:.0f}%",
        f":droplet: *Max humidity*\n{s['humidity_max_percent']:.0f}%",
    ])


def slack_blocks(manifest, publication=None, *, job_status="success", run_url=""):
    """Block Kit layout of the same content as `slack_message`."""
    publication = publication or {}
    valid = manifest.get("generation") == "success"
    published = bool(publication.get("url"))
    status = notification_status(manifest, publication, job_status)
    blocks = [header(f"LS4 night {manifest.get('evening_date', 'unknown')}")]
    if valid:
        blocks.append(section(f"{STATUS_EMOJI[status]} *{status.title()}* — "
                              f"{manifest['images']} exposures planned and validated"))
    else:
        text = (f"{STATUS_EMOJI[status]} *Failed* at `{manifest.get('stage', 'setup')}` — "
                "no new operational plan was published.")
        if manifest.get("error"):
            # Backticks would close the code block early.
            text += "\n```" + slack_escape(str(manifest["error"])[:500]).replace("`", "'") + "```"
        blocks.append(section(text))
    if valid and published and job_status != "success":
        blocks.append(section(":rotating_light: The plan is published, but another workflow step "
                              "failed; inspect the run logs."))
    if valid and not published:
        blocks.append(section(":no_entry: Publication failed or was disabled; "
                              "no new operational package link is available."))
    if manifest.get("weather"):
        blocks.extend([{"type": "divider"}, weather_block(manifest["weather"])])
    if manifest.get("warnings"):
        blocks.append(section(":warning: *Needs attention*\n"
                              + "\n".join("• " + slack_escape(w) for w in manifest["warnings"])))
    obsplan = obsplan_link(publication) if published else None
    if obsplan:
        blocks.append(section(f"*<{obsplan[1]}|Download {obsplan[0]}>*"))
    blocks.append({"type": "divider"})
    blocks.extend(github_links([(publication.get("url"), "Nightly package & weather report"),
                                (run_url, "Run logs")]))
    return blocks
