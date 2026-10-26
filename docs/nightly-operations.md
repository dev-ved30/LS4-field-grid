# Nightly observing plans

The production entry point is **Actions → Nightly observing plan**. The scheduled
run starts at **16:17 UTC daily**. Each run schedules the local evening at La Silla,
exports the telescope file, renders the GIF, retrieves weather, and publishes a
permanent GitHub Release. Slack channel **C0CB1J9LZPF** receives the status, GIF,
weather summary, and links. A separate deadline check runs at **19:47 UTC**.

## Activate

1. Review and merge the package code and workflows into the repository's default
   branch. The working `v2` branch and local science-code changes must be included
   in the version being reviewed. GitHub scheduled workflows use the default branch.
2. Create/install a Slack app with bot scopes `chat:write` and `files:write`. Invite
   the app to channel `C0CB1J9LZPF`. A plain incoming webhook cannot upload GIFs to
   a private channel.
3. In **GitHub → Settings → Secrets and variables → Actions**, add the bot token as
   repository secret **SLACK_BOT_TOKEN**. Enter it directly in GitHub, never in chat
   or a tracked file. Optionally set repository variable **SLACK_CHANNEL_ID** to
   override the default channel, and **SLACK_ALERT_USER** to the Slack member ID
   tagged on partial, failed, and deadline-check messages (default `U0888078546`).
4. Enable GitHub Actions and allow the workflow's `contents: write` permission so
   its built-in token can create releases. No personal GitHub token is needed.
5. Run the workflow manually with **publish unchecked** on `v2` or the default
   branch. Inspect the artifact's CSV, telescope file, GIF, weather report, and
   manifest. This mode uses existing production history but publishes nothing and
   sends no Slack message.
6. Run manually on the default branch with **publish checked**. Confirm the release
   and Slack delivery. Run the deadline-check workflow once to verify the monitor.

Standard GitHub-hosted Linux runners are sufficient. Dependencies and action
revisions are pinned. Science packages are installed from `requirements-nightly.lock`;
the source runs directly from the checked-out repository.

## Nightly package

Every published revision has a unique tag
`ls4-YYYYMMDD-<GitHub run ID>-<run attempt>`. Earlier revisions remain available.
Assets are uploaded to a draft release; the release becomes public to repository
readers only after every asset has uploaded successfully. A failed upload leaves
an unpublished draft for diagnosis. Never use draft releases operationally.
Publication also requires a clean checkout matching the manifest's exact code commit.

The release contains:

- `YYYYMMDD.csv`: the local evening's schedule.
- `<following morning YYYYMMDD>.obsplan`: telescope input, retaining the existing format.
- `night.gif`: pointings and dithered revisits through the night, if rendering succeeds.
- `REPORT.md`: observing window, metrics, hourly weather, history assumptions, and provenance.
- `weather.json`: forecast data, site/grid coordinates, units, source, and retrieval timestamp.
- `history.json`: cumulative per-field scheduling history, labelled as assumed execution.
- `manifest.json`: generation status, warnings, code/grid/lock identity, and file checksums.

Releases are durable; Actions artifacts are retained for 30 days for diagnostics.
Access to private-repository links requires GitHub access. The GIF is uploaded
directly to Slack, so channel members can view it without fetching a private GitHub URL.

## History and reruns

**History is based on generated plans, assuming execution. It is not confirmed
observing history.** Bad-weather losses and deviations are not automatically reconciled.

Generation creates a temporary SQLite database and never modifies
`assets/LS4_field_grid.db`. It restores the snapshot from the newest revision of
the most recent published night strictly before the requested evening. If no
previous release exists, it bootstraps complete pairs from repository
`plans/*.csv` before that date. Archived plans are excluded. Future and same-night
plans are excluded in both cases. A snapshot/grid mismatch, corrupt history, or
GitHub API failure aborts generation rather than resetting cadence silently.

Only successfully published nights advance production history. Replanning the
same date rebuilds from the preceding night's snapshot; retries cannot advance
history twice. Runs are serialized and each revision has its own URL. The newest
revision for a night is used for subsequent history; operators must communicate
which revision the telescope should use. A rerun does not retract a file already
downloaded by an observer. Backfilling an old night does not recompute snapshots
of subsequent nights; regenerate those chronologically if history needs correction.

To rerun: **Actions → Nightly observing plan → Run workflow**, select the default
branch, enter the local evening as `YYYY-MM-DD`, and leave publish checked.

## Weather and precision

Open-Meteo hourly forecasts use La Silla's coordinates and observatory elevation.
UTC samples bracket astronomical evening and morning twilight. Temperature, cloud
cover, precipitation probability, humidity, wind, and gusts are reported. Samples
just outside twilight are included to cover the boundary hours. The report includes
Open-Meteo attribution (CC BY 4.0) and the retrieval timestamp; the forecast API
does not supply the underlying model issue time. Missing/null/incomplete data are
reported as unavailable. Forecasts outside the provider's available dates will also
be unavailable on old backfills.

Weather is advisory. It does not cancel observing or change field selection.
Site weather limits and the telescope's operating procedures remain authoritative.
This initial provider choice suits noncommercial public research; review its service
terms if the usage changes.

The production entry point attempts to refresh IERS Earth-orientation tables.
Refresh warnings or inadequate coverage appear in the package and Slack status.
The regular CLI's offline configuration is unchanged.

## Status and recovery

- **SUCCESS:** validated plan published, weather and GIF available, no generation warnings.
- **PARTIAL SUCCESS:** validated plan published, but weather, GIF, precision, or another workflow step requires attention.
- **FAILED:** generation or publication failed. The Slack message links to logs and
  any available package. An unvalidated output is never published operationally.
- Slack/API upload failure makes the notification job fail; its error is visible
  in Actions. If GIF upload fails, the client attempts a text message with the plan link.

The separate notification job uses only Python's standard library. It can report
dependency-install and generation failures. Manual dry runs send no messages.
Canceled workflows and GitHub outages can prevent even the notification job from running.

The deadline checker alerts if the day's release is missing, cannot be verified,
or contains warnings. It runs separately from generation, but still depends on
GitHub Actions. **It cannot detect a GitHub-wide scheduling outage independently.**
For a strict operational deadline, add an externally hosted dead-man monitor and
a separately runnable fallback. GitHub schedules may be delayed/dropped, and
scheduled workflows in inactive public repositories can be disabled after 60 days.
Enable GitHub failure notifications as another way to spot Slack delivery failures.

## Website

**Actions → Observing-plan website** rebuilds the GitHub Pages site after every
nightly run on the default branch, after pushes that change `plans/` or `ls4/`, and
on demand. It reads repository plans and the newest revision of every published
release (`ls4 site --releases`). A failed site build never affects the nightly plan
or its Slack notification.

To activate it, open **GitHub → Settings → Pages** and set **Source** to **GitHub
Actions**, then run the workflow once manually. The site's address is
`https://<owner>.github.io/<repository>/`. Pages sites are public, as is this
repository, so plans, `.obsplan` files, and weather summaries become publicly viewable.

## Local generation and dependency updates

```bash
python -m pip install --require-hashes -r requirements-nightly.lock
python -m ls4.nightly generate --date 2026-10-21 --output-dir outputs/review-20261021
```

Local generation bootstraps from saved CSVs. Add `--github-history` with
`GITHUB_REPOSITORY` and `GITHUB_TOKEN` in the environment to use released history.
The output directory must not already exist; choose a new directory for a revision.
Local generation never sends Slack messages or publishes a release.

After deliberately updating dependencies, regenerate the lock and run the tests:

```bash
uv pip compile pyproject.toml --extra dev --python-version 3.11 --universal \
  --generate-hashes --output-file requirements-nightly.lock --no-header
python -m pytest -q
```

Review any scientific dependency changes and refresh the lock regularly so bundled
Earth-orientation data also remains current.

## Service documentation

- [GitHub scheduled workflows](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#schedule)
- [GitHub Releases API](https://docs.github.com/en/rest/releases/releases)
- [Slack messages](https://docs.slack.dev/reference/methods/chat.postMessage/)
- [Slack file upload](https://docs.slack.dev/reference/methods/files.getUploadURLExternal/)
- [Open-Meteo forecasts](https://open-meteo.com/en/docs) and [terms](https://open-meteo.com/en/terms)
