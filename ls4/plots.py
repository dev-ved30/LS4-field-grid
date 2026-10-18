"""Plotly figures: field grid, observing plans, coverage maps, and survey statistics."""

from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from astropy import units as u
from astropy.coordinates import SkyCoord
from m4opt.fov import footprint

from ls4.config import region

SKY_GEO = dict(showland=False, showcountries=False, showcoastlines=False,
               lonaxis=dict(showgrid=True, gridwidth=0.5), lataxis=dict(showgrid=True, gridwidth=0.5))


def _finish(fig, output=None, show=True):
    """Save a figure as standalone HTML and/or open it in the browser."""
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(output, include_plotlyjs="cdn")
        print(f"  wrote {output}")
    if show:
        fig.show(config={"scrollZoom": True, "displayModeBar": True})
    return fig


def footprint_outline(region):
    if hasattr(region, "to_polygon_sky"):
        vertices = region.to_polygon_sky().vertices
    elif hasattr(region, "vertices"):
        vertices = region.vertices
    else:
        return None, None
    lon = vertices.ra.wrap_at(180 * u.deg).deg
    lat = vertices.dec.deg
    return np.append(lon, lon[0]), np.append(lat, lat[0])


def _normalize_footprint_result(footprint_result):
    if hasattr(footprint_result, "to_polygon_sky") or hasattr(footprint_result, "vertices"):
        return [footprint_result]
    if isinstance(footprint_result, (list, tuple, np.ndarray)):
        if len(footprint_result) == 0:
            return []
        first_item = footprint_result[0]
        if isinstance(first_item, (list, tuple, np.ndarray)):
            return list(first_item)
        return list(footprint_result)
    return [footprint_result]


def field_grid(grid, output=None, show=True, example_fields=(1001, 1002, 1005, 1007)):
    """Grid centers plus a few example footprints, as a sanity check of the tiling."""
    vertices = SkyCoord(ra=grid["ra_deg"].to_numpy() * u.deg,
                        dec=grid["dec_deg"].to_numpy() * u.deg, frame="icrs")
    fig = go.Figure(go.Scattergeo(lon=vertices.ra.wrap_at(180 * u.deg).deg, lat=vertices.dec.deg,
                                  mode="markers", marker=dict(size=4, color="red"),
                                  name="Grid vertices"))
    examples = [index for index in example_fields if index < len(vertices)]
    for i, r in enumerate(footprint(region, vertices[examples])):
        lon, lat = footprint_outline(r)
        if lon is None:
            continue
        fig.add_trace(go.Scattergeo(lon=lon, lat=lat, mode="lines",
                                    line=dict(width=1.5, color="black"),
                                    name="FoV footprint" if i == 0 else None, showlegend=(i == 0)))
    fig.update_layout(title="LS4 Field Grid (interactive)", height=750,
                      geo=dict(projection_type="orthographic", **SKY_GEO))
    return _finish(fig, output, show)


def observing_plan(plan, output=None, show=True):
    """Animated footprint-by-footprint view of one plan (slider and play button)."""
    plan = plan.copy()
    required_columns = {"target", "start time (UTC)", "ra", "dec"}
    missing_columns = required_columns.difference(plan.columns)
    if missing_columns:
        raise ValueError(f"plan must include {sorted(required_columns)} columns")

    if "block_number" not in plan.columns:
        raise ValueError("plan must include 'block_number' for coloring")
    block_column = "block_number"

    plan = plan[plan["target"].astype(str) != "TransitionBlock"].copy()
    plan["start time (UTC)"] = pd.to_datetime(plan["start time (UTC)"])
    if "end time (UTC)" in plan.columns:
        plan["end time (UTC)"] = pd.to_datetime(plan["end time (UTC)"])
    plan["ra"] = pd.to_numeric(plan["ra"], errors="coerce")
    plan["dec"] = pd.to_numeric(plan["dec"], errors="coerce")
    plan.sort_values("start time (UTC)", inplace=True)
    plan.reset_index(drop=True, inplace=True)

    parity_colors = {
        "Even": "rgba(65,105,225,1.0)",
        "Odd": "rgba(220,20,60,1.0)",
    }
    parity_legend = {
        "Even": "first_pointing",
        "Odd": "dither",
    }

    fig = go.Figure()
    scheduled_rows = []
    row_trace_counts = []
    row_trace_starts = []

    legend_traces = [
        go.Scattergeo(
            lon=[None],
            lat=[None],
            mode="lines",
            line=dict(width=2.5, color=parity_colors["Even"]),
            name="First Pointing",
            visible="legendonly",
            legendgroup=parity_legend["Even"],
            hoverinfo="skip",
            showlegend=True,
        ),
        go.Scattergeo(
            lon=[None],
            lat=[None],
            mode="lines",
            line=dict(width=2.5, color=parity_colors["Odd"]),
            name="Revisit with dither",
            visible="legendonly",
            legendgroup=parity_legend["Odd"],
            hoverinfo="skip",
            showlegend=True,
        ),
    ]

    for trace in legend_traces:
        fig.add_trace(trace)

    for _, row in plan.iterrows():
        if pd.isna(row["ra"]) or pd.isna(row["dec"]):
            continue

        target_name = str(row["target"])
        target_coord = SkyCoord(
            ra=row["ra"] * u.deg,
            dec=row["dec"] * u.deg,
            frame="icrs",
        )
        footprint_regions = _normalize_footprint_result(footprint(region, target_coord))
        if not footprint_regions:
            continue

        start_time = row["start time (UTC)"]
        end_time = row["end time (UTC)"] if "end time (UTC)" in row and pd.notna(row["end time (UTC)"]) else None
        block_value = row[block_column]
        block_label = int(block_value) if pd.notna(block_value) and float(block_value).is_integer() else block_value
        block_parity = "Odd" if int(block_value) % 2 else "Even"
        block_color = parity_colors[block_parity]
        hover_text = f"Target: {target_name}<br>RA: {row['ra']:.3f}°<br>Dec: {row['dec']:.3f}°<br>Start: {start_time}"
        if end_time is not None:
            hover_text += f"<br>End: {end_time}"
        hover_text += f"<br>Block: {block_label} ({block_parity})"

        trace_count = 0
        for footprint_region in footprint_regions:
            lon, lat = footprint_outline(footprint_region)
            if lon is None:
                continue
            fig.add_trace(go.Scattergeo(
                lon=lon,
                lat=lat,
                mode="lines",
                line=dict(width=1.8, color=block_color),
                name=f"{block_parity} block",
                showlegend=False,
                legendgroup=parity_legend[block_parity],
                visible=False,
                hovertext=hover_text,
                hoverinfo="text",
            ))
            trace_count += 1

        if trace_count > 0:
            row_trace_starts.append(len(fig.data) - trace_count)
            scheduled_rows.append(row)
            row_trace_counts.append(trace_count)

    if scheduled_rows:
        animation_duration = 300
        steps = []
        total_traces = len(fig.data)
        scheduled_trace_total = total_traces - len(legend_traces)
        total_images = len(scheduled_rows)
        for i, row in enumerate(scheduled_rows):
            visible_trace_count = row_trace_starts[i] + row_trace_counts[i]
            visibility = [True] * len(legend_traces) + [False] * scheduled_trace_total
            for trace_index in range(len(legend_traces), visible_trace_count):
                visibility[trace_index] = True
            steps.append(dict(
                method="animate",
                args=[
                    [str(i)],
                    {
                        "mode": "immediate",
                        "frame": {"duration": animation_duration, "redraw": True},
                        "transition": {"duration": 0},
                    },
                ],
                label=f"{i + 1}: {row['target']}"
            ))

        initial_visibility = ["legendonly"] * len(legend_traces) + [False] * scheduled_trace_total
        for trace_index in range(len(legend_traces), len(legend_traces) + row_trace_counts[0]):
            initial_visibility[trace_index] = True
        for trace_index, visible in enumerate(initial_visibility):
            fig.data[trace_index].visible = visible

        frames = []
        for i, row in enumerate(scheduled_rows):
            visible_trace_count = row_trace_starts[i] + row_trace_counts[i]
            visibility = [True] * len(legend_traces) + [False] * scheduled_trace_total
            for trace_index in range(len(legend_traces), visible_trace_count):
                visibility[trace_index] = True
            frame_data = [go.Scattergeo(visible=visibility[trace_index]) for trace_index in range(len(legend_traces), total_traces)]
            frames.append(go.Frame(
                name=str(i),
                data=frame_data,
                traces=list(range(len(legend_traces), total_traces)),
                layout=go.Layout(title=f"LS4 Observing Plan through {row['target']}")
            ))

        fig.frames = frames

        fig.update_layout(
            sliders=[dict(
                active=0,
                currentvalue={"prefix": f"Schedule ({total_images} images): "},
                pad={"t": 50},
                steps=steps,
            )],
            updatemenus=[dict(
                type="buttons",
                direction="left",
                showactive=True,
                y=-0.20,
                x=0.5,
                xanchor="center",
                yanchor="top",
                buttons=[
                    dict(
                        label="Play / Pause",
                        method="animate",
                        args=[
                            None,
                            {
                                "fromcurrent": True,
                                "frame": {"duration": animation_duration, "redraw": True},
                                "transition": {"duration": 0},
                                "mode": "immediate",
                            },
                        ],
                        args2=[
                            [None],
                            {
                                "frame": {"duration": 0, "redraw": False},
                                "transition": {"duration": 0},
                                "mode": "immediate",
                            },
                        ],
                    ),
                ],
            )],
            legend=dict(title="Block parity", groupclick="togglegroup"),
            margin=dict(b=180)
        )


    fig.update_layout(
        title=f"LS4 Observing Plan ({len(scheduled_rows)} images)",
        height=1000,
        geo=dict(projection_type="orthographic", projection_rotation=dict(lon=0, lat=25), **SKY_GEO),
    )
    return _finish(fig, output, show)


def coverage_map(visits, title="Observed Sky Area", output=None, show=True):
    """HEALPix visit-count map; pixels with zero visits are hidden."""
    m = np.asarray(visits)
    observed_pixels = np.flatnonzero(m > 0)
    fig = go.Figure()
    if observed_pixels.size > 0:
        theta, phi = hp.pix2ang(hp.npix2nside(m.size), observed_pixels)
        fig.add_trace(go.Scattergeo(
            lon=((np.degrees(phi) + 180) % 360) - 180,
            lat=90.0 - np.degrees(theta),
            mode="markers",
            marker=dict(size=4, color=m[observed_pixels], colorscale="Viridis",
                        showscale=True, colorbar=dict(title="Visits")),
            hovertemplate="Lon: %{lon:.2f}°<br>Lat: %{lat:.2f}°<br>Visits: %{marker.color}<extra></extra>",
            showlegend=False,
        ))
    fig.update_layout(title=title, height=700, margin=dict(l=20, r=20, t=60, b=20),
                      hovermode="closest",
                      geo=dict(projection_type="mollweide", showframe=False, **SKY_GEO))
    return _finish(fig, output, show)


# ---- survey statistics ------------------------------------------------------

def visit_histogram(per_field, output_dir, show=False):
    counts = per_field["visits"]
    fig = go.Figure(go.Histogram(x=counts, nbinsx=int(counts.max()) + 1, marker_color="#4169ab"))
    fig.update_layout(title=f"Distribution of visits per field ({int((counts == 0).sum())} fields never visited)",
                      xaxis_title="Visits", yaxis_title="Number of fields", bargap=0.05)
    _finish(fig, Path(output_dir) / "visits_per_field_hist.html", show)


def cadence_histograms(cadence, per_field_cadence, output_dir, show=False):
    """Gaps between initial pointings to the same field, overall and per-field median."""
    fig = go.Figure(go.Histogram(x=cadence["gap_days"], nbinsx=60, marker_color="#8a63c4",
                                 name="initial-pointing gaps"))
    fig.update_layout(title="Cadence: time between initial pointings to the same field",
                      xaxis_title="Gap (days)", yaxis_title="Count (log)", yaxis_type="log", bargap=0.02)
    _finish(fig, Path(output_dir) / "cadence_gaps_hist.html", show)

    fig = go.Figure(go.Histogram(x=per_field_cadence["median_gap_days"], nbinsx=40, marker_color="#c48a63"))
    fig.update_layout(title="Per-field median gap between initial pointings",
                      xaxis_title="Median gap (days)", yaxis_title="Number of fields")
    _finish(fig, Path(output_dir) / "median_cadence_per_field_hist.html", show)


def cadence_cdf(cadence, output_dir, show=False):
    gaps = np.sort(cadence["gap_days"].to_numpy())
    percent = 100.0 * np.arange(1, len(gaps) + 1) / len(gaps)
    fig = go.Figure(go.Scatter(x=gaps, y=percent, mode="lines", line=dict(color="#4169ab", width=2),
                               hovertemplate="%{x:.2f} d: %{y:.1f}% of gaps<extra></extra>"))
    for value, level, label in ((np.median(gaps), 50, "median"), (np.percentile(gaps, 90), 90, "90th pct")):
        fig.add_hline(y=level, line=dict(width=1, dash="dot", color="#8a63c4"), opacity=0.5)
        fig.add_vline(x=value, line=dict(width=1, dash="dot", color="#8a63c4"),
                      annotation_text=f"{label}: {value:.2f} d", annotation_position="top", opacity=0.5)
    fig.update_layout(title=f"Cadence CDF: initial-pointing gaps shorter than x ({len(gaps)} gaps)",
                      xaxis_title="Gap between consecutive visits (days)", yaxis_title="% of gaps ≤ x",
                      xaxis_type="log", height=600)
    _finish(fig, Path(output_dir) / "cadence_gaps_cdf.html", show)


def nightly_summary(nightly, output_dir, show=False):
    """Images taken and time efficiency for each night."""
    nightly = nightly.sort_values("night").copy()
    nightly["time_efficiency_pct"] = 100 * nightly["used_minutes"] / nightly["night_minutes"]
    fig = go.Figure()
    fig.add_trace(go.Bar(x=nightly["night"], y=nightly["images"], name="images", marker_color="#44b8e8"))
    fig.add_trace(go.Scatter(x=nightly["night"], y=nightly["time_efficiency_pct"], name="time efficiency (%)",
                             yaxis="y2", mode="lines+markers", line=dict(color="#e652a0")))
    fig.update_layout(title="Nightly summary", xaxis_title="Plan", yaxis_title="Images",
                      yaxis2=dict(title="Time efficiency (%)", overlaying="y", side="right", range=[0, 105]),
                      legend=dict(x=0.01, y=0.99))
    _finish(fig, Path(output_dir) / "nightly_summary.html", show)
