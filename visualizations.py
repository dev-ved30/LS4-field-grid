import pandas as pd
import numpy as np
import healpy as hp
import plotly.graph_objects as go

from astropy import units as u
from astropy.coordinates import SkyCoord
from constants import *
from m4opt.fov import footprint

def footprint_outline(region):
    if hasattr(region, "to_polygon_sky"):
        vertices = region.to_polygon_sky().vertices
    elif hasattr(region, "vertices"):
        vertices = region.vertices
    else:
        return None, None
    lon = vertices.ra.wrap_at(180 * u.deg).deg
    lat = vertices.dec.deg
    lon = np.append(lon, lon[0])
    lat = np.append(lat, lat[0])
    return lon, lat

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

def plot_field_grid(field_grid_df):
    """Plot the field grid vertices and a few example footprints for testing.
    
    Arguments:
    field_grid_df -- a DataFrame containing the field grid vertices with columns "ra_deg", "dec_deg", and "Field Name"
    """

    vertices = SkyCoord(
        ra=field_grid_df["ra_deg"].to_numpy() * u.deg,
        dec=field_grid_df["dec_deg"].to_numpy() * u.deg,
        frame="icrs",
    )

    field_names = field_grid_df["Field Name"].astype(str)
    colors = []
    for f in field_names:
        ra_idx, dec_idx = f.split("_")
        ra_idx = int(ra_idx)
        dec_idx = int(dec_idx)
        if dec_idx % 2 == 0 and ra_idx % 2 == 0:
            colors.append("red")
        else:
            colors.append("gray")

    all_footprints = footprint(region, vertices)
    footprints = [all_footprints[1001], all_footprints[1002], all_footprints[1005], all_footprints[1007]]

    fig = go.Figure()

    fig.add_trace(go.Scattergeo(
        lon=vertices.ra.wrap_at(180 * u.deg).deg,
        lat=vertices.dec.deg,
        mode="markers",
        marker=dict(size=4, color=colors),
        name="Grid vertices"
    ))

    for i, r in enumerate(footprints):
        lon, lat = footprint_outline(r)
        if lon is None:
            continue
        fig.add_trace(go.Scattergeo(
            lon=lon,
            lat=lat,
            mode="lines",
            line=dict(width=1.5, color="black"),
            name="FoV footprint" if i == 0 else None,
            showlegend=(i == 0)
        ))

    fig.update_layout(
        title="LS4 Field Grid (interactive)",
        height=750,
        geo=dict(
            projection_type="orthographic",
            showland=False,
            showcountries=False,
            showcoastlines=False,
            lonaxis=dict(showgrid=True, gridwidth=0.5),
            lataxis=dict(showgrid=True, gridwidth=0.5)
        )
    )
    fig.show()

def plot_visible_fields(fields, visibility):

    names = [field.name for field in fields]
    ras = [field.coord.ra.wrap_at(180 * u.deg).deg for field in fields]
    decs = [field.coord.dec.deg for field in fields]
    colors = ['rgba(65,105,225,1.0)' if vis else 'rgba(128,128,128,0.3)' 
              for vis in visibility]

    
    fig = go.Figure()

    fig.add_trace(go.Scattergeo(
        lon=ras,
        lat=decs,
        mode="markers",
        marker=dict(size=5, color=colors),
        text=[f"RA: {ra:.2f}°, Dec: {dec:.2f}°<br>visible: {vis}<br>field: {name}" 
            for ra, dec, vis, name in zip(ras, decs, visibility, names)],
        hoverinfo="text",
        name="Grid vertices"
    ))

    fig.update_layout(
        title="LS4 Field Grid (180-day visibility)",
        height=1000,
        geo=dict(
            projection_type="orthographic",
            projection_rotation=dict(lon=0, lat=25),
            showland=False,
            showcountries=False,
            showcoastlines=False,
            lonaxis=dict(showgrid=True, gridwidth=0.5),
            lataxis=dict(showgrid=True, gridwidth=0.5)
        )
    )
    fig.show()

def plot_obs_plan(plan):
    plan = plan.copy()
    required_columns = {"target", "start time (UTC)", "ra", "dec"}
    missing_columns = required_columns.difference(plan.columns)
    if missing_columns:
        raise ValueError(f"plan must include {sorted(required_columns)} columns")

    if "block_number" in plan.columns:
        block_column = "block_number"
    elif "block_numer" in plan.columns:
        block_column = "block_numer"
    else:
        raise ValueError("plan must include 'block_number' or 'block_numer' for coloring")

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
        geo=dict(
            projection_type="orthographic",
            projection_rotation=dict(lon=0, lat=25),
            showland=False,
            showcountries=False,
            showcoastlines=False,
            lonaxis=dict(showgrid=True, gridwidth=0.5),
            lataxis=dict(showgrid=True, gridwidth=0.5),
        )
    )
    fig.show()


def plot_coverage_map(observed_mask, title="Observed Sky Area"):
    """Plot a HEALPix coverage map with Plotly.

    Parameters
    ----------
    observed_mask : array-like
        HEALPix mask or visit-count map. Values greater than zero are shown.
    """

    m = np.asarray(observed_mask)
    observed_pixels = np.flatnonzero(m > 0)

    fig = go.Figure()

    if observed_pixels.size > 0:
        nside = hp.npix2nside(m.size)
        theta, phi = hp.pix2ang(nside, observed_pixels)
        lon = np.degrees(phi)
        lon = ((lon + 180) % 360) - 180
        lat = 90.0 - np.degrees(theta)

        fig.add_trace(go.Scattergeo(
            lon=lon,
            lat=lat,
            mode="markers",
            marker=dict(
                size=4,
                color=m[observed_pixels],
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title="Visits"),
            ),
            hovertemplate="Lon: %{lon:.2f}°<br>Lat: %{lat:.2f}°<br>Visits: %{marker.color}<extra></extra>",
            showlegend=False,
        ))

    fig.update_layout(
        title=title,
        height=700,
        margin=dict(l=20, r=20, t=60, b=20),
        hovermode="closest",
        geo=dict(
            projection_type="mollweide",
            showland=False,
            showcountries=False,
            showcoastlines=False,
            showframe=False,
            lonaxis=dict(showgrid=True, gridwidth=0.5),
            lataxis=dict(showgrid=True, gridwidth=0.5),
        ),
    )

    fig.show(config={"scrollZoom": True, "displayModeBar": True})
    return fig