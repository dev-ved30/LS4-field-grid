import pandas as pd
import numpy as np
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

    fig = go.Figure()
    scheduled_rows = []
    row_trace_counts = []

    legend_traces = [
        go.Scattergeo(
            lon=[None],
            lat=[None],
            mode="lines",
            line=dict(width=2.5, color=parity_colors["Even"]),
            name="First Pointing",
            visible="legendonly",
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
                visible=False,
                hovertext=hover_text,
                hoverinfo="text",
            ))
            trace_count += 1

        if trace_count > 0:
            scheduled_rows.append(row)
            row_trace_counts.append(trace_count)

    if scheduled_rows:
        steps = []
        cumulative_traces = 0
        for i, row in enumerate(scheduled_rows):
            cumulative_traces += row_trace_counts[i]
            steps.append(dict(
                method="update",
                args=[
                    {"visible": ["legendonly", "legendonly"] + [j < cumulative_traces for j in range(sum(row_trace_counts))]},
                    {"title": f"LS4 Observing Plan through {row['target']}"},
                ],
                label=f"{i + 1}: {row['target']}"
            ))

        for trace_index in range(len(legend_traces), len(legend_traces) + row_trace_counts[0]):
            fig.data[trace_index].visible = True

        fig.update_layout(sliders=[dict(
            active=0,
            currentvalue={"prefix": "Schedule: "},
            pad={"t": 50},
            steps=steps,
        )], legend=dict(title="Block parity"))


    fig.update_layout(
        title="LS4 Observing Plan",
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


    