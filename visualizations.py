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

def plot_field_grid(vertices, footprints):
    fig = go.Figure()

    fig.add_trace(go.Scattergeo(
        lon=vertices.ra.wrap_at(180 * u.deg).deg,
        lat=vertices.dec.deg,
        mode="markers",
        marker=dict(size=4, color="royalblue"),
        name="Grid vertices"
    ))

    for i, region in enumerate(footprints):
        lon, lat = footprint_outline(region)
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

def plot_obs_plan(plan, region):

    field_grid = pd.read_csv("LS4_field_grid.csv")
    vertices = SkyCoord(
        ra=field_grid["ra_deg"].to_numpy() * u.deg,
        dec=field_grid["dec_deg"].to_numpy() * u.deg,
        frame="icrs",
    )
    names = field_grid["Field Name"].astype(str)

    plan = plan.copy()
    if "target" not in plan.columns or "start time (UTC)" not in plan.columns:
        raise ValueError("plan must include 'target' and 'start time (UTC)' columns")

    block_column = None
    if "block_numer" in plan.columns:
        block_column = "block_numer"
    elif "block_number" in plan.columns:
        block_column = "block_number"
    else:
        raise ValueError("plan must include 'block_numer' or 'block_number' for coloring")

    plan = plan[plan["target"].astype(str) != "TransitionBlock"].copy()
    plan["start time (UTC)"] = pd.to_datetime(plan["start time (UTC)"])
    if "end time (UTC)" in plan.columns:
        plan["end time (UTC)"] = pd.to_datetime(plan["end time (UTC)"])
    plan.sort_values("start time (UTC)", inplace=True)
    plan.reset_index(drop=True, inplace=True)

    parity_colors = {
        "Even": "rgba(65,105,225,1.0)",
        "Odd": "rgba(220,20,60,1.0)",
    }

    footprints = footprint(region, vertices)
    fig = go.Figure()
    scheduled_rows = []

    legend_traces = [
        go.Scattergeo(
            lon=[None],
            lat=[None],
            mode="lines",
            line=dict(width=2.5, color=parity_colors["Even"]),
            name="Even blocks",
            visible="legendonly",
            hoverinfo="skip",
            showlegend=True,
        ),
        go.Scattergeo(
            lon=[None],
            lat=[None],
            mode="lines",
            line=dict(width=2.5, color=parity_colors["Odd"]),
            name="Odd blocks",
            visible="legendonly",
            hoverinfo="skip",
            showlegend=True,
        ),
    ]

    for trace in legend_traces:
        fig.add_trace(trace)

    for _, row in plan.iterrows():
        target_name = str(row["target"])
        target_matches = names[names == target_name]
        if target_matches.empty:
            continue

        target_idx = target_matches.index[0]
        footprint_region = footprints[target_idx]
        lon, lat = footprint_outline(footprint_region)
        if lon is None:
            continue

        start_time = row["start time (UTC)"]
        end_time = row["end time (UTC)"] if "end time (UTC)" in row and pd.notna(row["end time (UTC)"]) else None
        block_value = row[block_column]
        block_label = int(block_value) if pd.notna(block_value) and float(block_value).is_integer() else block_value
        block_parity = "Odd" if int(block_value) % 2 else "Even"
        block_color = parity_colors[block_parity]
        hover_text = f"Target: {target_name}<br>Start: {start_time}"
        if end_time is not None:
            hover_text += f"<br>End: {end_time}"
        hover_text += f"<br>Block: {block_label} ({block_parity})"

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
        scheduled_rows.append(row)

    if scheduled_rows:
        steps = []
        for i, row in enumerate(scheduled_rows):
            steps.append(dict(
                method="update",
                args=[
                    {"visible": ["legendonly", "legendonly"] + [j <= i for j in range(len(scheduled_rows))]},
                    {"title": f"LS4 Observing Plan through {row['target']}"},
                ],
                label=f"{i + 1}: {row['target']}"
            ))

        fig.data[len(legend_traces)].visible = True

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


    