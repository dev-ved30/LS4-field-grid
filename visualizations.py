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
    print(field_grid.head())
    vertices = SkyCoord(ra=field_grid["ra_deg"].to_numpy()*u.deg, dec=field_grid["dec_deg"].to_numpy()*u.deg, frame="icrs")
    names = field_grid["Field Name"]

    fig = go.Figure()

    # plot all field vertices
    # fig.add_trace(go.Scattergeo(
    #     lon=vertices.ra.wrap_at(180 * u.deg).deg,
    #     lat=vertices.dec.deg,
    #     mode="markers",
    #     marker=dict(size=4, color="lightgray"),
    #     name="Grid vertices"
    # ))

    # sort the schedule by start time
    plan.observing_blocks.sort(key=lambda block: block.start_time)

    # plot scheduled fields by using the rectangular footprint and make a slider that adds and removes them by time
    footprints = footprint(region, vertices)
    for block in plan.observing_blocks:
        target_name = block.target.name
        target_idx = names[names == target_name].index[0]
        region = footprints[target_idx]
        lon, lat = footprint_outline(region)
        if lon is None:
            continue
        fig.add_trace(go.Scattergeo(
            lon=lon,
            lat=lat,
            mode="lines",
            line=dict(width=1.5, color="red"),
            name="Scheduled FoV" if block == plan.observing_blocks[0] else None,
            showlegend=(block == plan.observing_blocks[0])
        ))


    # Add a slider to show the schedule over time
    steps = []
    for i in range(len(plan.observing_blocks)):
        step = dict(
            method="update",
            args=[{"visible": [True] + [j <= i for j in range(len(plan.observing_blocks))]}],
            label=f"Block {i+1}: {plan.observing_blocks[i].target.name}"
        )
        steps.append(step)  
    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Schedule: "},
        pad={"t": 50},
        steps=steps
    )]
    fig.update_layout(sliders=sliders)


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
            lataxis=dict(showgrid=True, gridwidth=0.5)
        )
    )
    fig.show()


    