#####################################################
#####################################################

#####################################################
#####################################################

import geopandas as gpd
from ee import ImageCollection
from pyproj import Proj, Transformer
from shapely.geometry import Point
from shapely.ops import transform
from functools import partial
import json
import ee
from rasterio.transform import xy
from rasterio.transform import rowcol
from dotenv import load_dotenv
import pandas as pd
import wxee
import xarray
import numpy as np
from datetime import datetime
import os
import rasterio
from rasterio.transform import from_origin

load_dotenv()

ee.Authenticate()


def get_projections(anchor_point=None):
    if not anchor_point:  # Si aucun point d'ancrage, définir 0/0 (Mercator normal).
        anchor_point = Point(0, 0)
    crs = {
        'proj': 'omerc',
        'lat_0': anchor_point.y,
        'lonc': anchor_point.x,
        'alpha': 1e-6,
        'k': 1,
        'gamma': 0.0,
        'units': 'm',
        'ellps': 'WGS84',
        'no_defs': True,
    }
    proj = Proj(crs, preserve_units=True)
    proj_inverse = Transformer.from_proj(proj, 'EPSG:4326', always_xy=True).transform

    return partial(transform, proj), partial(transform, proj_inverse)





def load_file_apply_buffer_square(data_path, geometry="geometry", side=2560):
    if data_path is not None:
        data = gpd.read_file(data_path)


        center = data[geometry]

        buffers = []

        for centroid in center:
            proj, proj_inverse = get_projections(centroid)
            cartesian_center = proj(centroid)
            cartesian_square = cartesian_center.buffer(side / 2, cap_style=3)

            geodesic_square = proj_inverse(cartesian_square)

            # Convert to a list of lists, ensuring the polygon is closed
            coords = list(geodesic_square.exterior.coords)
            if coords[0] != coords[-1]:
                coords.append(coords[0])  # Ensure the polygon is closed

            buffers.append([coords])  # Ensure it's a list of lists of coordinates

        return buffers
    else:
        return "Error: File path is not provided."



if __name__ == "__main__":
    # with rasterio.open("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/gedi_output.nc/201904_042E_012N.tif") as src:
    # print(src.read().shape)
    # exit()
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))

    POLYGON_COORDS = load_file_apply_buffer_square(
        "/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/GeoJson/df_african.geojson")
    OUT_PATH_SENTINEL_1 = os.getenv("OUT_PATH_SENTINEL_1")
    OUT_PATH_SENTINEL_2 = os.getenv("OUT_PATH_SENTINEL_2")
    OUT_PATH_FUSED = os.getenv("OUT_PATH_FUSED")

    # Récupérer les données GEDI, Sentinel-1 et Sentinel-2
    get_gedi(POLYGON_COORDS,
             output_path="/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/gedi_output.nc")
    # get_sentinel1_monthly(POLYGON_COORDS, OUT_PATH_SENTINEL_1)
    # get_sentinel2_monthly(POLYGON_COORDS, OUT_PATH_SENTINEL_2)