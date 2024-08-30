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


import geopandas as gpd
import shapely.geometry
import json
from shapely.geometry import mapping


def load_file_apply_buffer_square(data_path, geometry="geometry", side=2560):
    if data_path is not None:
        data = gpd.read_file(data_path)

        data_africa = data[data['Continent_Code'] == 3]
        data_africa = data_africa[data_africa["End_Year"] >= 2019]

        center = data_africa[geometry].centroid

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


def get_gedi(polygons, output_path: str, start_date="2019-01-01", end_date="2019-12-31"):
    def quality_mask(image):
        quality_mask = image.select('quality_flag').eq(1)
        degrade_mask = image.select('degrade_flag').eq(0)
        return image.updateMask(quality_mask).updateMask(degrade_mask)

    def set_time_start(image):
        return image.set('system:time_start', ee.Date(image.get('system:time_start')).millis())

    all_datasets = []

    for i, polygon_coords in enumerate(polygons) :
        # Ensure the polygon coordinates are in the correct format
        try:
            roi = ee.Geometry.Polygon(polygon_coords)
        except Exception as e:
            print(f"Error creating polygon: {e}")
            continue

        dataset_collection = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY').filterBounds(roi) \
            .filterDate(start_date, end_date) \
            .map(quality_mask) \
            .select(['rh98', 'rh100']) \
            .map(set_time_start)

        da = dataset_collection.wx.to_xarray(region=roi, scale=10)
        df = da.to_dataframe().reset_index()

        if 'time' in df.columns:
            df['date'] = pd.to_datetime(df['time']).dt.date  # Extract only the date
            df = df.drop(columns=['time'])
            df = df.rename(columns={"date": "time"})


        all_datasets.append(df)

        print(f"Processing polygon {i + 1}/{len(polygons)}")
    final_df = pd.concat(all_datasets, ignore_index=True)
    final_df.to_csv("gedi_2019.csv", index=False)


def get_sentinel1_monthly(polygons, output_path, year=2019):
    def create_default_image(roi, band_names):

        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image


    all_datasets = []
    for i, polygon_coords in enumerate(polygons):
        # Ensure the polygon coordinates are in the correct format
        try:
            roi = ee.Geometry.Polygon(polygon_coords)
        except Exception as e:
            print(f"Error creating polygon {i + 1}: {e}")
            continue

        print(f"Processing polygon {i + 1}/{len(polygons)}")

        images = []
        for month in range(4, 13):
            # Define the date range for the month
            start_date = datetime(year, month, 1)
            if month == 12:
                end_date = datetime(year + 1, 1, 1)
            else:
                end_date = datetime(year, month + 1, 1)

            start_date_str = start_date.strftime('%Y-%m-%d')
            end_date_str = end_date.strftime('%Y-%m-%d')

            # Load the Sentinel-1 ImageCollection and apply filters
            sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD') \
                .filterBounds(roi) \
                .filterDate(start_date_str, end_date_str) \
                .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
                .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
                .filter(ee.Filter.eq('instrumentMode', 'IW'))

            if sentinel1.size().getInfo() > 0:

                composite = sentinel1.median().set('system:time_start', start_date.timestamp() * 1000)

            else:
                print(f"No valid images found for month {month} in polygon {i + 1}:")

                band_names =['VV','VH','angle']
                default_image = create_default_image(roi, band_names)
                composite = default_image.set('system:time_start', start_date.timestamp() * 1000)
            images.append(composite)
            # Convert images to xarray and append to dataset list
        imageCollection = ee.ImageCollection.fromImages(images)
        da = imageCollection.wx.to_xarray(region=roi, scale=10)
        da = da.to_dataframe().reset_index()
        all_datasets.append(da)
    final_df = pd.concat(all_datasets, ignore_index=True)
    final_df.to_csv(f"sentinel_1_{year}.csv", index=False)


def get_sentinel2_monthly(polygons, output_path, year=2019):
    def mask_s2_clouds(image):
        qa = image.select('QA60')
        cloud_bit_mask = 1 << 10
        cirrus_bit_mask = 1 << 11

        # Create a mask for clouds and cirrus
        mask = (
            qa.bitwiseAnd(cloud_bit_mask)
            .eq(0)
            .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0))
        )

        # Apply the mask and scale the image
        return image.updateMask(mask).divide(10000)

    def create_default_image(roi, band_names):

        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image
    all_datasets = []
    for i, polygon_coords in enumerate(polygons):
        # Ensure the polygon coordinates are in the correct format
        try:
            roi = ee.Geometry.Polygon(polygon_coords)
        except Exception as e:
            print(f"Error creating polygon {i}: {e}")
            continue

        print(f"Processing polygon {i + 1}/{len(polygons)}")

        images = []
        for month in range(4, 13):
            # Compute start and end dates for the given month and year
            start_date = datetime(year, month, 1)
            if month == 12:
                end_date = datetime(year + 1, 1, 1)
            else:
                end_date = datetime(year, month + 1, 1)

            start_date_str = start_date.strftime('%Y-%m-%d')
            end_date_str = end_date.strftime('%Y-%m-%d')

            # Load the Sentinel-2 ImageCollection and apply filters
            dataset = ee.ImageCollection('COPERNICUS/S2_SR') \
                .filterBounds(roi) \
                .filterDate(start_date_str, end_date_str) \
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)) \
                .map(mask_s2_clouds)

            if dataset.size().getInfo() > 0:

                composite = dataset.median().set('system:time_start', start_date.timestamp() * 1000)

            else:
                print(f"No valid images found for month {month} in polygon {i + 1}:")

                band_names =['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12','AOT','WVP','SCL','TCI_R','TCI_G','TCI_B','MSK_CLDPRB','MSK_SNWPRB','QA10','QA20','QA60']

                default_image = create_default_image(roi, band_names)
                composite = default_image.set('system:time_start', start_date.timestamp() * 1000)
            images.append(composite)

        imageCollection = ee.ImageCollection.fromImages(images)
        da = imageCollection.wx.to_xarray(region=roi, scale=10)
        da = da.to_dataframe().reset_index()
        all_datasets.append(da)

    final_df = pd.concat(all_datasets, ignore_index=True)
    final_df.to_csv(f"sentinel_2_{year}.csv", index=False)








if __name__ == "__main__":
    #with rasterio.open("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/gedi_output.nc/201904_042E_012N.tif") as src:
        #print(src.read().shape)
    #exit()
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))

    POLYGON_COORDS = load_file_apply_buffer_square("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/GeoJson/df_african.geojson")
    OUT_PATH_SENTINEL_1 = os.getenv("OUT_PATH_SENTINEL_1")
    OUT_PATH_SENTINEL_2 = os.getenv("OUT_PATH_SENTINEL_2")
    OUT_PATH_FUSED = os.getenv("OUT_PATH_FUSED")

    # Récupérer les données GEDI, Sentinel-1 et Sentinel-2
    get_gedi(POLYGON_COORDS, output_path="/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/gedi_output.nc")
    #get_sentinel1_monthly(POLYGON_COORDS, OUT_PATH_SENTINEL_1)
    #get_sentinel2_monthly(POLYGON_COORDS, OUT_PATH_SENTINEL_2)






