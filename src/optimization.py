import geopandas as gpd
from ee import ImageCollection
from pyproj import Proj, Transformer
from shapely.geometry import Point
from shapely.ops import transform
from functools import partial
import json
import ee
from rasterio.transform import xy, rowcol, from_origin
from rasterio.enums import Resampling
from dotenv import load_dotenv
import pandas as pd
import rasterio
import shapely
import geemap
import wxee
import xarray
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt

from common import load_file_apply_buffer_square

load_dotenv()

# Initialiser l'authentification Earth Engine
ee.Authenticate()


def get_projections(anchor_point=None):
    if not anchor_point:  # Si aucun point d'ancrage, définir 0/0 (Mercator normal)
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


def apply_gedi(polygon, year):
    def quality_mask(image):
        quality_mask = image.select('quality_flag').eq(1)
        degrade_mask = image.select('degrade_flag').eq(0)
        return image.updateMask(quality_mask).updateMask(degrade_mask)

    roi = ee.Geometry.Polygon(polygon)

    dataset_collection = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY').filterBounds(roi) \
        .filter(ee.Filter.calendarRange(year, year, 'year')) \
        .map(quality_mask) \
        .select(['rh98', 'rh100'])

    da = dataset_collection.wx.to_xarray(region=roi, scale=10)
    print(da.dims)

    return da


def apply_sentinel1(polygon, year, start_month, end_month):
    def create_default_image(roi, band_names):
        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image

    roi = ee.Geometry.Polygon(polygon)

    def process_month(month):
        sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD') \
            .filterBounds(roi) \
            .filter(ee.Filter.calendarRange(year, year, 'year')) \
            .filter(ee.Filter.calendarRange(month, month, 'month')) \
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
            .filter(ee.Filter.eq('instrumentMode', 'IW'))

        # Utiliser une condition côté serveur pour vérifier si la collection est vide
        composite = ee.Algorithms.If(
            sentinel1.size().gt(0),
            sentinel1.median().set('system:time_start', ee.Date.fromYMD(year, month, 1).millis()),
            create_default_image(roi, ['VV', 'VH', 'angle']).set('system:time_start',
                                                                 ee.Date.fromYMD(year, month, 1).millis())
        )

        return ee.Image(composite)

    # Appliquer la fonction de traitement pour chaque mois de l'année spécifiée
    months = ee.List.sequence(start_month, end_month)
    image_collection = ee.ImageCollection(months.map(process_month))
    da = image_collection.wx.to_xarray(region=roi, scale=10)

    return da


def apply_sentinel2(polygon, year, start_month, end_month):
    def mask_s2_clouds(image):
        qa = image.select('QA60')
        cloud_bit_mask = 1 << 10
        cirrus_bit_mask = 1 << 11

        # Créer un masque pour les nuages et les cirrus
        mask = (
            qa.bitwiseAnd(cloud_bit_mask)
            .eq(0)
            .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0))
        )

        return image.updateMask(mask).divide(10000)

    def create_default_image(roi, band_names):
        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image

    roi = ee.Geometry.Polygon(polygon)

    def process_month(month):
        dataset = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
            .filterBounds(roi) \
            .filter(ee.Filter.calendarRange(year, year, 'year')) \
            .filter(ee.Filter.calendarRange(month, month, 'month')) \
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)) \
            .map(mask_s2_clouds)

        composite = ee.Algorithms.If(
            dataset.size().gt(0),
            dataset.median().set('system:time_start', ee.Date.fromYMD(year, month, 1).millis()),
            create_default_image(roi, ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12',
                                       'AOT', 'WVP', 'SCL', 'TCI_R', 'TCI_G', 'TCI_B', 'MSK_CLDPRB', 'MSK_SNWPRB',
                                       'QA10', 'QA20', 'QA60']).set('system:time_start',
                                                                    ee.Date.fromYMD(year, month, 1).millis())
        )

        return ee.Image(composite)

    months = ee.List.sequence(start_month, end_month)
    image_collection = ee.ImageCollection(months.map(process_month))
    da = image_collection.wx.to_xarray(region=roi, scale=10)

    return da


# Cette fonction permet de renvoyer un dataframe propre sans les valeurs nulles dans les bandes rh98 ou rh100
def checking_dataframe(polygons, year, start_month=4, end_month=12):
    gedi = apply_gedi(polygons, year=year)

    non_null_months = []
    null_months = []

    for month in range(start_month, end_month + 1):
        gedi_month = gedi.sel(time=gedi['time'].dt.month == month)

        if not np.all(np.isnan(gedi_month.data_vars["rh100"].values)):
            non_null_months.append(month)
        else:
            null_months.append(month)

    if not non_null_months:
        print("Aucun mois n'a de données GEDI non nulles pour l'année spécifiée.")
        gedi_df = gedi.sel(time=gedi['time'].dt.month.isin(null_months))

        s1_month = apply_sentinel1(polygons, year=year, start_month=start_month, end_month=end_month).sel(
            time=gedi['time'].dt.month.isin(null_months))
        s2_month = apply_sentinel2(polygons, year=year, start_month=start_month, end_month=end_month).sel(
            time=gedi['time'].dt.month.isin(null_months))

        dataset = pd.concat([
            gedi_df.to_dataframe().reset_index(),
            s1_month.to_dataframe().reset_index(),
            s2_month.to_dataframe().reset_index()
        ], axis=1)
        dataset = dataset.dropna(subset=["rh100"])
        print(dataset.shape)
        return dataset
    else:
        print(f"Données disponibles pour les mois : {non_null_months}.")
        gedi_df = gedi.sel(time=gedi['time'].dt.month.isin(non_null_months))

        s1_month = apply_sentinel1(polygons, year=year, start_month=start_month, end_month=end_month).sel(
            time=gedi['time'].dt.month.isin(non_null_months))
        s2_month = apply_sentinel2(polygons, year=year, start_month=start_month, end_month=end_month).sel(
            time=gedi['time'].dt.month.isin(non_null_months))

        dataset = pd.concat([
            gedi_df.to_dataframe().reset_index(),
            s1_month.to_dataframe().reset_index(),
            s2_month.to_dataframe().reset_index()
        ], axis=1)
        dataset = dataset.dropna(subset=["rh98"])
        print(gedi_df.dims)
        print(s2_month.dims)
        return dataset


def process_multiple_polygons(polygons, year):
    results = {}
    all_data_frames = []
    for i, polygon in enumerate(polygons):
        print(f"Processing polygon {i + 1}/{len(polygons)}...")
        df = checking_dataframe(polygons=polygon, year=year)
        i = 305
        df['polygon'] = f'polygon_{i + 1}'
        all_data_frames.append(df)

    # Combiner tous les DataFrames en un seul
        combined_df = pd.concat(all_data_frames, ignore_index=True)
        combined_df.to_csv("Combined_Data.csv", index=False)

    results.update({f'polygon_{i + 1}': df})
    return results


if __name__ == "__main__":
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))

    POLYGON_COORDS = load_file_apply_buffer_square(os.getenv("DATASET_PATH_AFRICA"))[305:]
    process_multiple_polygons(polygons=POLYGON_COORDS, year=2019)