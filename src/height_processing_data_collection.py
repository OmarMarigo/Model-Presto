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


import ee


def add_indices(image, srtm, roi):
    # --- NDVI (Normalized Difference Vegetation Index) ---
    NDVI = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')

    # --- Elevation ,Slope and Aspect ---
    elevation = srtm.clip(roi).select('elevation').rename('elevation')

    slope = ee.Terrain.slope(srtm).rename('slope')

    aspect =ee.Terrain.aspect( srtm)

    # Adding all indices as bands to the image
    return image.addBands([NDVI,elevation, slope, aspect])




def add_nicfi(image, aoi, year, month):
    # Charger les mosaïques Planet NICFI pour la période donnée
    nicfi = ee.ImageCollection('projects/planet-nicfi/assets/basemaps/forest/planet_monthly_sr') \
        .filterBounds(aoi) \
        .filter(ee.Filter.calendarRange(year, year, 'year')) \
        .filter(ee.Filter.calendarRange(month, month, 'month')) \
        .mosaic()  # Créer une mosaïque pour le mois

    # Sélectionner les bandes RGB (Rouge, Vert, Bleu) de Planet NICFI
    nicfi_rgb = nicfi.select(['R', 'G', 'B']).clip(aoi).rename(['nicfi_red', 'nicfi_green', 'nicfi_blue'])

    # Ajouter les bandes Planet NICFI à l'image
    return image.addBands(nicfi_rgb)
# Cette fonction permet de renvoyer un dataframe propre sans les valeurs nulles dans les bandes rh98 ou rh100
def apply_landsat8(polygon, year, start_month, end_month):
    # Applique  scaling factors.
    def applyScaleFactors(image):

        opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
        thermalBands = image.select("ST_B.*").multiply(0.00341802).add(149.0)

        return image.addBands(opticalBands , overwrite=True)\
                    .addBands(thermalBands, overwrite=True)

    def create_default_image(roi, band_names):
        # Create a constant image with -1 value for all bands when no data is available
        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image

    srtm = ee.Image("USGS/SRTMGL1_003")
    roi = ee.Geometry.Polygon(polygon)

    def process_month(month):
        dataset = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
            .filterBounds(roi) \
            .filter(ee.Filter.calendarRange(year, year, 'year')) \
            .filter(ee.Filter.calendarRange(month, month, 'month')) \
            .filter(ee.Filter.lt('CLOUD_COVER', 20)) \
            .map(applyScaleFactors)

        composite = ee.Algorithms.If(
            dataset.size().gt(0),
            dataset.median().set('system:time_start', ee.Date.fromYMD(year, month, 1).millis()),
            create_default_image(roi, ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7","SR_QA_AEROSOL",
                                        "ST_B10", "ST_ATRAN", "ST_CDIST", "ST_DRAD", "ST_EMIS", "ST_EMSD", "ST_QA",
                                        "ST_TRAD", "ST_URAD", "QA_PIXEL", "QA_RADSAT"]).set('system:time_start',
                                                                    ee.Date.fromYMD(year, month, 1).millis()))

        return add_indices(ee.Image(composite), srtm, roi)

    months = ee.List.sequence(start_month, end_month)
    image_collection = ee.ImageCollection(months.map(process_month))

    # Convert the image collection to an xarray Dataset
    da = image_collection.wx.to_xarray(region=roi, scale=10)

    return da


def checking_dataframe(polygons, year, start_month=1, end_month=12):

    sentinel_1 = apply_sentinel1(polygons, year=year, start_month=start_month, end_month=end_month)

    sentinel_2= apply_sentinel2(polygons, year=year, start_month=start_month, end_month=end_month)

    landsat8_and_ncfi = apply_landsat8(polygons, year=year, start_month=start_month, end_month=end_month)

    all_contact = xarray.concat([sentinel_1, sentinel_2,landsat8_and_ncfi ], dim='source')

    gedi = apply_gedi(polygon=polygons, year=year)

    print(gedi)

    return all_contact


def process_multiple_polygons(polygons, year, output_dir="output"):
    # Créer un répertoire de sortie s'il n'existe pas encore
    os.makedirs(output_dir, exist_ok=True)

    for i, polygon in enumerate(polygons):
        try:
            # Vérifier si le format du polygone est valide
            if isinstance(polygon, list) and all(isinstance(coord, list) for coord in polygon[0]):
                print(f"Processing polygon {i + 1}/{len(polygons)}: {polygon}")

                # Appliquer les fonctions pour récupérer les données satellite pour l'année donnée
                da_s1 = apply_sentinel1(polygon, year, start_month=1, end_month=12)
                da_s2 = apply_sentinel2(polygon, year, start_month=1, end_month=12)
                da_gedi = apply_gedi(polygon, year)

                # Concaténer les datasets Xarray (Sentinel-1, Sentinel-2, GEDI) le long d'une nouvelle dimension 'source'
                all_data = xarray.concat([da_s1, da_s2, da_gedi], dim='source')

                # Définir le nom du fichier de sortie pour chaque polygone et année
                output_path = os.path.join(output_dir, f"polygon_{i + 1}_year_{year}.nc")

                # Sauvegarder l'objet Xarray en fichier NetCDF
                all_data.to_netcdf(output_path)

                print(f"Saved data for polygon {i + 1} year {year} to {output_path}")
            else:
                print(f"Invalid polygon format at index {i + 1}. It should be a list of coordinate pairs.")

        except ee.EEException as e:
            print(f"Error processing polygon {i + 1}: {e}")

    print("Processing completed.")

if __name__ == "__main__":
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))


    POLYGON_COORDS=load_file_apply_buffer_square("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CI/grid_ivory_coast.geojson")[:4]
    print(process_multiple_polygons(polygons=POLYGON_COORDS, year = 2023))