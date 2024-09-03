import geopandas as gpd

from pyproj import Proj, Transformer
from shapely.geometry import Polygon
from shapely.ops import transform
from functools import partial
import ee
from dotenv import load_dotenv

import wxee
import os


load_dotenv()
ee.Authenticate()


def get_projections(anchor_point=None):
    if not anchor_point:  # Si aucun point d'ancrage, définir 0/0 (Mercator normal).
        anchor_point =  Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        
    centroid = anchor_point.centroid
    crs = {
        'proj': 'omerc',
        'lat_0': centroid.y,
        'lonc': centroid.x,
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

            buffers.append([coords])  #

        return buffers
    else:
        return "Error: File path is not provided."


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
            .filter(ee.Filter.eq('instrumentMode', 'IW'))\
            .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))



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


if __name__ == "__main__":
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))
    
    POLYGONS_COORS = load_file_apply_buffer_square("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Processing_Collect_CI_culture/data/interne/data_culture_CI_exploded.geojson")[1]
    
    #print(POLYGONS_COORS)
    
    print(apply_sentinel1(polygon=POLYGONS_COORS, year=2024, start_month=1, end_month=12).dims)