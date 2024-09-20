import geopandas as gpd

import ee
from dotenv import load_dotenv

import wxee
import os


load_dotenv()
ee.Authenticate()







def load_file_apply_buffer_square(data_path, geometry="geometry"):
    
    
    if data_path is not None:
        data = gpd.read_file(data_path)
        center = data[geometry]
        buffers = []
        for centroid in center:
        
            coords = list(centroid.exterior.coords)
            if coords[0] != coords[-1]:
                coords.append(coords[0])  # Ensure the polygon is closed

            buffers.append([coords])  #

        return buffers
    else:
        return "Error: File path is not provided."


def apply_sentinel1(polygon, start_date, end_date):
    
    def create_default_image(roi, band_names):
        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image

    roi = ee.Geometry.Polygon(polygon)

    def process_month(month, year):
        start_date = ee.Date.fromYMD(year, month, 1)
        end_date = start_date.advance(1, 'month')

        sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD') \
            .filterBounds(roi) \
            .filter(ee.Filter.date(start_date, end_date)) \
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
            .filter(ee.Filter.eq('instrumentMode', 'IW')) \
            .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))

        composite = ee.Algorithms.If(
            sentinel1.size().gt(0),
            sentinel1.median().set('system:time_start', start_date.millis()),
            create_default_image(roi, ['VV', 'VH']).set('system:time_start', start_date.millis())
        )

        return ee.Image(composite)

    # Générer une liste de mois à partir de start_date et end_date
    start_date_obj = ee.Date(start_date)
    end_date_obj = ee.Date(end_date)
    

    months = ee.List.sequence(0, end_date_obj.difference(start_date_obj, 'month').subtract(1))

    
    monthly_composites = months.map(process_month)

    image_collection = ee.ImageCollection(monthly_composites)

    
    final_composite = image_collection.median()

    print(final_composite)
    return 

   
 
   
def calculate_indices(image):
    # Compute various indices
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    gndvi = image.normalizedDifference(['B8', 'B3']).rename('GNDVI')
    evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
        {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
            'BLUE': image.select('B2')
        }
    ).rename('EVI')
    evi2 = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 2.4 * RED + 1))',
        {
            'NIR': image.select('B8'),
            'RED': image.select('B4')
        }
    ).rename('EVI2')
    arvi = image.expression(
        '(NIR - (2 * RED - BLUE)) / (NIR + (2 * RED - BLUE))',
        {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
            'BLUE': image.select('B2')
        }
    ).rename('ARVI')
    ndre = image.normalizedDifference(['B8A', 'B4']).rename('NDRE')
    ndmi = image.normalizedDifference(['B8', 'B11']).rename('NDMI')
    msavi = image.expression(
        '(2 * NIR + 1 - sqrt((2 * NIR + 1) ** 2 - 8 * (NIR - RED))) / 2',
        {
            'NIR': image.select('B8'),
            'RED': image.select('B4')
        }
    ).rename('MSAVI')
    mndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI')
    ndwi = image.normalizedDifference(['B3', 'B11']).rename('NDWI')
    ndbi = image.normalizedDifference(['B11', 'B8']).rename('NDBI')
    sr = image.expression(
        'NIR / RED',
        {
            'NIR': image.select('B8'),
            'RED': image.select('B4')
        }
    ).rename('SR')
    bsi = image.expression(
        '(SWIR1 + RED - (NIR + BLUE)) / (SWIR1 + RED + NIR + BLUE)',
        {
            'SWIR1': image.select('B11'),
            'RED': image.select('B4'),
            'NIR': image.select('B8'),
            'BLUE': image.select('B2')
        }
    ).rename('BSI')
    nbwi = image.expression(
        '(NIR - SWIR) / (NIR + SWIR)',
        {
            'NIR': image.select('B8'),
            'SWIR': image.select('B11')
        }
    ).rename('NBWI')
    wetness = image.expression(
        '0.1509 * NIR + 0.1973 * SWIR1 + 0.3279 * SWIR2',
        {
            'NIR': image.select('B8'),
            'SWIR1': image.select('B11'),
            'SWIR2': image.select('B12')
        }
    ).rename('Wetness')
    brightness = image.expression(
        '0.2043 * BLUE + 0.4158 * GREEN + 0.5524 * RED',
        {
            'BLUE': image.select('B2'),
            'GREEN': image.select('B3'),
            'RED': image.select('B4')
        }
    ).rename('Brightness')
    greenness = image.expression(
        'NIR - (RED + GREEN)',
        {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
            'GREEN': image.select('B3')
        }
    ).rename('Greenness')
    
    # Add all indices to the image
    return image.addBands([
        ndvi, gndvi, evi, evi2, arvi, ndre, ndmi, msavi,
        mndwi, ndwi, ndbi, sr, bsi, nbwi, wetness, brightness, greenness
    ])

def apply_sentinel2(polygon, start_date, end_date):
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

    def process_month(month, year):
        start_date = ee.Date.fromYMD(year, month, 1)
        end_date = start_date.advance(1, 'month')
        dataset = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
            .filterBounds(roi) \
             .filter(ee.Filter.date(start_date, end_date)) \
            .filter(ee.Filter.calendarRange(month, month, 'month')) \
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)) \
            .map(mask_s2_clouds)\
            

        composite = ee.Algorithms.If(
            dataset.size().gt(0),
            dataset.median().set('system:time_start', ee.Date.fromYMD(year, month, 1).millis()),
            create_default_image(roi, ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12',
                                       'AOT', 'WVP', 'SCL', 'TCI_R', 'TCI_G', 'TCI_B', 'MSK_CLDPRB', 'MSK_SNWPRB',
                                       'QA10', 'QA20', 'QA60']).set('system:time_start',
                                                                    ee.Date.fromYMD(year, month, 1).millis())
        )

        return ee.Image(composite)

    start_date_obj = ee.Date(start_date)
    end_date_obj = ee.Date(end_date)
    

    months = ee.List.sequence(0, end_date_obj.difference(start_date_obj, 'month').subtract(1))

    
    monthly_composites = months.map(process_month)

    image_collection = ee.ImageCollection(monthly_composites)

    
    final_composite = image_collection.median()

    print(final_composite.bandNames())

    

    return 

if __name__ == "__main__":
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))
    
    POLYGONS_COORS = load_file_apply_buffer_square("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Processing_Collect_CI_culture/data/interne/data_culture_CI_exploded.geojson")[1]
    
    #print(POLYGONS_COORS)
    
    apply_sentinel2(polygon=POLYGONS_COORS,start_date="2023-05-01", end_date="2024-05-01")