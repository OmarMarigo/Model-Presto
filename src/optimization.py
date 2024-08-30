import ee
from dotenv import load_dotenv
import pandas as pd
import wxee
import numpy as np
import os 



load_dotenv()

ee.Authenticate()


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
   
    return da
    
    
def apply_sentinel1(polygon, year):
    
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
            create_default_image(roi, ['VV', 'VH', 'angle']).set('system:time_start', ee.Date.fromYMD(year, month, 1).millis())
        )

        return ee.Image(composite)

    # Appliquer la fonction de traitement pour chaque mois de l'année spécifiée
    months = ee.List.sequence(1, 12)
    image_collection = ee.ImageCollection(months.map(process_month))
    da = image_collection.wx.to_xarray(region=roi, scale=10)
   
    
    return da
    
   


def apply_sentinel2(polygon,year):
    
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
            create_default_image(roi,['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12','AOT','WVP','SCL','TCI_R','TCI_G','TCI_B','MSK_CLDPRB','MSK_SNWPRB','QA10','QA20','QA60']).set('system:time_start', ee.Date.fromYMD(year, month, 1).millis())
        )

        return ee.Image(composite)
    
    months = ee.List.sequence(1, 12)
    image_collection = ee.ImageCollection(months.map(process_month))
    da = image_collection.wx.to_xarray(region=roi, scale=10)
    
    return da


def checking_dataframe(polygons, year):
    
    # Chargement des données GEDI pour l'année entière
    gedi = apply_gedi(polygons, year=year)

 
    non_null_months = []
    for month in range(1, 13):  
        gedi_month = gedi.sel(time=gedi['time'].dt.month == month)
        
        if not np.all(np.isnan(gedi_month.data_vars["rh100"].values)):
            non_null_months.append(month)

    if not non_null_months:
        print("No month has non-zero GEDI data for the year specified.")
        return None
    else:
        print(f"Data available for the months: {non_null_months}.")
        # Sélectionner les données pour les mois où rh100 n'est pas entièrement nul
        gedi_df = gedi.sel(time=gedi['time'].dt.month.isin(non_null_months))

        s1_month = apply_sentinel1(polygons, year=year).sel(time=gedi['time'].dt.month.isin(non_null_months))
        s2_month = apply_sentinel2(polygons, year=year).sel(time=gedi['time'].dt.month.isin(non_null_months))
        dataset = pd.concat([
                gedi_df.to_dataframe().reset_index(),
                s1_month.to_dataframe().reset_index(),
                s2_month.to_dataframe().reset_index()
            ], axis=1)
        dataset = dataset.dropna(subset="rh98")
        print(dataset.shape)
        return dataset

if __name__ == "__main__":
    # Initialisation du projet Earth Engine
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))

    # Définition des coordonnées du polygone
    POLYGON_COORDS = [
        [43.15675751576177, 14.736616726173434],
        [43.156756261784274, 14.713479886543574],
        [43.1329848450281, 14.713479886947383],
        [43.132983591880425, 14.736616726577278],
        [43.15675751576177, 14.736616726173434]
    ]

    checking_dataframe(POLYGON_COORDS, year=2020)
    #apply_gedi(POLYGON_COORDS)
   
