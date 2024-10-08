import geopandas as gpd
from ee import ImageCollection
from pyproj import Proj, Transformer
from shapely.geometry import Point,Polygon
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
import xarray
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt
import geemap
from loguru import logger
from rasterio.features import rasterize
from google.cloud import storage
CLASSES_CODES={'mil': 8,
               'mais': 9,
                 'mil+mais':10, 
                 'arachide': 11,
                   'oseille': 12, 
                   'sorgho': 13,
                     'pasteque': 14,
                       'arachide+niebe': 15,
                         'verger': 16,
                           'niebe': 17,
                             'mil+niebe': 18,
                            'mil+sorgho': 19,
                              'riz': 20,
                            'arachide+mil': 21,
                            'eau':22}
def download_ee_image(
    image,
    filename,
    region=None,
    crs=None,
    crs_transform=None,
    scale=None,
    resampling="near",
    dtype=None,
    overwrite=True,
    num_threads=None,
    max_tile_size=None,
    max_tile_dim=None,
    shape=None,
    scale_offset=False,
    unmask_value=None,
    **kwargs,
):
    """Download an Earth Engine Image as a GeoTIFF. Images larger than the `Earth Engine size limit are split and downloaded as
        separate tiles, then re-assembled into a single GeoTIFF. See https://github.com/dugalh/geedim/blob/main/geedim/download.py#L574

    Args:
        image (ee.Image): The image to be downloaded.
        filename (str): Name of the destination file.
        region (ee.Geometry, optional): Region defined by geojson polygon in WGS84. Defaults to the entire image granule.
        crs (str, optional): Reproject image(s) to this EPSG or WKT CRS.  Where image bands have different CRSs, all are
            re-projected to this CRS. Defaults to the CRS of the minimum scale band.
        crs_transform (list, optional): tuple of float, list of float, rio.Affine, optional
            List of 6 numbers specifying an affine transform in the specified CRS.  In row-major order:
            [xScale, xShearing, xTranslation, yShearing, yScale, yTranslation].  All bands are re-projected to
            this transform.
        scale (float, optional): Resample image(s) to this pixel scale (size) (m).  Where image bands have different scales,
            all are resampled to this scale.  Defaults to the minimum scale of image bands.
        resampling (ResamplingMethod, optional): Resampling method, can be 'near', 'bilinear', 'bicubic', or 'average'. Defaults to None.
        dtype (str, optional): Convert to this data type (`uint8`, `int8`, `uint16`, `int16`, `uint32`, `int32`, `float32`
            or `float64`).  Defaults to auto select a minimum size type that can represent the range of pixel values.
        overwrite (bool, optional): Overwrite the destination file if it exists. Defaults to True.
        num_threads (int, optional): Number of tiles to download concurrently. Defaults to a sensible auto value.
        max_tile_size: int, optional
            Maximum tile size (MB).  If None, defaults to the Earth Engine download size limit (32 MB).
        max_tile_dim: int, optional
            Maximum tile width/height (pixels).  If None, defaults to Earth Engine download limit (10000).
        shape: tuple of int, optional
            (height, width) dimensions to export (pixels).
        scale_offset: bool, optional
            Whether to apply any EE band scales and offsets to the image.
        unmask_value (float, optional): The value to use for pixels that are masked in the input image. If the exported image contains
            zero values, you should set the unmask value to a  non-zero value so that the zero values are not treated as missing data. Defaults to None.

    """
    
    if os.environ.get("USE_MKDOCS") is not None:
        return

    try:
        import geedim as gd
    except ImportError:
        raise ImportError(
            "Please install geedim using `pip install geedim` or `conda install -c conda-forge geedim`"
        )

    if not isinstance(image, ee.Image):
        raise ValueError("image must be an ee.Image.")

    if unmask_value is not None:
        image = image.selfMask().unmask(unmask_value)
        if isinstance(region, ee.Geometry):
            image = image.clip(region)
        elif isinstance(region, ee.FeatureCollection):
            image = image.clipToCollection(region)

    if region is not None:
        kwargs["region"] = region

    if crs is not None:
        kwargs["crs"] = crs

    if crs_transform is not None:
        kwargs["crs_transform"] = crs_transform

    if scale is not None:
        kwargs["scale"] = scale

    if resampling is not None:
        kwargs["resampling"] = resampling

    if dtype is not None:
        kwargs["dtype"] = dtype

    if max_tile_size is not None:
        kwargs["max_tile_size"] = max_tile_size

    if max_tile_dim is not None:
        kwargs["max_tile_dim"] = max_tile_dim

    if shape is not None:
        kwargs["shape"] = shape

    if scale_offset:
        kwargs["scale_offset"] = scale_offset

    img=gd.download.BaseImage(image)

 
    img.download(filename, overwrite=overwrite, num_threads=num_threads, **kwargs)
def centroid_to_square(centroid, side=2560):

    centroid=Point(centroid)
    proj, proj_inverse = get_projections(centroid)
    cartesian_center = proj(centroid)
    cartesian_square = cartesian_center.buffer(side / 2, cap_style=3)

    geodesic_square = proj_inverse(cartesian_square)

    # Convert to a list of coordinates, ensuring the polygon is closed
    coords = list(geodesic_square.exterior.coords)
    if coords[0] != coords[-1]:
        coords.append(coords[0])  # Ensure the polygon is closed

    return coords  # Return the buffer coordinates





def load_file_apply_buffer_square(data_path, geometry="geometry", side=2560):
    if data_path is not None:
        data = gpd.read_file(data_path).head(5)


        center = data[geometry]
        ids=data['id']

        buffers = []

        for i,centroid in enumerate(center):
            proj, proj_inverse = get_projections(centroid)
            cartesian_center = proj(centroid)
            cartesian_square = cartesian_center.buffer(side / 2, cap_style=3)

            geodesic_square = proj_inverse(cartesian_square)

            # Convert to a list of lists, ensuring the polygon is closed
            coords = list(geodesic_square.exterior.coords)
            
            if coords[0] != coords[-1]:
                coords.append(coords[0])  # Ensure the polygon is closed

            buffers.append([ids[i],coords])  # Ensure it's a list of lists of coordinates

        return buffers
    else:
        return "Error: File path is not provided."
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




def apply_gedi_month(polygon, year, month,directory_path, filename,scale=10):
    def quality_mask(image):
        quality_mask = image.select('quality_flag').eq(1)
        degrade_mask = image.select('degrade_flag').eq(0)
        return image.updateMask(quality_mask).updateMask(degrade_mask)

    roi = ee.Geometry.Polygon(polygon)

    dataset_collection = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY').filterBounds(roi) \
        .filter(ee.Filter.calendarRange(year, year, 'year')) \
        .filter(ee.Filter.calendarRange(month, month, 'month')) \
        .map(quality_mask) \
        .select(['rh98'])

    # Get the single image for the specified month
    image = dataset_collection.first()
    download_ee_image(image,os.path.join(directory_path,filename), scale=10, region=roi, crs="EPSG:4326")
    return image

def get_S1_composite(polygon, start_date,end_date,directory_path, filename, scale=10):
    def create_default_image(roi, band_names):
        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image

    roi = ee.Geometry.Polygon(polygon)

    sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD') \
        .filterBounds(roi) \
        .filter(ee.Filter.date(start_date, end_date))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) \
        .filter(ee.Filter.eq('instrumentMode', 'IW')).select(['VV', 'VH'])

    # Use a server-side condition to check if the collection is empty
    image = ee.Algorithms.If(
        sentinel1.size().gt(0),
        sentinel1.median().set('system:time_start', ee.Date(end_date).millis()),
        create_default_image(roi, ['VV', 'VH']).set('system:time_start',
                                                             ee.Date(end_date).millis())
    )
    download_ee_image(ee.Image(image),os.path.join(directory_path,filename), scale=scale, region=roi, crs="EPSG:4326")

    blob_name=os.path.join(bucket_repository,directory_path.split("/")[-1],filename)
    upload_to_gcs(os.path.join(directory_path,filename),blob_name)

    return ee.Image(image)



def get_S2_composite(polygon, start_date,end_date,directory_path, filename,scale=10,cloud_pct=20):
    def create_default_image(roi, band_names):
        default_bands = [ee.Image.constant(-1).rename(band) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image
    roi = ee.Geometry.Polygon(polygon)
    dataset = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
        .filterBounds(roi) \
         .filter(ee.Filter.date(start_date, end_date)) \
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloud_pct))
    
    bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10','B11', 'B12']
    
    def scale_bands(image):
        scaled = image.select(bands[:-1]).divide(10000)
        return image.addBands(scaled, overwrite=True).select(bands)



    composite = ee.Algorithms.If(
        dataset.size().gt(0),
        dataset.median().set('system:time_start', ee.Date(end_date).millis()),
        create_default_image(roi, bands).set('system:time_start', ee.Date(end_date).millis())
    )

    image = ee.Image(composite)
    download_ee_image(image,os.path.join(directory_path,filename), scale=scale, region=roi, crs="EPSG:4326")
    blob_name=os.path.join(bucket_repository,directory_path.split("/")[-1],filename)
    upload_to_gcs(os.path.join(directory_path,filename),blob_name)
    return image







def get_srtm(polygon,directory_path,filename,scale=10):
    # --- NDVI (Normalized Difference Vegetation Index) ---
    roi = ee.Geometry.Polygon(polygon)
    srtm = ee.Image("USGS/SRTMGL1_003").clip(roi)
    # --- Elevation ,Slope and Aspect ---
    elevation = srtm.select('elevation').rename('elevation')

    slope = ee.Terrain.slope(srtm).rename('slope')

    aspect =ee.Terrain.aspect( srtm)
    image=elevation.addBands([slope, aspect])
    download_ee_image(image,os.path.join(directory_path,filename), scale=scale, region=roi, crs="EPSG:4326")
    blob_name=os.path.join(bucket_repository,directory_path.split("/")[-1],filename)
    upload_to_gcs(os.path.join(directory_path,filename),blob_name)
    # Adding all indices as bands to the image
    return image

def get_centroid(polygon):
    roi = Polygon(polygon)

    return roi.centroid




def get_nicfi_composite(polygon, start_date,end_date, directory_path,filename,scale=10):
    roi=ee.Geometry.Polygon(polygon)
    # Load Planet NICFI mosaic for the given month and year
    nicfi = ee.ImageCollection('projects/planet-nicfi/assets/basemaps/africa') \
        .filterBounds(roi) \
        .filter(ee.Filter.date(start_date, end_date)) \
        .median().clip(roi)# Get the first (and only) image for the month


    # Select the RGB and NIR bands from Planet NICFI
    nicfi_bands = nicfi.select(['R', 'G', 'B', 'N'])


    # Scale the NICFI bands (they are typically in the range 0-10000)
    #nicfi_scaled = nicfi_bands.divide(10000)
    download_ee_image(nicfi_bands,os.path.join(directory_path,filename), scale=scale, region=roi, crs="EPSG:4326")
    blob_name=os.path.join(bucket_repository,directory_path.split("/")[-1],filename)
    upload_to_gcs(os.path.join(directory_path,filename),blob_name)
    # Add the scaled NICFI bands to the original image
    return nicfi_bands

# Cette fonction permet de renvoyer un dataframe propre sans les valeurs nulles dans les bandes rh98 ou rh100
def get_landsat8_composite(polygon, start_date,end_date, directory_path, filename, scale=10,cloud_pct=20):
    # Apply scaling factors
    def applyScaleFactors(image):
        opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
        return opticalBands.addBands(image.select('QA_PIXEL'))

    roi = ee.Geometry.Polygon(polygon)

    dataset = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .filterBounds(roi) \
        .filter(ee.Filter.date(start_date, end_date)) \
        .filter(ee.Filter.lt("CLOUD_COVER", cloud_pct))\
        



    if dataset.size().getInfo() > 0:
        composite = dataset.median()
    else:
        # Create a default image if no data is available
        band_names = ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL"]
        default_bands = [ee.Image.constant(-1).rename(band) for band in band_names]
        composite = ee.Image.cat(default_bands)

    composite = composite.set('system:time_start', ee.Date(end_date).millis())

    download_ee_image(composite, os.path.join(directory_path, filename), scale=scale, region=roi, crs="EPSG:4326")
    blob_name=os.path.join(bucket_repository,directory_path.split("/")[-1],filename)
    upload_to_gcs(os.path.join(directory_path,filename),blob_name)
    return composite




def generate_metadata(id, year, month, directory_path_input, directory_path_output):
    metadata = {
        "dataset_description": "This dataset contains various satellite imagery and elevation data for forest height estimation.",
        "id": id,
        "year": year,
        "month": month,
        "images": {
            f"{id}_gedi_{year}_{month}.tif": {
                "description": "GEDI (Global Ecosystem Dynamics Investigation) data",
                "bands": ["rh98"],
                "source": "NASA's GEDI mission"
            },
            f"{id}_S1_{year}_{month}.tif": {
                "description": "Sentinel-1 SAR data",
                "bands": ["VV", "VH"],
                "source": "Copernicus Sentinel-1"
            },
            f"{id}_S2_{year}_{month}.tif": {
                "description": "Sentinel-2 multispectral data",
                "bands": ["B2", "B3", "B4", "B8", "B11", "B12"],
                "source": "Copernicus Sentinel-2"
            },
            f"{id}_srtm.tif": {
                "description": "SRTM (Shuttle Radar Topography Mission) elevation data",
                "bands": ["elevation" , "slope", "aspect"],
                "source": "NASA SRTM"
            },
            f"{id}_nicfi_{year}_{month}.tif": {
                "description": "NICFI (Norway's International Climate and Forest Initiative) high-resolution satellite imagery",
                "bands": ["red", "green", "blue", "nir"],
                "source": "Planet NICFI"
            },
            f"{id}_LC8_{year}_{month}.tif": {
                "description": "Landsat 8 multispectral data",
                "bands": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL"],
                "source": "USGS Landsat 8"
            }
        }
    }
    
    # Save metadata as JSON
    import json
    with open(os.path.join(directory_path_output, f"{id}_metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=4)

def dowload_for_one_polygon(id,polygon,start_date,end_date,directory_path_input,scale=10):
    os.makedirs(directory_path_input, exist_ok=True)
    
    try:
        get_S1_composite(polygon,start_date,end_date,directory_path_input,f"{id}_S1_composite.tif",scale=scale)
    except:
        print("no_sentinel1_data")
    try:
        get_S2_composite(polygon,start_date,end_date,directory_path_input,f"{id}_S2_composite.tif",scale=scale)
    except:
        print("no_sentinel2_data")
    try:
        get_nicfi_composite(polygon,start_date,end_date,directory_path_input,f"{id}_nicf8_composite.tif",scale=scale)
    except:
        print("no_nicfi_data")
    try:    
        get_landsat8_composite(polygon,start_date,end_date,directory_path_input,f"{id}_LC8_composite.tif",scale=scale)
    except:
        print("no_landsat8_data")
    try:    
        get_srtm(polygon,directory_path_input,f"{id}_srtm.tif",scale=scale)
    except:
        print("no_srtm")

    


    
    
    # Generate metadata after downloading all files

def create_mask(raster_path,directory_path,filename,shapes,landcover_path):
 
    with rasterio.open(raster_path) as src:
        raster = src.read(1)  # Lire la première bande du raster
        out_meta = src.meta
        raster_shape = raster.shape
        print(raster_shape)
      
        mask = np.zeros(raster_shape, dtype=np.uint8)

        # Rasterizer les polygones dans le masque
        mask = rasterize(shapes=shapes, out_shape=raster_shape, fill=0, out=mask, transform=out_meta['transform'], dtype='uint8',)

        # Sauvegarder le masque en tant que raster
        with rasterio.open(landcover_path) as lc_src:
            landcover = lc_src.read(1)
            lc_meta = lc_src.meta

            # Vérifier si les dimensions des rasters correspondent
            if mask.shape != landcover.shape:
                raise ValueError("Les dimensions du masque et du raster de classification ne correspondent pas")

            # Créer une copie du masque pour le mettre à jour
            updated_mask = np.copy(mask)
            condition_landcover_zero = (landcover == 0)
            landcover[condition_landcover_zero] = 22
            # Mettre à jour le masque selon la condition 
            condition = (mask == 0) & (landcover != 4)
            updated_mask[condition] = landcover[condition]
        
        out_meta.update({"count": 1, "dtype": 'uint8','nodata':100})

        with rasterio.open(os.path.join(directory_path,filename), "w", **out_meta) as dest:
            dest.write(updated_mask, 1)
        blob_name=os.path.join(bucket_repository,directory_path.split("/")[-1],filename)
        upload_to_gcs(os.path.join(directory_path,filename),blob_name)

        print(f"Le masque a été créé et sauvegardé dans {os.path.join(directory_path,filename)}")




    
def dynamic_world(
    region=None,
    start_date="2020-01-01",
    end_date="2021-01-01",
    clip=False,
    reducer=None,
    projection="EPSG:3857",
    scale=10,
    return_type="hillshade",
):
    """Create 10-m land cover composite based on Dynamic World. The source code is adapted from the following tutorial by Spatial Thoughts:
    https://developers.google.com/earth-engine/tutorials/community/introduction-to-dynamic-world-pt-1

    Args:
        region (ee.Geometry | ee.FeatureCollection): The region of interest.
        start_date (str | ee.Date): The start date of the query. Default to "2020-01-01".
        end_date (str | ee.Date): The end date of the query. Default to "2021-01-01".
        clip (bool, optional): Whether to clip the image to the region. Default to False.
        reducer (ee.Reducer, optional): The reducer to be used. Default to None.
        projection (str, optional): The projection to be used for creating hillshade. Default to "EPSG:3857".
        scale (int, optional): The scale to be used for creating hillshade. Default to 10.
        return_type (str, optional): The type of image to be returned. Can be one of 'hillshade', 'visualize', 'class', or 'probability'. Default to "hillshade".

    Returns:
        ee.Image: The image with the specified return_type.
    """

    if return_type not in ["hillshade", "visualize", "class", "probability"]:
        raise ValueError(
            f"{return_type} must be one of 'hillshade', 'visualize', 'class', or 'probability'."
        )

    if reducer is None:
        reducer = ee.Reducer.mode()

    dw = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").filter(
        ee.Filter.date(start_date, end_date)
    )

    if isinstance(region, ee.FeatureCollection) or isinstance(region, ee.Geometry):
        dw = dw.filterBounds(region)
    else:
        raise ValueError("region must be an ee.FeatureCollection or ee.Geometry.")

    # Create a Mode Composite
    classification = dw.select("label")
    #dwComposite = classification.sort('system:time_start', False).first()
    dwComposite = classification.mode()
    if clip and (region is not None):
        if isinstance(region, ee.Geometry):
            dwComposite = dwComposite.clip(region)
        elif isinstance(region, ee.FeatureCollection):
            dwComposite = dwComposite.clipToCollection(region)
        elif isinstance(region, ee.Feature):
            dwComposite = dwComposite.clip(region.geometry())

    dwVisParams = {
        "min": 0,
        "max": 8,
        "palette": [
            "#419BDF",
            "#397D49",
            "#88B053",
            "#7A87C6",
            "#E49635",
            "#DFC35A",
            "#C4281B",
            "#A59B8F",
            "#B39FE1",
        ],
    }

    if return_type == "class":
        return dwComposite
    elif return_type == "visualize":
        return dwComposite.visualize(**dwVisParams)
    else:
        # Create a Top-1 Probability Hillshade Visualization
        probabilityBands = [
            "water",
            "trees",
            "grass",
            "flooded_vegetation",
            "crops",
            "shrub_and_scrub",
            "built",
            "bare",
            "snow_and_ice",
        ]

        # Select probability bands
        probabilityCol = dw.select(probabilityBands)

        # Create a multi-band image with the average pixel-wise probability
        # for each band across the time-period
        meanProbability = probabilityCol.reduce(ee.Reducer.mean())

        # Composites have a default projection that is not suitable
        # for hillshade computation.
        # Set a EPSG:3857 projection with 10m scale
        proj = ee.Projection(projection).atScale(scale)
        meanProbability = meanProbability.setDefaultProjection(proj)

        # Create the Top1 Probability Hillshade
        top1Probability = meanProbability.reduce(ee.Reducer.max())

        if clip and (region is not None):
            if isinstance(region, ee.Geometry):
                top1Probability = top1Probability.clip(region)
            elif isinstance(region, ee.FeatureCollection):
                top1Probability = top1Probability.clipToCollection(region)
            elif isinstance(region, ee.Feature):
                top1Probability = top1Probability.clip(region.geometry())

        if return_type == "probability":
            return top1Probability
        else:
            top1Confidence = top1Probability.multiply(100).int()
            hillshade = ee.Terrain.hillshade(top1Confidence).divide(255)
            rgbImage = dwComposite.visualize(**dwVisParams).divide(255)
            probabilityHillshade = rgbImage.multiply(hillshade)

            return probabilityHillshade

def get_dynamic_world_mode(polygon,start_date="2020-12-01",end_date="2020-12-31",temp_directory="../data/temp",filename="dw.tif",scale=10):
    os.makedirs(temp_directory,exist_ok=True)
    region = ee.Geometry.Polygon(polygon)    
    image = dynamic_world(region, start_date, end_date,projection="EPSG:4326",return_type="class").clip(region)
    download_ee_image(image,os.path.join(temp_directory,filename), scale=scale, region=region, crs="EPSG:4326")
    blob_name=os.path.join(bucket_repository,temp_directory.split("/")[-1],filename)
    upload_to_gcs(os.path.join(temp_directory,filename),blob_name)
    return os.path.join(temp_directory,filename)
    
def main(data_path,start_date,end_date,scale=10,side=2560,class_name='Crop'):
    gdf=gpd.read_file(data_path)
    dataset_path=data_path.split(".geojson")[0]+"_dataset"
    os.makedirs(dataset_path,exist_ok=True)
    geoms=gdf.geometry
    gdf[class_name]=gdf[class_name].apply(lambda x:x.lower().replace(' - ','+').replace("é","e").replace("ï","i"))
    gdf['label'] = gdf[class_name].map(CLASSES_CODES) 
    shapes = [(geom, value) for geom, value in zip(gdf.geometry, gdf["label"])]
    for i,geom in enumerate(geoms):
        centroid=geom.centroid
        polygon=centroid_to_square(centroid,side=side)
        dowload_for_one_polygon(i,polygon,start_date,end_date,dataset_path,scale)
        filename_dw=get_dynamic_world_mode(polygon,start_date,end_date,dataset_path,'{}_dynamic_world.tif'.format(i),scale)
        create_mask(os.path.join(dataset_path,'{}_S1_composite.tif'.format(i)),dataset_path,'{}_label.tif'.format(i),shapes,filename_dw)
        break
def upload_to_gcs(source_file_path, destination_blob_name):
    
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_filename(source_file_path)
    print(f"Le fichier {source_file_path} a été téléchargé dans le bucket {bucket_name} avec le nom {destination_blob_name}")
if __name__ == "__main__":
    load_dotenv()
    ee.Initialize()
    storage_client = storage.Client()

    centroids_path=os.getenv("CENTROIDS_PATH")
    start_year=int(os.getenv("START_YEAR"))
    end_year=int(os.getenv("END_YEAR"))
    directory_path_input=os.getenv("TRAINING_DATA_PATH")
    gedi_path_output=os.getenv("TARGET_TRAINING_DATA_PATH")
    scale=int(os.getenv("SCALE"))
    side=int(os.getenv("SIDE"))
    bucket_name=os.getenv("BUCKET")
    bucket_repository=os.getenv("BUCKET_REPOSITORY")
    bucket = storage_client.bucket(bucket_name)

    main("../data/geometry/arachide_sample_4326.geojson","2024-07-01","2024-12-01")

    
   