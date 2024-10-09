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
import xarray
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt
import geemap
from loguru import logger


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
        sentinel1.median().set('system:time_start', ee.Date.fromYMD(end_date).millis()),
        create_default_image(roi, ['VV', 'VH']).set('system:time_start',
                                                             ee.Date.fromYMD(end_date).millis())
    )
    download_ee_image(ee.Image(image),os.path.join(directory_path,filename), scale=scale, region=roi, crs="EPSG:4326")

    return ee.Image(image)



def get_S2_composite(polygon, start_date,end_date,directory_path, filename,scale=10,cloud_pct=20):
    def create_default_image(roi, band_names):
        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        default_image = ee.Image.cat(default_bands)
        return default_image
    roi = ee.Geometry.Polygon(polygon)
    dataset = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
        .filterBounds(roi) \
         .filter(ee.Filter.date(start_date, end_date)) \
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloud_pct)).clip(roi)
    
    bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10','B11', 'B12']
    
    def scale_bands(image):
        scaled = image.select(bands[:-1]).divide(10000)
        return image.addBands(scaled, overwrite=True).select(bands)



    composite = ee.Algorithms.If(
        dataset.size().gt(0),
        dataset.median().set('system:time_start', start_date.millis()),
        create_default_image(roi, bands).set('system:time_start', start_date.millis())
    )

    image = ee.Image(composite)
    download_ee_image(image,os.path.join(directory_path,filename), scale=scale, region=roi, crs="EPSG:4326")

    return image







def apply_srtm(polygon,directory_path,filename,scale=10):
    # --- NDVI (Normalized Difference Vegetation Index) ---
    roi = ee.Geometry.Polygon(polygon)
    srtm = ee.Image("USGS/SRTMGL1_003").clip(roi)
    # --- Elevation ,Slope and Aspect ---
    elevation = srtm.select('elevation').rename('elevation')

    slope = ee.Terrain.slope(srtm).rename('slope')

    aspect =ee.Terrain.aspect( srtm)
    image=elevation.addBands([slope, aspect])
    download_ee_image(image,os.path.join(directory_path,filename), scale=scale, region=roi, crs="EPSG:4326")
    # Adding all indices as bands to the image
    return image

def get_centroid(polygon):
    roi = ee.Geometry.Polygon(polygon)
    centroid = roi.centroid()
    return centroid




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
        .clip(roi)



    if dataset.size().getInfo() > 0:
        composite = dataset.median()
    else:
        # Create a default image if no data is available
        band_names = ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL"]
        default_bands = [ee.Image.constant(-1).rename(band).clip(roi) for band in band_names]
        composite = ee.Image.cat(default_bands)

    composite = composite.set('system:time_start', ee.Date.fromYMD(year, month, 1).millis())

    download_ee_image(composite, os.path.join(directory_path, filename), scale=scale, region=roi, crs="EPSG:4326")

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
        get_S2_composite(polygon,start_date,end_date,directory_path_input,f"{id}_S2.tif",scale=scale)
    except:
        print("no_sentinel2_data")
    try:
        get_nicfi_composite(polygon,start_date,end_date,directory_path_input,f"{id}_nicfi.tif",scale=scale)
    except:
        print("no_nicfi_data")
    try:    
        get_landsat8_composite(polygon,start_date,end_date,directory_path_input,f"{id}_LC8.tif",scale=scale)
    except:
        print("no_landsat8_data")
    
    # Generate metadata after downloading all files

def create_mask(polygon,S2_img):
    pass
    

def download_all_month_for_one_polygon(id_polygon,polygon,year,directory_path_input,gedi_path_output,scale=10):
    os.makedirs(directory_path_input, exist_ok=True)
    os.makedirs(gedi_path_output, exist_ok=True)
    roi=ee.Geometry.Polygon(polygon)
    def vqualityMask(im):
        return im.updateMask(im.select('quality_flag').eq(1)).updateMask(im.select('degrade_flag').eq(0));
   
    dataset = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY').map(vqualityMask).select('rh98').filterDate(f"{year-1}-10-01",f"{year}-03-31");
    # Select the height layer (assuming height data is under 'rh100' band, adjust if necessary)
    gedi_image = dataset.median()

    # Mask NoData values (e.g., assuming no-data values are represented by 0)

    # Reduce over the region to check for valid data
    valid_pixel_count = gedi_image.reduceRegion(
        reducer=ee.Reducer.count(),
        geometry=roi,
        scale=30,  # Adjust scale according to the resolution you need
        maxPixels=1e13
    )

    # Get the result

    valid_pixel_count_value = valid_pixel_count.getInfo()
  

    # Check if there are any valid pixels
    
    if valid_pixel_count_value['rh98'] > 0:
        for month in range(1, 13):
            dowload_one_polygon_for_one_month(id_polygon,polygon,year-1,month,directory_path_input,scale=scale)
        apply_srtm(polygon,directory_path_input,f"{id_polygon}_srtm.tif",scale=scale)
        download_ee_image(gedi_image,os.path.join(gedi_path_output,f"{id_polygon}_gedi_{year}_01.tif"), scale=10, region=roi, crs="EPSG:4326")
    else:
        pass


def generate_data_for_one_polygon(feature,start_year,end_year,directory_path_input,gedi_path_output,scale=10,side=2560):
    for year in range(start_year,end_year+1):
        id_polygon=feature['properties']['id']
        logger.success(f"Processing polygon {id_polygon} for year {year}")
        polygon=centroid_to_square(feature['geometry']['coordinates'],side)

        download_all_month_for_one_polygon(id_polygon,polygon,year,directory_path_input,gedi_path_output,scale=scale)

def main(centroids_path,start_year,end_year,directory_path_input,gedi_path_output,scale=10,side=2560):
    data=json.load(open(centroids_path) )
    features=data['features']
    for feat in features:
        generate_data_for_one_polygon(feat,start_year,end_year,directory_path_input,gedi_path_output,scale=scale,side=side)


if __name__ == "__main__":
    load_dotenv()
    ee.Initialize()
    centroids_path=os.getenv("CENTROIDS_PATH")
    start_year=int(os.getenv("START_YEAR"))
    end_year=int(os.getenv("END_YEAR"))
    directory_path_input=os.getenv("TRAINING_DATA_PATH")
    gedi_path_output=os.getenv("TARGET_TRAINING_DATA_PATH")
    scale=int(os.getenv("SCALE"))
    side=int(os.getenv("SIDE"))
    print(centroids_path,start_year,end_year,directory_path_input,gedi_path_output,scale,side)
    main(centroids_path,start_year,end_year,directory_path_input,gedi_path_output,scale=scale,side=side)
   