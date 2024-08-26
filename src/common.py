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
import rasterio
import shapely
import geemap
import wxee
import xarray 
import numpy as np
from datetime import datetime
import os 
import rasterio
from rasterio.transform import from_origin
from rasterio.enums import Resampling
import  matplotlib.pyplot as plt


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

def load_file_apply_buffer_square(data_path, geometry="geometry", side =2560):

    if data_path is not None :
        
        data = gpd.read_file(data_path)
        

        data_africa = data[data['Continent_Code'] == 3]
        
        center = data_africa[geometry].centroid
        

        proj, proj_inverse = get_projections(center)
        cartesian_center = proj(center)
        cartesian_square = cartesian_center.buffer(side/2 , cap_style=3)
        
        geodesic_square = proj_inverse(cartesian_square)
        print(json.dumps(shapely.geometry.mapping(geodesic_square)))
        
        return json.dumps(shapely.geometry.mapping(geodesic_square)) 
    else:
        return "Error: File path is not provided."


def get_gedi(polygon,output_path: str,monthly_data=list(),filename="data_buffer_square.tif", monthly_times=list(), start_date= "2019-01-01", end_date="2019-12-31"):

        def quality_mask(image):
            quality_mask = image.select('quality_flag').eq(1)
            degrade_mask = image.select('degrade_flag').eq(0)
            return image.updateMask(quality_mask).updateMask(degrade_mask)
        def set_time_start(image):
            return image.set('system:time_start', ee.Date(image.get('system:time_start')).millis())

        roi = ee.Geometry.Polygon(polygon)
        dataset_collection = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY').filterBounds(roi) \
            .filterDate(start_date, end_date) \
            .map(quality_mask)\
            .select(['rh98', 'rh100'])\
            .map(set_time_start)

        image = dataset_collection.mean()


        #geemap.download_ee_image_collection(dataset_collection, region=roi, out_dir=output_path,crs="EPSG:4326")

        da = dataset_collection.wx.to_xarray(region = roi,  scale=10)


        df = da.to_dataframe().reset_index()
        if 'time' in df.columns:
            df['date'] = pd.to_datetime(df['time']).dt.date  # Extract only the date
        # Optionally, you can drop the 'time' column if not needed
            df = df.drop(columns=['time'])
            df = df.rename(columns={"date":"time"})


        # Save DataFrame to CSV
        df.to_csv("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/gedi.csv", index=False)



        """datasets = []
        for image_file in os.listdir(output_path):
            if image_file.endswith('.tif'):
                file_path = os.path.join(output_path, image_file)
                with rasterio.open(file_path) as src:
                    data = src.read()
                    transform = src.transform
                    height, width = src.height, src.width

                    # Générer les coordonnées
                    rows, cols = np.meshgrid(np.arange(height), np.arange(width), indexing='ij')
                    lon, lat = xy(transform, rows, cols)
                    lon = np.array(lon)
                    lat = np.array(lat)




                    # Convertir en xarray.Dataset
                    band_names = ['rh98', 'rh100']
                    ds = xarray.Dataset(
                        data_vars={
                            band_names[i]: (['time', 'lat', 'lon'], data[i, np.newaxis, :, :])
                            for i in range(data.shape[0])
                        },
                        coords={
                            'time': (['time'], [datetime.strptime(dt, '%Y-%m-%d') for dt in monthly_times]),
                            'lon': (['lat', 'lon'], lon),
                            'lat': (['lat', 'lon'], lat)
                        },
                        attrs=dict(description="GEDI data")
                    )

                    # Ajouter le Dataset à la liste
                    datasets.append(ds)


                print(datasets) """







def get_sentinel1_monthly(polygon, output_path, year=2019, monthly_data =[],  monthly_times = []):

    roi = ee.Geometry.Polygon(polygon)

    i = 0
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
        
        # Calculate the median composite for the selected month
        composite = sentinel1.median().set('system:time_start', start_date.timestamp() * 1000)
        images.append(composite)


        # Define the filename
        filename = f"sentinel1_composite_{year}_{month:02d}.tif"
        full_output_path = os.path.join(output_path, filename)

        # Download the composite image as a GeoTIFF
        """geemap.download_ee_image(
            image=composite,
            filename=full_output_path,
            scale=10,
            region=roi,
            crs="EPSG:4326"
        )"""

       # print(f"Downloaded {filename} to {output_path}")

        i+=1
        

    # Convert the list to a 4D numpy array (time, band, y, x)
    imageCollection= ee.ImageCollection.fromImages(images)
    da = imageCollection.wx.to_xarray(region=roi, scale=10)

    df = da.to_dataframe().reset_index()

    # Save DataFrame to CSV
    df.to_csv("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/sentinel_1.csv",index=False)

    """
    data_array = np.stack(monthly_data, axis=0)
    num_bands = data.shape[0]
    
    ### verifciation de nombre de bandes 
    if num_bands == 2:
        band_names = ["VV", "VH"]
    elif num_bands == 3:
        band_names = ["VV_Asc", "VH_Asc", "VV_Desc"]
    else:
        band_names



    print(data_array.shape)

    data_vars = {}
    for i in range(data_array.shape[1]):
        data_vars[band_names[i]] = (['time', 'lat', 'lon'], data_array[:, i, :])

    ds = xarray.Dataset(
        data_vars={
            band_names[i]: (['time', 'lat', 'lon'], data_array[:, i, :, :])
            for i in range(data_array.shape[1])
        },
        coords={
            'lon': (['lat', 'lon'], lon),
            'lat': (['lat', 'lon'], lat),
            'time': (['time'], monthly_times)
        },
        attrs=dict(description="Sentinel 1 bands")
    )

    ad = ds.to_dataframe()
    ad.to_csv('/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/sentinel_1.csv')
"""






def get_sentinel2_monthly(polygon, output_path, year=2019,monthly_data = list(),monthly_times = list()):
    
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

    # Define the region of interest (ROI)
    roi = ee.Geometry.Polygon(polygon)
    images=[]
    i=0
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
        dataset = ee.ImageCollection('COPERNICUS/S2_SR')\
                    .filterBounds(roi)\
                    .filterDate(start_date_str, end_date_str)\
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))\
                    .map(mask_s2_clouds)
        
        # Create a composite image for the month
        composite = dataset.median().set('system:time_start', start_date.timestamp() * 1000)
        images.append(composite)
        # Define the filename
        filename = f"Sentinel_2_{year}_{month:02d}.tif"
        
        full_output_path = os.path.join(output_path, filename)
        
        """geemap.download_ee_image(
            image=composite,
            filename=full_output_path,
            scale=10,
            region=roi,
            crs="EPSG:4326"
        )

        print(f"Downloaded {filename} to {output_path}")"""

        
        i+=1


    imageCollection: ImageCollection=ee.ImageCollection.fromImages(images)
    da = imageCollection.wx.to_xarray(region = roi,  scale=10)


    df = da.to_dataframe().reset_index()

        # Save DataFrame to CSV
    df.to_csv("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/sentinel_2.csv", index=False)
    
    # Print the coordinates (or process them as needed)
    """
    data_array = np.stack(monthly_data, axis=0)
    print(data_array.shape)
    band_names = ['B{}'.format(i+1) for i in range(data_array.shape[1])]
    data_vars={}
    for i in range(data_array.shape[1]):
        data_vars[band_names[i]]= (['time', 'lat', 'lon'], data_array[:, i, :])

    ds = xarray.Dataset(
        data_vars={
            band_names[i]: (['time', 'lat', 'lon'], data_array[:, i, :, :])
            for i in range(data_array.shape[1])
        },
        coords={
            'lon': (['lat', 'lon'], lon),
            'lat': (['lat', 'lon'], lat),
            'time': (['time'], monthly_times)
        },
        attrs=dict(description="Sentinel 2 bands")
    )
    print(ds["lon"])
    ad = ds.to_dataframe()
    ad.to_csv('/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/sentinel_2.csv')
  

   
    da = xarray.DataArray(
    data_array,
    dims=("time", "band", "y", "x"),
    coords={
        "time": monthly_times,
        "band": band_names,
        "y": np.arange(data_array.shape[2]),
        "x": np.arange(data_array.shape[3])
    },
    attrs={
        "transform": transform,
        "crs": crs
    }
)
    print(da)

    ds = xarray.Dataset({"data": da})
    print(ds)
    """
    






if __name__ == "__main__":
    #with rasterio.open("/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/gedi_output.nc/201904_042E_012N.tif") as src:
        #print(src.read().shape)
    #exit()
    wxee.Initialize(project=os.getenv("ID_NAME_PROJECT_EE"))

    POLYGON_COORDS = json.loads(os.getenv("POLYGON_COORDS"))
    OUT_PATH_SENTINEL_1 = os.getenv("OUT_PATH_SENTINEL_1")
    OUT_PATH_SENTINEL_2 = os.getenv("OUT_PATH_SENTINEL_2")
    OUT_PATH_FUSED = os.getenv("OUT_PATH_FUSED")

    # Récupérer les données GEDI, Sentinel-1 et Sentinel-2
    #get_gedi(POLYGON_COORDS, output_path="/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/gedi_output.nc")
    #get_sentinel1_monthly(POLYGON_COORDS, OUT_PATH_SENTINEL_1)
    #get_sentinel2_monthly(POLYGON_COORDS, OUT_PATH_SENTINEL_2)
    load_file_apply_buffer_square()


    #
    # Save to TIFF
    #save_to_tiff(fused_da, OUT_PATH_FUSED)
    #with rasterio.open('/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/gedi_output.nc/201904_042E_012N.tif','r') as src:
        #data = src.read(1)
        #plt.imshow(data)
        #plt.show()


def get_gedi_bimass():
    return