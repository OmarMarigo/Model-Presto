import ee
import numpy as np
from datetime import datetime
import ccd
import json
ee.Initialize()
point = ee.Geometry.Point( [ -71.083640209143056, 1.41497978268366 ])
import matplotlib.pyplot as plt
import pandas as pd


# Define the time range
start_date = '2020-01-01'
end_date = '2022-12-31'

# Load the Landsat image collection for the specified time range and location
# collection = ee.ImageCollection('LANDSAT/LC08/C01/T1') \
#                 .filterDate(start_date, end_date) \
#                 .filterBounds(point)

#Function to extract the required band values for a point
def extract_bands(image):
    # Replace 'B4', 'B3', etc., with the correct band identifiers for your Landsat product
    bands = image.select(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP',"QA_PIXEL"]).reduceRegion(
        geometry=point,
        reducer=ee.Reducer.toList(),
        scale=30
    )
    date = image.date().format('YYYY-MM-dd')
    # Combine date and band values into a dictionary
    feature = ee.Feature(None, {'date': date, 'bands': bands})
    return feature




# print(len(extracted_data_list.get('features')))

# # Process the data as needed to create your final array structure
# # This will depend on the specifics of your data and requirements

# # Print or use the structured data
# print(extracted_data_list)

# Define the functions to prepare Landsat images with masked noisy pixels.
def prepare_l4l5(image):
    band_list = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6',"QA_PIXEL"]
    name_list = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP',"QA_PIXEL"]
    scale = 2.75e-05

    scale_temp=0.00341802
    scaling = [scale, scale, scale, scale, scale, scale, scale_temp,1]
    scaled = image.select(band_list).rename(name_list).multiply(scaling)

    valid_qa = [66, 130, 68, 132]
 
    mask1 = image.select(['QA_PIXEL']).remap(valid_qa, ee.List.repeat(1, len(valid_qa)), 0)
    mask2 = image.select('radsat_qa').eq(0)
    mask3 = image.select(band_list).reduce(ee.Reducer.min()).gt(0)
    mask4 = image.select("sr_atmos_opacity").unmask().lt(300)
    return image.addBands(scaled)#.updateMask(mask1)#.And(mask2).And(mask3).And(mask4))

def prepare_l7(image):
    # Prepare L7 images similarly to L4/L5, with an additional erosion step to remove scan line artifacts.
    prepared_image = prepare_l4l5(image)
    mask5 = prepared_image.mask().reduce(ee.Reducer.min()).focal_min(2.5)
    return prepared_image#.updateMask(mask5)

def prepare_l8(image):
    # Adjusted band_list and validation for Landsat 8
    band_list = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10',"QA_PIXEL"]
    name_list = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP',"QA_PIXEL"]
    scale = 2.75e-05

  

    scale_temp=0.00341802
    scaling = [scale, scale, scale, scale, scale, scale, scale_temp,1]

    scaled = image.select(band_list).rename(name_list).multiply(scaling)
    mask1 = image.select('QA_PIXEL').bitwiseAnd(int('11111', 2)).eq(0)
    mask2 = image.select('radsat_qa').eq(0)
    mask3 = image.select(band_list).reduce(ee.Reducer.min()).gt(0)
    #mask4 = image.select(['sr_aerosol']).remap(valid_toa, ee.List.repeat(1, len(valid_toa)), 0)
    return image.addBands(scaled)#.updateMask(mask1)#.And(mask2).And(mask3)#.And(mask4))

# Example usage for a given region and date range.
options = {
    'start': '1985-01-01',
    'end': '2023-01-01',
    'useMask': True,
    'sensors': {'l4': True, 'l5': True, 'l7': True, 'l8': True}
}

# Adjust the filtering and collection merging process according to options.

collection4 = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2')#.filterDate(options['start'], options['end'])
collection5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')#.filterDate(options['start'], options['end'])
collection7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')#.filterDate(options['start'], options['end'])
collection8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')#.filterDate(options['start'], options['end'])
collection9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')#.filterDate(options['start'], options['end'])

if options['useMask']:
    collection9 = collection9.map(prepare_l8)
    collection8 = collection8.map(prepare_l8)
    collection7 = collection7.map(prepare_l7)
    collection5 = collection5.map(prepare_l4l5)
    collection4 = collection4.map(prepare_l4l5)

#all_collections = collection4.merge(collection5).merge(collection7)#.merge(collection8).merge(collection9)
all_collections = collection8#.merge(collection9)
params = {'QA_BITPACKED': False,
              'QA_FILL': 255,
              'QA_CLEAR': 0,
              'QA_WATER': 1,
              'QA_SHADOW': 2,
              'QA_SNOW': 3,
              'QA_CLOUD': 4}
# Optional: Filter by geometry if a region is specified.
# geometry = ee.Geometry.Point([longitude, latitude])  # Define your region's geometry here.
all_dict=[]
print("yes")
data=json.load(open('/Users/maika/Desktop/Earth-Phenotype/data/processed/representative_sample.geojson'))
#compute the time for each feature iteration
# Create a dictionary for Level 1 land cover values
level1_land_cover = {
    1: "Water",
    2: "Ice/snow",
    3: "Developed",
    4: "Barren/sparsely vegetated",
    5: "Trees",
    6: "Shrub",
    7: "Herbaceous"
}
class1_vals=[]
all_swir1s=[]
all_dates=[]
all_qas=[]
for num, feat in enumerate(data['features']):


    point = ee.Geometry.Point(feat['geometry']['coordinates'])
    start_date = feat['properties']['Start_Year']
    end_date = feat['properties']['End_Year']
    id=feat['properties']["Glance_ID"]
    class1_val=feat['properties']["Glance_Class_ID_level1"]
    class2_val=feat['properties']["Glance_Class_ID_level2"]
    lon=feat["properties"]["Lon"]
    lat=feat["properties"]["Lat"]
    col = all_collections.filterBounds(point).filterDate(f'{start_date}-01-01', f'{end_date}-12-31')    
  
    # Map the function over the collection
    extracted_data = col.map(extract_bands)

    
  
    
    # # Get the data out of the GEE environment
    try:

        extracted_data_list = extracted_data.getInfo()
        features=extracted_data_list.get('features')
    except:
    
        raise ValueError
        continue
    # Initialize a list to hold our data rows
    data_rows = []

    # Process each feature

    for feature in features:
        bands = feature['properties']['bands']
        # Skip the feature if any band is an empty list
        if any(len(bands[band]) == 0 for band in bands):
            continue
        # Convert date to ordinal
        date_ordinal = int(datetime.strptime(feature['properties']['date'], "%Y-%m-%d").toordinal())
        # Extract band values in the order: Date (ordinal), RED, GREEN, NIR, SWIR1, SWIR2, TEMP
        row = [
            date_ordinal,
            bands['BLUE'][0] if 'BLUE' in bands and bands['BLUE'] else np.nan,
            bands['GREEN'][0] if 'GREEN' in bands and bands['GREEN'] else np.nan,
            bands['RED'][0] if 'RED' in bands and bands['RED'] else np.nan,
            bands['NIR'][0] if 'NIR' in bands and bands['NIR'] else np.nan,
            bands['SWIR1'][0] if 'SWIR1' in bands and bands['SWIR1'] else np.nan,
            bands['SWIR2'][0] if 'SWIR2' in bands and bands['SWIR2'] else np.nan,
            bands['TEMP'][0] if 'TEMP' in bands and bands['TEMP'] else np.nan,
            0,
            bands['QA_PIXEL'][0] if 'QA_PIXEL' in bands and bands['QA_PIXEL'] else np.nan,
        ]
        data_rows.append(row)

    # Convert the list of rows into a NumPy array
    data_rows.sort(key=lambda x: x[0])
    data_array = np.array(data_rows)
    
    if data_array.size == 0:
        continue
    dates = data_array[:, 0]  # First column
    blues = data_array[:, 1]  # Second column
    greens = data_array[:, 2]  # and so on...
    reds = data_array[:, 3]
    nirs = data_array[:, 4]
    swir1s = data_array[:, 5]
    swir2s = data_array[:, 6]
    thermals = data_array[:, 7]
    qas = data_array[:, 9]  # Ninth column
    dict_vectors={
        "id":id,
        "Glance_Class_ID_level1":class1_val,
        "Glance_Class_ID_level2":class2_val,
        "longitude":lon,
        "latitude":lat,
        "start_date":start_date,
        "end_date":end_date,
        "dates":dates,
        "blues":blues,
        "greens":greens,
        "reds":reds,
        "nirs":nirs,
        "swir1s":swir1s,
        "swir2s":swir2s,
        "thermals":thermals,
        "qa_pixel":qas
    }
    

    all_swir1s.append(swir1s)
    all_dates.append(dates)
    class1_vals.append(class1_val)
    all_qas.append(qas)
    

grouped_data = {}
for i, class1_val in enumerate(class1_vals):
    if class1_val not in grouped_data:
        grouped_data[class1_val] = []
    grouped_data[class1_val].append((all_dates[i], all_swir1s[i],all_qas[i]))


# Get all unique classes
all_classes = sorted(grouped_data.keys())

# Create subplots for each class separately
for class_val in all_classes:

    class_data = grouped_data[class_val]
    
    # Create a new figure for each class
    fig, axs = plt.subplots(5, 2, figsize=(20, 25))
    fig.suptitle(f'SWIR1 Values Over Time for Class: {level1_land_cover[class_val]}', fontsize=16)
    
    for i, (dates, swir1,qas) in enumerate(class_data):
        dates_datetime = [datetime.fromordinal(int(date)) for date in dates]
        print(len(dates_datetime))
    
        row = i // 2
        col = i % 2
        ax = axs[row, col]


        dates_datetime = [datetime.fromordinal(int(date)) for date in dates]
        # Create a mask for cloud pixels
        qas_int = np.array(qas).astype(np.int64)
        print(qas_int)
        cloud_mask = [((int(qa) & 0b0000000000001000) != 0) or  # Cloud (bit 3)
              ((int(qa) & 0b0000000000000010) != 0)     # Dilated cloud (bit 1)
              for qa in qas_int]

        
        # Plot SWIR1 values
        ax.scatter(dates_datetime, swir1, label=f'Record {i+1}', alpha=0.7,color="green")
        
        # Plot red dots for cloud pixels
        cloud_dates = [date for date, is_cloud in zip(dates_datetime, cloud_mask) if is_cloud]
        cloud_swir1 = [value for value, is_cloud in zip(swir1, cloud_mask) if is_cloud]
        ax.scatter(cloud_dates, cloud_swir1, color='red', s=20, label='Cloud Pixels')
        ax.plot(dates_datetime, swir1, label=f'Record {i+1}', alpha=0.7)
        
        ax.set_xlabel('Date')
        ax.set_ylabel('SWIR1 Value')
        ax.set_title(f'Record {i+1}')
        ax.legend()
        ax.grid(True)
        

    
    # Remove any unused subplots
    for i in range(len(class_data), 10):
        row = i // 2
        col = i % 2
        fig.delaxes(axs[row, col])
    
    plt.tight_layout()
    plt.savefig(f'../../data/visualizations/with_cloud_tag/swir1_subplots_class_{class_val}.png')
    plt.close()


