import os
import rasterio
import numpy as np

def print_raster_shapes(directory):
    for filename in os.listdir(directory):
        if filename.endswith(('.tif', '.tiff')):
            filepath = os.path.join(directory, filename)
            with rasterio.open(filepath) as src:
                shape = src.read().shape
                print(f"{filename}: {shape}")
                print(src.read())
                #count na for each raster
                na_count = np.sum(src.read() == src.nodata)
                print(f"Number of NA values: {na_count}")
                #print the no data value
                print(f"No data value: {src.nodata}")

# Directory path
directory = '/Users/maika/Desktop/preprocessing-geospatial-data/src/target_train_features'

# Call the function to print shapes
print_raster_shapes(directory)
