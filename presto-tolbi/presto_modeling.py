r"""
Chesapeake CVPR Data Processing Script
======================================

This script processes GeoTIFF files from the Chesapeake CVPR dataset to create
image chips for segmentation tasks.

Dataset Source:
---------------
Chesapeake CVPR data from LILA:
https://lila.science/datasets/chesapeakelandcover

For this experiment, we will use images from NY.

Notes:
------
1. Only copy *_lc.tif & *_naip-new.tif files that we will use for our
segmentation downstream task.
   Using s5cmd for this: https://github.com/peak/s5cmd
   - Train:
   s5cmd cp \
        --no-sign-request \
        --include "*_lc.tif" \
        --include "*_naip-new.tif" \
        "s3://us-west-2.opendata.source.coop/agentmorris/lila-wildlife/lcmcvpr2019/cvpr_chesapeake_landcover/ny_1m_2013_extended-debuffered-train_tiles/*" \
        data/cvpr/files/train/
   - Val:
   s5cmd cp \
        --no-sign-request \
        --include "*_lc.tif" \
        --include "*_naip-new.tif" \
        "s3://us-west-2.opendata.source.coop/agentmorris/lila-wildlife/lcmcvpr2019/cvpr_chesapeake_landcover/ny_1m_2013_extended-debuffered-val_tiles/*" \
        data/cvpr/files/val/

2. We will create chips of size `224 x 224` to feed them to the model, feel
free to experiment with other chip sizes as well.
   Run the script as follows:
   python preprocess_data.py <data_dir> <output_dir> <chip_size>

   Example:
   python preprocess_data.py data/cvpr/files data/cvpr/ny 224
"""  # noqa E501
import rioxarray
from pyproj import Transformer
import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score

from tqdm import tqdm

import torch
#from torch.utils.data import DataLoader, TensorDataset

import presto

# this is to silence the xarray deprecation warning.
# Our version of xarray is pinned, but we'll need to fix this
# when we upgrade
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
import os
import sys
from pathlib import Path

import numpy as np
import rasterio as rio


from pathlib import Path

def process_pixel(id,pos_i,pos_j, data_dir,sensor_prefix="_S2_*.tif"):
    # Get the list of files and sort them
    globs=list((data_dir).glob(f"{id}{sensor_prefix}"))
    if len(globs)>1:
        
        filenames = sorted(list((data_dir).glob(f"{id}{sensor_prefix}")), key=lambda x: int(x.stem.split('_')[-1]))
    else:
        filenames = list((data_dir).glob(f"{id}{sensor_prefix}"))

    
    # Rest of your code...def process_images(id,data_dir):
    arrays = []
    for filename in (filenames):
        tif_file = rioxarray.open_rasterio(filename)
        crs = tif_file.rio.crs
        transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)

        # firstly, get the latitudes and longitudes
        x, y = tif_file.x[pos_j], tif_file.y[pos_i]
        lon, lat = transformer.transform(x, y) 
        # then, get the eo_data, mask and dynamic world
        s2_data_for_pixel = tif_file.values[:, pos_i, pos_j].astype(int)
        arrays.append(s2_data_for_pixel)

    return np.stack(arrays,axis=0) ,np.array(lat,lon)
      


def process_pixel_for_all_sensors(id,pos_i,pos_j, data_dir):
    # Get the list of files and sort them


    s2_data,_=process_pixel(id,pos_i,pos_j, data_dir,sensor_prefix="_S2_*.tif")
    s1_data,_=process_pixel(id,pos_i,pos_j, data_dir,sensor_prefix="_S1_*.tif")
    dynamic_world_data,_=process_pixel(id,pos_i,pos_j, data_dir,sensor_prefix="_dynamic_world_*.tif")
    srtm_data,_=process_pixel(id,pos_i,pos_j, data_dir,sensor_prefix="_srtm*.tif")
    label_data,latlon=process_pixel(id,pos_i,pos_j, data_dir,sensor_prefix="_label.tif")
    x, mask, dynamic_world = presto.construct_single_presto_input(
                    s2=s2_data, s2_bands=['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10','B11', 'B12'],
                    s1=s1_data, s1_bands=['VV', 'VH'],
                    dynamic_world=dynamic_world_data,
                    srtm=srtm_data[:,:2],
                    srtm_bands= ["elevation", "slope"]
                )
    
    return x, mask, dynamic_world,latlon, label_data.rehsape(1)
def process_id(id,data_dir):
    ff_path=list((data_dir).glob(f"{id}_S2_*.tif"))[0]
    tif_file = rioxarray.open_rasterio(ff_path)
    _,h,w=tif_file.shape
    arrays, masks, dynamic_worlds, latlons, labels= [], [], [], [], []
    for i in range(h):
        for j in range(w):
            x, mask, dynamic_world,latlon,label_data=process_pixel_for_all_sensors(id,i,j,data_dir)
            arrays.append(x)
            masks.append(mask)
            dynamic_worlds.append(dynamic_world)
            latlons.append(latlon)
            labels.append(label_data.item())
    return arrays, masks, dynamic_worlds, latlons, labels

def process_all_ids(ids,data_dir):
    arrays, masks, dynamic_worlds, latlons, labels= [], [], [], [], []
    for id in ids:
        xs, _masks, _dynamic_worlds, _latlons,_labels=process_id(id,data_dir)
        #not append it is a list ,  want to append the elements

        arrays.extend(xs)
        masks.extend(_masks)
        dynamic_worlds.extend(_dynamic_worlds)
        latlons.extend(_latlons)
        labels.extend(_labels.item())
    
    return (torch.stack(arrays, axis=0),
            torch.stack(masks, axis=0),
            torch.stack(dynamic_worlds, axis=0),
            torch.stack(latlons, axis=0),
            torch.tensor(labels),
        )

    



    










       

                 
    # for filename in tqdm(filenames):
    #     tif_file = rioxarray.open_rasterio(filename)
    #     crs = tif_file.rio.crs
    #     transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)

    #     for x_idx in range(tif_file.shape[2]):
    #         for y_idx in range(tif_file.shape[1]):
                
      
    #             # firstly, get the latitudes and longitudes
    #             x, y = tif_file.x[x_idx], tif_file.y[y_idx]
    #             lon, lat = transformer.transform(x, y) 
    #             latlons.append(torch.tensor([lat, lon]))
                
    #             # then, get the eo_data, mask and dynamic world
    #             s2_data_for_pixel = torch.from_numpy(tif_file.values[:, x_idx, y_idx].astype(int)).float()
    #             s2_data_with_time_dimension = s2_data_for_pixel.unsqueeze(0)

    #             print(s2_data_with_time_dimension.shape)





    # x, mask, dynamic_world = presto.construct_single_presto_input(
    #     s2=s2_data_with_time_dimension, s2_bands=TREESATAI_S2_BANDS
    # )
    # arrays.append(x)
    # masks.append(mask)
    # dynamic_worlds.append(dynamic_world)
    
    # labels.append(0 if filename.startswith("Abies") else 1)
    # image_names.append(filename)

#     return (torch.stack(arrays, axis=0),
#             torch.stack(masks, axis=0),
#             torch.stack(dynamic_worlds, axis=0),
#             torch.stack(latlons, axis=0),
#             torch.tensor(labels),
#             image_names,
#         )
def main():
    """
    Main function to process files and create chips.
    Expects three command line arguments:
        - data_dir: Directory containing the input GeoTIFF files.
        - output_dir: Directory to save the output chips.
        - chip_size: Size of the square chips.
    """
    if len(sys.argv) != 4:  # noqa: PLR2004
        print("Usage: python script.py <data_dir> <output_dir> <chip_size>")
        sys.exit(1)

    data_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    chip_size = int(sys.argv[3])


    ids = list(set([ str(filename).split("/")[-1].split("_")[0]
             for filename in list((data_dir).glob("*_S2_*"))]))
    label_paths = list((data_dir).glob("*_label.tif"))

   
    res,_=process_pixel_for_all_sensors(ids[0],2,2,data_dir)
    print(res)
          

    # process_files(input_image_paths, output_dir / "S2/chips", chip_size)
    
    # process_files(label_paths, output_dir / "labels", chip_size)


if __name__ == "__main__":
    main()


