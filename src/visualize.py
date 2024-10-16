import rasterio
from matplotlib import pyplot as plt

# Replace these with the paths to your files
label_file_path = '/Users/maika/Desktop/preprocessing-geospatial-data/data/cartagraphie_des_cultures_dataset_2024-07-01_2024-10-15_monthly_data_hivernal_2024/5066_label.tif'
dynamic_world_file_path = '/Users/maika/Desktop/preprocessing-geospatial-data/data/cartagraphie_des_cultures_dataset_2024-07-01_2024-10-15_monthly_data_hivernal_2024/5066_dynamic_world_2024_7.tif'

# Open the label file
with rasterio.open(label_file_path) as src:
    # Read the first band
    label_img = src.read(1)

# Open the dynamic world file
with rasterio.open(dynamic_world_file_path) as src:
    # Read the first band
    dynamic_world_img = src.read(1)

# Create a new figure with 2 subplots
fig, axs = plt.subplots(1, 2, figsize=(20,10))

# Display the label image in the first subplot
axs[0].imshow(label_img, cmap='viridis')
axs[0].set_title('Label')

# Display the dynamic world image in the second subplot
axs[1].imshow(dynamic_world_img, cmap='viridis')
axs[1].set_title('Dynamic World')

# Show the plot
plt.show()