import ee
cloud_project = 'ee-tolbico2024'

try:
        ee.Initialize(project=cloud_project)
except:
        ee.Authenticate()
        ee.Initialize(project=cloud_project)
import ee
cloud_project = 'ee-tolbico2024'

try:
        ee.Initialize(project=cloud_project)
except:
        ee.Authenticate()
        ee.Initialize(project=cloud_project)


def dowload_one_polygon(id,polygon,year,month,directory_path_input,directory_path_output):
    apply_gedi_month(polygon,year,month,directory_path_output,f"{id}_gedi_{year}_{month}.tif")
    apply_sentinel1_month(polygon,year,month,directory_path_input,f"{id}_S1_{year}_{month}.tif")
    apply_sentinel2_monthly(polygon,year,month,directory_path_input,f"{id}_S2_{year}_{month}.tif")
    apply_srtm(polygon,directory_path,f"{id}_srtm.tif")
    add_nicfi_monthly(id,polygon,year,month,directory_path,f"{id}_nicfi_{year}_{month}.tif")
    apply_landsat8_month(polygon,year,month,directory_path,f"{id}_LC8_{year}_{month}.tif")

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