import rasterio
import geopandas as gpd
from rasterio.features import rasterize
import numpy as np
from matplotlib import pyplot as plt
polygons_path = "/Users/maika/Desktop/preprocessing-geospatial-data/data/geometry/Cartagraphie des cultures.geojson"
gdf = gpd.read_file(polygons_path)

CLASSES=gdf["Renseigner la culture"].unique()
print(CLASSES)
CODES= {CLASSES[i].lower().replace(' - ','+').replace("é","e").replace("ï","i"):i+8 for i in range(len(CLASSES)) if CLASSES[i]!=None }

polygons_path = "/Users/maika/Desktop/preprocessing-geospatial-data/data/geometry/arachide_sample_4326.geojson"
gdf = gpd.read_file(polygons_path)

CLASSES=gdf["Crop"].unique()
print(CLASSES)
CODES2= {CLASSES[i].lower().replace(' - ','+').replace("é","e").replace("ï","i"):i+8 for i in range(len(CLASSES)) if CLASSES[i]!=None }
print(len(CODES2))
print(CODES2)
CODES2.update({k: len(CODES2) for k in CODES if k not in CODES2})
print(len(CODES2))
print(CODES2)


# Charger le raster avec rasterio
raster_path = "/Users/maika/Desktop/preprocessing-geospatial-data/data/geometry/arachide_sample_4326_dataset/0_S1_composite.tif"
raster_path="/Users/maika/Desktop/preprocessing-geospatial-data/data/geometry/arachide_sample_4326_dataset/0_dynamic_world.tif"
with rasterio.open(raster_path) as src:
    raster = src.read(1)  # Lire la première bande du raster
    out_meta = src.meta
    raster_shape = raster.shape
    print(raster_shape)
    plt.imshow(raster)
    plt.show()
    
exit()
# Charger les polygones avec Geopandas
polygons_path = "/Users/maika/Desktop/preprocessing-geospatial-data/data/geometry/arachide_sample_4326.geojson"

gdf = gpd.read_file(polygons_path)
print(gdf)


# Assigner les classes à chaque polygone
 # Ajuster les noms de culture selon vos données

# Créer une liste de tuples (geometry, value) pour rasterize
shapes = [(geom, value) for geom, value in zip(gdf.geometry, gdf['Code_crop'])]

# Créer un tableau numpy initialisé avec des valeurs de 0 (pour les zones hors des polygones)
mask = np.zeros(raster_shape, dtype=np.uint8)

# Rasterizer les polygones dans le masque
mask = rasterize(shapes=shapes, out_shape=raster_shape, fill=0, out=mask, transform=out_meta['transform'], dtype='uint8',)

# Sauvegarder le masque en tant que raster
output_path = "masque_classes.tif"
out_meta.update({"count": 1, "dtype": 'uint8','nodata':100})
print(out_meta)
with rasterio.open(output_path, "w", **out_meta) as dest:
    dest.write(mask, 1)

print(f"Le masque a été créé et sauvegardé dans {output_path}")


