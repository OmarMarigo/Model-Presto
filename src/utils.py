import os
import rasterio
import numpy as np

# Dossier contenant vos fichiers TIF
input_directory = "/Users/mac/Desktop/docko/tolbi-next-gen/preprocessing-geospatial-data/data/cartagraphie_des_cultures_dataset_2024-07-01_2024-12-01"

# Liste des bandes Sentinel-2 L1C
ALL_BANDS_S2_L1C = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']

# Initialiser des dictionnaires pour stocker les sommes, les comptes et les sommes des carrés des pixels pour chaque bande
band_sums = {band: 0 for band in ALL_BANDS_S2_L1C}
band_counts = {band: 0 for band in ALL_BANDS_S2_L1C}
band_sums_squared = {band: 0 for band in ALL_BANDS_S2_L1C}



# Parcourir chaque fichier dans le dossier
for filename in os.listdir(input_directory):
    if filename.endswith("S2_composite.tif"):
        filepath = os.path.join(input_directory, filename)
        
        # Ouvrir le fichier TIF avec rasterio
        with rasterio.open(filepath) as src:
            # Parcourir chaque bande
            for band_id in range(1, src.count + 1):
                band_name = ALL_BANDS_S2_L1C[band_id - 1]
                
                # Lire les valeurs de la bande
                band = src.read(band_id)
                
                # Masquer les valeurs no-data
                band = np.ma.masked_equal(band, src.nodata)
                
                # Ajouter la somme, le nombre de pixels valides et la somme des carrés
                band_sums[band_name] += band.sum()
                band_counts[band_name] += band.count()
                band_sums_squared[band_name] += (band ** 2).sum()

# Calculer la moyenne et l'écart-type pour chaque bande
band_means = {}
band_stds = {}
for band_name in ALL_BANDS_S2_L1C:
    if band_counts[band_name] > 0:
        band_mean = band_sums[band_name] / band_counts[band_name]
        band_variance = (band_sums_squared[band_name] / band_counts[band_name]) - (band_mean ** 2)
        band_std = np.sqrt(band_variance)
        band_means[band_name] = band_mean
        band_stds[band_name] = band_std
    else:
        band_means[band_name] = None
        band_stds[band_name] = None

# Afficher le dictionnaire des moyennes et l'écart-type moyen
print("Moyennes par bande:", band_means)
print("Écart-type moyen par bande:", band_stds)

# Calculer et afficher l'écart-type des moyennes
mean_values = [value for value in band_means.values() if value is not None]
print(f"mean des moyennes: {mean_values}")
if mean_values:
    mean_std = np.std(mean_values)
    print(f"Écart-type des moyennes: {mean_std}")
else:
    print("Écart-type des moyennes: No valid data")


