
from shapely.geometry import Point
from joblib import Parallel, delayed

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Point, shape
import geopandas as gpd
from spatialdata_io import visium
from spatialdata_io import visium_hd
import pandas as pd
import spatialdata as sd
from pathlib import Path

# Define the function to check if a point is within the mask
def check_point_within(row, mask_geom):
  return mask_geom.contains(Point(row['x'], row['y']))


sdata = visium_hd('/scratch/01481067/Colorectal/outs')

shapes = sdata.shapes["Visium_HD_Human_Colon_Cancer_square_008um"]



geojson_path = r'/scratch/01481067/KODAMA-Analysis/data/Annotations/Visium_HD_Human_Colon_Cancer_290325.geojson'  # Replace with your actual path to the GeoJSON file
geojson_data = gpd.read_file(geojson_path)

mask_geometry = geojson_data.union_all()  # Assuming the mask is a single unified geometry

# Load the low-resolution image
lowres_image = sdata.images["Visium_HD_Human_Colon_Cancer_lowres_image"]
image = lowres_image.values.transpose((1, 2, 0))
x_min, y_min, x_max, y_max = shapes.total_bounds





points_df = shapes[['geometry']].copy()

points_df['x'] = points_df.geometry.centroid.x
points_df['y'] = points_df.geometry.centroid.y

points_df['spot_id'] = points_df.index


# Determine which spots are within the mask
#points_within_mask = points_df.apply(lambda row: mask_geometry.contains(Point(row['x'], row['y'])), axis=1)
# Use joblib to parallelize
spots_in_mask = Parallel(n_jobs=-1, prefer="threads")(
  delayed(check_point_within)(row, mask_geometry) for _, row in points_df.iterrows()
)

# Add result back to DataFrame
points_df['within_mask'] = spots_in_mask

expression_data = sdata.tables["square_008um"]



# Map the spot indices in mask to gene expression data
#spot_indices_in_mask = spots_in_mask.index
spot_indices_in_mask = [i for i, val in enumerate(spots_in_mask) if val]

expression_data_in_mask = expression_data[spot_indices_in_mask, :]

# Create DataFrame for gene expression data
expression_df = pd.DataFrame(expression_data_in_mask.X.toarray(), index=expression_data_in_mask.obs_names, columns=expression_data_in_mask.var_names)

# Add sequences (barcodes) to spots_in_mask
#spots_in_mask.loc[:, 'sequence'] = expression_data.obs_names[spot_indices_in_mask]
# Get the index (e.g., spot names) from AnnData
sequence_names = expression_data.obs_names[spots_in_mask]

# Assign to the DataFrame (only the rows within the mask)
points_df.loc[spots_in_mask, 'sequence'] = sequence_names


# Extract classifications from GeoJSON data
classifications = geojson_data[['geometry', 'classification']]  # Assuming 'classification' is the relevant field

# Merge classifications with spots data
#spots_in_mask_gdf = gpd.GeoDataFrame(spots_in_mask, geometry=gpd.points_from_xy(spots_in_mask.x, spots_in_mask.y))
# Subset the DataFrame to spots within the mask
spots_in_mask_df = points_df[spots_in_mask].copy()

# Now create a GeoDataFrame from that
spots_in_mask_gdf = gpd.GeoDataFrame(
  spots_in_mask_df,
  geometry=gpd.points_from_xy(spots_in_mask_df.x, spots_in_mask_df.y)
)




# Ensure there are no conflicting column names before the join
spots_in_mask_gdf = spots_in_mask_gdf.rename(columns={'spot_id': 'spot_id_left'})
classifications = classifications.rename(columns={'spot_id': 'spot_id_right'})

# Perform spatial join
spots_with_classification = gpd.sjoin(spots_in_mask_gdf, classifications, how="inner", predicate='intersects')

# Select relevant columns and rename spot_id to sequence
spots_with_classification = spots_with_classification[['sequence', 'classification']]

# Save to CSV
output_csv_path = r'/scratch/01481067/KODAMA-Analysis/data/Annotations/spots_classification_VisiumHD.csv'  # Replace with your desired output path
spots_with_classification.to_csv(output_csv_path, index=False)


