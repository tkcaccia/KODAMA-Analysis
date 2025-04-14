import os
import numpy as np
from PIL import Image
from spatialdata_io import visium
from spatialdata.dataloader.datasets import ImageTilesDataset
# Define regions and paths
regions = {
    "151672": "/media/user/Data/KODAMA/DLPFC/151672/",
    "151671": "/media/user/Data/KODAMA/DLPFC/151671/",
    "151670": "/media/user/Data/KODAMA/DLPFC/151670/",
    "151669": "/media/user/Data/KODAMA/DLPFC/151669/"
}
output_dir = "/media/user/Data/KODAMA/DLPFC/Tiles2"
os.makedirs(output_dir, exist_ok=True)

# Loop through each region
for region_name, region_path in regions.items():
    print(f"Processing region: {region_name}")
    # Load the SpatialData object
    sdata = visium(
        region_path,
        fullres_image_file=os.path.join(region_path, "spatial", f"{region_name}_full_image.tif")
    )
    # Create the dataset
    dataset = ImageTilesDataset(
        sdata=sdata,
        regions_to_images={region_name: f"{region_name}_full_image"},
        regions_to_coordinate_systems={region_name: "global"},
        table_name="table",
        tile_dim_in_units=500,
        rasterize=True,
        rasterize_kwargs={"target_width": 128},
    )
    # Get row names
    row_names = sdata["table"].obs.index
    # Save tiles for this region
    for idx, element in enumerate(dataset):
        try:
            # Access the image data
            image_dataarray = element.images[f"{region_name}_full_image"]
            image_np = image_dataarray.values  # Convert to NumPy array
            # Rearrange channels for saving as an RGB image
            image_np = np.transpose(image_np, (1, 2, 0))  # Shape: (height, width, 3)
            # Convert to uint8 if needed
            if image_np.dtype != np.uint8:
                image_np = (255 * (image_np - np.min(image_np)) / (np.max(image_np) - np.min(image_np))).astype(np.uint8)
            # Construct filename and save
            row_name = row_names[idx].replace("-1", "")
            output_filename = f"{row_name}-{region_name}.tiff"
            output_path = os.path.join(output_dir, output_filename)
            Image.fromarray(image_np).save(output_path, format="TIFF", compression="tiff_adobe_deflate")
            print(f"Saved image {idx + 1}/{len(dataset)} for {region_name} to {output_path}")
        except Exception as e:
            print(f"Error processing image {idx + 1} for {region_name}: {e}")
