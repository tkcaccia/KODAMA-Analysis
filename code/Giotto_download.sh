mkdir Giotto_Mouse_brain
cd Giotto_Mouse_brain

# Input Files
curl -O https://cg.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_image.tif

# Output Files
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_raw_feature_bc_matrix.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_spatial.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_cloupe.cloupe

tar -xvzf V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz
tar -xvzf V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz
tar -xvzf V1_Adult_Mouse_Brain_spatial.tar.gz

rm V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz
rm V1_Adult_Mouse_Brain_raw_feature_bc_matrix.tar.gz
rm V1_Adult_Mouse_Brain_spatial.tar.gz
