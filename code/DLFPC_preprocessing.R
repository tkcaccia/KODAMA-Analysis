
seuList=list()

labels_list=NULL
xy_list=NULL
samples_list=NULL
spe_list=NULL
subjects_list=NULL
data_list=NULL
pca_list=NULL
for(i in 1:length(sample_names)){
  print(i)
  spe_list[[i]] <- spe[top[1:gene_number], colData(spe)$sample_id ==  sample_names[i]]
  spe_list[[i]] <- scater::runPCA(spe_list[[i]], 50,subset_row = top[1:gene_number], scale=TRUE)
  pca_list[[i]]=reducedDim(spe_list[[i]],type = "PCA")[,1:50]
  rownames(pca_list[[i]])=rownames(colData(spe_list[[i]]))
  labels_list[[i]]=as.factor(colData(spe_list[[i]])$layer_guess_reordered)
  names(labels_list[[i]])=rownames(colData(spe_list[[i]]))
  xy_list[[i]]=as.matrix(spatialCoords(spe_list[[i]]))
  rownames(xy_list[[i]])=rownames(colData(spe_list[[i]]))
  samples_list[[i]]=colData(spe_list[[i]])$sample_id
  names(samples_list[[i]])=rownames(colData(spe_list[[i]]))
  subjects_list[[i]]=colData(spe_list[[i]])$subject
  names(subjects_list[[i]])=rownames(colData(spe_list[[i]]))
  data_list[[i]]=t(logcounts(spe_list[[i]]))
  rownames(data_list[[i]])=rownames(colData(spe_list[[i]]))
  colnames(xy_list[[i]])=c("row","col")
  
  colnames(spe_list[[i]])=gsub("-1",paste("-",i,sep=""),colnames((spe_list[[i]])))
  colnames(xy_list[[i]])=gsub("-1",paste("-",i,sep=""),colnames((xy_list[[i]])))
  colnames(data_list[[i]])=gsub("-1",paste("-",i,sep=""),colnames((data_list[[i]])))
  
  
  
  seuList[[i]]  <- CreateSeuratObject(counts = spe_list[[i]]@assays@data$counts,
                                      meta.data = as.data.frame(colData(spe_list[[i]])))
  
  
  segmentations.data <- list(
    "centroids" = CreateCentroids(xy_list[[i]])
  )
  
  coords <- CreateFOV(
    coords = segmentations.data,
    type = "centroids",
    assay="RNA"
  )
  
  seuList[[i]] [['RNA']] <- AddMetaData(seuList[[i]] [['RNA']], metadata=as.data.frame( rowData(spe_list[[i]])))
  seuList[[i]] @tools$platform <- "Visium"
  
  colnames(seuList[[i]]@meta.data) [colnames(seuList[[i]]@meta.data) %in% c("array_row","array_col")]=c("row","col")
  
  coords <- subset(x = coords, cells = Cells(x = seuList[[i]]))
  seuList[[i]][["fov"]] <- coords
  
  
  PCA = CreateDimReducObject(embeddings = pca_list[[i]],  # should we choose larger number of dims
                             key = "Dimensions_", assay = "RNA")
  seuList[[i]]@reductions$pca = PCA
}
save(seuList,data_list,xy_list,labels_list,pca_list,top,samples_list,subject_names,subjects_list,file="data/DLPFC-general.RData")
