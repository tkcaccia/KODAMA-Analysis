pca_subject=list()
pca_subject[[1]]=pca_Br5292
pca_subject[[2]]=pca_Br5595
pca_subject[[3]]=pca_Br8100
names(pca_subject)=c("Br5292","Br5595","Br8100")

labels_subject=list()
labels_subject[[1]]=labels_Br5292
labels_subject[[2]]=labels_Br5595
labels_subject[[3]]=labels_Br8100
names(labels_subject)=c("Br5292","Br5595","Br8100")

samples_subject=list()
samples_subject[[1]]=samples_Br5292
samples_subject[[2]]=samples_Br5595
samples_subject[[3]]=samples_Br8100
names(samples_subject)=c("Br5292","Br5595","Br8100")

xy_subject=list()
xy_subject[[1]]=xy_Br5292
xy_subject[[2]]=xy_Br5595
xy_subject[[3]]=xy_Br8100
names(xy_subject)=c("Br5292","Br5595","Br8100")

data_subject=list()
data_subject[[1]]=data_Br5292
data_subject[[2]]=data_Br5595
data_subject[[3]]=data_Br8100
names(data_subject)=c("Br5292","Br5595","Br8100")



seuList=list()

labels_list=NULL
xy_list=NULL
samples_list=NULL
spe_list=NULL
subjects_list=NULL
data_list=NULL
pca_list=NULL
ma=0
for(i in 1:length(sample_names)){
  print(i)
  spe_list[[i]] <- spe[top[1:gene_number], colData(spe)$sample_id ==  sample_names[i]]
  spe_list[[i]] <- scater::runPCA(spe_list[[i]], 50,subset_row = top[1:gene_number], scale=TRUE)
  pca_list[[i]]=reducedDim(spe_list[[i]],type = "PCA")[,1:50]
  rownames(pca_list[[i]])=rownames(colData(spe_list[[i]]))
  labels_list[[i]]=as.factor(colData(spe_list[[i]])$layer_guess_reordered)
  names(labels_list[[i]])=rownames(colData(spe_list[[i]]))

#Orizontalization of the slides
  xy_temp=as.matrix(spatialCoords(spe_list[[i]]))
  xy_temp[,1]=xy_temp[,1]+ma
  ran=range(xy_temp[, 1])
  ma=ran[2]+ dist(ran)[1]*0.5
  xy_list[[i]]=xy_temp
  
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
save(pca_subject,labels_subject,samples_subject,xy_subject,data_subject,spe_list,
     seuList,data_list,xy_list,labels_list,pca_list,top,samples_list,subject_names,subjects_list,file="data/DLPFC-general.RData")
