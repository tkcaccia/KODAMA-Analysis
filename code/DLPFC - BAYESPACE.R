
library(BayesSpace)

load("data/DLPFC-general.RData")

result_BayesSpace=list(clusters=list(),
                        tissue_segments=list(),
                        feature_extraction=list(),
                        xy=xy_list,
                        samples=samples_list,
                        subjects=subjects_list)


for(i in 1:3){
  ncluster=sum(table(labels_subject[[i]])>0)


  info <- (xy_subject[[i]])
  colnames(info) <- c("row", "col")
  data=t(data_subject[[i]])
  colnames(data)=NULL
  sce <- SingleCellExperiment(assays=list(logcounts=data),
                              reducedDims=SimpleList(PCA=pca_subject[[i]]), colData = info)
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

  sce <- spatialCluster(sce, q = ncluster, d = 30,
                        init.method = "mclust", model = "t",
                        nrep = 10000, burn.in = 1000)

  umap_result <- umap::umap(reducedDims(sce)$PCA)$layout
  clu <- colData(sce)$spatial.cluster
  result_BayesSpace$clusters=c(result_BayesSpace$clusters,tapply(clu,samples_subject[[i]],function(x) x))
  result_BayesSpace$tissue_segments=c(result_BayesSpace$tissue_segments,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))
  result_BayesSpace$feature_extraction<-c(result_BayesSpace$feature_extraction,by(umap_result[,1:2],samples_subject[[i]],function(x) x))
}




save(result_BayesSpace,file="output/DLPFC-BayesSpace-results.RData")

