
library(BayesSpace)

load("data/DLPFC-general.RData")

results_BayesSpace=list()
results_BayesSpace$clusters=list()
results_BayesSpace$labels=list()
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
  
  
  
  clu <- colData(sce)$spatial.cluster
  
  results_BayesSpace$clusters=c(results_BayesSpace$clusters,tapply(clu,samples_subject[[i]],function(x) x))
  results_BayesSpace$labels=c(results_BayesSpace$labels,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))

}
save(results_BayesSpace,file="output/BayesSpace-results.RData")

