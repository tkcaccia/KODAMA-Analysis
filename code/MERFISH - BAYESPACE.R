
# Load the preprocessed MERFISH data
load("data/MERFISH-input.RData")

# Load the necessary libraries
library(SingleCellExperiment)
library(BayesSpace)

result_BayesSpace=list(clusters=list(),feature_extraction=list(),tissue_segments=list(),xy=list())


# Perform BayesSpace clustering for each subject

  # Create SingleCellExperiment object
  info <- as.data.frame(xyz)
  colnames(info) <- c("row", "col", "z")
  x <- t(RNA)
  colnames(x) <- NULL
  sce <- SingleCellExperiment(assays = list(logcounts = x), colData = info)

  # Preprocess and cluster
  sce <- spatialPreprocess(sce, n.PCs = 20, n.HVGs = 100, log.normalize = FALSE)
  sce <- spatialCluster(sce, q = 8, d = 20,
                        init.method = "mclust", model = "t",
                        nrep = 10000, burn.in = 1000)



  umap_result <- umap::umap(reducedDims(sce)$PCA)$layout
  zlabels=colData(sce)$spatial.cluster

  z=xyz[,3]
  uz=unique(z)

  for(i in 1:5){
    sel=z==uz[i]
    result_BayesSpace$clusters[[i]]=zlabels[sel]
    result_BayesSpace$tissue_segments[[i]]=tissue_segments[sel]
    result_BayesSpace$xy[[i]]=xyz[sel,]
    result_BayesSpace$feature_extraction[[i]]=umap_result[sel,]
  }

  result_BayesSpace$xyz=xyz
  result_BayesSpace$clusters_allslide=zlabels
  result_BayesSpace$feature_extraction_allslide=umap_result
  result_BayesSpace$tissue_segments_allslide=tissue_segments


# Save BayesSpace results
save(result_BayesSpace, file = "output/MERFISH-BayesSpace-results.RData")


