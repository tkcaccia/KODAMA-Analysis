
load("data/MERFISH-input.RData")

library(Banksy)
library(SpatialExperiment)



result_BANKSY=list(clusters=list(),feature_extraction=list(),tissue_segments=list(),xy=list())


spe <- SpatialExperiment(
  assay = list(counts=t(RNA)),
  colData = data.frame(tissue_segments),
  spatialCoords = xyz)


compute_agf <- FALSE
k_geom <- 6
spe_banksy <-  computeBanksy(spe, assay_name = "counts",
                             compute_agf = compute_agf, k_geom = k_geom)

use_agf=FALSE
lambda=0.2
res <- 0.7

spe_banksy <- runBanksyPCA(spe_banksy,seed = 543210,npcs = 20)
spe_banksy <- runBanksyUMAP(spe_banksy, use_agf = use_agf, lambda = lambda, seed = 543210)

cnm_umap= sprintf("UMAP_M%s_lam%s", as.numeric(use_agf), lambda)
dimred=reducedDim(spe_banksy,cnm_umap)

ncluster=8

t=Inf
res=0.05
while(ncluster!=t){
  res=res+0.005

  spe_banksy_2 <- clusterBanksy(spe_banksy, use_agf = use_agf, lambda = lambda, resolution = res, seed = 543210)
  cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
  cnm_smooth= sprintf("clust_M%s_lam%s_k50_res%s_smooth", as.numeric(use_agf), lambda, res)

  spe_banksy_2 <- connectClusters(spe_banksy_2,map_to=cnm)
  spe_banksy_2 = smoothLabels(spe_banksy_2,coord_names =xyz,cluster_names = cnm, k = 6L, verbose = FALSE)
  cluster=spe_banksy_2[[cnm_smooth]]

  t=length(levels(cluster))
  print(t)
}



zlabels=cluster

z=xyz[,3]
uz=unique(z)



for(i in 1:5){
  sel=z==uz[i]
  result_BANKSY$clusters[[i]]=zlabels[sel]
  result_BANKSY$tissue_segments[[i]]=tissue_segments[sel]
  result_BANKSY$xy[[i]]=xyz[sel,-3]
  result_BANKSY$feature_extraction[[i]]=dimred[sel,]
}

result_BANKSY$xyz=xyz
result_BANKSY$clusters_allslide=zlabels
result_BANKSY$feature_extraction_allslide=dimred
result_BANKSY$tissue_segments_allslide=tissue_segments



save(result_BANKSY,file="output/MERFISH-BANSKY-results.RData")
############################################################################################3



