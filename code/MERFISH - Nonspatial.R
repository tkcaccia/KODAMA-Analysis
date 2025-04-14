library(igraph)
library(bluster)

load("data/MERFISH-input.RData")

# Clustering with Walktrap, Leiden, and Louvain
Walktrap <- list()
Leiden <- list()
Louvain <- list()
result_Walktrap<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Leiden<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Louvain<-list(clusters = list(),tissue_segments = list(),xy = list())

result_Walktrap.PM<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Leiden.PM<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Louvain.PM<-list(clusters = list(),tissue_segments = list(),xy = list())

result_Walktrap.PM.ref<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Leiden.PM.ref<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Louvain.PM.ref<-list(clusters = list(),tissue_segments = list(),xy = list())

result_UMAP.PM<-list(tissue_segments = list(),feature_extraction = list())
result_tSNE.PM<-list(tissue_segments = list(),feature_extraction = list())
result_PCA.PM<-list(tissue_segments = list(),feature_extraction = list())

result_UMAP<-list(tissue_segments = list(),feature_extraction = list())
result_tSNE<-list(tissue_segments = list(),feature_extraction = list())
result_PCA<-list(tissue_segments = list(),feature_extraction = list())



u1=umap::umap(pca[,1:20])$layout
u2=umap::umap(pca.PM[,1:20])$layout
t1=Rtsne(pca[,1:20],pca=FALSE)$Y
t2=Rtsne(pca.PM[,1:20],pca=FALSE)$Y

ncluster = sum(table(tissue_segments) > 0)

# K-nearest neighbors graph
g <- makeSNNGraph(pca[,1:20], k = 20)
g.PM <- makeSNNGraph(pca.PM[,1:20], k = 20)

# Clustering with Walktrap
g_walk <- cluster_walktrap(g)
clu_walktrap = as.character(igraph::cut_at(g_walk, no = ncluster))
Walktrap$clusters = clu_walktrap


g_walk <- cluster_walktrap(g.PM)
clu_walktrap = as.character(igraph::cut_at(g_walk, no = ncluster))
Walktrap$clusters.PM = clu_walktrap
Walktrap$ref=refine_SVM(xyz,clu_walktrap,cost=1000)


#Clustering with Leiden
t=Inf
res=0.012
while(ncluster!=t){
  clu_leiden = cluster_leiden(g, resolution_parameter = res)
  res = res + 0.001
  t = clu_leiden$nb_clusters
  print(t)
}
clu_leiden = clu_leiden$membership
Leiden$clusters = clu_leiden
#Clustering with Leiden
t=Inf
res=0.012
while(ncluster!=t){
  clu_leiden = cluster_leiden(g.PM, resolution_parameter = res)
  res = res + 0.001
  t = clu_leiden$nb_clusters
  print(t)
}
clu_leiden = clu_leiden$membership
Leiden$clusters.PM = clu_leiden
Leiden$ref=refine_SVM(xyz,clu_leiden,cost=1000)


# Clustering with Louvain
res = 0.05
while (ncluster != length(unique(clu_louvain <- cluster_louvain(g, resolution = res)$membership))) {
  res = res + 0.07
}
Louvain$clusters = clu_louvain

res = 0.05
while (ncluster != length(unique(clu_louvain <- cluster_louvain(g.PM, resolution = res)$membership))) {
  res = res + 0.07
}
Louvain$clusters.PM = clu_louvain
Louvain$ref=refine_SVM(xyz,clu_louvain,cost=1000)



z=xyz[,3]
uz=unique(z)
for(i in 1:5){
  sel=z==uz[i]
  result_Walktrap$tissue_segments[[i]]=tissue_segments[sel]
  result_Leiden$tissue_segments[[i]]=tissue_segments[sel]
  result_Louvain$tissue_segments[[i]]=tissue_segments[sel]
  result_Walktrap.PM$tissue_segments[[i]]=tissue_segments[sel]
  result_Leiden.PM$tissue_segments[[i]]=tissue_segments[sel]
  result_Louvain.PM$tissue_segments[[i]]=tissue_segments[sel]
  result_Walktrap.PM.ref$tissue_segments[[i]]=tissue_segments[sel]
  result_Leiden.PM.ref$tissue_segments[[i]]=tissue_segments[sel]
  result_Louvain.PM.ref$tissue_segments[[i]]=tissue_segments[sel]

  result_UMAP$tissue_segments[[i]]=tissue_segments[sel]
  result_tSNE$tissue_segments[[i]]=tissue_segments[sel]
  result_PCA$tissue_segments[[i]]=tissue_segments[sel]
  result_UMAP.PM$tissue_segments[[i]]=tissue_segments[sel]
  result_tSNE.PM$tissue_segments[[i]]=tissue_segments[sel]
  result_PCA.PM$tissue_segments[[i]]=tissue_segments[sel]


  result_Walktrap$xy[[i]]=xyz[sel,-3]
  result_Leiden$xy[[i]]=xyz[sel,-3]
  result_Louvain$xy[[i]]=xyz[sel,-3]
  result_Walktrap.PM$xy[[i]]=xyz[sel,-3]
  result_Leiden.PM$xy[[i]]=xyz[sel,-3]
  result_Louvain.PM$xy[[i]]=xyz[sel,-3]
  result_Walktrap.PM.ref$xy[[i]]=xyz[sel,-3]
  result_Leiden.PM.ref$xy[[i]]=xyz[sel,-3]
  result_Louvain.PM.ref$xy[[i]]=xyz[sel,-3]
#############################################
  result_Walktrap$clusters[[i]] <- Walktrap$clusters[sel]
  result_Leiden$clusters[[i]] <- Leiden$clusters[sel]
  result_Louvain$clusters[[i]] <- Louvain$clusters[sel]

  result_Walktrap.PM$clusters[[i]] <- Walktrap$clusters.PM[sel]
  result_Leiden.PM$clusters[[i]] <- Leiden$clusters.PM[sel]
  result_Louvain.PM$clusters[[i]] <- Louvain$clusters.PM[sel]

  result_Walktrap.PM.ref$clusters[[i]] <- Walktrap$ref[sel]
  result_Leiden.PM.ref$clusters[[i]] <- Leiden$ref[sel]
  result_Louvain.PM.ref$clusters[[i]] <- Louvain$ref[sel]

  #############################################

  result_UMAP$feature_extraction[[i]]=u1[sel,]
  result_UMAP.PM$feature_extraction[[i]]=u2[sel,]
  result_tSNE$feature_extraction[[i]]=t1[sel,]
  result_tSNE.PM$feature_extraction[[i]]=t2[sel,]
  result_PCA$feature_extraction[[i]]=pca[sel,1:2]
  result_PCA.PM$feature_extraction[[i]]=pca.PM[sel,1:2]

}
result_Louvain$xyz=xyz
result_Louvain$clusters_allslide=Louvain$clusters

result_Leiden$xyz=xyz
result_Leiden$clusters_allslide=Leiden$clusters

result_Walktrap$xyz=xyz
result_Walktrap$clusters_allslide=Walktrap$clusters

result_Louvain.PM$xyz=xyz
result_Louvain.PM$clusters_allslide=Louvain$clusters.PM

result_Leiden.PM$xyz=xyz
result_Leiden.PM$clusters_allslide=Leiden$clusters.PM

result_Walktrap.PM$xyz=xyz
result_Walktrap.PM$clusters_allslide=Walktrap$clusters.PM

result_Louvain.PM.ref$xyz=xyz
result_Louvain.PM.ref$clusters_allslide=Louvain$ref

result_Leiden.PM.ref$xyz=xyz
result_Leiden.PM.ref$clusters_allslide=Leiden$ref

result_Walktrap.PM.ref$xyz=xyz
result_Walktrap.PM.ref$clusters_allslide=Walktrap$ref

result_UMAP$feature_extraction_allslide=u1
result_UMAP.PM$feature_extraction_allslide=u2
result_tSNE$feature_extraction_allslide=t1
result_tSNE.PM$feature_extraction_allslide=t2
result_PCA$feature_extraction_allslide=pca[,1:2]
result_PCA.PM$feature_extraction_allslide=pca.PM[,1:2]

result_UMAP$tissue_segments_allslide=tissue_segments
result_UMAP.PM$tissue_segments_allslide=tissue_segments
result_tSNE$tissue_segments_allslide=tissue_segments
result_tSNE.PM$tissue_segments_allslide=tissue_segments
result_PCA$tissue_segments_allslide=tissue_segments
result_PCA.PM$tissue_segments_allslide=tissue_segments

###############################################################################


# Saving the results
save(result_Walktrap, result_Leiden, result_Louvain,
     result_Walktrap.PM, result_Leiden.PM, result_Louvain.PM,
     result_Walktrap.PM.ref, result_Leiden.PM.ref, result_Louvain.PM.ref,
     result_UMAP,result_UMAP.PM,result_tSNE,result_tSNE.PM,result_PCA,result_PCA.PM,
     file = "output/MERFISH-NONSPATIAL-results.RData")
