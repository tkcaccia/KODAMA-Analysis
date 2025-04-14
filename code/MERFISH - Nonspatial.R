
load("data/MERFISH-input.RData")

library(KODAMAextra)
library(igraph)
library(bluster)
library(Rtsne)


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
result_Walktrap$clusters_allslide = clu_walktrap


g_walk <- cluster_walktrap(g.PM)
clu_walktrap = as.character(igraph::cut_at(g_walk, no = ncluster))
result_Walktrap.PM$clusters_allslide = clu_walktrap
result_Walktrap.PM.ref$clusters_allslide=refine_SVM(xyz,clu_walktrap,cost=1000)


#Clustering with Leiden
result_Leiden$clusters_allslide = leiden(g,ncluster)$membership
result_Leiden.PM$clusters_allslide = leiden(g.PM,ncluster)$membership
result_Leiden.PM.ref$clusters_allslide=refine_SVM(xyz,result_Leiden.PM$clusters_allslide,cost=1000)


# Clustering with Louvain
result_Louvain$clusters_allslide = louvain(g,ncluster)$membership
result_Louvain.PM$clusters_allslide = louvain(g.PM,ncluster)$membership
result_Louvain.PM.ref$clusters_allslide=refine_SVM(xyz,result_Louvain.PM$clusters_allslide,cost=1000)




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
  result_Walktrap$clusters[[i]] <- result_Walktrap$clusters_allslide[sel]
  result_Leiden$clusters[[i]] <- result_Leiden$clusters_allslide[sel]
  result_Louvain$clusters[[i]] <- result_Louvain$clusters_allslide[sel]

  result_Walktrap.PM$clusters[[i]] <- result_Walktrap.PM$clusters_allslide[sel]
  result_Leiden.PM$clusters[[i]] <- result_Leiden.PM$clusters_allslide[sel]
  result_Louvain.PM$clusters[[i]] <- result_Louvain.PM$clusters_allslide[sel]

  result_Walktrap.PM.ref$clusters[[i]] <- result_Walktrap.PM.ref$clusters_allslide[sel]
  result_Leiden.PM.ref$clusters[[i]] <- result_Leiden.PM.ref$clusters_allslide[sel]
  result_Louvain.PM.ref$clusters[[i]] <- result_Louvain.PM.ref$clusters_allslide[sel]

  #############################################

  result_UMAP$feature_extraction[[i]]=u1[sel,]
  result_UMAP.PM$feature_extraction[[i]]=u2[sel,]
  result_tSNE$feature_extraction[[i]]=t1[sel,]
  result_tSNE.PM$feature_extraction[[i]]=t2[sel,]
  result_PCA$feature_extraction[[i]]=pca[sel,1:2]
  result_PCA.PM$feature_extraction[[i]]=pca.PM[sel,1:2]

}
result_Louvain$xyz=xyz
result_Leiden$xyz=xyz
result_Walktrap$xyz=xyz

result_Louvain.PM$xyz=xyz
result_Leiden.PM$xyz=xyz
result_Walktrap.PM$xyz=xyz

result_Louvain.PM.ref$xyz=xyz
result_Leiden.PM.ref$xyz=xyz
result_Walktrap.PM.ref$xyz=xyz

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
     file = "output/MERFISH-Nonspatial-results.RData")
