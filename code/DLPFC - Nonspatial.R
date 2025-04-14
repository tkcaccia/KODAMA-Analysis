load("data/DLPFC-general.RData")


library(igraph)
library(bluster)
library(Rtsne)
library(KODAMA)
library(KODAMAextra)


result_Walktrap<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Leiden<-list(clusters = list(),tissue_segments = list(),xy = list())
result_Louvain<-list(clusters = list(),tissue_segments = list(),xy = list())


result_UMAP<-list(tissue_segments = list(),feature_extraction = list())
result_tSNE<-list(tissue_segments = list(),feature_extraction = list())
result_PCA<-list(tissue_segments = list(),feature_extraction = list())



for(i in 1:3){


  ncluster=sum(table(labels_subject[[i]])>0)

  g <- makeSNNGraph(pca_subject[[i]],k = 20)

  # Wlaktrap
  g_walk <- cluster_walktrap(g)
  clu = as.character(igraph::cut_at(g_walk, no=ncluster))
  result_Walktrap$clusters=c(result_Walktrap$clusters,tapply(clu,samples_subject[[i]],function(x) x))

  print("Leiden algorithm")
  ##### Leiden algorithm
  clu=leiden(g,ncluster=ncluster)
  clu=clu$membership
  result_Leiden$clusters=c(result_Leiden$clusters,tapply(clu,samples_subject[[i]],function(x) x))
  print("Louvain algorithm")
  ##### Louvain algorithm

  clu=louvain(g,ncluster=ncluster)
  clu=clu$membership
  result_Louvain$clusters=c(result_Louvain$clusters,tapply(clu,samples_subject[[i]],function(x) x))


  result_Walktrap$tissue_segments=c(result_Walktrap$tissue_segments,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))
  result_Leiden$tissue_segments=c(result_Leiden$tissue_segments,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))
  result_Louvain$tissue_segments=c(result_Louvain$tissue_segments,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))

  result_Walktrap$xy=c(result_Walktrap$xy,by(xy_subject[[i]],samples_subject[[i]],function(x) x))
  result_Leiden$xy=c(result_Leiden$xy,by(xy_subject[[i]],samples_subject[[i]],function(x) x))
  result_Louvain$xy=c(result_Louvain$xy,by(xy_subject[[i]],samples_subject[[i]],function(x) x))



u1=umap::umap(pca_subject[[i]][,1:50])$layout
t1=Rtsne(pca_subject[[i]][,1:50],pca=FALSE)$Y

  result_UMAP$feature_extraction<-c(result_UMAP$feature_extraction,by(u1[,1:2],samples_subject[[i]],function(x) x))
  result_tSNE$feature_extraction<-c(result_tSNE$feature_extraction,by(t1[,1:2],samples_subject[[i]],function(x) x))
  result_PCA$feature_extraction<-c(result_PCA$feature_extraction,by(pca_subject[[i]][,1:2],samples_subject[[i]],function(x) x))

  result_UMAP$tissue_segments=c(result_UMAP$tissue_segments,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))
  result_tSNE$tissue_segments=c(result_tSNE$tissue_segments,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))
  result_PCA$tissue_segments=c(result_PCA$tissue_segments,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))

}


result_Walktrap$samples=samples_list
result_Walktrap$subjects=subjects_list


result_Leiden$samples=samples_list
result_Leiden$subjects=subjects_list

result_Louvain$samples=samples_list
result_Louvain$subjects=subjects_list

result_UMAP$samples=samples_list
result_tSNE$samples=samples_list
result_PCA$samples=samples_list

result_UMAP$subjects=subjects_list
result_tSNE$subjects=subjects_list
result_PCA$subjects=subjects_list


save(result_Walktrap,
     result_Leiden,
     result_Louvain,
     result_UMAP,
     result_tSNE,
     result_PCA,
     file="output/DLPFC-Nonspatial-results.RData")



