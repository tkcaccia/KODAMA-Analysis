load("data/DLPFC-general.RData")


library(igraph)
library(bluster)


results_Walktrap=list()
results_Walktrap$clusters=list()
results_Walktrap$labels=list()
results_Leiden=list()
results_Leiden$clusters=list()
results_Leiden$labels=list()
results_Louvain=list()
results_Louvain$clusters=list()
results_Louvain$labels=list()

for(i in 1:3){
  ncluster=sum(table(labels_subject[[i]])>0)

  g <- makeSNNGraph(pca_subject[[i]],k = 20)

  # Wlaktrap
  g_walk <- cluster_walktrap(g)
  clu = as.character(igraph::cut_at(g_walk, no=ncluster))
  results_Walktrap$clusters=c(results_Walktrap$clusters,tapply(clu,samples_subject[[i]],function(x) x))


  ##### Leiden algorithm
  t=Inf
  res=0.16
  while(ncluster!=t){
    clu=cluster_leiden(g,resolution_parameter=res)
    res=res+0.001
    t=clu$nb_clusters
    print(t)
  }
  clu=clu$membership
  results_Leiden$clusters=c(results_Leiden$clusters,tapply(clu,samples_subject[[i]],function(x) x))

  ##### Louvain algorithm
  t=Inf
  res=0.05
  while(ncluster!=t){
    clu=cluster_louvain(g,resolution=res)
    res=res+0.05
    t=length(unique(clu$membership))
    print(t)
  }
  clu=clu$membership
  results_Louvain$clusters=c(results_Louvain$clusters,tapply(clu,samples_subject[[i]],function(x) x))


  results_Walktrap$labels=c(results_Walktrap$labels,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))
  results_Leiden$labels=c(results_Leiden$labels,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))
  results_Louvain$labels=c(results_Louvain$labels,tapply(labels_subject[[i]],samples_subject[[i]],function(x) x))

}
save(results_Walktrap,results_Leiden,results_Louvain,file="output/Nonspatial-results.RData")



