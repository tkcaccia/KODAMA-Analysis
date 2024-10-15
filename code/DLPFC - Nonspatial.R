

library(igraph)
library(bluster)


results_Walktrap=list()
results_Walktrap$cluster=list()
results_Leiden=list()
results_Leiden$cluster=list()
results_Louvain=list()
results_Louvain$cluster=list()

for(i in 1:3){
  ncluster=sum(table(labels_subject[[i]])>0)
  
  g <- makeSNNGraph(pca_subject[[i]],k = 20)
  
  # Wlaktrap
  g_walk <- cluster_walktrap(g)
  clu = as.character(igraph::cut_at(g_walk, no=ncluster))
  results_Walktrap$cluster=c(results_Walktrap$cluster,tapply(clu,samples[[i]],function(x) x))
  
  
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
  results_Leiden$cluster=c(results_Leiden$cluster,tapply(clu,samples[[i]],function(x) x))
  
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
  results_Louvain$cluster=c(results_Louvain$cluster,tapply(clu,samples[[i]],function(x) x))

}
save(results_Walktrap,results_Leiden,results_Louvain,file="output/Nonspatial-results.RData")



