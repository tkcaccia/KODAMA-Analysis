
library(PRECAST)
library(mclust)

load("data/DLPFC-general.RData")

results_PRECAST=list()
results_PRECAST$clusters=list()
results_PRECAST$feature_extraction=list()
results_PRECAST$labels=list()

for(j in 1:3){

  ncluster=sum(table(labels_subject[[j]])>0)

  PRECASTObj <- CreatePRECASTObject(seuList[1:4+4*(j-1)], project = "DLPFC", gene.number = 2000,customGenelist=top[1:2000],
                                    premin.spots = 0, premin.features = 0, postmin.spots = 0, postmin.features = 0)
  PRECASTObj <- AddParSetting(PRECASTObj, maxIter=30, Sigma_equal = FALSE,
                              verbose = TRUE, int.model = "EEE", seed = 1)

  PRECASTObj <- AddAdjList(PRECASTObj, platform = "Visium")
  PRECASTObj <- PRECAST(PRECASTObj, K = ncluster)

  PRECASTObj <- SelectModel(PRECASTObj)

  seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
  seuInt <- AddTSNE(seuInt, n_comp = 2)

  slide=seuInt@meta.data$batch
  for(i in 1:4){
    results_PRECAST$clusters[[i+4*(j-1)]]=as.numeric(PRECASTObj@resList$cluster[[i]])
    results_PRECAST$labels[[i+4*(j-1)]]=PRECASTObj@seulist[[i]]$layer_guess_reordered
    results_PRECAST$feature_extraction[[i+4*(j-1)]]=seuInt@reductions$tSNE@cell.embeddings[slide==i,]
  }
}
#########################################################################################3

save(results_PRECAST,file="output/PRECAST-results.RData")
