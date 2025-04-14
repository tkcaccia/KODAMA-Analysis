
library(PRECAST)
library(mclust)

load("data/DLPFC-general.RData")

result_PRECAST=list(clusters=list(),
                   tissue_segments=list(),
                   feature_extraction=list(),
                   xy=list(),
                   samples=list(),
                   subjects=list())
for(j in 1:3){

  ncluster=sum(table(labels_subject[[j]])>0)

  PRECASTObj <- CreatePRECASTObject(seuList[1:4+4*(j-1)], project = "DLPFC", gene.number = 2000,customGenelist=top[1:2000],
                                    premin.spots = 0, premin.features = 0, postmin.spots = 0, postmin.features = 0)

  identifiers <- rownames(PRECASTObj@seulist[[1]]@meta.data)

  PRECASTObj <- AddParSetting(PRECASTObj, maxIter=30, Sigma_equal = FALSE,
                              verbose = TRUE, int.model = "EEE", seed = 1)

  PRECASTObj <- AddAdjList(PRECASTObj, platform = "Visium")
  PRECASTObj <- PRECAST(PRECASTObj, K = ncluster)

  PRECASTObj <- SelectModel(PRECASTObj)

  seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
  seuInt <- AddTSNE(seuInt, n_comp = 2)

  slide=seuInt@meta.data$batch
  for(i in 1:4){
    result_PRECAST$clusters[[i+4*(j-1)]]=as.numeric(PRECASTObj@resList$cluster[[i]])
    result_PRECAST$tissue_segments[[i+4*(j-1)]]=PRECASTObj@seulist[[i]]$layer_guess_reordered
    result_PRECAST$feature_extraction[[i+4*(j-1)]]=seuInt@reductions$tSNE@cell.embeddings[slide==i,]
    result_PRECAST$xy[[i+4*(j-1)]]=PRECASTObj@seulist[[i]]@meta.data[,c("col","row")]

    result_PRECAST$samples[[i+4*(j-1)]]=PRECASTObj@seulist[[i]]$sample_id
    result_PRECAST$subjects[[i+4*(j-1)]]=PRECASTObj@seulist[[i]]$subject

  }
}
#########################################################################################3




save(result_PRECAST,file="output/DLPFC-PRECAST-results.RData")
