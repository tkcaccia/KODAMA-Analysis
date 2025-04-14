load("data/DLPFC-general.RData")

library(SpatialExperiment)
library(S4Vectors)
library(Banksy)



result_BANKSY=list(clusters=list(),
                   tissue_segments=list(),
                   feature_extraction=list(),
                   xy=xy_list,
                   samples=samples_list,
                   subjects=subjects_list)

compute_agf <- FALSE
k_geom <- 6

lambda <- 0.2
use_agf <- FALSE

for(j in 1:3){
  ncluster=sum(table(labels_subject[[j]])>0)
  spe_list_banksy <- lapply(spe_list[1:4+4*(j-1)], computeBanksy, assay_name = "logcounts",
                            compute_agf = compute_agf, k_geom = k_geom)


  spe_joint <- do.call(cbind, spe_list_banksy)

  invisible(gc())
  spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000,npcs = 50)
  spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)



  t=Inf
  res=0.1
  while(ncluster!=t){
    res=res+0.01
    spe_joint2=spe_joint
    spe_joint2 <- clusterBanksy(spe_joint2, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
    cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)
    #spe_joint <- connectClusters(spe_joint,map_to="clust_M0_lam0.2_k50_res0.7")
    spe_joint2 <- connectClusters(spe_joint2,map_to=cnm)

    t=length(unique(spe_joint2[[cnm]]))
    print(t)
  }




  slide=spe_joint2$sample_id
  slide_names=unique(slide)

  spe_list_banksy <- lapply(slide_names, function(x) spe_joint2[, spe_joint2$sample_id == x])

  invisible(gc())



  smooth <- sprintf("clust_M%s_lam%s_k50_res%s_smooth", as.numeric(use_agf), lambda, res)

  for(i in 1:4){
    spe_list_banksy[[i]]=smoothLabels(spe_list_banksy[[i]],coord_names = xy_list[[i+4*(j-1)]],cluster_names = cnm, k = 6L, verbose = FALSE)
    result_BANKSY$clusters[[i+4*(j-1)]]=spe_list_banksy[[i]][[smooth]]
    result_BANKSY$tissue_segments[[i+4*(j-1)]]=spe_list_banksy[[i]]$layer_guess_reordered
    result_BANKSY$feature_extraction[[i+4*(j-1)]]=reducedDim(spe_joint2,type = "UMAP_M0_lam0.2")[slide==slide_names[i],]
  }

}
############################################################################################3



 save(result_BANKSY,file="output/DLPFC-BANSKY-results.RData")

