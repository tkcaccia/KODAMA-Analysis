
load("data/MERFISH-input.RData")
#save(ref,tissue_segments,vis,xyz,pca,pca.PM,file=

rownames(xyz)=rownames(RNA)
colnames(xyz)=c("x","y","z")
xyz.list=list(xyz)

library("BASS")
# Initialize the KODAMA list
result_BASS =list(clusters=list(),tissue_segments=list(),xy=list())


pca_harmony <- RunHarmony(pca, data.frame(z=xyz[,3]),"z",lambda=NULL)
rownames(pca_harmony)=rownames(RNA)



# hyper-parameters
# We set the number of cell types to a relatively large
# number (20) to capture the expression heterogeneity.
C <- 20
# number of spatial domains
R <- 8

set.seed(0)
# Set up BASS object
RNA.list=list(t(RNA))

res <- createBASSObject(RNA.list, xyz.list, C = C, R = R, beta_method = "SW",)


res <- BASS.preprocess(res, doLogNormalize = FALSE, doPCA = FALSE, scaleFeature = FALSE, nPC = 30)

res@X_run=t(pca_harmony[,1:30])

res <- BASS.run(res)
res <- BASS.postprocess(res)
zlabels <- res@results$z[[1]]

result_BASS=list()

z=xyz[,3]
uz=unique(z)
result_BASS$clusters=list()
result_BASS$tissue_segments=list()
result_BASS$xy=list()
for(i in 1:5){
  sel=z==uz[i]
  result_BASS$clusters[[i]]=zlabels[sel]
  result_BASS$tissue_segments[[i]]=tissue_segments[sel]
  result_BASS$xy[[i]]=xyz[sel,-3]
}

result_BASS$xyz=xyz
result_BASS$clusters_allslide=zlabels
result_BASS$tissue_segments_allslide=tissue_segments


save(result_BASS, file = "output/MERFISH-BASS-results.RData")

