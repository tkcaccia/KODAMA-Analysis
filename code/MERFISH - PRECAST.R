
# Load necessary packages
library(ProFAST)
library(PRECAST)

# Initialize PRECAST clusters list
result_PRECAST=list(clusters=list(),feature_extraction=list(),tissue_segments=list(),xy=list())

# Load the required data
load("data/MERFISH-input.RData")


# Combine data from all slides
data <- NULL
for (i in as.character(slides)) {
  data <- cbind(data, cnts_mult[[i]])
}


names(tissue_segments) <-  rownames(RNA)
colnames(xyz) <- c("col", "row", "z")
rownames(xyz) <- rownames(RNA)

# Create a Seurat object
library(Seurat)
seu <- CreateSeuratObject(counts = t(RNA),meta.data = data.frame(xyz))

# Set a seed for reproducibility
set.seed(543210)



# Convert the Seurat object into a list
seuList <- list(seu)


# Create a PRECAST object
PRECASTObj <- CreatePRECASTObject(seuList, project = "MERFISH", gene.number = 100,
                                  selectGenesMethod = "SPARK-X",
                                  premin.spots = 20, premin.features = 20,
                                  postmin.spots = 1, postmin.features = 10)


identifiers <- rownames(PRECASTObj@seulist[[1]]@meta.data)
colnames(xyz) <- c("x", "y", "z")
posList <- list(xyz[identifiers,])

# Extract 3D coordinates and compute the adjacency matrix

PRECASTObj@AdjList <- pbapply::pblapply(posList, AddAdj, neighbors = 15)

# Set model parameters
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, verbose = TRUE, maxIter = 30)

# Run the PRECAST clustering algorithm
PRECASTObj <- PRECAST(PRECASTObj, K = 8)

# Select the optimal model
PRECASTObj <- SelectModel(PRECASTObj)


seuInt <- IntegrateSpaData(PRECASTObj, species = "Mouse")
seuInt <- AddTSNE(seuInt, n_comp = 2)

PRECASTtSNE=seuInt@reductions$tSNE@cell.embeddings

zlabels=as.numeric(PRECASTObj@resList$cluster[[1]])
names(zlabels)=rownames(PRECASTtSNE)
result_PRECAST=list()

z=xyz[,3]
uz=unique(z)
result_PRECAST$clusters=list()
result_PRECAST$tissue_segments=list()
result_PRECAST$xy=list()
result_PRECAST$feature_extraction=list()
for(i in 1:5){
  sel=rownames(xyz)[z==uz[i] & (rownames(xyz) %in% identifiers)]
  result_PRECAST$clusters[[i]]=zlabels[sel]
  result_PRECAST$tissue_segments[[i]]=tissue_segments[sel]
  result_PRECAST$xy[[i]]=xyz[sel,]
  result_PRECAST$feature_extraction[[i]]=PRECASTtSNE[sel,]
}

result_PRECAST$xyz=xyz[identifiers,]
result_PRECAST$clusters_allslide=zlabels
result_PRECAST$tissue_segments_allslide=tissue_segments[identifiers]
result_PRECAST$feature_extraction_allslide=PRECASTtSNE[identifiers,]

save(result_PRECAST, file = "output/MERFISH-PRECAST-results.RData")















