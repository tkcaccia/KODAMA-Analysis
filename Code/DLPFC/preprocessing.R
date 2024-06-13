library("spatialLIBD")
library("nnSVG")
library("scater")
library("scran")
library("scry")
library("SPARK")
library("harmony")
library("Seurat")

address="/Users/stefano/HPC-scratch/KODAMA/DLPFC/"
setwd(address)
dir.create("Data")

sample_names=c("151507",
               "151508",
               "151509",
               "151510",
               "151669",
               "151670",
               "151671",
               "151672",
               "151673",  
               "151674",
               "151675",
               "151676")

subject_names= c("Br5292","Br5595", "Br8100")
################################################################################
# preprocessing
spe <- fetch_data(type = 'spe')


metaData = SingleCellExperiment::colData(spe)
expr = SingleCellExperiment::counts(spe)
sample_names <- paste0("sample_", unique(colData(spe)$sample_id))
sample_names <-  unique(colData(spe)$sample_id)
dim(spe)
# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# calculate per-spot QC metrics
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

# select QC thresholds
qc_lib_size <- colData(spe)$sum < 500
qc_detected <- colData(spe)$detected < 250
qc_mito <- colData(spe)$subsets_mito_percent > 30
qc_cell_count <- colData(spe)$cell_count > 12

# spots to discard
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
table(discard)
colData(spe)$discard <- discard
# filter low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)



################################################################################

# readjust xy
xy=spatialCoords(spe)
samples=unique(colData(spe)$sample_id)
for(j in 1:length(samples)){
  sel=samples[j]==colData(spe)$sample_id
  xy[sel,1]=spatialCoords(spe)[sel,1]+12000*(j-1)
}
spatialCoords(spe)=xy

# filter genes
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 2,   #ncounts
  filter_genes_pcspots = 0.5, 
  filter_mito = TRUE
)

dim(spe)

sel= !is.na(colData(spe)$layer_guess_reordered)
spe = spe[,sel]
dim(spe)


dim(spe)

# normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)


gene_i=NULL
pvalue_i=NULL
for(i in 1:length(sample_names)){
  sel=colData(spe)$sample_id==sample_names[i]
  spe_sub= spe[,sel]
  
  sparkX <- sparkx(logcounts(spe_sub),spatialCoords(spe_sub),numCores=1,option="mixture")
  
  gene_i=c(gene_i,rowData(spe)$gene_id)
  pvalue_i=c(pvalue_i,sparkX$res_mtest$combinedPval)
  
  print(sample_names[i])
}
oo=order(pvalue_i)
top_genes=gene_i[oo]
n=ave(1:length(top_genes), top_genes, FUN = seq_along)
top_genes=top_genes[n==1]

top=top_genes

################################################################################
#OPTION2: PCA without harmony, analysis all samples from the same subject together
################################################################################
seuList=list()

labels=NULL
xy=NULL
samples=NULL
spe_list=NULL
subjects=NULL
data=NULL
pca=NULL
for(i in 1:length(sample_names)){
  print(i)
  spe_list[[i]] <- spe[top[1:2000], colData(spe)$sample_id ==  sample_names[i]]
  spe_list[[i]] <- runPCA(spe_list[[i]], 50,subset_row = top[1:2000], scale=TRUE)
  pca[[i]]=reducedDim(spe_list[[i]],type = "PCA")[,1:50]
  labels[[i]]=as.factor(colData(spe_list[[i]])$layer_guess_reordered)
  xy[[i]]=as.matrix(spatialCoords(spe_list[[i]]))
  samples[[i]]=colData(spe_list[[i]])$sample_id
  subjects[[i]]=colData(spe)$subject
  data[[i]]=t(logcounts(spe_list[[i]]))
  colnames(xy[[i]])=c("row","col")

  seuList[[i]]  <- CreateSeuratObject(counts = spe_list[[i]]@assays@data$counts, 
                             meta.data = as.data.frame(colData(spe_list[[i]])))
  
  seuList[[i]] [['RNA']] <- AddMetaData(seuList[[i]] [['RNA']], metadata=as.data.frame( rowData(spe_list[[i]])))
  seuList[[i]] @tools$platform <- "Visium"
  
  colnames(seuList[[i]]@meta.data) [colnames(seuList[[i]]@meta.data) %in% c("array_row","array_col")]=c("row","col")
 
#  colnames(seuList[[i]]@meta.data)=c("orig.ident","row","col")
#  seuList[[i]]@meta.data=seuList[[i]]@meta.data[,c("row","col")]
  
}


 save(data,xy,labels,pca,top,samples,subjects,file=paste("Data/",subject_names[i],".RData",sep=""))





for(i in 1:length(subject_names)){
  print(i)
  spe_sub <- spe[, colData(spe)$subject ==  subject_names[i]]
  subjects=colData(spe_sub)$subject
  dim(spe_sub)
  spe_sub <- runPCA(spe_sub, 50,subset_row = top[1:2000], scale=TRUE)
  plot(reducedDim(spe_sub,type = "PCA"), col=as.factor(colData(spe_sub)$sample_id))
  
  pca=reducedDim(spe_sub,type = "PCA")[,1:50]
  labels=as.factor(colData(spe_sub)$layer_guess_reordered)
  xy=as.matrix(spatialCoords(spe_sub))
  samples=colData(spe_sub)$sample_id
  data=t(logcounts(spe_sub)[top[1:2000],])
  
  save(data,xy,labels,pca,top,samples,subjects,file=paste("Data/",subject_names[i],".RData",sep=""))
}
################################################################################

for(i in 1:length(subject_names)){
  print(i)
  spe_sub <- spe[, colData(spe)$subject ==  subject_names[i]]
  subjects=colData(spe_sub)$subject
  dim(spe_sub)
  spe_sub <- runPCA(spe_sub, 50,subset_row = top[1:2000], scale=TRUE)
  plot(reducedDim(spe_sub,type = "PCA"), col=as.factor(colData(spe_sub)$sample_id))
  
  spe_sub <- RunHarmony(spe_sub, group.by.vars = "sample_id",lambda=NULL)

  pca=reducedDim(spe_sub,type = "HARMONY")[,1:50]
  
  labels=as.factor(colData(spe_sub)$layer_guess_reordered)
  xy=as.matrix(spatialCoords(spe_sub))
  samples=colData(spe_sub)$sample_id
  data=t(logcounts(spe_sub)[top[1:2000],])
  
  save(data,xy,labels,pca,top,samples,subjects,file=paste("Data/",subject_names[i],"_harmony.RData",sep=""))
}
################################################################################

 
 data=t(logcounts(spe)[top[1:2000],])
 
 subjects=colData(spe)$subject
 dim(spe_sub)
 spe <- runPCA(spe, 50,subset_row = top[1:2000], scale=TRUE)
 
 labels=as.factor(colData(spe)$layer_guess_reordered)
 xy=as.matrix(spatialCoords(spe))
 samples=colData(spe)$sample_id
 
 spe <- RunHarmony(spe, subjects,lambda=NULL)
 pca=reducedDim(spe,type = "HARMONY")[,1:50]
 
 save(data,xy,labels,pca,top,samples,subjects,file=paste("Data/All.RData",sep=""))
 
 
 
 data=t(logcounts(spe)[top[1:2000],])
 
 subjects=colData(spe)$subject
 dim(spe_sub)
 spe <- runPCA(spe, 50,subset_row = top[1:2000], scale=TRUE)
 
 labels=as.factor(colData(spe)$layer_guess_reordered)
 xy=as.matrix(spatialCoords(spe))
 samples=colData(spe)$sample_id
 
 spe <- RunHarmony(spe, group.by.vars = c("sample_id","subject"),lambda=NULL)
 pca=reducedDim(spe,type = "HARMONY")[,1:50]
 
 save(data,xy,labels,pca,top,samples,subjects,file=paste("All-H2.RData",sep=""))
 
 








