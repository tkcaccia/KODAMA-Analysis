
library(SpatialExperiment)
library(scater)
library(nnSVG)
library(SPARK)
library(harmony)
library(scuttle)
library(BiocSingular)

tissues=c("Normal_prostate",
          "Acinar_Cell_Carcinoma",
          "Adjacent_normal_section",
          "Adenocarcinoma")

dir="/Users/stefano/Team/Data/Visium/Prostate/"
address <- file.path(dir, tissues, "")

(spe <- read10xVisium(address, tissues, 
                      type = "sparse", data = "raw", 
                      images = "lowres", load = FALSE))



metaData = SingleCellExperiment::colData(spe)
expr = SingleCellExperiment::counts(spe)
sample_names <- paste0("sample_", unique(colData(spe)$sample_id))
sample_names <-  unique(colData(spe)$sample_id)



spe=spe[,colData(spe)$in_tissue]

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
if(length(discard)>0){
  table(discard)
   colData(spe)$discard <- discard


# filter low-quality spots
   spe <- spe[, !colData(spe)$discard]
}
dim(spe)

colnames(rowData(spe)) = "gene_name"

spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 2,   #ncounts
  filter_genes_pcspots = 0.5, 
  filter_mito = TRUE
)

dim(spe)

xy=spatialCoords(spe)
samples=unique(colData(spe)$sample_id)




for(j in 1:length(samples)){
  sel=samples[j]==colData(spe)$sample_id
  xy[sel,1]=spatialCoords(spe)[sel,1]+25000*(j-1)
}
spatialCoords(spe)=xy



# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)






sparkX <- sparkx(logcounts(spe),spatialCoords(spe),numCores=1,option="mixture")



ix=order(sparkX$res_mtest$combinedPval)[1:3000]
rowData(spe)$gene_id=rownames(rowData(spe))
top <- rowData(spe)$gene_id[ix]





sample_id=colData(spe)$sample_id


spe <- runPCA(spe, subset_row = top,scale=TRUE)

spe <- RunHarmony(spe, group.by.vars = "sample_id",lambda=NULL)

plot(reducedDim(spe,type = "PCA"),col=as.factor(colData(spe)$sample_id))
plot(reducedDim(spe,type = "HARMONY"),col=as.factor(colData(spe)$sample_id))

pca=reducedDim(spe,type = "HARMONY")[,1:50]
samples=as.factor(colData(spe)$sample_id)
xy=as.matrix(spatialCoords(spe))
data=t(logcounts(spe))

###############################

patho=read.csv("/Users/stefano/Team/Data/Visium/Prostate/Adenocarcinoma/outs/Pathology.csv")
rownames(patho)=patho[,1]
pathology=rep(NA,ncol(spe))
sel=colData(spe)$sample_id=="Adenocarcinoma"
pathology[sel]=patho[rownames(colData(spe))[sel],"Pathology"]
pathology[pathology==""]=NA
pathology=factor(pathology,levels=c("Invasive carcinoma",
                                    "Blood vessel",
                                    "Fibro-muscular tissue",
                                    "Fibrous tissue",
                                    "Immune Cells" ,
                                    "Nerve",
                                    "Normal gland"))

col_pathology=c("#0000ff","#e41a1c","#006400","#000000","#ffd700","#00ff00","#b2dfee")
#plot(reducedDim(spe,type = "PCA"),col=pathology)



library(KODAMAextra)



#save(spe,file="KODAMA_prostate.RData")
save(pca,
     samples,
     xy,
     col_pathology,
     pathology,file="Prostate_data.RData")



spe <- RunKODAMAmatrix(spe,
                        reduction = "HARMONY",
                        FUN="PLS",
                        landmarks = 100000,
                        splitting = 300,
                        f.par.pls =  50,
                        spatial.resolution = 0.4,
                        n.cores=4)



config=umap.defaults
config$n_threads = 4
config$n_sgd_threads = "auto"


spe=RunKODAMAvisualization(spe,method="UMAP",config=config)

plot(reducedDim(spe,type = "KODAMA"),col=as.factor(colData(spe)$sample_id))
plot(reducedDim(spe,type = "KODAMA"),col=pathology)



save(spe,file="/Users/stefano/Desktop/KODAMA_prostate.RData")

load("/Users/stefano/Desktop/KODAMA_prostate.RData")


kk_UMAP1=reducedDim(spe,type = "KODAMA")

plot(reducedDim(spe,type = "KODAMA"),col=pathology)




pca=reducedDim(spe,type = "HARMONY")


save(xy,pca,file="/home/user/HPC-Ebtesam/KODAMA/4-Prostate/Data/All_harmony.RData")


###################################################################################

plot(kk_UMAP1,col=pathology)

plot(kk_UMAP1,col=as.factor(colData(spe)$sample_id))



plot(kk_UMAP1[pathology=="Normal gland",],col=as.factor(colData(spe)$sample_id)[pathology=="Normal gland"])

da=kk_UMAP1[which(pathology=="Normal gland"),]
sel=which(pathology=="Normal gland")
da=kk_UMAP1[sel,]
ncluster=4

g <- bluster::makeSNNGraph(da,k = 20)
g_walk <- igraph::cluster_walktrap(g)
cluster_normal = as.character(igraph::cut_at(g_walk, no=ncluster))
cluster_normal=as.numeric(cluster_normal)
plot(da,col=cluster_normal)
plot(xy[sel,],col=cluster_normal)
plot(xy[sel,],col=pathology[sel])

plot(da,col=as.factor(colData(spe)$sample_id)[which(pathology=="Normal gland")])
names(cluster_normal)=rownames(da)

pathology2=as.vector(pathology)
pathology2[sel]=paste("Normal gland",cluster_normal)
pathology2=as.factor(pathology2)


clu=rep(NA,ncol(spe))
sel=which(pathology=="Normal gland")
clu[sel]=cluster_normal


library(Seurat)
library(SeuratData)
library(ggplot2)


colData(spe)$groups=as.factor(clu)

markers <- FindMarkers(object = spe, ident.1 = "1", ident.2 = "2")



data("pbmc_small")
# Find markers for cluster 2
markers <- FindMarkers(object = pbmc_small, ident.1 = 2)
head(x = markers)



logC=logcounts(spe)[top,which(pathology=="Normal gland")]


plot(xy[which(pathology=="Normal gland"),],col=clu)

library(KODAMA)
m=multi_analysis(t(logC),clu)
m=cbind(rowData(spe)[top,1],m)
m[order(as.numeric(m$`p-value`)),]



library("DESeq2")
clu=cluster_normal
clu[clu>1]=2
sel1=which(pathology2=="Normal gland 2" | pathology2=="Normal gland 3")
countdata <- as.matrix(logcounts(spe)[,sel1])
condition=factor(pathology2[sel1])
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,DataFrame(condition),~condition)
  ma=cbind(
ma=multi_analysis(t(countdata),condition)


ma=cbind(rowData(spe)$gene_name,ma)




library("GSVA")
library("GSA")

sel1=which(pathology2=="Normal gland 2" | pathology2=="Normal gland 3")
countdata <- as.matrix(logcounts(spe)[,sel1])
genes=countdata
t=rowData(spe)$gene_name
selt=ave(1:length(t), t, FUN = length)
genes=genes[selt==1,]
rownames(genes)=t[selt==1]

geneset=GSA.read.gmt("/Users/stefano/OneDrive/Projects/KODAMA spatial transcriptomics/Liver analysis/Genesets/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.symbols.gmt")
names(geneset$genesets)=geneset$geneset.names



gsva_TCGA <- gsva(genes, geneset$genesets,min.sz = 5)


ll=condition
ma=multi_analysis(t(gsva_TCGA),ll)

ma=ma[order(as.numeric(ma$`p-value`)),]
