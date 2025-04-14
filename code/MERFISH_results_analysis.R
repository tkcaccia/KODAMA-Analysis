library(ggpubr)
library(mclust)
library(cluster)
#BiocManager::install("RDRToolbox")
library("RDRToolbox")
library(nnet)
library(DescTools)
#setwd("/Users/stefano/HPC-scratch/KODAMA-Analysis/")

#####################################################################################
# Define method colors

method_colors_all=c("Real" = "#aa7ddddb",
              "KODAMA" = "#0073c2bb",
              "BASS"= "#efc000bb",
              "BayesSpace"= "#868686bb",
              "BANKSY" = "#cd534cbb",
              "PRECAST"= "#7aabdcbb",
              "Louvain"= "#003c67bb",
              "Leiden"= "#7b6f5ebb",
              "Walktrap" = "#3c3c3cbb",
              "Louvain.PM"= "#8e44adbb",
              "Leiden.PM"= "#3498dbbb",
              "Walktrap.PM"= "#3498dbbb",
              "Louvain.PM.ref"= "#3498dbbb",
              "Leiden.PM.ref"= "#3498dbbb",
              "Walktrap.PM.ref"= "#3498dbbb",
              "UMAP"= "#3498dbbb",
              "tSNE"= "#3498dbbb",
              "PCA"= "#3498dbbb",
              "UMAP.PM"= "#3498dbbb",
              "tSNE.PM"= "#3498aabb",
              "PCA.PM"= "#3498dbbb")

all_methods=names(method_colors_all)

load("output/MERFISH-KODAMA-results.RData")
load("output/MERFISH-BANSKY-results.RData")
load("output/MERFISH-BASS-results.RData")
load("output/MERFISH-BayesSpace-results.RData")
load("output/MERFISH-KODAMA-results.RData")
load("output/MERFISH-NONSPATIAL-results.RData")
load("output/MERFISH-PRECAST-results.RData")


# Define Bregma levels and clustering methods
bregma_levels <- c("0.24", "0.19", "0.14", "0.09", "0.04")



all_results=list(Real=ground_true,
                 KODAMA=result_KODAMA,
                 BASS=result_BASS,
                 BayesSpace=result_BayesSpace,
                 BANKSY=result_BANKSY,
                 PRECAST=result_PRECAST,
                 Louvain=result_Louvain,
                 Leiden=result_Leiden,
                 Walktrap=result_Walktrap,
                 Louvain.PM=result_Louvain.PM,
                 Leiden.PM=result_Leiden.PM,
                 Walktrap.PM=result_Walktrap.PM,
                 Louvain.PM.ref=result_Louvain.PM.ref,
                 Leiden.PM.ref=result_Leiden.PM.ref,
                 Walktrap.PM.ref=result_Walktrap.PM.ref,
                 UMAP=result_UMAP,
                 tSNE=result_tSNE,
                 PCA=result_PCA,
                 UMAP.PM=result_UMAP.PM,
                 tSNE.PM=result_tSNE.PM,
                 PCA.PM=result_PCA.PM)



methods_tested_for_clustering_1 <- c("KODAMA", "BASS", "BayesSpace", "BANKSY", "PRECAST", "Louvain", "Leiden", "Walktrap")

lARI=lapply(all_results[methods_tested_for_clustering_1],function(x) c(adjustedRandIndex(x$tissue_segments[[1]],x$clusters[[1]]),
                                 adjustedRandIndex(x$tissue_segments[[2]],x$clusters[[2]]),
                                 adjustedRandIndex(x$tissue_segments[[3]],x$clusters[[3]]),
                                 adjustedRandIndex(x$tissue_segments[[4]],x$clusters[[4]]),
                                 adjustedRandIndex(x$tissue_segments[[5]],x$clusters[[5]])))

ARI=matrix(unlist(lARI),ncol=5,byrow = TRUE)
rownames(ARI)=methods_tested_for_clustering_1
colnames(ARI)=bregma_levels


ari_data <- data.frame(
  Method = factor(rep(methods_tested_for_clustering_1,each=5),levels=all_methods),
  ARI = unlist(lARI))


################################################################



methods_tested <- c("KODAMA",  "BayesSpace","BANKSY", "PRECAST","UMAP","tSNE","PCA")


lsilhouette=lapply(all_results[methods_tested],function(x) c(summary(silhouette(as.numeric(x$tissue_segments[[1]]),
                                                                     dist(x$feature_extraction[[1]])))$si.summary["Mean"],
                                                  summary(silhouette(as.numeric(x$tissue_segments[[2]]),
                                                                     dist(x$feature_extraction[[2]])))$si.summary["Mean"],
                                                  summary(silhouette(as.numeric(x$tissue_segments[[3]]),
                                                                     dist(x$feature_extraction[[3]])))$si.summary["Mean"],
                                                  summary(silhouette(as.numeric(x$tissue_segments[[4]]),
                                                                     dist(x$feature_extraction[[4]])))$si.summary["Mean"],
                                                  summary(silhouette(as.numeric(x$tissue_segments[[5]]),
                                                                     dist(x$feature_extraction[[5]])))$si.summary["Mean"]
                                                  ))

SILHOUETTE=matrix(unlist(lsilhouette),ncol=5,byrow = TRUE)
rownames(SILHOUETTE)=methods_tested
colnames(SILHOUETTE)=bregma_levels


silhouette_data <- data.frame(
  Method = factor(rep(methods_tested,each=5),levels=all_methods),
  silhouette = unlist(lsilhouette))



lDBIndex=lapply(all_results[methods_tested],function(x) c(1-DBIndex(x$feature_extraction[[1]],as.numeric(x$tissue_segments[[1]])),
                                                      1-DBIndex(x$feature_extraction[[2]],as.numeric(x$tissue_segments[[2]])),
                                                      1-DBIndex(x$feature_extraction[[3]],as.numeric(x$tissue_segments[[3]])),
                                                      1-DBIndex(x$feature_extraction[[4]],as.numeric(x$tissue_segments[[4]])),
                                                      1-DBIndex(x$feature_extraction[[5]],as.numeric(x$tissue_segments[[5]]))
))


DBINDEX=matrix(unlist(lDBIndex),ncol=5,byrow = TRUE)
rownames(DBINDEX)=methods_tested
colnames(DBINDEX)=bregma_levels


DBIndex_data <- data.frame(
  Method = factor(rep(methods_tested,each=5),levels=all_methods),
  DBIndex = unlist(lDBIndex))


pseudoR2_McFaddenAdj = function(data,lab){
  fit <- multinom(lab ~ data, family = multinomial(), data = as.data.frame(data), model = TRUE)
  PseudoR2(fit, c("McFaddenAdj"))
}
lpseudoR2=lapply(all_results[methods_tested],function(x) c(pseudoR2_McFaddenAdj(x$feature_extraction[[1]],as.numeric(x$tissue_segments[[1]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[2]],as.numeric(x$tissue_segments[[2]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[3]],as.numeric(x$tissue_segments[[3]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[4]],as.numeric(x$tissue_segments[[4]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[5]],as.numeric(x$tissue_segments[[5]]))
))


PSEUDOR2=matrix(unlist(lpseudoR2),ncol=5,byrow = TRUE)
rownames(PSEUDOR2)=methods_tested
colnames(PSEUDOR2)=bregma_levels


pseudoR2_data <- data.frame(
  Method = factor(rep(methods_tested,each=5),levels=all_methods),
  pseudoR2 = unlist(lpseudoR2))


FE_data <- data.frame(
  Method = factor(rep(methods_tested,each=5),levels=all_methods),
  DBIndex = unlist(lDBIndex),
  pseudoR2 = unlist(lpseudoR2))





# Create the box plot
Nplot2 <- ggboxplot(ari_data, x = "Method", y = "ARI", palette = method_colors_all,
                    fill = "Method", add = "jitter",
                    add.params = list(size = 2, jitter = 0.2, fill = 3, shape = 21)) +
  ylab("Adjusted Rand Index") +
  labs(title = "",
       x = "Clustering Method") +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt x-axis labels for better readability
  ylim(0, NA) # Set the start of the y-axis at 0


svg("output/Figures/Figure 2 - boxplot_ARI_merfish.svg")
# Display the plot
print(Nplot2)
dev.off()




# Create the scatter plot with well-defined axes
Nplot1 <- ggplot(FE_data, aes(x = pseudoR2, y = DBIndex, color = Method, shape = Method)) +
  geom_point(size = 3) +
  scale_color_manual(values = method_colors_all) +
  labs(title = "Scatter Plot of Pseudo R^2 vs DB Index by Method",
       x = "Pseudo R^2 (Coefficient of Determination)",
       y = "DB Index (Davies-Bouldin Index)") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  theme_minimal(base_size = 15) +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks = element_line(size = 0.8),
    axis.line = element_line(size = 0.8, color = "black")
  )

pdf("output/Figure/Figure 1 - Scatter Plot of Pseudo R^2 vs DB Index by Method MERFISH.pdf")
print(Nplot1)
dev.off()



###########################################################################


#####################################################################

########PLOT 1 SLIDE

##################################
# Load necessary libraries
library(grid)
library(gridExtra)

# Define Bregma levels and clustering methods
methods <- c("Real Labels", "KODAMA", "BASS", "BayesSpace", "BANKSY", "PRECAST", "Louvain", "Leiden", "Walktrap")

# Colors for each cluster with a softer palette
cols_cluster <- c("#0000b6", "#81b29a", "#f2cc8f", "#e07a5f", "#81ccff", "#cc00b6", "#33b233", "#ff5733", "#FFD111")




c("Real Labels", "KODAMA", "BASS", "BayesSpace", "BANKSY", "PRECAST", "Louvain", "Leiden", "Walktrap")

slide_number=4
xy_slide_comparison=rbind(ground_true$xy[[slide_number]][,1:2],
                          result_KODAMA$xy[[slide_number]][,1:2],
                          result_BASS$xy[[slide_number]][,1:2],
                          result_BayesSpace$xy[[slide_number]][,1:2],
                          result_BANKSY$xy[[slide_number]][,1:2],
                          result_PRECAST$xy[[slide_number]][,1:2],
                          result_Louvain$xy[[slide_number]][,1:2],
                          result_Leiden$xy[[slide_number]][,1:2],
                          result_Walktrap$xy[[slide_number]][,1:2])
cluster_slide_comparison=c(ground_true$tissue_segments[[slide_number]],
                          result_KODAMA$clusters[[slide_number]],
                          result_BASS$clusters[[slide_number]],
                          result_BayesSpace$clusters[[slide_number]],
                          result_BANKSY$clusters[[slide_number]],
                          result_PRECAST$clusters[[slide_number]],
                          result_Louvain$clusters[[slide_number]],
                          result_Leiden$clusters[[slide_number]],
                          result_Walktrap$clusters[[slide_number]])

method_slide_comparison=factor(rep(methods,c(length(ground_true$tissue_segments[[slide_number]]),
                          length(result_KODAMA$clusters[[slide_number]]),
                          length(result_BASS$clusters[[slide_number]]),
                          length(result_BayesSpace$clusters[[slide_number]]),
                          length(result_BANKSY$clusters[[slide_number]]),
                          length(result_PRECAST$clusters[[slide_number]]),
                          length(result_Louvain$clusters[[slide_number]]),
                          length(result_Leiden$clusters[[slide_number]]),
                          length(result_Walktrap$clusters[[slide_number]]))),levels=methods)

Ngplot3 <- plot_slide(xy_slide_comparison,method_slide_comparison,cluster_slide_comparison,col = cols_cluster, nrow = 3) +
  ggtitle("Bregma -0.09") + theme(plot.title = element_text(size = 15, face = "bold"))



# Create and save the PDF with customized layout
svg("output/Figures/Figure 3_Comparative_Methods.svg", width = 10, height = 10)

Ngplot3

# Close the PDF device
dev.off()



#################################################

######################################################################

#####################################

########PLOT 12 SLIDE

# Load necessary libraries
library(grid)
library(gridExtra)



method_tested_all=c("Real",
                    "KODAMA",
                    "BASS",
                    "BayesSpace",
                    "BANKSY",
                    "PRECAST",
                    "Louvain",
                    "Leiden",
                    "Walktrap" ,
                    "Louvain.PM",
                    "Leiden.PM",
                    "Walktrap.PM",
                    "Louvain.PM.ref",
                    "Leiden.PM.ref",
                    "Walktrap.PM.ref")


plotall=lapply(all_results[method_tested_all],function(x) plot_slide(x$xyz[,-3],x$xyz[,3],x$clusters_allslide,  col = cols_cluster))



# Define grid dimensions
ncol_value <- 1  # 5 columns for Bregma levels
nrow_value <- length(method_tested_all)  # 9 rows for Methods

# Create and save the PDF with customized layout
png("output/Figures/Figure_4_plots_12_slide.png", width = 2000, height = 5000)

# Arrange the plots with a bold title and customized layout
Nplot4=grid.arrange(grobs = plotall, ncol = ncol_value, nrow = nrow_value,
             padding = unit(c(2), "cm"),
             top = textGrob("Spatial Domain Detection in MERFISH Dataset",
                            gp = gpar(fontsize = 24, fontface = "bold", col = "black"))
)

# Add Bregma level titles at the top of each column in larger font
for (k in seq_along(bregma_levels)) {
  grid.text(paste("Bregma", bregma_levels[k]),
            x = 0.11 + (k - 1) * 0.2,  # Adjust x position to align with the first plot's title
            y = 0.963,  # Adjust y to match the position of the first method's title
            gp = gpar(fontsize = 18, fontface = "bold", col = "black"))
}


# Close the PDF device
dev.off()

############################################################################################





result_KODAMA$feature_extraction_allslide=rbind(result_KODAMA$feature_extraction[[1]],result_KODAMA$feature_extraction[[2]],result_KODAMA$feature_extraction[[3]],result_KODAMA$feature_extraction[[4]],result_KODAMA$feature_extraction[[5]])

methods=c("KODAMA", "BayesSpace", "BANKSY", "PRECAST")

xy_slide_comparison=rbind(result_KODAMA$feature_extraction_allslide,
                          result_BayesSpace$feature_extraction_allslide,
                          result_BANKSY$feature_extraction_allslide,
                          result_PRECAST$feature_extraction_allslide)
cluster_slide_comparison=c(unlist(result_KODAMA$tissue_segments),
                           unlist(result_BayesSpace$tissue_segments),
                           unlist(result_BANKSY$tissue_segments),
                           unlist(result_PRECAST$tissue_segments))

method_slide_comparison=factor(rep(methods,c(length(unlist(result_KODAMA$tissue_segments)),
                                             length(unlist(result_BayesSpace$tissue_segments)),
                                             length(unlist(result_BANKSY$tissue_segments)),
                                             length(unlist(result_PRECAST$tissue_segments)))),levels=methods)

Ngplot3 <- plot_slide(xy_slide_comparison,method_slide_comparison,cluster_slide_comparison,col = cols_cluster, nrow = 2) +
  ggtitle("Bregma -0.09") + theme(plot.title = element_text(size = 15, face = "bold"))



# Create and save the PDF with customized layout
svg("output/Figures/Figure 3_Comparative_Methods.svg", width = 10, height = 10)

Ngplot3

# Close the PDF device
dev.off()


































#############################################################################################


###################################

##########################################

color_palette <-c("#0000b6", "#81b29a", "#f2cc8f", "#e07a5f", "#cc00b6",
  "#81ccff", "#33b233", "#ff5733")

########plot 12 slide feature extraction
# Load the necessary libraries
library(ggplot2)
library(dplyr)

# Initialize combined lists for each method
plot_data_combined_kodama <- list()
plot_data_combined_precast <- list()
plot_data_combined_banksy <- list()
plot_data_combined_pca <- list()
plot_data_combined_umap <- list()
plot_data_combined_t_SNE <- list()

# # Create a "slide" list for each method
 sld = split(slide, slide)

# Prepare the data for each method
# KODAMA
for (i in 1:5) {
  data_kodama <- as.data.frame(results_KODAMA$feature_extraction[[i]][, 1:2])
  colnames(data_kodama) <- c("x", "y")  # matching with plot_slide
  data_kodama$Label <- factor(tissue_segments_list[[i]])
  data_kodama$Method <- "KODAMA"
  data_kodama$slide <- sld[[i]]
  plot_data_combined_kodama[[i]] <- data_kodama
}

# PRECAST
for (i in 1:5) {
  data_precast <- as.data.frame(result_PRECAST$feature_extraction[[i]][, 1:2])
  colnames(data_precast) <- c("x", "y")
  data_precast$Label <- factor(result_PRECAST$tissue_segments[[i]])
  data_precast$Method <- "PRECAST"
  data_precast$slide <- sld[[i]][1:nrow(data_precast)]
  plot_data_combined_precast[[i]] <- data_precast
}

# BANKSY
for (i in 1:5) {
  data_banksy <- as.data.frame(result_BANKSY$feature_extraction[[i]][, 1:2])
  colnames(data_banksy) <- c("x", "y")
  data_banksy$Label <- factor(tissue_segments_list[[i]])
  data_banksy$Method <- "BANKSY"
  data_banksy$slide <- sld[[i]]
  plot_data_combined_banksy[[i]] <- data_banksy
}

# PCA
for (i in 1:5) {
  pca_data_matrix <- matrix(pca_list[[i]], ncol=2, byrow=TRUE)
  pca_data <- as.data.frame(pca_data_matrix)
  colnames(pca_data) <- c("x", "y")
  pca_data$Label <- factor(tissue_segments_list[[i]])
  pca_data$Method <- "PCA"
  pca_data$slide <- sld[[i]]
  plot_data_combined_pca[[i]] <- pca_data
}

# UMAP
for (i in 1:5) {
  umap_data <- as.data.frame(umap_results_list[[i]])
  colnames(umap_data) <- c("x", "y")
  umap_data$Label <- factor(tissue_segments_list[[i]])
  umap_data$Method <- "UMAP"
  umap_data$slide <- sld[[i]]
  plot_data_combined_umap[[i]] <- umap_data
}

# t-SNE
for (i in 1:5) {
  tsne_data <- as.data.frame(tsne_results_list[[i]])
  colnames(tsne_data) <- c("x", "y")
  tsne_data$Label <- factor(tissue_segments_list[[i]])
  tsne_data$Method <- "t-SNE"
  tsne_data$slide <- sld[[i]]
  plot_data_combined_t_SNE[[i]] <- tsne_data
}

# Combine all data for the final plot
final_plot_data <- bind_rows(
  do.call(rbind, plot_data_combined_kodama),
  do.call(rbind, plot_data_combined_precast),
  do.call(rbind, plot_data_combined_banksy),
  do.call(rbind, plot_data_combined_pca),
  do.call(rbind, plot_data_combined_umap),
  do.call(rbind, plot_data_combined_t_SNE)
)

# Create a list to store ggplot objects for each method
plot_list <- list()

# Create plots for each method and store them in plot_list
methods <- unique(final_plot_data$Method)

for (method in methods) {
  # Filter data for the current method
  method_data <- final_plot_data %>% filter(Method == method)

  # Check the structure of method_data to confirm the correct columns
  if ("slide" %in% colnames(method_data)) {
    slide <- method_data$slide
  } else {
    stop("The 'slide' column is missing from the data.")
  }

  # Prepare xy and labels for the plot_slide function
  xy <- method_data[, c("x", "y")]
  labels <- method_data$Label

  # Call the plot_slide function with a title for each method
  plot <- plot_slide(xy, slide, labels, color_palette, nrow = 1)
  plot <- plot + ggtitle(bquote(bold(.(method))))
  plot_list[[method]] <- plot
}

# Create a single long PDF page
pdf("output/Figure/Figure 5 - methods_feature.pdf", width = 20, height = 18)  # Adjust height

# Arrange and print all plots
Nplot5 = grid.arrange(grobs = plot_list, ncol = 1)  # Arrange plots in one column

dev.off()  # Close the PDF device



ggsave("output/svg/Figure 5 - methods_feature.svg", plot = Nplot5, width = 50, height = 50, units = "cm", device = "svg")

###############################################################
####refine comparaison

library(gridExtra)

####### refine performing
refreallabels=refine_SVM(xyz,unlist(tissue_segments_list),cost=10000)
refkodama=refine_SVM(xyz,unlist(results_KODAMA$cluster),cost=10000)
refBASS=refine_SVM(xyz,unlist(result_BASS$clusters),cost=10000)
refBayesSpace=refine_SVM(xyz,unlist(result_BayesSpace$clusters),cost=10000)
refBANKSY=refine_SVM(xyz,unlist(result_BANKSY$clusters),cost=10000)
refPRECAST=refine_SVM(result_PRECAST$xyz,unlist(result_PRECAST$clusters),cost=10000)
reflouvain=refine_SVM(xyz,unlist(result_Louvain$clusters),cost=10000)
refLeiden=refine_SVM(xyz,unlist(result_Leiden$clusters),cost=10000)
refWalktrap=refine_SVM(xyz,unlist(result_Walktrap$clusters),cost=10000)

################refine plot

plot_refreallabels <- plot_slide(xyz, refreallabels, rep("real labels", length(refreallabels)), col = cols_cluster,nrow = 1)
plot_refkodama <- plot_slide(xyz, refkodama, rep("KODAMA", length(refkodama)), col = cols_cluster,nrow = 1)
plot_refBASS <- plot_slide(xyz, refBASS, rep("BASS", length(refBASS)), col = cols_cluster, nrow = 1)
plot_refBayesSpace <- plot_slide(xyz, refBayesSpace, rep("BayesSpace", length(refBayesSpace)), col = cols_cluster, nrow = 1)
plot_refBANKSY <- plot_slide(xyz, refBANKSY, rep("BANKSY", length(refBANKSY)), col = cols_cluster, nrow = 1)
plot_refPRECAST <- plot_slide(result_PRECAST$xyz, refPRECAST, rep("PRECAST", length( refPRECAST)), col = cols_cluster, nrow = 1)
plot_Louvainref <- plot_slide(xyz, reflouvain, rep("Louvain", length(reflouvain)), col = cols_cluster,  nrow = 1)
plot_refLeiden <- plot_slide(xyz, refLeiden, rep("Leiden", length( refLeiden)), col = cols_cluster, nrow = 1)
plot_refWalktrap <- plot_slide(xyz, refWalktrap, rep("Walktrap", length( refWalktrap)), col = cols_cluster, nrow = 1)

combined_plots_ref <- list(
  plot_refreallabels, plot_refkodama, plot_refBASS, plot_refBayesSpace,
  plot_refBANKSY, plot_refPRECAST, plot_Louvainref, plot_refLeiden,
  plot_refWalktrap
)


# Create the grid layout with 9 rows and 5 columns
final_grid <- gridExtra::arrangeGrob(
  grobs = combined_plots_ref,
  ncol = 1,
  top = "Comparative Analysis of Methods"
)

# Save the final grid to a PDF
pdf("output/Figure/figure 6 refine svm.pdf", width = 35, height = 50)
grid.draw(final_grid)
dev.off()
###########################
