#setwd("/home/user/HEX-scratch/KODAMA-Analysis/")

library(ggpubr)
library(ggplot2)

load("data/DLPFC-general.RData")
load("output/BANSKY-results.RData")
load("output/BASS-results.RData")
load("output/BayesSpace-results.RData")
load("output/KODAMA-results.RData")
load("output/Nonspatial-results.RData")
load("output/PRECAST-results.RData")


library(nnet) # multinom function
library(DescTools)

PCA_R2=NULL
for(i in 1:12){
  data=as.matrix(pca_list[[i]][,1:2])
  lab=as.numeric(labels_list[[i]])
  fit=multinom(lab ~ data, maxit = 1000, MaxNWts = 2000,model = TRUE)
  PCA_R2[i]= PseudoR2(fit,c("McFaddenAdj"))
}

KODAMA_R2=NULL
KODAMA_ARI=NULL
for(i in 1:12){
  KODAMA_ARI[i]=mclust::adjustedRandIndex(results_KODAMA$labels[[i]],results_KODAMA$clusters[[i]])
  data=as.matrix(results_KODAMA$feature_extraction[[i]][,1:2])
  lab=as.numeric(results_KODAMA$labels[[i]])
  fit=multinom(lab ~ data, maxit = 1000, MaxNWts = 2000,model = TRUE)
  KODAMA_R2[i]= PseudoR2(fit,c("McFaddenAdj"))
}

BASS_ARI=NULL
for(i in 1:12){
  BASS_ARI[i]=mclust::adjustedRandIndex(results_BASS$labels[[i]],results_BASS$clusters[[i]])
}

BayesSpace_ARI=NULL
for(i in 1:12){
  BayesSpace_ARI[i]=mclust::adjustedRandIndex(results_BayesSpace$labels[[i]],results_BayesSpace$clusters[[i]])
}

Leiden_ARI=NULL
for(i in 1:12){
  Leiden_ARI[i]=mclust::adjustedRandIndex(results_Leiden$labels[[i]],results_Leiden$clusters[[i]])
}

Louvain_ARI=NULL
for(i in 1:12){
  Louvain_ARI[i]=mclust::adjustedRandIndex(results_Louvain$labels[[i]],results_Louvain$clusters[[i]])
}

Walktrap_ARI=NULL
for(i in 1:12){
  Walktrap_ARI[i]=mclust::adjustedRandIndex(results_Walktrap$labels[[i]],results_Walktrap$clusters[[i]])
}

PRECAST_ARI=NULL
for(i in 1:12){
  PRECAST_ARI[i]=mclust::adjustedRandIndex(results_PRECAST$labels[[i]],results_PRECAST$clusters[[i]])
}

BANKSY_ARI=NULL
for(i in 1:12){
  BANKSY_ARI[i]=mclust::adjustedRandIndex(results_BANKSY$labels[[i]],results_BANKSY$clusters[[i]])
}









# Convert results into a data frame
ari_data <- data.frame(
  Method = c(rep("KODAMA", length(KODAMA_ARI)),
             rep("BANKSY", length(BANKSY_ARI)),
             rep("PRECAST", length(PRECAST_ARI)),
             rep("BayesSpace", length(BayesSpace_ARI)),
             rep("BASS", length(BASS_ARI)),
             rep("Walktrap", length(Walktrap_ARI)),
             rep("Leiden", length(Leiden_ARI)),
             rep("Louvain", length(Louvain_ARI))),
  ARI = c(KODAMA_ARI, BANKSY_ARI, PRECAST_ARI, BayesSpace_ARI,BASS_ARI, Walktrap_ARI, Leiden_ARI, Louvain_ARI)
)

# Colors for the box plot
colors <- c("#0073c2bb", "#efc000bb", "#868686bb", "#cd534cbb", "#7aabdcbb", "#003c67bb", "#7b6f5e", "#3c3c3c")

# Create the box plot
Nplot1 <- ggboxplot(ari_data, x = "Method", y = "ARI", palette = colors,
                    fill = "Method", add = "jitter",
                    add.params = list(size = 2, jitter = 0.2, fill = 3, shape = 21)) +
  ylab("Adjusted Rand Index") +
  labs(title = "",
       x = "Clustering Method") +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt x-axis labels for better readability
  ylim(0, NA) # Set the start of the y-axis at 0


pdf("output/Figure 1 - boxplot.pdf")
# Display the plot
print(Nplot1)
dev.off()


















##################################
# Define cluster colors
cols_cluster <- c("#0000b6", "#81b29a", "#f2cc8f", "#e07a5f", "#cc00b6", "#81ccff", "#33b233")

##################################
# Set up the plot layout for 3 rows and 3 columns
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))  # margins: bottom, left, top, right

# Define the list of methods and their associated clusters
methods <- c("Real Labels", "KODAMA", "BayesSpace", "PRECAST", "BASS", "BANKSY", "Louvain", "Leiden", "Walktrap")

# List of clusters (first element will be the real labels)
clusters <- list(
  as.vector(as.factor(labels_list[[10]])),            # Real labels
  as.vector(as.factor(results_KODAMA$clusters[[10]])),         # KODAMA clusters
  as.vector(as.factor(results_BayesSpace$clusters[[10]])), # BayesSpace clusters
  as.vector(as.factor(results_PRECAST$clusters[[10]])),   # PRECAST clusters
  as.vector(as.factor(results_BASS$clusters[[10]])),   # BASS clusters
  as.vector(as.factor(results_BANKSY$clusters[[10]])),    # BANKSY clusters
  as.vector(as.factor(results_Louvain$clusters[[10]])),    # Louvain clusters
  as.vector(as.factor(results_Leiden$clusters[[10]])),     # Leiden clusters
  as.vector(as.factor(results_Walktrap$clusters[[10]]))    # Walktrap clusters
)

# Initialize variables to combine data for all methods
methods_multi <- NULL
xy_multi <- NULL
cluster_multi <- NULL

# Combine data for each method
for (i in 1:9) {
  xy_multi <- rbind(xy_multi, xy_list[[10]])
  methods_multi <- c(methods_multi, rep(methods[i], nrow(xy_list[[10]])))
  cluster_multi <- c(cluster_multi, as.factor(clusters[[i]]))
}

# Convert methods to factor
methods_multi <- factor(methods_multi, levels = methods)


#pdf("output/Figure 2 - DLPFC 10.pdf")
# Plot the clusters with a slide per method
#plot_slide(xy_multi, methods_multi, cluster_multi, col = cols_cluster, nrow = 3)

#dev.off()

# Reset the plot layout
par(mfrow = c(1, 1))

