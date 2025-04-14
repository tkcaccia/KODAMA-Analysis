library(ggpubr)
library(mclust)
library(cluster)
library("RDRToolbox")
library(nnet)
library(DescTools)
#setwd("/Users/stefano/HPC-scratch/KODAMA-Analysis/")
#setwd("/home/user/HEX-scratch/KODAMA-Analysis/")


library(KODAMAextra)


load("output/DLPFC-BANSKY-results.RData")
load("output/DLPFC-BASS-results.RData")
load("output/DLPFC-BayesSpace-results.RData")
load("output/KODAMA-results.RData")
load("output/DLPFC-Nonspatial-results.RData")
load("output/DLPFC-PRECAST-results.RData")


cols_cluster <- c("#0000b6", "#81b29a", "#f2cc8f","#e07a5f",
                  "#cc00b6", "#81ccff", "#33b233")

xy=do.call(rbind, ground_true$xy)
samples=do.call(c,result_PCA$samples)
labels=do.call(c,ground_true$tissue_segments)

svg("output/Figures/DLPFC/Figure 1 - ground true.svg", width = 36, height = 4.5)
plot_slide(xy,samples,labels,col=cols_cluster)
dev.off()


load("output/DLFPC-Br5292.RData")
load("output/DLFPC-Br5595.RData")
load("output/DLFPC-Br8100.RData")

KODAMA=rbind(kk_UMAP_Br5292,kk_UMAP_Br5595,kk_UMAP_Br8100)
subjects=do.call(c,result_PCA$subjects)
labels=do.call(c,ground_true$tissue_segments)



svg("output/Figures/DLPFC/Figure 2 - KODAMA.svg", width = 12, height = 4.5)
plot_slide(KODAMA,subjects,labels,col=cols_cluster)
dev.off()



############################################################################################################
samples_name=unique(ground_true$samples)


# Define method colors

method_colors_all <- c(
  "Real" = "#FF5733aa",      # Rouge vif
  "KODAMA" = "#99CC37aa",    # Vert vif
  "BASS" = "#3357FFaa",      # Bleu vif
  "BayesSpace" = "#FF33A1aa",# Rose vif
  "BANKSY" = "#CCBB33aa",    # Jaune doré
  "PRECAST" = "#33aaa0aa",   # Turquoise
  "Louvain" = "#B833FFaa",   # Violet vif
  "Leiden" = "#8D8833aa",    # Vert lime
  "Walktrap" = "#FF8C33aa",  # Orange vif
  "UMAP" = "#5733FFaa",       # Bleu marine
  "tSNE" = "#aa7733aa",       # Vert menthe
  "PCA" = "#ccaaffaa"
)


all_methods=names(method_colors_all)



all_results=list(Real=ground_true,
                 KODAMA=result_KODAMA,
                 BASS=result_BASS,
                 BayesSpace=result_BayesSpace,
                 BANKSY=result_BANKSY,
                 PRECAST=result_PRECAST,
                 Louvain=result_Louvain,
                 Leiden=result_Leiden,
                 Walktrap=result_Walktrap,
                 UMAP=result_UMAP,
                 tSNE=result_tSNE,
                 PCA=result_PCA)



methods_tested_for_clustering_1 <- c("KODAMA", "BASS", "BayesSpace", "BANKSY", "PRECAST", "Louvain", "Leiden", "Walktrap")

lARI=lapply(all_results[methods_tested_for_clustering_1],function(x) c(adjustedRandIndex(x$tissue_segments[[1]],x$clusters[[1]]),
                                                                       adjustedRandIndex(x$tissue_segments[[2]],x$clusters[[2]]),
                                                                       adjustedRandIndex(x$tissue_segments[[3]],x$clusters[[3]]),
                                                                       adjustedRandIndex(x$tissue_segments[[4]],x$clusters[[4]]),
                                                                       adjustedRandIndex(x$tissue_segments[[5]],x$clusters[[5]]),
                                                                       adjustedRandIndex(x$tissue_segments[[6]],x$clusters[[6]]),
                                                                       adjustedRandIndex(x$tissue_segments[[7]],x$clusters[[7]]),
                                                                       adjustedRandIndex(x$tissue_segments[[8]],x$clusters[[8]]),
                                                                       adjustedRandIndex(x$tissue_segments[[9]],x$clusters[[9]]),
                                                                       adjustedRandIndex(x$tissue_segments[[10]],x$clusters[[10]]),
                                                                       adjustedRandIndex(x$tissue_segments[[11]],x$clusters[[11]]),
                                                                       adjustedRandIndex(x$tissue_segments[[12]],x$clusters[[12]])))

ARI=matrix(unlist(lARI),ncol=12,byrow = TRUE)
rownames(ARI)=methods_tested_for_clustering_1
colnames(ARI)=samples_name


ari_data <- data.frame(
  Method = factor(rep(methods_tested_for_clustering_1,each=12),levels=all_methods),
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
                                                                                dist(x$feature_extraction[[5]])))$si.summary["Mean"],
                                                             summary(silhouette(as.numeric(x$tissue_segments[[6]]),
                                                                                dist(x$feature_extraction[[6]])))$si.summary["Mean"],
                                                             summary(silhouette(as.numeric(x$tissue_segments[[7]]),
                                                                                dist(x$feature_extraction[[7]])))$si.summary["Mean"],
                                                             summary(silhouette(as.numeric(x$tissue_segments[[8]]),
                                                                                dist(x$feature_extraction[[8]])))$si.summary["Mean"],
                                                             summary(silhouette(as.numeric(x$tissue_segments[[9]]),
                                                                                dist(x$feature_extraction[[9]])))$si.summary["Mean"],
                                                             summary(silhouette(as.numeric(x$tissue_segments[[10]]),
                                                                                dist(x$feature_extraction[[10]])))$si.summary["Mean"],
                                                             summary(silhouette(as.numeric(x$tissue_segments[[11]]),
                                                                                dist(x$feature_extraction[[11]])))$si.summary["Mean"],
                                                             summary(silhouette(as.numeric(x$tissue_segments[[12]]),
                                                                                dist(x$feature_extraction[[12]])))$si.summary["Mean"]
))

SILHOUETTE=matrix(unlist(lsilhouette),ncol=12,byrow = TRUE)
rownames(SILHOUETTE)=methods_tested
colnames(SILHOUETTE)=samples_name


silhouette_data <- data.frame(
  Method = factor(rep(methods_tested,each=12),levels=all_methods),
  silhouette = unlist(lsilhouette))





lDBIndex=lapply(all_results[methods_tested],function(x) c(1-DBIndex(x$feature_extraction[[1]],as.numeric(x$tissue_segments[[1]])),
                                                          1-DBIndex(x$feature_extraction[[2]],as.numeric(x$tissue_segments[[2]])),
                                                          1-DBIndex(x$feature_extraction[[3]],as.numeric(x$tissue_segments[[3]])),
                                                          1-DBIndex(x$feature_extraction[[4]],as.numeric(x$tissue_segments[[4]])),
                                                          1-DBIndex(x$feature_extraction[[5]],as.numeric(x$tissue_segments[[5]])),
                                                          1-DBIndex(x$feature_extraction[[6]],as.numeric(x$tissue_segments[[6]])),
                                                          1-DBIndex(x$feature_extraction[[7]],as.numeric(x$tissue_segments[[7]])),
                                                          1-DBIndex(x$feature_extraction[[8]],as.numeric(x$tissue_segments[[8]])),
                                                          1-DBIndex(x$feature_extraction[[9]],as.numeric(x$tissue_segments[[9]])),
                                                          1-DBIndex(x$feature_extraction[[10]],as.numeric(x$tissue_segments[[10]])),
                                                          1-DBIndex(x$feature_extraction[[11]],as.numeric(x$tissue_segments[[11]])),
                                                          1-DBIndex(x$feature_extraction[[12]],as.numeric(x$tissue_segments[[12]]))
))


DBINDEX=matrix(unlist(lDBIndex),ncol=12,byrow = TRUE)
rownames(DBINDEX)=methods_tested
colnames(DBINDEX)=samples_name


DBIndex_data <- data.frame(
  Method = factor(rep(methods_tested,each=12),levels=all_methods),
  DBIndex = unlist(lDBIndex))


#library("VGAM")
#####
#pseudoR2_McFaddenAdj = function(data,lab){
#  data=as.matrix(data)
##  fit <- multinom(lab ~ . , family = multinomial(), data = data, model = TRUE)
##  fit <- vglm(lab ~ data_scaled, family = acat(reverse = TRUE),model=TRUE)
#
#  fit <- vglm(lab ~ data,
#              family = cumulative(link = "logitlink", parallel = TRUE),
#              model = TRUE)
#  PseudoR2(fit, c("McFaddenAdj"))
#}

pseudoR2_McFaddenAdj <- function(data, lab) {
  require("ordinal")
  df <- as.data.frame(data)  # important!
  names(df) <- paste0("V", seq_len(ncol(df)))  # rename numeric columns safely
  df$lab <- factor(lab, ordered = TRUE)

  fit <- clm(lab ~ ., data = df)
  ll_full <- logLik(fit)
  ll_null <- logLik(update(fit, . ~ 1))  # intercept-only model

  k <- length(coef(fit)) + length(fit$Theta)  # coefficients + thresholds
  r2_mcfadden_adj <- 1 - ((as.numeric(ll_full) - k) / as.numeric(ll_null))

  return(r2_mcfadden_adj)
}

lpseudoR2 <- lapply(all_results[methods_tested], function(x) {
  sapply(1:12, function(i) {
    pseudoR2_McFaddenAdj(x$feature_extraction[[i]], (x$tissue_segments[[i]]))
  })
})




PSEUDOR2=matrix(unlist(lpseudoR2),ncol=12,byrow = TRUE)
rownames(PSEUDOR2)=methods_tested
colnames(PSEUDOR2)=samples_name


pseudoR2_data <- data.frame(
  Method = factor(rep(methods_tested,each=12),levels=all_methods),
  pseudoR2 = unlist(lpseudoR2))


FE_data <- data.frame(
  Method = factor(rep(methods_tested,each=12),levels=all_methods),
  DBIndex = unlist(lDBIndex),
  pseudoR2 = unlist(lpseudoR2))


###########################################################################





# Load necessary libraries
library(grid)
library(gridExtra)

method_tested_all <- c("KODAMA",  "BayesSpace","BANKSY", "PRECAST","UMAP","tSNE","PCA")


plotall=lapply(all_results[method_tested_all],function(x) plot_slide(do.call(rbind,x$feature_extraction),
                                                                     do.call(c,x$subjects),
                                                                     do.call(c,x$tissue_segments),
                                                                     col = cols_cluster,
                                                                     nrow = 1,
                                                                     size.dot = 3,
                                                                     size.legend.text=50,
                                                                     size.legend.title=50,
                                                                     size.legend.dot=20,
                                                                     size.strip.text=50)+
                 labs(y = "Your Vertical Title"))



# Define grid dimensions
ncol_value <- 1  # 5 columns for Bregma levels
nrow_value <- length(method_tested_all)  # 9 rows for Methods

# Create and save the PDF with customized layout

# Arrange the plots with a bold title and customized layout

png("output/Figures/DLPFC/Fig S7.png", width = 3000, height = 600*length(method_tested_all))

#grid.arrange(grobs = plotall, ncol = ncol_value, nrow = nrow_value,
#             padding = unit(c(2), "cm"),
#          #   top = textGrob("Spatial Domain Detection in MERFISH Dataset",
#                            gp = gpar(fontsize = 50, fontface = "bold", col = "black"))

row_titles <-c("KODAMA",  "BayesSpace","BANKSY","PRECAST","UMAP","tSNE","PCA")

# Adapt as needed
stopifnot(length(row_titles) == length(plotall))

title_width <- unit(3, "cm")  # ← increase this value for more space

# Combine each vertical title with its corresponding plot row
rows_with_titles <- mapply(function(p, title) {
  title_grob <- textGrob(title, rot = 90, gp = gpar(fontsize = 50, fontface = "bold"))
  arrangeGrob(
    title_grob, p,
    ncol = 2,
    widths = unit.c(title_width, unit(1, "npc") - title_width)
  )
}, plotall, row_titles, SIMPLIFY = FALSE)

# Stack all titled rows vertically
final_grid <- do.call(grid.arrange, c(rows_with_titles, ncol = 1))


# Draw the final plot
grid.newpage()
grid.draw(final_grid)

# Close the PDF device
dev.off()

######################################################


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
                    "Walktrap" )


plotall=lapply(all_results[method_tested_all],function(x) plot_slide(do.call(rbind,x$xy),
                                                                     do.call(c,x$samples),
                                                                     do.call(c,x$clusters),
                                                                     col = cols_cluster,
                                                                     nrow = 1,
                                                                     size.dot = 3,
                                                                     size.legend.text=50,
                                                                     size.legend.title=50,
                                                                     size.legend.dot=20,
                                                                     size.strip.text=50)+
                 labs(y = "Your Vertical Title"))



# Define grid dimensions
ncol_value <- 1  # 5 columns for Bregma levels
nrow_value <- length(method_tested_all)  # 9 rows for Methods

# Create and save the PDF with customized layout

# Arrange the plots with a bold title and customized layout

png("output/Figures/DLPFC/Fig S8.png", width = 3000, height = 600*length(method_tested_all))

#grid.arrange(grobs = plotall, ncol = ncol_value, nrow = nrow_value,
#             padding = unit(c(2), "cm"),
#          #   top = textGrob("Spatial Domain Detection in MERFISH Dataset",
#                            gp = gpar(fontsize = 50, fontface = "bold", col = "black"))

row_titles <- c("ground truth",
                "KODAMA",
                "BASS",
                "BayesSpace",
                "BANKSY",
                "PRECAST",
                "Louvain",
                "Leiden",
                "Walktrap" )  # Adapt as needed
stopifnot(length(row_titles) == length(plotall))

title_width <- unit(3, "cm")  # ← increase this value for more space

# Combine each vertical title with its corresponding plot row
rows_with_titles <- mapply(function(p, title) {
  title_grob <- textGrob(title, rot = 90, gp = gpar(fontsize = 50, fontface = "bold"))
  arrangeGrob(
    title_grob, p,
    ncol = 2,
    widths = unit.c(title_width, unit(1, "npc") - title_width)
  )
}, plotall, row_titles, SIMPLIFY = FALSE)

# Stack all titled rows vertically
final_grid <- do.call(grid.arrange, c(rows_with_titles, ncol = 1))


# Draw the final plot
grid.newpage()
grid.draw(final_grid)

# Close the PDF device
dev.off()






















####################################################################################################

# Create the box plot
Nplot1 <- ggboxplot(ari_data, x = "Method", y = "ARI", palette = method_colors_all,
                    fill = "Method", add = "jitter",
                    add.params = list(size = 2, jitter = 0.2, fill = 3, shape = 21)) +
  ylab("Adjusted Rand Index") +
  labs(title = "",
       x = "Clustering Method") +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt x-axis labels for better readability
  ylim(0, NA) # Set the start of the y-axis at 0


svg("output/Figures/DLPFC/Figure 3 ARI.svg")
# Display the plot
print(Nplot1)
dev.off()





# Create the scatter plot with well-defined axes
Nplot2 <- ggplot(FE_data, aes(x = pseudoR2, y = DBIndex, color = Method, shape = Method)) +
  geom_point(size = 5,stroke = 1.5) +
  scale_shape_manual(values = c(16, 17, 18, 19, 8, 15, 3)) +
  scale_color_manual(values = method_colors_all) +

  scale_x_continuous(
    name = "Pseudo R^2 (Coefficient of Determination)",  # Custom axis label
    limits = c(0.2, 0.85),
    breaks = seq(0.2, 0.8, 0.2)
  ) +
  # Set continuous scale for y-axis
  scale_y_continuous(
    name = "1- DB Index",
    limits = c(0.2, 0.85),
    breaks = seq(0.2, 0.8, 0.2)
  ) +
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

svg("output/Figures/DLPFC/Figure 4 - Scatter Plot of Pseudo R^2 vs DB Index .svg",height = 4)
print(Nplot2)
dev.off()


############################################################



methods <- c("Real Labels", "KODAMA", "BASS", "BayesSpace", "BANKSY", "PRECAST", "Louvain", "Leiden", "Walktrap")

slide_number=11
colnames(result_BayesSpace$xy[[slide_number]])=colnames(result_KODAMA$xy[[slide_number]])
colnames(result_BANKSY$xy[[slide_number]])=colnames(result_KODAMA$xy[[slide_number]])
colnames(result_PRECAST$xy[[slide_number]])=colnames(result_KODAMA$xy[[slide_number]])
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
svg("output/Figures/DLPFC/Figure 5_Comparative_Methods.svg", width = 11, height = 13)

Ngplot3

# Close the PDF device
dev.off()
























######################################################


plot_slide <- function(xy, slide, labels, col = NULL, nrow = 1, scales = "free", size.dot = 3,
                       size.legend.text = 10, size.legend.title = 10, size.legend.dot = 5, size.strip.text = 10) {

  if (is.null(col)) {
    labels <- as.factor(labels)
    nn <- length(levels(labels))
    col <- rainbow(nn)
  }

  df <- data.frame(xy, slide, labels)
  colnames(df) <- c("x", "y", "slide", "labels")
  df$slide <- as.factor(df$slide)
  df$labels <- as.factor(df$labels)

  ggplot(df, aes(x, y, color = labels)) +
    geom_point(size = size.dot) +
    facet_wrap(~slide, nrow = nrow, scales = scales) +
    theme_bw() +
    theme(
      legend.position = "none",  # Remove legend
      axis.title = element_blank(),
      strip.text = element_text(size = size.strip.text, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_color_manual(values = col, drop = FALSE)
}




########PLOT 12 SLIDE

load("output/DLFPC-All.RData")

MALAT1=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "MALAT1"]
PPP3CA=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "PPP3CA"]
CAMK2N1=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "CAMK2N1"]

CCK=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "CCK"]
RGS4=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "RGS4"]
CALM2=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "CALM2"]

MBP=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "MBP"]
PLP1=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "PLP1"]
CNP=spe@rowRanges@elementMetadata$gene_id[spe@rowRanges@elementMetadata$gene_name %in% "CNP"]


MALAT1_values=data[subjects=="Br5595",MALAT1]
PPP3CA_values=data[subjects=="Br5595",PPP3CA]
CAMK2N1_values=data[subjects=="Br5595",CAMK2N1]

CCK_values=data[subjects=="Br5595",CCK]
RGS4_values=data[subjects=="Br5595",RGS4]
CALM2_values=data[subjects=="Br5595",CALM2]


MBP_values=data[subjects=="Br5595",MBP]
PLP1_values=data[subjects=="Br5595",PLP1]
CNP_values=data[subjects=="Br5595",CNP]


maxx=max(c(MALAT1_values,PPP3CA_values,CAMK2N1_values,
           CCK_values,RGS4_values,CALM2_values,
           MBP_values,PLP1_values,CNP_values
           ))

plotall=list(MALAT1=list(),
             PPP3CA=list(),
             CAMK2N1=list(),
             CCK=list(),
             RGS4=list(),
             CALM2=list(),
             MBP=list(),
             PLP1=list(),
             CNP=list())

MALAT1_values=round(255*MALAT1_values/maxx)+1
PPP3CA_values=round(255*PPP3CA_values/maxx)+1
CAMK2N1_values=round(255*CAMK2N1_values/maxx)+1
CCK_values=round(255*CCK_values/maxx)+1
RGS4_values=round(255*RGS4_values/maxx)+1
CALM2_values=round(255*CALM2_values/maxx)+1
MBP_values=round(255*MBP_values/maxx)+1
PLP1_values=round(255*PLP1_values/maxx)+1
CNP_values=round(255*CNP_values/maxx)+1


plotall$MALAT1=plot_slide(xy_Br5595,
                       samples_Br5595,
                       MALAT1_values,
                       col =  viridis(256),
                       nrow = 1,
                       size.dot = 3,
                       size.legend.text=50,
                       size.legend.title=50,
                       size.legend.dot=20,
                       size.strip.text=50)

plotall$PPP3CA=plot_slide(xy_Br5595,
                       samples_Br5595,
                       PPP3CA_values,
                       col =  viridis(256),
                       nrow = 1,
                       size.dot = 3,
                       size.legend.text=50,
                       size.legend.title=50,
                       size.legend.dot=20,
                       size.strip.text=50)



plotall$CAMK2N1=plot_slide(xy_Br5595,
                           samples_Br5595,
                           CAMK2N1_values,
                           col =  viridis(256),
                           nrow = 1,
                           size.dot = 3,
                           size.legend.text=50,
                           size.legend.title=50,
                           size.legend.dot=20,
                           size.strip.text=50)


plotall$CCK=plot_slide(xy_Br5595,
                           samples_Br5595,
                       CCK_values,
                           col =  viridis(256),
                           nrow = 1,
                           size.dot = 3,
                           size.legend.text=50,
                           size.legend.title=50,
                           size.legend.dot=20,
                           size.strip.text=50)


plotall$RGS4=plot_slide(xy_Br5595,
                           samples_Br5595,
                        RGS4_values,
                           col =  viridis(256),
                           nrow = 1,
                           size.dot = 3,
                           size.legend.text=50,
                           size.legend.title=50,
                           size.legend.dot=20,
                           size.strip.text=50)



plotall$CALM2=plot_slide(xy_Br5595,
                         samples_Br5595,
                         CALM2_values,
                         col =  viridis(256),
                         nrow = 1,
                         size.dot = 3,
                         size.legend.text=50,
                         size.legend.title=50,
                         size.legend.dot=20,
                         size.strip.text=50)


plotall$MBP=plot_slide(xy_Br5595,
                         samples_Br5595,
                         MBP_values,
                         col =  viridis(256),
                         nrow = 1,
                         size.dot = 3,
                         size.legend.text=50,
                         size.legend.title=50,
                         size.legend.dot=20,
                         size.strip.text=50)


plotall$PLP1=plot_slide(xy_Br5595,
                         samples_Br5595,
                         PLP1_values,
                         col =  viridis(256),
                         nrow = 1,
                         size.dot = 3,
                         size.legend.text=50,
                         size.legend.title=50,
                         size.legend.dot=20,
                         size.strip.text=50)


plotall$CNP=plot_slide(xy_Br5595,
                         samples_Br5595,
                       CNP_values,
                         col =  viridis(256),
                         nrow = 1,
                         size.dot = 3,
                         size.legend.text=50,
                         size.legend.title=50,
                         size.legend.dot=20,
                         size.strip.text=50)



# Load necessary libraries
library(grid)
library(gridExtra)



row_titles=c("MALAT1",
                    "PPP3CA",
                    "CAMK2N1",
                    "CCK",
                    "RGS4",
                    "CALM2",
                    "MBP",
                    "PLP1",
                    "CNP")



# Define grid dimensions
ncol_value <- 1  # 5 columns for Bregma levels
nrow_value <- length(method_tested_all)  # 9 rows for Methods

# Create and save the PDF with customized layout

# Arrange the plots with a bold title and customized layout

png("output/Figures/DLPFC/Fig Sgene.png", width = 2000, height = 600*length(method_tested_all))

#grid.arrange(grobs = plotall, ncol = ncol_value, nrow = nrow_value,
#             padding = unit(c(2), "cm"),
#          #   top = textGrob("Spatial Domain Detection in MERFISH Dataset",
#                            gp = gpar(fontsize = 50, fontface = "bold", col = "black"))


stopifnot(length(row_titles) == length(plotall))

title_width <- unit(3, "cm")  # ← increase this value for more space

# Combine each vertical title with its corresponding plot row
rows_with_titles <- mapply(function(p, title) {
  title_grob <- textGrob(title, rot = 90, gp = gpar(fontsize = 50, fontface = "bold"))
  arrangeGrob(
    title_grob, p,
    ncol = 2,
    widths = unit.c(title_width, unit(1, "npc") - title_width)
  )
}, plotall, row_titles, SIMPLIFY = FALSE)

# Stack all titled rows vertically
final_grid <- do.call(grid.arrange, c(rows_with_titles, ncol = 1))


# Draw the final plot
grid.newpage()
grid.draw(final_grid)

# Close the PDF device
dev.off()


