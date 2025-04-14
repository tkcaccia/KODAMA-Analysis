library(ggpubr)
library(mclust)
library(cluster)
#BiocManager::install("RDRToolbox")
library("RDRToolbox")
library(nnet)
library(DescTools)
#setwd("/Users/stefano/HPC-scratch/KODAMA-Analysis/")



plot_slide =  function (xy, slide, labels, col = NULL, nrow = 1, scales = "free",size.dot = 3,
                        size.legend.text=10,size.legend.title=10,size.legend.dot=5,size.strip.text=10)
{
  if (is.null(col)) {
    labels = as.factor(labels)
    nn = length(levels(labels))
    col = rainbow(nn)
  }
  df <- data.frame(xy, slide, labels)
  colnames(df) = c("x", "y", "slide", "labels")
  df$slide = as.factor(df$slide)
  df$labels = as.factor(df$labels)
  ggplot(df, aes(x, y, color = labels)) + geom_point(size = size.dot) +
    facet_wrap(~slide, nrow = nrow, scales = scales) + theme_bw() +
    theme(legend.position = "bottom", axis.title = element_blank(),
          legend.text = element_text(size = size.legend.text),
          legend.title = element_text(size = size.legend.title),
          strip.text = element_text(size = size.strip.text, face = "bold"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    scale_color_manual("Domain",  values = col) + guides(color = guide_legend(nrow = 1,  override.aes = list(size = size.legend.dot)))  # Legend
}





library(KODAMAextra)
#####################################################################################

load("output/MERFISH-BANSKY-results.RData")
load("output/MERFISH-BASS-results.RData")
load("output/MERFISH-BayesSpace-results.RData")
load("output/MERFISH-KODAMA-results.RData")
load("output/MERFISH-Nonspatial-results.RData")
load("output/MERFISH-PRECAST-results.RData")



all_results=list(Real=ground_true,
                 KODAMA=result_KODAMA,
                 BASS=result_BASS,
                 BayesSpace=result_BayesSpace,
                 BANKSY=result_BANKSY,
                 PRECAST=result_PRECAST,
                 Louvain=result_Louvain,
                 Leiden=result_Leiden,
                 Walktrap=result_Walktrap,
                 "Louvain-APM"=result_Louvain.PM,
                 "Leiden-APM"=result_Leiden.PM,
                 "Walktrap-APM"=result_Walktrap.PM,
                 "Louvain-APM-SVM"=result_Louvain.PM.ref,
                 "Leiden-APM-SVM"=result_Leiden.PM.ref,
                 "Walktrap-APM-SVM"=result_Walktrap.PM.ref,
                 UMAP=result_UMAP,
                 "tSNE"=result_tSNE,
                 PCA=result_PCA,
                 "UMAP-APM"=result_UMAP.PM,
                 "tSNE-APM"=result_tSNE.PM,
                 "PCA-APM"=result_PCA.PM)


cols <- c("#669bbc", "#81b29a", "#f2cc8f", "#adc178",
          "#dde5b6", "#a8dadc", "#e5989b", "#e07a5f",
          "#aae5b6", "#a8aadc", "#e59811", "#aa7900")

####FIGURE 1 B
sele=c("PCA","tSNE","UMAP","PCA-APM","tSNE-APM","UMAP-APM")
xx=all_results[sele]
met=NULL
labels=NULL
slide=NULL
for(i in 1:length(xx)){
  met=rbind(met,xx[[i]]$feature_extraction_allslide)
  labels=c(labels,as.vector(xx[[i]]$tissue_segments_allslide))
  slide=c(slide,rep(names(xx)[i],length(xx[[i]]$tissue_segments_allslide)))
}
slide=factor(slide,levels=sele)




p1=plot_slide(met,slide,labels,col=cols,nrow = 2,size.dot = 3,
              size.legend.text=50,size.legend.title=50,size.legend.dot=20,size.strip.text=50)

png("output/Figures/MERFISH/Fig S2.png",width = 1800,height = 1200)
p1
dev.off()


######################################

bregma_levels=c("-0.04 mm", "-0.09 mm", "-0.14 mm", "-0.19 mm", "-0.24 mm")


#sele=c("KODAMA","PCA","tSNE","UMAP","BayesSpace","BANKSY","PRECAST")
#sele2=paste(rep(sele,each=5),rep(bregma_levels,length(sele)),sep=" ")

#met=NULL
#labels=NULL
#slide=NULL
#for(i in 1:length(sele)){
#  for(j in 1:5){
#    met=rbind(met,all_results[[sele[i]]]$feature_extraction[[j]])
#    temp=as.vector(all_results[[sele[i]]]$tissue_segments[[j]])
#    labels=c(labels,temp)
#    slide=c(slide,rep(sele2[j+(i-1)*5],length(temp)))
#  }
#}
#slide=factor(slide,levels=sele2)

#p2=plot_slide(met,slide,labels,col=cols,nrow = length(sele),size.dot = 3,
#              size.legend.text=50,size.legend.title=50,size.legend.dot=20,size.strip.text=50)

#png("output/Figures/MERFISH/Fig S1.png",width = 3000,height = 600*length(sele))
#p2
#dev.off()


sele=c("KODAMA","PCA","tSNE","UMAP","BayesSpace","BANKSY","PRECAST")
xx=all_results[sele]
met=NULL
labels=NULL
slide=NULL
for(i in 1:length(xx)){
  met=rbind(met,xx[[i]]$feature_extraction_allslide)
  labels=c(labels,as.vector(xx[[i]]$tissue_segments_allslide))
  slide=c(slide,rep(names(xx)[i],length(xx[[i]]$tissue_segments_allslide)))
}
slide=factor(slide,levels=sele)




p1=plot_slide(met,slide,labels,col=cols,nrow = 3,size.dot = 3,
              size.legend.text=50,size.legend.title=50,size.legend.dot=20,size.strip.text=50)

png("output/Figures/MERFISH/Fig S1.png",width = 1800,height = 1800)
p1
dev.off()










#######################################

####FIGURE C AND D



# Define method colors

method_colors_all <- c(
  "Real" = "#FF5733aa",
  "KODAMA" = "#99CC37aa",
  "BASS" = "#3357FFaa",
  "BayesSpace" = "#FF33A1aa",
  "BANKSY" = "#CCBB33aa",
  "PRECAST" = "#33aaa0aa",
  "Louvain" = "#B833FFaa",
  "Leiden" = "#8D8833aa",
  "Walktrap" = "#FF8C33aa",
  "UMAP" = "#5733FFaa",
  "tSNE" = "#aa7733aa",
  "PCA" = "#ccaaffaa",


  "Louvain-APM" = "#33A1FFaa",# Bleu clair
  "Leiden-APM" = "#FF33E6aa", # Rose magenta
  "Walktrap-APM" = "#A133FFaa",# Pourpre
  "Louvain-APM-SVM" = "#FF3333aa",# Rouge intense
  "Leiden-APM-SVM" = "#33FF8Caa", # Vert émeraude
  "Walktrap-APM-SVM" = "#FF5733aa",# Rouge mandarine
  "UMAP-APM" = "#33D1FFaa",    # Bleu azur
  "tSNE-APM" = "#C433FFaa",    # Violet profond
  "PCA-APM" = "#FF7133aa"      # Orange rougeâtre
)


all_methods=names(method_colors_all)




methods_tested_for_clustering <- c("KODAMA", "BASS", "BayesSpace", "BANKSY", "PRECAST", "Louvain", "Leiden", "Walktrap",
                                   "Louvain-APM","Leiden-APM" ,"Walktrap-APM","Louvain-APM-SVM" , "Leiden-APM-SVM","Walktrap-APM-SVM")


lARI=lapply(all_results[methods_tested_for_clustering],function(x) c(adjustedRandIndex(x$tissue_segments[[1]],x$clusters[[1]]),
                                 adjustedRandIndex(x$tissue_segments[[2]],x$clusters[[2]]),
                                 adjustedRandIndex(x$tissue_segments[[3]],x$clusters[[3]]),
                                 adjustedRandIndex(x$tissue_segments[[4]],x$clusters[[4]]),
                                 adjustedRandIndex(x$tissue_segments[[5]],x$clusters[[5]])))

ARI=matrix(unlist(lARI),ncol=5,byrow = TRUE)
rownames(ARI)=methods_tested_for_clustering
colnames(ARI)=bregma_levels




################################################################



methods_tested <- c("KODAMA",  "BayesSpace","BANKSY", "PRECAST","UMAP","tSNE","PCA","UMAP-APM","tSNE-APM","PCA-APM")


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
lpseudoR2=lapply(all_results[methods_tested],function(x) c(pseudoR2_McFaddenAdj(x$feature_extraction[[1]],as.factor(x$tissue_segments[[1]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[2]],as.factor(x$tissue_segments[[2]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[3]],as.factor(x$tissue_segments[[3]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[4]],as.factor(x$tissue_segments[[4]])),
                                                       pseudoR2_McFaddenAdj(x$feature_extraction[[5]],as.factor(x$tissue_segments[[5]]))
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




ari_data <- data.frame(
  Method = factor(rep(methods_tested_for_clustering,each=5),levels=all_methods),
  ARI = unlist(lARI))

selARI=ari_data$Method %in% c( "KODAMA","BASS","BayesSpace","BANKSY","PRECAST","Louvain","Leiden"   ,"Walktrap"    )

ari_data_sel=ari_data[selARI,]

# Create the box plot
Nplot1 <- ggboxplot(ari_data_sel, x = "Method", y = "ARI", palette = method_colors_all,
                    fill = "Method", add = "jitter",
                    add.params = list(size = 2, jitter = 0.2, fill = 3, shape = 21)) +
  ylab("ARI") +
  labs(title = "",
       x = "Clustering Method") +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt x-axis labels for better readability
  ylim(0, NA) # Set the start of the y-axis at 0


pdf("output/Figures/MERFISH/Fig 2-1.pdf")
# Display the plot
print(Nplot1)
dev.off()




selARI=ari_data$Method %in% c("Louvain","Leiden","Walktrap","Louvain-APM","Leiden-APM","Walktrap-APM","Louvain-APM-SVM","Leiden-APM-SVM","Walktrap-APM-SVM")

ari_data_sel=ari_data[selARI,]

# Create the box plot
Nplot1 <- ggboxplot(ari_data_sel, x = "Method", y = "ARI", palette = method_colors_all,
                    fill = "Method", add = "jitter",
                    add.params = list(size = 2, jitter = 0.2, fill = 3, shape = 21)) +
  ylab("ARI") +
  labs(title = "",
       x = "Clustering Method") +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt x-axis labels for better readability
  ylim(0, NA) # Set the start of the y-axis at 0


pdf("output/Figures/MERFISH/Fig S5.pdf")
# Display the plot
print(Nplot1)
dev.off()






sel_FE=FE_data$Method %in% c("KODAMA",  "BayesSpace","BANKSY", "PRECAST","UMAP","tSNE","PCA")

FE_data_sel=FE_data[sel_FE,]

# Create the scatter plot with well-defined axes
Nplot2 <- ggplot(FE_data_sel, aes(x = pseudoR2, y = DBIndex, color = Method, shape = Method)) +
  geom_point(size = 5,stroke = 1.5) +
  scale_shape_manual(values = c(16, 17, 18, 19, 8, 15, 3)) +
  scale_color_manual(values = method_colors_all) +

  scale_x_continuous(
    name = "Pseudo R^2",  # Custom axis label
    limits = c(0.1, 0.75),
    breaks = seq(0.2, 0.8, 0.2)
  ) +
  # Set continuous scale for y-axis
  scale_y_continuous(
    name = "1- DB Index",
    limits = c(0.1, 0.75),
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

pdf("output/Figures/MERFISH/Fig 2-2.pdf",height = 4)
print(Nplot2)
dev.off()








###########################################################################

### FIGURE E

########PLOT 1 SLIDE

##################################
# Load necessary libraries
library(grid)
library(gridExtra)

# Define Bregma levels and clustering methods
methods <- c("Real Labels", "KODAMA", "BASS", "BayesSpace", "BANKSY", "PRECAST", "Louvain", "Leiden", "Walktrap")

# Colors for each cluster with a softer palette
#cols_cluster <- c("#0000b6", "#81b29a", "#f2cc8f", "#e07a5f", "#81ccff", "#cc00b6", "#33b233", "#ff5733", "#FFD111")

cols_cluster <- c("#81b29a","#adc178", "#f2cc8f", "#e5989b",
          "#dde5b6", "#a8dadc", "#669bbc", "#e07a5f","#0000b6")

cols_cluster <- c("#669bbc", "#81b29a", "#f2cc8f", "#adc178",
          "#dde5b6", "#a8dadc", "#e5989b", "#e07a5f",
          "#aae5b6", "#a8aadc", "#e59811", "#aa7900")

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


cols <- c("#669bbc", "#81b29a", "#f2cc8f", "#adc178",
          "#dde5b6", "#a8dadc", "#e5989b", "#e07a5f",
          "#aae5b6", "#a8aadc", "#e59811", "#aa7900")

Ngplot3 <- plot_slide(xy_slide_comparison,method_slide_comparison,
                      as.factor(cluster_slide_comparison),col = cols, nrow = 3,size.dot = 1) +
  ggtitle("Bregma -0.09") + theme(plot.title = element_text(size = 15, face = "bold"))



# Create and save the PDF with customized layout
pdf("output/Figures/MERFISH/Fig 2e.pdf", width = 10, height = 10)

Ngplot3

# Close the PDF device
dev.off()


######################################

#######################################################


### FIGURE 2

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


plotall=lapply(all_results[method_tested_all],function(x) plot_slide(x$xyz[,-3],
                                                                     x$xyz[,3],
                                                                     x$clusters_allslide,
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

png("output/Figures/MERFISH/Fig S4.png", width = 3000, height = 600*length(method_tested_all))

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

#################################################





######################################################################

#####################################

### FIGURE 4

########PLOT 12 SLIDE



# Load necessary libraries
library(grid)
library(gridExtra)


method_tested_all <- c(
  "Louvain",
  "Leiden",
  "Walktrap",
  "Louvain-APM",
  "Leiden-APM",
  "Walktrap-APM",
  "Louvain-APM-SVM",
  "Leiden-APM-SVM",
  "Walktrap-APM-SVM"
)


plotall=lapply(all_results[method_tested_all],function(x) plot_slide(x$xyz[,-3],
                                                                     x$xyz[,3],
                                                                     x$clusters_allslide,
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

png("output/Figures/MERFISH/Fig S6.png", width = 3000, height = 600*length(method_tested_all))

#grid.arrange(grobs = plotall, ncol = ncol_value, nrow = nrow_value,
#             padding = unit(c(2), "cm"),
#          #   top = textGrob("Spatial Domain Detection in MERFISH Dataset",
#                            gp = gpar(fontsize = 50, fontface = "bold", col = "black"))

row_titles <- c(
  "Louvain",
  "Leiden",
  "Walktrap",
  "Louvain-APM",
  "Leiden-APM",
  "Walktrap-APM",
  "Louvain-APM-SVM",
  "Leiden-APM-SVM",
  "Walktrap-APM-SVM")  # Adapt as needed
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

###################
##################


####FIGURE B


######################################














sel_FE=FE_data$Method %in% c("UMAP","tSNE","PCA","UMAP-APM","tSNE-APM","PCA-APM")

FE_data_sel=FE_data[sel_FE,]

# Create the scatter plot with well-defined axes
Nplot2 <- ggplot(FE_data_sel, aes(x = pseudoR2, y = DBIndex, color = Method, shape = Method)) +
  geom_point(size = 5,stroke = 1.5) +
  scale_shape_manual(values = c(16, 17, 18, 19, 8, 15, 3)) +
  scale_color_manual(values = method_colors_all) +

  scale_x_continuous(
    name = "Pseudo R^2",  # Custom axis label
    limits = c(0.1, 0.75),
    breaks = seq(0.2, 0.8, 0.2)
  ) +
  # Set continuous scale for y-axis
  scale_y_continuous(
    name = "1- DB Index",
    limits = c(0.1, 0.75),
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

pdf("output/Figures/MERFISH/Fig S3.pdf",height = 4)
print(Nplot2)
dev.off()






