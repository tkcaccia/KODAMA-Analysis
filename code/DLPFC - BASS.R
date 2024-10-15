

# Load the necessary BASS library
library(BASS)

load("data/DLPFC-general.RData")

# Initialize the subject names
subject_names <- c("Br5292", "Br5595", "Br8100")

# Initialize the BASS list
results_BASS <- list()
results_BASS$cluster <- list()

# Loop over each subject

for(i in 1:3){
  ncluster=sum(table(labels[[i]])>0)
  
  # Create lists of xy coordinates and data

  data <- data_list[(1:4)+4*(i-1)]
  data=lapply(data, function(x) t(as.matrix(x)))
  xy <- xy_list[(1:4)+4*(i-1)]
  
  # Hyperparameters
  C <- 20  # Number of cell types
  R <- ncluster   # Number of spatial domains
  
  set.seed(0)
  
  # Create the BASS object
  res <- createBASSObject(data , xy, C = C, R = R,
                          beta_method = "SW", init_method = "mclust", 
                          nsample = 10000)
  
  # Add PCA results
  res@X_run <- t(pca_subject[[i]][, 1:30])
  
  # Run BASS
  res <- BASS.run(res)
  
  # Post-process BASS
  res <- BASS.postprocess(res)
  
  # Get the cluster labels
  zlabels <- res@results$z
  
  results_BASS$cluster = c(results_BASS$cluster, zlabels)
  

}


# Save the results
save(results_BASS, file = "output/BASS-results.RData")
