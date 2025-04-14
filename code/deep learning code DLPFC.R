A set of histological image tiles from subject “Br5595” was selected, focusing on layer 3 with reference classes 1, 4, and 5, and balanced by randomly drawing 700 tiles per class for training while the remaining tiles formed the test set. Each tile was stored as a 128×128 RGB TIFF image, and a custom R function read these files into a four-dimensional tensor suitable for convolutional input. The network architecture comprised two convolutional blocks, each with ReLU activation and 16 or 32 filters, followed by 2×2 max pooling to reduce spatial dimensions. Rather than using a flattening layer, global average pooling aggregated feature maps into a single value per channel, and a final dense layer with a softmax activation classified tiles into three possible reference classes. The model was compiled using the Adam optimizer and a sparse categorical crossentropy loss, and was trained for 20 epochs. Performance was evaluated on a held-out test set, where accuracy and loss were measured, and class predictions were compared to ground-truth labels via a confusion matrix.

setwd("/Users/stefano/Downloads/Tiles2/")


load("/Users/stefano/Downloads/DL.RData")


library(keras)

library(tiff)


# Load TIFF images from Layer 3

load_tiles <- function(samples) {
  tiff_data <- array(dim  =c(length(samples),128,128,3))
  for (i in 1:length(samples)) {
    file_name <- paste0(samples[i], ".tiff",sep="")
    if (file.exists(file_name)) {
      temp=readTIFF(file_name)[,,]
      tiff_data[i,,,] <- temp # Use rast() instead of raster()
    } else {
      cat("Missing file:", file_name, "\n")
    }
  }
  return(tiff_data)
}



setwd("/Users/stefano/Downloads/Tiles2/")





sel_t=names(ref[sel_subcluster_layer3])

ref2=ref
ref2[ref2==5]=3

acc=NULL
for (i in 1:100){
ss=c(sample(which(ref2[sel_t]==1),2566),sample(which(ref2[sel_t]==3),2566))



tiles_training=names(labels_clear[sel_t][ss])
training_labels=as.numeric(as.factor(as.vector(ref2[sel_t][ss])))-1

tiles_test=names(labels_clear[sel_t][-ss])
test_labels=as.numeric(as.factor(as.vector(ref2[sel_t][-ss])))-1
test_labels2=as.vector(ref2[sel_t][-ss])

trainingset=load_tiles(tiles_training)
testset=load_tiles(tiles_test)




model <- keras_model_sequential() %>%
  # First Convolutional Block with fewer filters
  layer_conv_2d(filters = 16, kernel_size = c(4,4), activation = 'relu',
                input_shape = c(128, 128, 3)) %>%
  layer_max_pooling_2d(pool_size = c(2,2)) %>%

  # Second Convolutional Block with fewer filters
  layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu') %>%
  layer_max_pooling_2d(pool_size = c(2,2)) %>%

  # Global Average Pooling replaces Flatten + Dense
  layer_global_average_pooling_2d() %>%

  # Final Classification Layer
  layer_dense(units = 2, activation = 'softmax')



model %>% compile(
  optimizer = 'adam',
  loss = 'sparse_categorical_crossentropy',
  metrics = c('accuracy')
)



model %>% fit(trainingset, training_labels, epochs = 20, verbose = 2)

score <- model %>% evaluate(testset, test_labels, verbose = 0)
cat('Test loss:', score["loss"], "\n")
cat('Test accuracy:', score["accuracy"], "\n")

predictions <- model %>% predict(testset)
predicted_classes <- apply(predictions, 1, which.max) - 1

table(predicted_classes)
table(test_labels)

table(predicted_classes,test_labels2)

table(predicted_classes,test_labels)

ta=table(predicted_classes>0,test_labels>0)
ta
acc[i]=sum(diag(ta))/sum(ta)
print(acc[i])
}


