library(readxl)
library(glmnet)
library(mclust)
library(caret)
library(ncvreg)

# Load the data from the .mat file
data <- read_excel("C:/Users/dell/Desktop/FCR/Empirical application/data.xlsx", sheet = "sbp852")
# Extract the first column as the dependent variable
y <- data[, 1]
y <- as.matrix(y)
# Extract all columns except the first one as independent variables
X <- data[, -1] 
X <- as.matrix(X)

# Define parameters for cross-validation
folds <- 20  # Number of folds
repeats <- 25  # Number of repeats

# Define parameters for k-means clustering
num_clusters <- 4  # Number of clusters

# Define the minimum sample threshold
min_samples_threshold <- 10

# Store the linear regression models for each cluster
models <- list()

# Store the MSE for each fold of each repeat
k_mse <- numeric(folds * repeats)

# Repeat cross-validation
for (i in 1:repeats) {
  # Split the dataset
  folds_indices <- createMultiFolds(y, k = folds, times = repeats)
  lasso_mse <- numeric(folds)  # Store MSE for each fold
  
  # Iterate over each fold for cross-validation
  for (j in 1:folds) {
    # Get training and testing indices
    train_indices <- unlist(folds_indices[j])
    test_indices <- setdiff(1:length(y), train_indices)
    X_train <- X[train_indices, ]
    y_train <- y[train_indices]
    X_test <- X[test_indices, ]
    y_test <- y[test_indices]
    
    # Perform k-means clustering on the training set
    kmeans_result <- kmeans(X_train, centers = num_clusters)
    cluster_labels <- kmeans_result$cluster
    
    # Fit linear regression models for each cluster
    for (k in 1:num_clusters) {
      # Get indices of data points in the current cluster
      cluster_indices <- which(cluster_labels == k)
      
      # Check if the number of samples in the cluster is below the threshold
      if (length(cluster_indices) < min_samples_threshold) {
        next  # Skip this cluster if below threshold
      }
      
      # Extract data for the current cluster
      X_cluster <- X_train[cluster_indices, ]
      y_cluster <- y_train[cluster_indices]
      
      # Fit linear regression model
      fit_success <- TRUE
      tryCatch({
        lasso_fit <- suppressWarnings(ncvreg::cv.ncvreg(X_cluster, y_cluster, penalty = "lasso", nfolds = 5, maxit = 2000))
      }, warning = function(w) {
        fit_success <<- FALSE
      })
      
      if (!fit_success) {
        next  # Skip to next iteration if fitting failed
      }
      
      # 在测试集上进行预测
      predictions <- predict(lasso_fit$fit, X_test, s = lasso_fit$lambda.min)
      
      # 计算 MSE
      lmse <- mean((y_test - predictions)^2)
      lasso_mse[j] <- lmse
    }
  }
  # Store the MSE for this repeat
  start_index <- (i - 1) * folds + 1
  end_index <- start_index + folds - 1
  k_mse[start_index:end_index] <- lasso_mse
}

# Print or use k_mse as needed
k_mse
mean(k_mse)

# Save the k_mse vector
save(k_mse, file = "kmeans_CV_MSE.RData")


    


  

