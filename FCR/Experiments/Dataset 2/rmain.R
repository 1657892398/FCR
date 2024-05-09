# Load required libraries
library(R.matlab)
library(WCLR)
library(mclust)
library(caret)

# Set working directory
setwd("D:/OneDrive - mail.sdu.edu.cn/机器学习实验/r语言代码非线性版本")

# Load the data from the .mat file
da <- readMat("data1.mat")
x <- da$x.data
y <- da$y.data
label <- da$label.data

# Set the values of k, alpha, m, and gamma
k <- 3 # 分类种数
m <- 1 # 模糊值（这里一直设置成1就行）

# Function to compute RMSE (Root Mean Squared Error)
compute_rmse <- function(x, y, w, b) {
  N <- nrow(x)
  C <- ncol(b)
  mse <- 0
  
  # Create a matrix of ones with the same number of rows as x
  ones_column <- matrix(1, nrow = N, ncol = 1)
  
  # Concatenate the ones_column with the original matrix x
  x <- cbind(ones_column, x)
  
  for (i in 1:N) {
    yy <- 0
    for (j in 1:C) {
      yy <- yy + w[i, j] * t(x[i,]) %*% b[, j]
    }
    mse <- mse + (y[i] - yy)^2
  }
  
  mse <- sqrt(mse / N)
  return(mse)
}


# Perform 3*5 cross-validation using caret package
n_folds <- 5
n_repeats <- 3
n_iterations <- 100


# Create a custom cross-validation function
custom_cv <- function (x, y, gamma, clustering_function) {
  all_ari <- c()
  all_rmse <- c()
  all_pre <- c()
  for (repeats in 1:n_repeats) {
    indices <- createFolds(1:length(label), k = n_folds, list = TRUE, returnTrain = FALSE)
    for (fold in 1:n_folds) {
      # Split data into training and validation sets
      validation_indices <- indices[[fold]]
      training_indices <- setdiff(1:length(label), validation_indices)
      
      x_train <- x[training_indices, ]
      y_train <- y[training_indices, ]
      label_train <- label[training_indices,]
      x_valid <- x[validation_indices, ]
      y_valid <- y[validation_indices, ]
      
      # Perform 100 runs of clr and find the best one based on loss
      best_clr <- NULL
      best_loss <- Inf
      for (iter in 1:n_iterations) {
        clr_result <- wclr(x_train, y_train, k, alpha = gamma,wnorm = clustering_function, m)
        if (clr_result$loss < best_loss) {
          best_clr <- clr_result
          best_loss <- clr_result$loss
        }
      }
      # Compute ARI and RMSE
      ari <- adjustedRandIndex(label_train, best_clr$cluster)
      rmse <- compute_rmse(x_train,y_train,best_clr$membership ,best_clr$coefficients)
      rmse_pre <- compute_rmse(x_valid,y_valid,best_clr$membership ,best_clr$coefficients)
      
      # Store the results
      all_ari <- c(all_ari,ari)
      all_rmse <- c(all_rmse,rmse)
      all_pre <- c(all_pre,rmse_pre)
    }
  }
  return(list(ARI = all_ari, RMSE = all_rmse, PRE = all_pre))
}

custom_cv_kplane <- function (x, y, gamma) {
  all_ari <- c()
  all_rmse <- c()
  all_pre <- c()
  for (repeats in 1:n_repeats) {
    indices <- createFolds(1:length(label), k = n_folds, list = TRUE, returnTrain = FALSE)
    for (fold in 1:n_folds) {
      # Split data into training and validation sets
      validation_indices <- indices[[fold]]
      training_indices <- setdiff(1:length(label), validation_indices)
      
      x_train <- x[training_indices, ]
      y_train <- y[training_indices, ]
      label_train <- label[training_indices,]
      x_valid <- x[validation_indices, ]
      y_valid <- y[validation_indices, ]
      
      # Perform 100 runs of clr and find the best one based on loss
      best_clr <- NULL
      best_loss <- Inf
      for (iter in 1:n_iterations) {
        clr_result <- kplane(x_train, y_train, k, gamma = gamma, m)
        if (clr_result$loss < best_loss) {
          best_clr <- clr_result
          best_loss <- clr_result$loss
        }
      }
      # Compute ARI and RMSE
      ari <- adjustedRandIndex(label_train, best_clr$cluster)
      rmse <- compute_rmse(x_train,y_train,best_clr$membership ,best_clr$coefficients)
      rmse_pre <- compute_rmse(x_valid,y_valid,best_clr$membership ,best_clr$coefficients)
      
      # Store the results
      all_ari <- c(all_ari,ari)
      all_rmse <- c(all_rmse,rmse)
      all_pre <- c(all_pre,rmse_pre)
    }
  }
  return(list(ARI = all_ari, RMSE = all_rmse, PRE = all_pre))
}

custom_cv_clr <- function (x, y, gamma) {
  all_ari <- c()
  all_rmse <- c()
  all_pre <- c()
  for (repeats in 1:n_repeats) {
    indices <- createFolds(1:length(label), k = n_folds, list = TRUE, returnTrain = FALSE)
    for (fold in 1:n_folds) {
      # Split data into training and validation sets
      validation_indices <- indices[[fold]]
      training_indices <- setdiff(1:length(label), validation_indices)
      
      x_train <- x[training_indices, ]
      y_train <- y[training_indices, ]
      label_train <- label[training_indices,]
      x_valid <- x[validation_indices, ]
      y_valid <- y[validation_indices, ]
      
      # Perform 100 runs of clr and find the best one based on loss
      best_clr <- NULL
      best_loss <- Inf
      for (iter in 1:n_iterations) {
        clr_result <- clr(x_train, y_train, k, m)
        if (clr_result$loss < best_loss) {
          best_clr <- clr_result
          best_loss <- clr_result$loss
        }
      }
      # Compute ARI and RMSE
      ari <- adjustedRandIndex(label_train, best_clr$cluster)
      rmse <- compute_rmse(x_train,y_train,best_clr$membership ,best_clr$coefficients)
      rmse_pre <- compute_rmse(x_valid,y_valid,best_clr$membership ,best_clr$coefficients)
      
      # Store the results
      all_ari <- c(all_ari,ari)
      all_rmse <- c(all_rmse,rmse)
      all_pre <- c(all_pre,rmse_pre)
    }
  }
  return(list(ARI = all_ari, RMSE = all_rmse, PRE = all_pre))
}

custom_cv_kmeans <- function (x, y, gamma) {
  all_ari <- c()
  all_rmse <- c()
  all_pre <- c()
  for (repeats in 1:n_repeats) {
    indices <- createFolds(1:length(label), k = n_folds, list = TRUE, returnTrain = FALSE)
    for (fold in 1:n_folds) {
      # Split data into training and validation sets
      validation_indices <- indices[[fold]]
      training_indices <- setdiff(1:length(label), validation_indices)
      
      x_train <- x[training_indices, ]
      y_train <- y[training_indices, ]
      label_train <- label[training_indices,]
      x_valid <- x[validation_indices, ]
      y_valid <- y[validation_indices, ]
      
      # Perform 100 runs of clr and find the best one based on loss
      best_clr <- NULL
      best_loss <- Inf
      for (iter in 1:n_iterations) {
        clr_result <- kmeans(x_train, y_train, k, m)
        if (clr_result$loss < best_loss) {
          best_clr <- clr_result
          best_loss <- clr_result$loss
        }
      }
      # Compute ARI and RMSE
      ari <- adjustedRandIndex(label_train, best_clr$cluster)
      rmse <- compute_rmse(x_train,y_train,best_clr$membership ,best_clr$coefficients)
      rmse_pre <- compute_rmse(x_valid,y_valid,best_clr$membership ,best_clr$coefficients)
      
      # Store the results
      all_ari <- c(all_ari,ari)
      all_rmse <- c(all_rmse,rmse)
      all_pre <- c(all_pre,rmse_pre)
    }
  }
  return(list(ARI = all_ari, RMSE = all_rmse, PRE = all_pre))
}
#####################################################
#clr
cv_results <- custom_cv_clr(x, y,gamma)

# Calculate mean and variance
mean_ari <- mean(cv_results$ARI)
mean_rmse <- mean(cv_results$RMSE)
mean_pre <- mean(cv_results$PRE)
var_ari <- var(cv_results$ARI)
var_rmse <- var(cv_results$RMSE)
var_pre <- var(cv_results$PRE)

# Display the results
print("clr")
cat("Mean ARI:", mean_ari, "\n")
cat("Variance ARI:", var_ari, "\n")
cat("Mean RMSE:", mean_rmse, "\n")
cat("Variance RMSE:", var_rmse, "\n")
cat("Mean PRE:", mean_pre, "\n")
cat("Variance PRE:", var_pre, "\n")

#####################################################
#kmeans
cv_results <- custom_cv_kmeans(x, y,gamma)

# Calculate mean and variance
mean_ari <- mean(cv_results$ARI)
mean_rmse <- mean(cv_results$RMSE)
mean_pre <- mean(cv_results$PRE)
var_ari <- var(cv_results$ARI)
var_rmse <- var(cv_results$RMSE)
var_pre <- var(cv_results$PRE)

# Display the results
print("kmeans")
cat("Mean ARI:", mean_ari, "\n")
cat("Variance ARI:", var_ari, "\n")
cat("Mean RMSE:", mean_rmse, "\n")
cat("Variance RMSE:", var_rmse, "\n")
cat("Mean PRE:", mean_pre, "\n")
cat("Variance PRE:", var_pre, "\n")

#####################################################
#kplane
# Run the cross-validation using the custom_cv function
for(gamma in c(0.00001,0.001,1,100,100000)){
  #kplane
  cv_results <- custom_cv_kplane(x, y,gamma)
  
  # Calculate mean and variance
  mean_ari <- mean(cv_results$ARI)
  mean_rmse <- mean(cv_results$RMSE)
  mean_pre <- mean(cv_results$PRE)
  var_ari <- var(cv_results$ARI)
  var_rmse <- var(cv_results$RMSE)
  var_pre <- var(cv_results$PRE)
  
  # Display the results
  #for kplane cat("gamma:", gamma, "\n")
  print("kplane")
  cat("gamma:", gamma, "\n")
  cat("Mean ARI:", mean_ari, "\n")
  cat("Variance ARI:", var_ari, "\n")
  cat("Mean RMSE:", mean_rmse, "\n")
  cat("Variance RMSE:", var_rmse, "\n")
  cat("Mean PRE:", mean_pre, "\n")
  cat("Variance PRE:", var_pre, "\n")
}

##############################################################
#wclr

clustering_functions <- list("epg","epl","qpg","qpl")

for(gamma in c(0.00001,0.001,1,100,100000)){
  for (clustering_function in clustering_functions){
    start_time <- Sys.time()
    cv_results <- custom_cv(x, y,gamma,clustering_function)
    # Stop the timer
    end_time <- Sys.time()
    # Calculate the elapsed time
    elapsed_time <- end_time - start_time
    # Print the elapsed time
    print(elapsed_time)
    
    # Calculate mean and variance
    mean_ari <- mean(cv_results$ARI)
    mean_rmse <- mean(cv_results$RMSE)
    var_ari <- var(cv_results$ARI)
    var_rmse <- var(cv_results$RMSE)
    
    # Display the results
    #for kplane cat("gamma:", gamma, "\n")
    print(clustering_function)
    cat("gamma:", gamma, "\n")
    cat("Mean ARI:", mean_ari, "\n")
    cat("Variance ARI:", var_ari, "\n")
    cat("Mean RMSE:", mean_rmse, "\n")
    cat("Variance RMSE:", var_rmse, "\n")
  }
}

