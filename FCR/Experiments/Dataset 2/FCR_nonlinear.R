library(R.matlab)
library(glmnet)
library(mclust)
library(caret)


# 主函数
setwd("C:/Users/dell/Desktop/写论文/r语言代码非线性版本")
da <- readMat("data1.mat")
x <- da$x.data
y <- da$y.data
label <- da$label.data
x_non <- cbind(x,x^2,x^3)
#x_non <- x
n_folds <- 5
n_repeats <- 3
n_iterations <- 10
C <- 3
m <- 1.86
lamda <- 1 * rep(1, C)
gamma <- c(0.62, 1)
eps = 0.001

#函数定义部分
iterative_mu <- function(w, x, m) {
  # w:权重矩阵 N*C
  # x:数据集 N*p
  # m:隶属度
  p <- ncol(x)
  C <- ncol(w)
  mu <- matrix(0, nrow = p, ncol = C)
  
  for (j in 1:C) {
    mu[, j] <- t(x) %*% (w[, j]^m) / sum(w[, j]^m)
  }
  
  return(mu)
}

iterative_b <- function(x, y, w, b, m, lamda) {
  # x: 自变量 N*p
  # y: 因变量 N*1
  # w: 成员矩阵 N*C
  # b: 模型参数变量 p*C
  # m: 成员指数
  # lambda: Lasso 参数
  
  w <- w^m
  b1 <- b
  b2 <- matrix(100, nrow = nrow(b), ncol = ncol(b))
  k <- 0
  n <- ncol(x)  # n = p+1
  C <- ncol(b)
  
  while (sum(abs(b1 - b2)) > 0.001) {
    #print(sum(abs(b1 - b2)))
    b2 <- b1
    
    for (p in 1:C) {
      for (t in 1:n) {
        b1[t, p] <- 0
        fenmu <- sum(w[, p] * x[, t]^2)
        fenzi1 <- sum(w[, p] * y * x[, t])
        xx <- x %*% b1
        fenzi2 <- sum(xx[, p] * w[, p] * x[, t])
        fenzi <- fenzi1 - fenzi2
        
        if (fenzi - lamda[p]/2 > 0) {
          b1[t, p] <- (fenzi - lamda[p]/2) / fenmu
        } else if (fenzi + lamda[p]/2 < 0) {
          b1[t, p] <- (fenzi + lamda[p]/2) / fenmu
        } else {
          b1[t, p] <- 0
        }
      }
    }
    
    k <- k + 1
    if (k >= 1000) {
      break
    }
  }
  return(b1)
}

iterative_w <- function(x_tilde, x, y, b, mu, m, ga) {
  # x_tilde:用于回归的变量 N*p+1
  # x:自变量 N*p
  # y:因变量 N*1
  # b:模型的参数变量 N*C
  # mu:簇心矩阵 p*C
  # m:隶属度的次数
  # ga:gamma参数
  N <- nrow(x)
  C <- ncol(b)
  w <- matrix(0, nrow = N, ncol = C)
  
  for (i in 1:N) {
    ww <- numeric(C)
    
    for (k in 1:C) {
      ww[k] <- ga[1] * (y[i] - x_tilde[i,] %*% b[, k])^2 + ga[2] * sum((x[i,] - mu[, k])^2)
    }
    
    ww <- 1 / ww^(1/(m-1))
    ww <- sum(ww)
    for (j in 1:C) {
      w[i, j] <- 1 / (ga[1] * (y[i] - x_tilde[i,] %*% b[, j])^2 + ga[2] * sum((x[i,] - mu[, j])^2))^(1/(m-1)) / ww
    }
  }
  return(w)
}

compute_cost <- function(x, x_t, y, w, b, mu, m, gamma, lambda) {
  # 计算Lasso模型的损失函数值
  # 输入：
  # x: 自变量
  # y: 因变量
  # w: 隶属度
  # b: 模型的参数变量
  # m: 隶属度的次数
  # lamda: 惩罚参数
  # 输出：
  # loss: 损失函数值
  N <- nrow(x)
  C <- ncol(b)
  cost <- 0
  
  for (j in 1:C) {
    for (i in 1:N) {
      L <- gamma[1] * (y[i] - x_t[i,] %*% b[, j])^2 + gamma[2] * sum((x[i,] - mu[, j])^2)
      cost <- cost + w[i, j]^m * L
    }
    
    ma <- max(abs(b[, j]))
    cost <- cost + lambda[j] * ma
  }
  
  cost <- cost / N
  return(cost)
}
w_pre <- function(x,mu) {#预测w
  M <- nrow(x)
  C <- ncol(mu)
  w <- matrix(0, nrow = M, ncol = C)
  
  # 计算新点到每个簇心的距离
  for (k in 1:M) {
    x_dis <- numeric(C)
    for (j in 1:C) {
      x_dis[j] <- sqrt(sum((x[k,] - mu[,j])^2))
    }
    mindis <- which.min(x_dis)
    w[k, mindis] <- 1
  }
  return(w)
}

compute_mse <- function(x, y, w, b) {
  # 计算Lasso模型的损失函数值
  # 输入：
  # x: 自变量
  # y: 因变量
  # w: 隶属度
  # b: 模型的参数变量
  # m: 隶属度的次数
  # lamda: 惩罚参数
  # 输出：
  # loss: 损失函数值
  N <- nrow(x)
  C <- ncol(b)
  mse <- 0
  
  for (i in 1:N) {
    yy <- 0
    for (j in 1:C) {
      yy <- yy + w[i, j] * x[i,] %*% b[, j]
    }
    mse <- mse + (y[i] - yy)^2
  }
  
  mse <- mse / N
  return(mse)
}

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

# Define the function
FCR <- function(x_data, x_non, y_data, C, m, lamda, gamma=c(1,1), nstart = 1, eps = 0.001) {
  # Prepare data
  N <- nrow(x_data)
  p <- ncol(x_non)
  lamda <- lamda / gamma[1]
  x_tilde <- cbind(1, x_non)
  
  min_cost <- Inf
  best_w <- NULL
  best_b <- NULL
  best_mu <- NULL
  best_label <- NULL
  best_ari <- NULL
  
  for (start in 1:nstart) {
    cat(paste0(round(start/nstart*100), '% completed'))
    if(start==nstart)cat('\nDONE!\n')
    else cat('\r')
    # Initialize values
    b1 <- matrix(rnorm((p + 1) * C), nrow = p + 1, ncol = C)
    w <- matrix(0, nrow = N, ncol = C)
    mu1 <- matrix(0, nrow = p, ncol = C)
    
    b2 <- matrix(100000, nrow = p + 1, ncol = C)
    mu2 <- matrix(100000, nrow = p, ncol = C)
    
    times <- 0
    
    while (sum(abs(b1 - b2)) > eps || sum(abs(mu1 - mu2)) > eps) {
      b2 <- b1
      mu2 <- mu1
      w <- iterative_w(x_tilde, x_data, y_data, b1, mu1, m, gamma)
      mu1 <- iterative_mu(w, x_data, m)
      b1 = iterative_b(x_tilde, y_data, w, b1, m, lamda)
      
      times <- times + 1
      if (times > 10000) {
        stop("循环次数过多")
      }
    }
    #print(b1)
    #print(mu1)
    w <- iterative_w(x_tilde, x_data, y_data, b1, mu1, m, gamma)
    
    # Compute cost
    cost <- compute_cost(x_data, x_tilde, y_data, w, b1, mu1, m, gamma, lamda)
    
    if (cost < min_cost) {
      # Update the best results
      #label <- apply(w, 1, function(row) which(row == max(row))[1])
      #best_ari <- adjustedRandIndex(label, label_real)
      best_label <- apply(w, 1, function(row) which(row == max(row))[1])
      best_w <- w
      best_b <- b1
      best_mu <- mu1
      min_cost <- cost
    }
  }
  
  return(list(membership = best_w, coefficients = best_b, mu = best_mu, cluster = best_label, loss = min_cost))
}

custom_cv_fcr <- function (x,x_non,y, C, m, lamda, gamma=c(1,1)) {
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
      x_non_train <- x_non[training_indices, ]
      label_train <- label[training_indices,]
      x_valid <- x[validation_indices, ]
      y_valid <- y[validation_indices, ]
      x_non_valid <- x_non[validation_indices, ]
      
      # Perform 100 runs of clr and find the best one based on loss
      best_clr <- NULL
      best_loss <- Inf
      for (iter in 1:n_iterations) {
        clr_result <- FCR(x_train,x_non_train, y_train, C, m, lamda, gamma=gamma)
        if (clr_result$loss < best_loss) {
          best_clr <- clr_result
          best_loss <- clr_result$loss
        }
      }
      # Compute ARI and RMSE
      ari <- adjustedRandIndex(label_train, best_clr$cluster)
      rmse <- compute_rmse(x_non_train,y_train,best_clr$membership ,best_clr$coefficients)
      rmse_pre <- compute_rmse(x_non_valid,y_valid,w_pre(x_non_valid,best_clr$mu) ,best_clr$coefficients)
      print(rmse_pre)
      
      # Store the results
      all_ari <- c(all_ari,ari)
      all_rmse <- c(all_rmse,rmse)
      all_pre <- c(all_pre,rmse_pre)
    }
  }
  w = iterative_w(cbind(1, x_non), x, y, best_clr$coefficients, best_clr$mu, m, gamma)
  return(list(ARI = all_ari, RMSE = all_rmse, PRE = all_pre, best_w = w))
}

cv_results <- custom_cv_fcr(x,x_non, y, C, m, lamda, gamma=gamma)

# Calculate mean and variance
mean_ari <- mean(cv_results$ARI)
mean_rmse <- mean(cv_results$RMSE)
mean_pre <- mean(cv_results$PRE)
var_ari <- var(cv_results$ARI)
var_rmse <- var(cv_results$RMSE)
var_pre <- var(cv_results$PRE)

# Display the results
print("FCR")
cat("Mean ARI:", mean_ari, "\n")
cat("Variance ARI:", var_ari, "\n")
cat("Mean RMSE:", mean_rmse, "\n")
cat("Variance RMSE:", var_rmse, "\n")
cat("Mean PRE:", mean_pre, "\n")
cat("Variance PRE:", var_pre, "\n")

