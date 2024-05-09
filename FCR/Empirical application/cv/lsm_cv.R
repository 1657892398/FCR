library(glmnet)
library(caret)
library(readxl)
library(ncvreg)
library(openxlsx)
library(readxl)
library(ggplot2)

setwd("C:/Users/dell/Desktop/FCR/Empirical application")
 # 读取Excel文件
revise  <- read_excel("C:/Users/dell/Desktop/FCR/4.Empirical application/data.xlsx", sheet = "text")
revise <- as.matrix(revise)
# 提取第一列作为因变量
y <- revise [, 1]
# 提取除第一列之外的所有列作为自变量
x <- revise [, -1]  
f=20  #flod
r=25 #repate

# # 将数据框转换为向量
load("CV.RData")#这里改成要读取的data文件
fcr_mse <- all_mse
# 去除NA值
fcr_mse <- na.omit(fcr_mse)
mean(fcr_mse)

load("kmeans_CV_MSE.RData")#这里改成要读取的data文件

# 存储所有 MSE 的向量
all_lasso_mse <- numeric(r*f)
all_scad_mse <- numeric(r*f)
all_mcp_mse <- numeric(r*f)

#lasso
# 重复20次交叉验证
for (j in 1:r) {
  # 执行交叉验证并计算 MSE
  lasso_mse <- numeric(f)  # 存储每次交叉验证的 MSE 值
  for (i in 1:f) {
    # 将数据划分为训练集和测试集
    # set.seed(3*j + i)  # 设置随机种子以确保结果的可重复性
    set.seed(as.integer(Sys.time()))  # 使用系统时间作为种子
    indices <- createDataPartition(y, p = 0.75, list = FALSE)
    x_train <- x[indices, ]
    y_train <- y[indices]
    x_test <- x[-indices, ]
    y_test <- y[-indices]
    
    # 使用 lasso 拟合线性模型
    lasso_fit <- ncvreg::cv.ncvreg(X = x_train , y = y_train, penalty = "lasso")
    
    # 在测试集上进行预测
    predictions1 <- predict(lasso_fit$fit, x_test, s = lasso_fit$lambda.min)
    
    # 计算 MSE
    lmse <- mean((y_test - predictions1)^2)
    lasso_mse[i] <- lmse
  }
  # 将每次交叉验证的 MSE 添加到所有 MSE 向量中
  start_index <- (j - 1) * f + 1
  end_index <- start_index + f - 1
  all_lasso_mse[start_index:end_index] <- lasso_mse
}


#scad
# 重复10次交叉验证
for (j in 1:r) {
  # 执行交叉验证并计算 MSE
  scad_mse <- numeric(f)  # 存储每次交叉验证的 MSE 值
  
  for (i in 1:f) {
    # 将数据划分为训练集和测试集
    # set.seed(3*j + i)  # 设置随机种子以确保结果的可重复性
    set.seed(as.integer(Sys.time()))  # 使用系统时间作为种子
    indices <- createDataPartition(y, p = 0.75, list = FALSE)
    # indices <- sample(1:nrow(revise), size = 0.8 * nrow(revise), replace = FALSE)
    x_train <- x[indices, ]
    y_train <- y[indices]
    x_test <- x[-indices, ]
    y_test <- y[-indices]
    
    # 使用 SCAD 拟合线性模型
    scad_fit <- ncvreg::cv.ncvreg(X = x_train , y = y_train, penalty = "SCAD")
    
    # 在测试集上进行预测
    predictions2 <- predict(scad_fit$fit, x_test, s = scad_fit$lambda.min)
    
    # 计算 MSE
    smse <- mean((y_test - predictions2)^2)
    scad_mse[i] <- smse
  }
  
  # 将每次交叉验证的 MSE 添加到所有 MSE 向量中
  start_index <- (j - 1) * f + 1
  end_index <- start_index + f - 1
  all_scad_mse[start_index:end_index] <- scad_mse
}


#mcp
# 重复10次交叉验证
for (j in 1:r) {
  # 执行交叉验证并计算 MSE
  mcp_mse <- numeric(f)  # 存储每次交叉验证的 MSE 值
  
  for (i in 1:f) {
    # 将数据划分为训练集和测试集
    # set.seed(3*j + i)  # 设置随机种子以确保结果的可重复性
    set.seed(as.integer(Sys.time()))  # 使用系统时间作为种子
    indices <- createDataPartition(y, p = 0.75, list = FALSE)
    x_train <- x[indices, ]
    y_train <- y[indices]
    x_test <- x[-indices, ]
    y_test <- y[-indices]
    
    # 使用 MCP 拟合线性模型
    mcp_fit <- ncvreg::cv.ncvreg(X = x_train , y = y_train, penalty = "MCP")
    
    # 在测试集上进行预测
    predictions3 <- predict(mcp_fit$fit, x_test, s = mcp_fit$lambda.min)
    
    # 计算 MSE
    mmse <- mean((y_test - predictions3)^2)
    mcp_mse[i] <- mmse
  }
  
  # 将每次交叉验证的 MSE 添加到所有 MSE 向量中
  start_index <- (j - 1) * f + 1
  end_index <- start_index + f - 1
  all_mcp_mse[start_index:end_index] <- mcp_mse
}

# 输出lasso_MSE值
all_lasso_mse
# 输出所有 scad_mse 值
all_scad_mse
# 输出所有 mcp_mse 值
all_mcp_mse
# 输出所有 fcr_mse 值
fcr_mse
# 输出所有 kmeans_mse 值
k_mse

mean(k_mse)
mean(fcr_mse)
mean(all_lasso_mse)
mean(all_scad_mse)
mean(all_mcp_mse)

#画图BOXplot
# 创建一个示例数据集
data <- list(
  fcr = fcr_mse,
  kmeans = k_mse,
  lasso = all_lasso_mse,
  scad = all_scad_mse,
  mcp = all_mcp_mse
  
)
# Combine the MSE values into a data frame
df <- data.frame(
  Method = rep(c("FCR", "K-means", "LASSO", "SCAD", "MCP"), 
               sapply(data, length)),
  MSE = unlist(data)
)

# 设置绘图设备为PNG文件
png("cv.png", width = 800, height = 600, res = 120)


# Create the boxplot
ggplot(df, aes(x = Method, y = MSE)) +
  geom_boxplot(aes(fill = Method), alpha = 0.7) +
  labs(x = "Methods", y = "MSE") +
  ggtitle("Boxplot of MSE calculated by various methods")



# 关闭绘图设备
dev.off()
