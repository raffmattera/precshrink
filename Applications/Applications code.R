############### Empirical comparisons chapter 1 ################################

rm(list=ls())

library(rio)
library(RiskPortfolios)
library(nlshrink)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)

#----------- Prepare the data

folder_path <- getwd()
datalist <- list()
zip_files <- list.files(path = folder_path, pattern = "\\.zip$", full.names = TRUE)
for (zip_file in zip_files) {
  # Extract the file names from the .zip (assume single CSV per .zip)
  csv_name <- unzip(zip_file, list = TRUE)$Name[1]
  
  # Extract the CSV file to a temporary directory
  temp_dir <- tempdir()
  unzip(zip_file, files = csv_name, exdir = temp_dir)
  
  # Read the CSV file
  csv_path <- file.path(temp_dir, csv_name)
  data <- as.matrix(fread(csv_path,skip=20))
  
  # Append the data to the list
  datalist[[zip_file]] <- data
}
names(datalist) <- paste0("Data",1:length(datalist))

#----------- Parallelize for all datasets

num_cores <- detectCores() - 1  # Use one less than total cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- list()

# Change M below to get 120, 240 and 360.

for(p in 1:length(datalist)) {
  
  Dataset <- apply(as.matrix(datalist[[p]][-1,-1]),2,as.numeric)
  Dataset <-Dataset[!apply(Dataset == -99.99, 1, any), ]
  n <- ncol(Dataset) # no. Assets
  t <- nrow(Dataset) # no. Observations
  
  M <- 360 # Sample size
  g <- 1 # Risk aversion coefficient

# Equally weighted:

ew_r<-vector(mode="numeric", length = t-M) # ew portfolio returns

for (i in 1:(t-M)) {
  ew_r[i]<-mean(Dataset[M+i,])
}

# Identity matrix:

I<-matrix(0, n, n)
diag(I)<-1

w_I <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "I"))
}

I_r<-vector(mode="numeric", length = t-M) # portfolio returns

for (i in 1:(t-M)) {
  I_r[i]<-weighted.mean(Dataset[M+i,], t(w_I)[i,])
}

# Sample covariance estimator 

w_sc <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "sc"))
}

sc_r<-vector(mode="numeric", length = t-M) 

for (i in 1:(t-M)) {
  sc_r[i]<-weighted.mean(Dataset[M+i,], t(w_sc)[i,])
}

# Unbiased scaling covariance estimator 

w_us <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "us"))
}

us_r<-vector(mode="numeric", length = t-M) 

for (i in 1:(t-M)) {
  us_r[i]<-weighted.mean(Dataset[M+i,], t(w_us)[i,])
}

# Unbiased quasi optimal scaling "c3" covariance estimator 

w_c3 <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "c3"))
}

c3_r<-vector(mode="numeric", length = t-M) 

for (i in 1:(t-M)) {
  c3_r[i]<-weighted.mean(Dataset[M+i,], t(w_c3)[i,])
}

# Optimal scaling c* estimator

w_c <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "c"))
}

c_r<-vector(mode="numeric", length = t-M) 

for (i in 1:(t-M)) {
  c_r[i]<-weighted.mean(Dataset[M+i,], t(w_c)[i,])
}

# Ledoit & Wolf (2004) linear shrinkage toward Identity

w_lw2 <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "lw_I"))
}

lw2_r<-vector(mode="numeric", length = t-M)

for (i in 1:(t-M)) {
  lw2_r[i]<-weighted.mean(Dataset[M+i,], t(w_lw2)[i,])
}

# Ledoit & Wolf (2003) linear shrinkage estimator

w_lw <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "lw_mkt"))
}

lw_r<-vector(mode="numeric", length = t-M)

for (i in 1:(t-M)) {
  lw_r[i]<-weighted.mean(Dataset[M+i,], t(w_lw)[i,])
}

# NLS shrinkage Ledoit & Wolf (2017)

w_NLS <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "nls"))
}

NLS_r<-vector(mode="numeric", length = t-M)

for (i in 1:(t-M)) {
  NLS_r[i]<-weighted.mean(Dataset[M+i,], t(w_NLS)[i,])
}

# Our shrinkage estimator US towards Identity

w_usI <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "mp_usI"))
}

usI_r<-vector(mode="numeric", length = t-M) 

for (i in 1:(t-M)) {
  usI_r[i]<-weighted.mean(Dataset[M+i,], t(w_usI)[i,])
}

# Our shrinkage estimator US towards Equally Weighted

w_usEW <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "mp_usEW"))
}

usEW_r<-vector(mode="numeric", length = t-M)

for (i in 1:(t-M)) {
  usEW_r[i]<-weighted.mean(Dataset[M+i,], t(w_usEW)[i,])
}

# Our shrinkage estimator c3 towards Identity

w_c3I <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "mp_c3I"))
}

c3I_r<-vector(mode="numeric", length = t-M) 

for (i in 1:(t-M)) {
  c3I_r[i]<-weighted.mean(Dataset[M+i,], t(w_usI)[i,])
}

# Our shrinkage estimator c3 towards Equally Weighted

w_c3EW <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "mp_c3EW"))
}

c3EW_r<-vector(mode="numeric", length = t-M)

for (i in 1:(t-M)) {
  c3EW_r[i]<-weighted.mean(Dataset[M+i,], t(w_c3EW)[i,])
}

# Our shrinkage estimator c* towards Identity

w_cI <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "mp_cI"))
}

cI_r<-vector(mode="numeric", length = t-M) 

for (i in 1:(t-M)) {
  cI_r[i]<-weighted.mean(Dataset[M+i,], t(w_cI)[i,])
}

# Our shrinkage estimator c* towards Equally Weighted

w_cEW <- foreach(i = 1:(t - M), .combine = cbind, .packages = c("RiskPortfolios")) %dopar% {
  mvp(Dataset[i:(M + i - 1), ], g = g, control = list(type = "mp_cEW"))
}

cEW_r<-vector(mode="numeric", length = t-M)

for (i in 1:(t-M)) {
  cEW_r[i]<-weighted.mean(Dataset[M+i,], t(w_cEW)[i,])
}

results[[p]] <- (cbind(ew_r,I_r,sc_r,us_r,c3_r,c_r,lw_r,lw2_r, NLS_r,usI_r,usEW_r,c3I_r,c3EW_r,cI_r,cEW_r))

print(paste0("Dataset ",p," completed"))
}
stopCluster(cl)

# Evaluation

res <- matrix(NA, nrow=15, ncol=12)
for(i in 1:12){
  res[,i] <- sharpe(results[[i]]) 
}
rownames(res) <-  c("EW","I","SC","PM","c3","c*","LW1",
                    "LW2","LW3","PM+I","PM+EW","c3+I","c3+EW",
                    "c*+I", "c*+EW")
xtable::xtable(round(res,6)*100)
apply(res, 2, which.max)

# p-values

respval <- matrix(NA, nrow=15, ncol=12)
for(i in 1:12){
  for (j in 1:14) {
    respval[j,i] <- sharpeTesting(results[[i]][,j], results[[i]][,15])$pval 
  }
}
rownames(respval) <-  c("EW","I","SC","PM","c3","c*","LW1",
                    "LW2","LW3","PM+I","PM+EW","c3+I","c3+EW",
                    "c*+I", "c*+EW")
xtable::xtable(respval)
