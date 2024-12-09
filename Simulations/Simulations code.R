#-----------------------------------------
#
#      Precision shrinkage portfolio
#
#-----------------------------------------

rm(list=ls())

setwd("H:/Il mio Drive/Precision shrinkage/Simulations")

library(rio)
library(RiskPortfolios)
library(nlshrink)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)

# Suppose to consider an i.i.d. economy where the true expected mean vector and covariance matrix are:

M <- 1000
SRp10 <- list()
SRp30 <- list()
SRp48 <- list()
Err10 <- list()
Err30 <- list()
Err48 <- list()
leng <- c(50, 150, 300, 600)

num_cores <- detectCores() - 1  # Use one less than total cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

for (j in 1:3){
  
  p <- j
  
  if (p==1){
    
    Dataset<-import("10industry-US.csv") 
    Dataset <- na.omit(apply(as.matrix(Dataset[,-1]),2,as.numeric))
    
  } else if (p==2){
    
    Dataset<-import("30industry-US.csv") 
    Dataset <- na.omit(apply(as.matrix(Dataset[,-1]),2,as.numeric))
    
  } else {
    
    Dataset<-import("48industry-US.csv") 
    Dataset <- na.omit(apply(as.matrix(Dataset[,-1]),2,as.numeric))
    
  }
  
  P <- ncol(Dataset)
  
  true_mu<-meanEstimation(Dataset)/100
  true_cov<-covEstimation(Dataset)/N
  gamma <- 1
  
# True weights and expected utility are:

w_true<-(1/gamma)*solve(true_cov)%*%true_mu
EX_true<-t(w_true)%*%true_mu-(gamma/2)*(t(w_true)%*%true_cov%*%w_true)
SR_true<-t(true_mu)%*%solve(true_cov)%*%true_mu

for (n in 1:length(leng)){

  N<-leng[n] # number of observations 
  
results <- foreach(i = 1:M, .packages = c("data.table","RiskPortfolios","nlshrink","MASS")) %dopar% {
  
  set.seed(i)
  Data<-mvrnorm(N, mu=true_mu, Sigma=true_cov)
  mu_ML<-meanEstimation(Data)
  D<-Data-mu_ML  # Scatter matrix: important to compute ML Cov, Sample Cov, Scaled Cov, etc.
  
  # Equally weighted portfolio:
  
  EX_ew<-t(rep(1/P, P))%*%true_mu-(gamma/2)*(t(rep(1/P, P))%*%true_cov%*%rep(1/P, P))
  SR_ew<-EX_true-EX_ew
  SR_ew<-(rep(1/P, P)%*%true_mu)/(t(rep(1/P, P))%*%(true_cov)%*%rep(1/P, P))
  
  # Equally weighted implied covariance:
  
  EW<-matrix(0, P, P)
  diag(EW)<-P*mu_ML
  w_EW<-(1/gamma)*solve(EW)%*%mu_ML
  EX_impEW<-t(w_EW)%*%true_mu-(gamma/2)*(t(w_EW)%*%true_cov%*%w_EW)
  SR_impEW<-EX_true-EX_impEW
  SR_impEW<-(t(w_EW)%*%true_mu)/t(w_EW)%*%(true_cov)%*%w_EW
  
  # Identity:
  
  I<-matrix(0, P, P)
  diag(I)<-1
  w_I<-(1/gamma)*I%*%mu_ML
  EX_I<-t(w_I)%*%true_mu-(gamma/2)*(t(w_I)%*%true_cov%*%w_I)
  SR_I<-EX_true-EX_I
  SR_I<-(t(w_I)%*%true_mu)/t(w_I)%*%(true_cov)%*%w_I
  
  # Maximum Likelihood Estimator:
  
  cov_ML<-(1/N)*(t(D)%*%D)
  w_ML<-(1/gamma)*solve(cov_ML)%*%mu_ML
  EX_ML<-t(w_ML)%*%true_mu-(gamma/2)*(t(w_ML)%*%true_cov%*%w_ML)
  SR_ML<-EX_true-EX_ML
  SR_ML<-(t(w_ML)%*%true_mu)/t(w_ML)%*%(true_cov)%*%w_ML
  
  # Sample Estimator:
  
  cov_SC<-(1/(N-1))*(t(D)%*%D)
  w_SC<-(1/gamma)*solve(cov_SC)%*%mu_ML
  EX_SC<-t(w_SC)%*%true_mu-(gamma/2)*(t(w_SC)%*%true_cov%*%w_SC)
  SR_SC<-EX_true-EX_SC
  SR_SC<-(t(w_SC)%*%true_mu)/t(w_SC)%*%(true_cov)%*%w_SC
  
  # Unbiased scaling Estimator:
  
  cov_SE<-(1/(N-P-2))*(t(D)%*%D)
  w_SE<-(1/gamma)*solve(cov_SE)%*%mu_ML
  EX_SE<-t(w_SE)%*%true_mu-(gamma/2)*(t(w_SE)%*%true_cov%*%w_SE)
  SR_SE<-EX_true-EX_SE
  SR_SE<-(t(w_SE)%*%true_mu)/t(w_SE)%*%(true_cov)%*%w_SE
  
  # Biased quasi-optimal scaling Estimator:
  
  c3<-(((N-P-1)*(N-P-4))/(N*(N-2)))
  cov_c3<-cov_ML/c3
  w_c3<-(1/gamma)*solve(cov_c3)%*%mu_ML
  EX_c3<-t(w_c3)%*%true_mu-(gamma/2)*(t(w_c3)%*%true_cov%*%w_c3)
  SR_c3<-EX_true-EX_c3
  SR_c3<-(t(w_c3)%*%true_mu)/t(w_c3)%*%(true_cov)%*%w_c3
  
  # Biased optimal scaling Estimator:
  
  SR_u<-as.numeric(((N-P-2)*(t(mu_ML)%*%solve(cov_ML)%*%mu_ML)+P)/N)
  c<-as.numeric((((N-P-1)*(N-P-4))/(N*(N-2)))*(SR_u/(SR_u+P/N)))
  cov_c<-cov_ML/c
  w_c<-(1/gamma)*solve(cov_c)%*%mu_ML
  EX_c<-t(w_c)%*%true_mu-(gamma/2)*(t(w_c)%*%true_cov%*%w_c)
  SR_c<-EX_true-EX_c
  SR_c<-(t(w_c)%*%true_mu)/t(w_c)%*%(true_cov)%*%w_c
  
  # Ledoit & Wolf - shrinkage tw identity (LW1):
  
  cov_LW1<-covEstimation(Data,control = list(type = 'diag'))
  w_LW1<-(1/gamma)*solve(cov_LW1)%*%mu_ML
  EX_LW1<-t(w_LW1)%*%true_mu-(gamma/2)*(t(w_LW1)%*%true_cov%*%w_LW1)
  SR_LW1<-EX_true-EX_LW1
  SR_LW1<-(t(w_LW1)%*%true_mu)/t(w_LW1)%*%(true_cov)%*%w_LW1
  
  # Ledoit & Wolf - shrinkage tw market (LW2):
  
  cov_LW2<-covEstimation(Data,control = list(type = 'lw'))
  w_LW2<-(1/gamma)*solve(cov_LW2)%*%mu_ML
  EX_LW2<-t(w_LW2)%*%true_mu-(gamma/2)*(t(w_LW2)%*%true_cov%*%w_LW2)
  SR_LW2<-EX_true-EX_LW2
  SR_LW2<-(t(w_LW2)%*%true_mu)/t(w_LW2)%*%(true_cov)%*%w_LW2
  
  # LW nls
  
  cov_LW3<-nlshrink::nlshrink_cov(Data)
  w_LW3<-(1/gamma)*solve(cov_LW3)%*%mu_ML
  EX_LW3<-t(w_LW3)%*%true_mu-(gamma/2)*(t(w_LW3)%*%true_cov%*%w_LW3)
  SR_LW3<-EX_true-EX_LW3
  SR_LW3<-(t(w_LW3)%*%true_mu)/t(w_LW3)%*%(true_cov)%*%w_LW3
  
  
  # Shrinkage Unbiased Scaling Estimator (SE) towards equally weighted (equivalent to Tu and Zhou, 2011):
  
  EW<-matrix(0, P, P)
  diag(EW)<-P*mu_ML
  c1<-((N-2)*(N-P-2))/((N-P-1)*(N-P-4))
  num<-(t(solve(EW)%*%mu_ML)%*%covEstimation(Data)%*%solve(EW)%*%mu_ML)-2*t(solve(EW)%*%mu_ML)%*%mu_ML+SR_u
  den<-c1*SR_u+c1*(P/N)+(t(solve(EW)%*%mu_ML)%*%covEstimation(Data)%*%solve(EW)%*%mu_ML)-2*t(solve(EW)%*%mu_ML)%*%mu_ML
  a<-as.numeric(num/den)
  
  Sigma_s<-a*solve(cov_SE)+(1-a)*solve(EW)
  w_sEW<-(1/gamma)*Sigma_s%*%mu_ML
  EX_sEW<-t(w_sEW)%*%true_mu-(gamma/2)*(t(w_sEW)%*%true_cov%*%w_sEW)
  SR_sEW<-EX_true-EX_sEW
  SR_sEW<-(t(w_sEW)%*%true_mu)/t(w_sEW)%*%(true_cov)%*%w_sEW
  
  # Shrinkage US towards I:
  
  c1<-((N-2)*(N-P-2))/((N-P-1)*(N-P-4))
  I<-matrix(0, P, P)
  diag(I)<-1
  lambda<-t(mu_ML)%*%mu_ML
  Qhat<-(N/(N-1)) * (t(mu_ML)%*%cov_ML%*%mu_ML)
  num<-Qhat-2*lambda+SR_u
  den<-c1*(SR_u+P/N)+Qhat-2*lambda
  a<-as.numeric(num/den)
  Sigma_sI<-a*solve(cov_SE)+(1-a)*I
  w_sI<-(1/gamma)*Sigma_sI%*%mu_ML
  EX_sI<-t(w_sI)%*%true_mu-(gamma/2)*(t(w_sI)%*%true_cov%*%w_sI)
  SR_sI<-EX_true-EX_sI
  SR_sI<-(t(w_sI)%*%true_mu)/t(w_sI)%*%(true_cov)%*%w_sI

  # Shrinkage US towards biased quasi optimal scaling c3:
  
  c1<-((N-2)*(N-P-2))/((N-P-1)*(N-P-4))
  num<-((1-c1)/c1)*(P/N)
  den<-((1-c1)^2/c1)*SR_u+((1-c1)^2/c1)*(P/N)
  a<-as.numeric(num/den)
  Sigma_sc3<-a*solve(cov_SE)+(1-a)*solve(cov_c3)
  w_sc3<-(1/gamma)*Sigma_sc3%*%mu_ML
  EX_sc3<-t(w_sc3)%*%true_mu-(gamma/2)*(t(w_sc3)%*%true_cov%*%w_sc3)
  SR_sc3<-EX_true-EX_sc3
  SR_sc3<-(t(w_sc3)%*%true_mu)/t(w_sc3)%*%(true_cov)%*%w_sc3
  
  # Shrinkage c3 towards equally weighted:
  
  EW<-matrix(0, P, P)
  diag(EW)<-P*mu_ML
  c1<-((N-2)*(N-P-2))/((N-P-1)*(N-P-4))
  num<-(t(solve(EW)%*%mu_ML)%*%covEstimation(Data)%*%solve(EW)%*%mu_ML)-(1+(1/c1))*t(solve(EW)%*%mu_ML)%*%mu_ML+(1/c1)*SR_u
  den<-(1/c1)*SR_u+(1/c1)*(P/N)+(t(solve(EW)%*%mu_ML)%*%covEstimation(Data)%*%solve(EW)%*%mu_ML)-2*t(solve(EW)%*%mu_ML)%*%mu_ML
  a<-as.numeric(num/den)
  Sigma_s<-a*solve(cov_c3)+(1-a)*solve(EW)
  w_sc3EW<-(1/gamma)*Sigma_s%*%mu_ML
  EX_sc3EW<-t(w_sc3EW)%*%true_mu-(gamma/2)*(t(w_sc3EW)%*%true_cov%*%w_sc3EW)
  SR_sc3EW<-EX_true-EX_sc3EW
  SR_sc3EW<-(t(w_sc3EW)%*%true_mu)/t(w_sc3EW)%*%(true_cov)%*%w_sc3EW
  
  # Shrinkage c3 towards I:
  
  c1<-((N-2)*(N-P-2))/((N-P-1)*(N-P-4))
  I<-matrix(0, P, P)
  diag(I)<-1
  lambda<-t(mu_ML)%*%mu_ML
  Qhat<-(N/(N-1)) * (t(mu_ML)%*%cov_ML%*%mu_ML)
  num<-Qhat-(1+(1/c1))*lambda+(1/c1)*SR_u
  den<-(1/c1)*(SR_u+P/N-2*lambda)+Qhat
  a<-as.numeric(num/den)
  Sigma_sI<-a*solve(cov_c3)+(1-a)*I
  w_sc3I<-(1/gamma)*Sigma_sI%*%mu_ML
  EX_sc3I<-t(w_sc3I)%*%true_mu-(gamma/2)*(t(w_sc3I)%*%true_cov%*%w_sc3I)
  SR_sc3I<-EX_true-EX_sc3I
  SR_sc3I<-(t(w_sc3I)%*%true_mu)/t(w_sc3I)%*%(true_cov)%*%w_sc3I
  
  # Shrinkage c* towards equally weighted:
  
  EW<-matrix(0, P, P)
  diag(EW)<-P*mu_ML
  c1<-((N-2)*(N-P-2))/((N-P-1)*(N-P-4))
  num<-(1/c1)*(SR_u^2/(SR_u+(P/N)))+(t(solve(EW)%*%mu_ML)%*%covEstimation(Data)%*%solve(EW)%*%mu_ML)-(1/c1)*(SR_u/(SR_u+(P/N)))*t(solve(EW)%*%mu_ML)%*%mu_ML-t(solve(EW)%*%mu_ML)%*%mu_ML
  den<-(1/c1)*(P/N+SR_u)+(t(solve(EW)%*%mu_ML)%*%covEstimation(Data)%*%solve(EW)%*%mu_ML)-2*(1/c1)*(SR_u/(SR_u+(P/N)))*t(solve(EW)%*%mu_ML)%*%mu_ML
  a<-as.numeric(num/den)
  Sigma_s<-a*solve(cov_c)+(1-a)*solve(EW)
  w_scEW<-(1/gamma)*Sigma_s%*%mu_ML
  EX_scEW<-t(w_scEW)%*%true_mu-(gamma/2)*(t(w_scEW)%*%true_cov%*%w_scEW)
  SR_scEW<-EX_true-EX_scEW
  SR_scEW<-(t(w_scEW)%*%true_mu)/t(w_scEW)%*%(true_cov)%*%w_scEW
  
  # Shrinkage c* towards I:
  
  c1<-((N-2)*(N-P-2))/((N-P-1)*(N-P-4))
  I<-matrix(0, P, P)
  diag(I)<-1
  lambda<-t(mu_ML)%*%mu_ML
  Qhat<-(N/(N-1)) * (t(mu_ML)%*%cov_ML%*%mu_ML)
  num<-(1/c1)*(SR_u^2/(SR_u+(P/N)))-lambda+Qhat-(1/c1)*(SR_u/(SR_u+(P/N)))*lambda
  den<-(1/c1)*(P/N+SR_u)+Qhat-2*(1/c1)*(SR_u/(SR_u+(P/N)))*lambda
  a<-as.numeric(num/den)
  Sigma_sI<-a*solve(cov_c)+(1-a)*I
  w_scI<-(1/gamma)*Sigma_sI%*%mu_ML
  EX_scI<-t(w_scI)%*%true_mu-(gamma/2)*(t(w_scI)%*%true_cov%*%w_scI)
  SR_scI<-EX_true-EX_scI
  SR_scI<-(t(w_scI)%*%true_mu)/t(w_scI)%*%(true_cov)%*%w_scI

  list(cbind("ML"=(EX_ML), "SC"=(EX_SC), "SE"=(EX_SE), "c3"=(EX_c3), 
             "c"=(EX_c), "I"=(EX_I), "Implied EW"=(EX_impEW), 
             "EW"=(EX_ew),"SE tw EW"=(EX_sEW), "SE tw I"=(EX_sI), 
             "SE tw c3"=(EX_sc3), "c3 tw EW"=(EX_sc3EW), "c3 tw I"=(EX_sc3I), 
             "c* tw EW"=(EX_scEW), "c* tw I"=(EX_scI), "LW1"=(EX_LW1), "LW2"=(EX_LW2), "LW3"=(EX_LW3)),
       cbind("ML"=(SR_ML), "SC"=(SR_SC), "SE"=(SR_SE), "c3"=(SR_c3), 
        "c"=(SR_c), "I"=(SR_I), "Implied EW"=(SR_impEW), 
        "EW"=(SR_ew),"SE tw EW"=(SR_sEW), "SE tw I"=(SR_sI), 
        "SE tw c3"=(SR_sc3), "c3 tw EW"=(SR_sc3EW), "c3 tw I"=(SR_sc3I), 
        "c* tw EW"=(SR_scEW), "c* tw I"=(SR_scI), "LW1"=(SR_LW1), "LW2"=(SR_LW2), "LW3"=(SR_LW3))
  )
  
}

if (p==1){
  
  SRmat <- matrix(unlist(results[[2]]),nrow=18, ncol=M)
  Errmat <- matrix(unlist(results[[1]]),nrow=18, ncol=M)
  SRp10[[n]]<- rowMeans(SRmat)
  Err10[[n]] <- rowMeans(Errmat)
  
} else if (p==2){
  
  SRmat <- matrix(unlist(results[[2]]),nrow=18, ncol=M)
  Errmat <- matrix(unlist(results[[1]]),nrow=18, ncol=M)
  SRp30[[n]]<- rowMeans(SRmat)
  Err30[[n]] <- rowMeans(Errmat)
  
} else {
  
  SRmat <- matrix(unlist(results[[2]]),nrow=18, ncol=M)
  Errmat <- matrix(unlist(results[[1]]),nrow=18, ncol=M)
  SRp48[[n]]<- rowMeans(SRmat)
  Err48[[n]] <- rowMeans(Errmat)
  
}

}

print(paste0("Simulation with P=",j," finished."))

}
stopCluster(cl)

#-------- Tables

PanelA <- matrix(unlist(lapply(SRp5, colMeans)), nrow=18, ncol=5)
rownames(PanelA) <- c("ML","SC","UP","c3","c*","I","Naive","EW","UP+EW",
                      "UP+I","UP+c3","c3+EW","c3+I","c*+EW","c*+I",
                      "LW1", "LW1", "NLS")

PanelB <- matrix(unlist(lapply(SRp30, function(x){lapply(x, mean)})), nrow=18, ncol=5)
rownames(PanelB) <- c("ML","SC","UP","c3","c*","I","Naive","EW","UP+EW",
                      "UP+I","UP+c3","c3+EW","c3+I","c*+EW","c*+I",
                      "LW1", "LW1", "NLS")

xtable::xtable(PanelA, digits=4)
xtable::xtable(PanelB, digits=4)

#------------------ Estimators Theta, lambda and Q

rm(list=ls())

nn <- c(50, 150, 300, 600)
M <- 1000

BIAS1 <- matrix(NA, nrow=3, ncol=length(nn))
BIAS2 <- matrix(NA, nrow=3, ncol=length(nn))
BIAS3 <- matrix(NA, nrow=3, ncol=length(nn))

for (j in 1:3){

p <- j

if (p==1){
  
  Dataset<-import("10industry-US.csv") 
  Dataset <- na.omit(apply(as.matrix(Dataset[,-1]),2,as.numeric))
  
} else if (p==2){
  
  Dataset<-import("30industry-US.csv") 
  Dataset <- na.omit(apply(as.matrix(Dataset[,-1]),2,as.numeric))
  
} else {

  Dataset<-import("48industry-US.csv") 
  Dataset <- na.omit(apply(as.matrix(Dataset[,-1]),2,as.numeric))
  
}

P <- ncol(Dataset)

for (n in 1:length(nn)){

N <- nn[n]
Qest <- numeric(M) 
Lest <- numeric(M) 
SRest <- numeric(M) 
true_mu<-meanEstimation(Dataset)/100
true_cov<-covEstimation(Dataset)/N

for (i in 1:M){

set.seed(i)
Data<-mvrnorm(N, mu=true_mu, Sigma=true_cov)
mu_ML<-colMeans(Data)
cov_SC<-cov(Data)
D<-Data-mu_ML  # Scatter matrix: important to compute ML Cov, Sample Cov, Scaled Cov, etc.
cov_ML<-(1/N)*(t(D)%*%D)

Qest[i]<- ((N)/(N-1))*( (t(mu_ML)%*%cov_ML%*%mu_ML) )# ((N)/(N-1)) ,((N-1)/(N-p-2))*
Lest[i]<- (t(mu_ML)%*%mu_ML)
SRest[i]<- as.numeric(((N-P-2)*(t(mu_ML)%*%solve(cov_ML)%*%mu_ML)+P)/N)

print((i/M)*100)
}

BIAS1[p,n] <- mean(Qest)- (sum(diag(true_cov %*% true_cov)) / N + t(true_mu) %*% true_cov %*% true_mu)

BIAS2[p,n] <-mean(Lest)- ( sum(diag(true_cov)) / N + t(true_mu) %*% true_mu )

BIAS3[p,n] <- mean(SRest)- ( t(true_mu) %*% solve(true_cov) %*% true_mu )

}

print(paste0("Simulation with P=",j," finished."))

}
B1 <- abs(BIAS1)
B2 <- abs(BIAS2)
B3 <- abs(BIAS3)

matplot(nn, t(B1), type = "b", pch = 1:3, lty = 1, col = 1:4,
         xlab = "Sample size (T)", ylab = "Bias", main = "Bias for the Q estimator")
legend("topright", legend = paste("N =", c(10, 30, 48)), col = 1:3, pch = 1:4, lty = 1, cex=0.8)

matplot(nn, t(B2), type = "b", pch = 1:3, lty = 1, col = 1:3,
        xlab = "Sample size (T)", ylab = "Bias", main = "Bias for the Lambda estimator")
legend("topright", legend = paste("N =", c(10, 30, 48)), col = 1:3, pch = 1:3, lty = 1, cex=0.8)

matplot(nn, t(B3), type = "b", pch = 1:3, lty = 1, col = 1:3,
        xlab = "Sample size (T)", ylab = "Bias", main = "Bias for the Theta estimator")
legend("topright", legend = paste("N =", c(10, 30, 48)), col = 1:3, pch = 1:3, lty = 1, cex=0.8)

# Reset plot layout
par(mfrow = c(1, 1))
