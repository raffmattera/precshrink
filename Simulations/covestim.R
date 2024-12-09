cov_us <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  covus<-((t-1)/(t-n-2))*cov(rets)

  return(covus)
}

cov_mkt<-function(rets) {

  if (missing(rets)) d
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  D<-rets-m
  mkt<-rowMeans(rets)
  sample<-cov(cbind(rets, mkt))
  covmkt<-sample[1:n, n + 1]
  varmkt<-as.numeric(sample[n + 1, n + 1])
  sample<-sample[-(n + 1), -(n + 1)]

  target<-outer(covmkt, covmkt)/varmkt
  diag(target) <- diag(sample)

  return(target)
}

cov_factor<-function(rets) {

  std <- apply(rets, 2, sd)
  sigma <- cov_us(rets)

  K<-ifelse(berryFunctions::is.error(POETKhat(rets)), 1, POETKhat(rets)$K1HL)

  loading <- factanal(rets, K)$loadings
  uniquenesses <- factanal(rets, K)$uniquenesses

  R <- tcrossprod(loading) + diag(uniquenesses)
  diagStd <- diag(std)
  Sigma <- diagStd %*% R %*% diagStd

  return(Sigma)
}

cov_lwI<-function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  D<-rets-m
  sample<-cov(rets)

  I<-matrix(0, n, n)
  diag(I) <- 1

  y<-D^2

  # Phi hat
  phiMat <-  crossprod(y)/t - 2 * crossprod(D) * sample/t + (sample^2)/t
  phi <- sum(apply(phiMat, 2, sum))

  # Gamma hat
  gamma <- norm(sample - I, "F")^2

  # Kappa hat
  kappa <- phi/gamma

  # shrinkage value
  shrinkage <- pmax(0, pmin(1, kappa/t))

  # Sigma hat
  Sigma <- shrinkage * I + (1 - shrinkage) * sample

  return(list("cov"=Sigma, "shrinkage"=shrinkage))
}

mpcov_usI <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  I<-matrix(0, n, n)
  diag(I)<-1
  SR_u<-as.numeric(((t-n-2)*(t(m)%*%solve(cov(rets))%*%m)+n)/t)
  lambda<-t(m)%*%m
  Qhat<-(t(m)%*%cov(rets)%*%m)
  c1<-((n-2)*(t-n-2))/((t-n-1)*(t-n-4))
  num<-Qhat-2*lambda+SR_u
  den<-c1*(SR_u+(n/t))+Qhat-2*lambda
  a<-as.numeric(num/den)
  covus<-((t-1)/(t-n-2))*cov(rets)

  shrinkage<-pmax(0, pmin(1, a))

  Sigma_usI<-shrinkage*solve(covus)+(1-shrinkage)*I

  return(list("cov"=Sigma_usI, "shrinkage"=shrinkage))
}

mpcov_usEW <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  covus<-((t-1)/(t-n-2))*cov(rets)
  EW<-matrix(0, n, n)
  diag(EW)<-n*m
  wew<-rep(1/n, n)
  c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
  SR_u<-as.numeric(((t-n-2)*(t(m)%*%solve(cov(rets))%*%m)+n)/t)
  num<-(t(wew)%*%cov(rets)%*%wew)-(1+(1/c1))*t(wew)%*%m+(1/c1)*SR_u
  den<-(1/c1)*(SR_u+(n/t))+(t(wew)%*%cov(rets)%*%wew)-2*t(wew)%*%m
  a<-as.numeric(num/den)

  shrinkage<-pmax(0, pmin(1, a))

  Sigma_usEW<-shrinkage*solve(covus)+(1-shrinkage)*solve(EW)

  return(list("cov"=Sigma_usEW, "shrinkage"=shrinkage))
}

cov_c3 <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
  covc3<-(cov(rets)*((t-1)/t))/c3

  return(covc3)
}

mpcov_c3I <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
  covc3<-(cov(rets)*((t-1)/t))/c3
  c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
  I<-matrix(0, n, n)
  diag(I)<-1
  SR_u<-as.numeric(((t-n-2)*(t(m)%*%solve(cov(rets))%*%m)+n)/t)
  lambda<-t(m)%*%m
  Qhat<-(t(m)%*%cov(rets)%*%m)
  num<-Qhat-(1+(1/c1))*lambda+(1/c1)*SR_u
  den<-(1/c1)*(SR_u+(n/t)-2*lambda)+Qhat
  a<-as.numeric(num/den)

  shrinkage<-pmax(0, pmin(1, a))

  Sigma_c3I<-shrinkage*solve(covc3)+(1-shrinkage)*I

  return(list("cov"=Sigma_c3I, "shrinkage"=shrinkage))
}


mpcov_c3EW <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
  covc3<-(cov(rets)*((t-1)/t))/c3
  EW<-matrix(0, n, n)
  diag(EW)<-n*m
  wew<-rep(1/n, n)
  c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
  SR_u<-as.numeric(((t-n-2)*(t(m)%*%solve(cov(rets))%*%m)+n)/t)
  num<-(t(wew)%*%cov(rets)%*%wew)-(1+(1/c1))*t(wew)%*%m+(1/c1)*SR_u
  den<-(1/c1)*(SR_u+(n/t))+(t(wew)%*%cov(rets)%*%wew)-2*t(wew)%*%m
  a<-as.numeric(num/den)

  shrinkage<-pmax(0, pmin(1, a))

  Sigma_c3EW<-shrinkage*solve(covc3)+(1-shrinkage)*solve(EW)

  return(list("cov"=Sigma_c3EW, "shrinkage"=shrinkage))
}

cov_c <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  SR_u<-as.numeric(((t-n-2)*(t(m)%*%solve(cov(rets))%*%m)+n)/t)
  c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
  covc<-(cov(rets)*((t-1)/t))/c

  return(covc)
}

mpcov_cI <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  SR_u<-as.numeric(((t-n-2)*(t(m)%*%solve(cov(rets))%*%m)+n)/t)
  c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
  covc<-(cov(rets)*((t-1)/t))/c
  c1<-((n-2)*(t-n-2))/((t-n-1)*(t-n-4))
  I<-matrix(0, n, n)
  diag(I)<-1
  lambda<-t(m)%*%m
  Qhat<-(t(m)%*%cov(rets)%*%m)
  num<-(1/c1)*(SR_u^2/(SR_u+(n/t)))-lambda+Qhat-(1/c1)*(SR_u/(SR_u+(n/t)))*lambda
  den<-(1/c1)*((n/t)+SR_u)+Qhat-2*(1/c1)*(SR_u/(SR_u+(n/t)))*lambda
  a<-as.numeric(num/den)

  shrinkage<-pmax(0, pmin(1, a))

  Sigma_cI<-shrinkage*solve(covc)+(1-shrinkage)*I

  return(list("cov"=Sigma_cI, "shrinkage"=shrinkage))
}

mpcov_cEW <- function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  SR_u<-as.numeric(((t-n-2)*(t(m)%*%solve(cov(rets))%*%m)+n)/t)
  c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
  covc<-(cov(rets)*((t-1)/t))/c
  EW<-matrix(0, n, n)
  diag(EW)<-n*m
  wew<-rep(1/n, n)
  c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
  num<-(1/c1)*(SR_u^2/(SR_u+(n/t)))+(t(wew)%*%cov(rets)%*%wew)-(1/c1)*(SR_u/(SR_u+(n/t)))*t(wew)%*%m-t(wew)%*%m
  den<-(1/c1)*(n/t+SR_u)+(t(wew)%*%cov(rets)%*%wew)-2*(1/c1)*(SR_u/(SR_u+(n/t)))*t(wew)%*%m
  a<-as.numeric(num/den)

  shrinkage<-pmax(0, pmin(1, a))

  Sigma_cEW<-shrinkage*solve(covc)+(1-shrinkage)*solve(EW)

  return(list("cov"=Sigma_cEW, "shrinkage"=shrinkage))
}

cov_lwmkt<-function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  D<-rets-m
  mkt<-rowMeans(rets)
  sample<-cov(cbind(rets, mkt))
  covmkt<-sample[1:n, n + 1]
  varmkt<-as.numeric(sample[n + 1, n + 1])
  sample<-sample[-(n + 1), -(n + 1)]

  target<-outer(covmkt, covmkt)/varmkt
  diag(target) <- diag(sample)

  y<-D^2
  z <- sweep(D, 1, mkt, "*")


  # Phi hat
  phiMat <-  crossprod(y)/t - 2 * crossprod(D) * sample/t + (sample^2)
  phi <- sum(apply(phiMat, 2, sum))

  # Rho hat
  rhoMat1 <- 1/t * sweep(crossprod(y, z), 2, covmkt, "*")/varmkt
  rhoMat3 <- 1/t * crossprod(z) * outer(covmkt, covmkt)/varmkt^2
  rhoMat <- 2 * rhoMat1 - rhoMat3 - target * sample
  diag(rhoMat) <- diag(phiMat)
  rho <- sum(apply(rhoMat, 2, sum))

  # Gamma hat
  gamma <- norm(sample - target, "F")^2

  # Kappa hat
  kappa <- (phi - rho)/gamma

  # shrinkage value
  shrinkage <- pmax(0, pmin(1, kappa/t))

  # Sigma hat
  Sigma <- shrinkage * target + (1 - shrinkage) * sample

  return(list("cov"=Sigma, "shrinkage"=shrinkage))
}

cov_lwI<-function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  D<-rets-m
  sample<-cov(rets)

  I<-matrix(0, n, n)
  diag(I) <- 1

  y<-D^2

  # Phi hat
  phiMat <-  crossprod(y)/t - 2 * crossprod(D) * sample/t + (sample^2)/t
  phi <- sum(apply(phiMat, 2, sum))

  # Gamma hat
  gamma <- norm(sample - I, "F")^2

  # Kappa hat
  kappa <- phi/gamma

  # shrinkage value
  shrinkage <- pmax(0, pmin(1, kappa/t))

  # Sigma hat
  Sigma <- shrinkage * I + (1 - shrinkage) * sample

  return(list("cov"=Sigma, "shrinkage"=shrinkage))
}

cov_lwsI<-function(rets) {

  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }

  t<-nrow(rets)
  n<-ncol(rets)
  m<-colMeans(rets)
  D<-rets-m
  sample<-cov(rets)

  I<-matrix(0, n, n)
  diag(I) <- mean(diag(sample))

  y<-D^2

  # Phi hat
  phiMat <-  crossprod(y)/t - 2 * crossprod(D) * sample/t + (sample^2)/t
  phi <- sum(apply(phiMat, 2, sum))

  # Gamma hat
  gamma <- norm(sample - I, "F")^2

  # Kappa hat
  kappa <- phi/gamma

  # shrinkage value
  shrinkage <- pmax(0, pmin(1, kappa/t))

  # Sigma hat
  Sigma <- shrinkage * I + (1 - shrinkage) * sample

  return(list("cov"=Sigma, "shrinkage"=shrinkage))
}
