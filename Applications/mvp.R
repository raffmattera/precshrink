mvp <- function(rets, g, control = list()) {
  
  if (missing(rets)) {
    stop("rets is missing")
  }
  if (!is.matrix(rets)) {
    stop("rets must be a (T x n) matrix")
  }
  if (is.null(g)) {
    stop("Relative risk aversion gamma is required to compute the mean-variance portfolio")
  }
  
  ctr <- .ctrCov(control)
  
  if (ctr$type[1] == "sc") {
    SigmaInv <- .SigmaInv_sc(rets)
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else if (ctr$type[1] == "us") {
    SigmaInv <- .SigmaInv_us(rets)
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else if (ctr$type[1] == "c3") {
    SigmaInv <- .SigmaInv_c3(rets)
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else if (ctr$type[1] == "c") {
    SigmaInv <- .SigmaInv_c(rets)
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else if (ctr$type[1] == "I") {
    
    n<-ncol(rets)
    I<-matrix(0, n, n)
    diag(I)<-1
    SigmaInv <- I
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else if (ctr$type[1] == "mp_usI") {
    
    n<-ncol(rets)
    I<-matrix(0, n, n)
    diag(I)<-1
    SigmaInv <- .SigmaInv_us(rets)
    a  <- .mp_usI(rets = rets)
    w<-(a/g)*SigmaInv%*%colMeans(rets)+((1-a)/g)*I%*%colMeans(rets)
    
  } else if (ctr$type[1] == "mp_usEW") {
    SigmaInv <- .SigmaInv_us(rets)
    a <- .mp_usEW(rets)
    w<-(a/g)*SigmaInv%*%colMeans(rets)+((1-a)/g)*rep(1/ncol(rets), ncol(rets))
    
  } else if (ctr$type[1] == "mp_c3I") {
    n<-ncol(rets)
    I<-matrix(0, n, n)
    diag(I)<-1
    SigmaInv <- .SigmaInv_c3(rets)
    a <- .mp_c3I(rets)
    w<-(a/g)*SigmaInv%*%colMeans(rets)+((1-a)/g)*I%*%colMeans(rets)
    
  } else if (ctr$type[1] == "mp_c3EW") {
    SigmaInv <- .SigmaInv_c3(rets)
    a <- .mp_c3EW(rets)
    w<-(a/g)*SigmaInv%*%colMeans(rets)+((1-a)/g)*rep(1/ncol(rets), ncol(rets))
    
  } else if (ctr$type[1] == "mp_cI") {
    n<-ncol(rets)
    I<-matrix(0, n, n)
    diag(I)<-1
    SigmaInv <- .SigmaInv_c(rets)
    a <- .mp_cI(rets)
    w<-(a/g)*SigmaInv%*%colMeans(rets)+((1-a)/g)*I%*%colMeans(rets)
    
  } else if (ctr$type[1] == "mp_cEW") {
    SigmaInv <- .SigmaInv_c(rets)
    a <- .mp_cEW(rets)
    w<-(a/g)*SigmaInv%*%colMeans(rets)+((1-a)/g)*rep(1/ncol(rets), ncol(rets))
    
  } else if (ctr$type[1] == "lw_mkt") {
    SigmaInv <- .lwSigmaInv_mkt(rets)
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else if (ctr$type[1] == "lw_I") {
    SigmaInv <- .lwSigmaInv_I(rets)
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else if (ctr$type[1] == "nls") {
    SigmaInv <- .nlsSigmaInv(rets)
    w<-(1/g)*SigmaInv%*%colMeans(rets)
    
  } else {
    
    stop("control$type is not well defined")
  }
  
  return(w)
}

  .ctrCov <- function(control = list()) {

  if (!is.list(control)) {
    stop("control must be a list")
  }
  if (length(control) == 0) {
    control <- list(type = "sc")
  }
  nam <- names(control)
  if (!("type" %in% nam) || is.null(control$type)) {
    control$type <- c("sc", "us", "c3", "c", "I", "mp_usI", 
                      "mp_usEW", "mp_c3I", "mp_c3EW", "mp_cI", "mp_cEW",
                      "lw_mkt", "lw_I")
  }
  
  return(control)
}
 
.SigmaInv_sc <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
    t<-nrow(rets)
    n<-ncol(rets)
    
    if (t>n+4) {
      
      covsc<-cov(rets)
      
      return(solve(covsc))
      
    } else {
      
      covsc<-MASS::ginv(cov(rets))
      
      return(covsc)
  }
}
  
.SigmaInv_us <- function(rets) {
  
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
   D <- rets-colMeans(rets)
   t<-nrow(rets)
   n<-ncol(rets)
  
    if (t>n+4) {
    
    covus<- (1/(t-n-2))*t(D)%*%D
    
    return(solve(covus))
    
    } else {
    
      covus<- (1/(t-n-2))*t(D)%*%D
      covus<-MASS::ginv(covus)
      
      return(covus)
    }
  }
  
.SigmaInv_c3 <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
  
    D <- rets-colMeans(rets)
    t<-nrow(rets)
    n<-ncol(rets)
    
    if (t>n+4) {
      
    c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
    covc3 <- ((1/t)*(t(D)%*%D))/c3
    
    return(solve(covc3))
    
    } else {
    
      c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
      covc3 <- ((1/t)*(t(D)%*%D))/c3
      covc3<-MASS::ginv(covc3)
      
      return(covc3)
    }
  }
  
.SigmaInv_c <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
  
  D <- rets-colMeans(rets)
  t<-nrow(rets)
  n<-ncol(rets)
  
    if (t>n+4) {
    
    SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
    c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
    covc <- ((1/t)*(t(D)%*%D))/c
    
    return(solve(covc))
    
    } else {
    
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
      covc <- ((1/t)*(t(D)%*%D))/c
      covc<-MASS::ginv(covc)
      
      return(covc)
    }
  }
  
.mp_usI <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
  
  D <- rets-colMeans(rets)
  t<-nrow(rets)
  n<-ncol(rets)
    
    if (t>n+4) {
     
      I<-matrix(0, n, n)
      diag(I)<-1
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      lambda<-t(colMeans(rets))%*%colMeans(rets)
      Qhat<-(t/(t-1))*(t(colMeans(rets))%*%cov(rets)%*%colMeans(rets))
      c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
      num<-Qhat-2*lambda+SR_u
      den<-c1*(SR_u+(n/t))+Qhat-2*lambda
      a<-as.numeric(num/den)
      covus<-(1/(t-n-2))*t(D)%*%D
      
      shrinkage<-pmax(0, pmin(1, a))
      
      return(shrinkage)
      
    } else {
    
        I<-matrix(0, n, n)
        diag(I)<-1
        SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*% MASS::ginv(cov(rets))%*%colMeans(rets))+n)/t)
        lambda<-t(colMeans(rets))%*%colMeans(rets)
        Qhat<-(t/(t-1))*(t(colMeans(retw))%*%cov(rets)%*%colMeans(rets))
        c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
        num<-Qhat-2*lambda+SR_u
        den<-c1*(SR_u+(n/t))+Qhat-2*lambda
        a<-as.numeric(num/den)
        covus<-(1/(t-n-2))*t(D)%*%D
    
    return(shrinkage)
    
    }
}
  
.mp_usEW <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
  
  D <- rets-colMeans(rets)
  t<-nrow(rets)
  n<-ncol(rets)
    
    if (t>n+4) {
      
      covus<-(1/(t-n-2))*t(D)%*%D
      
      EW<-matrix(0, n, n)
      diag(EW)<-n*colMeans(rets)
      wew<-rep(1/n, n)
      c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      num<-(t(wew)%*%cov(rets)%*%wew)-2*(t(wew)%*%colMeans(rets))+SR_u
      den<-(c1*SR_u+c1*(n/t))+(t(wew)%*%cov(rets)%*%wew)-2*t(wew)%*%colMeans(rets)
      a<-as.numeric(num/den)
      
      shrinkage<-pmax(0, pmin(1, a))

      return(shrinkage) 
    } else {
    
    mEW<-matrix(0, n, n)
    diag(EW)<-n*colMeans(rets)
    wew<-rep(1/n, n)
    c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
    SR_u<-as.numeric(((t-n-2)*(t(m)%*%MASS::ginv(cov(rets))%*%m)+n)/t)
    num<-(t(wew)%*%cov(rets)%*%wew)-2*(t(wew)%*%colMeans(rets))+SR_u
    den<-(c1*SR_u+c1*(n/t))+(t(wew)%*%cov(rets)%*%wew)-2*t(wew)%*%colMeans(rets)
    a<-as.numeric(num/den)
    
    shrinkage<-pmax(0, pmin(1, a))

    return(shrinkage)
  }
}
  
  
.mp_c3I <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
  
    D <- rets-colMeans(rets)
    t<-nrow(rets)
    n<-ncol(rets)
    
    if (t>n+4) {
      
    c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
    covc3 <- ((1/t)*(t(D)%*%D))/c3
    
    c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
    I<-matrix(0, n, n)
    diag(I)<-1
    SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
    lambda<-t(colMeans(rets))%*%colMeans(rets)
    Qhat<-(t/(t-1))*(t(colMeans(rets))%*%cov(rets)%*%colMeans(rets))
    num<-Qhat-(1+(1/c1))*lambda+(1/c1)*SR_u
    den<-(1/c1)*(SR_u+(n/t)-2*lambda)+Qhat
    a<-as.numeric(num/den)
    
    shrinkage<-pmax(0, pmin(1, a))

    return(shrinkage)
    } else {
      
      c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
      covc3 <- ((1/t)*(t(D)%*%D))/c3
      
      c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
      I<-matrix(0, n, n)
      diag(I)<-1
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%MASS::ginv(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      lambda<-t(colMeans(rets))%*%colMeans(rets)
      Qhat<-(t/(t-1))*(t(colMeans(rets))%*%cov(rets)%*%colMeans(rets))
      num<-Qhat-(1+(1/c1))*lambda+(1/c1)*SR_u
      den<-(1/c1)*(SR_u+(n/t)-2*lambda)+Qhat
      a<-as.numeric(num/den)
      
      shrinkage<-pmax(0, pmin(1, a))

      return(shrinkage)
  }
}
  
  .mp_c3EW <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
    D <- rets-colMeans(rets)
    t<-nrow(rets)
    n<-ncol(rets)
    
    if (t>n+4){
      
    c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
    covc3 <- ((1/t)*(t(D)%*%D))/c3
      
      
    EW<-matrix(0, n, n)
    diag(EW)<-n*colMeans(rets)
    wew<-rep(1/n, n)
    c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
    SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
    num<-(t(wew)%*%cov(rets)%*%wew)-(1+(1/c1))*t(wew)%*%colMeans(rets)+(1/c1)*SR_u
    den<-(1/c1)*(SR_u+(n/t))+(t(wew)%*%cov(rets)%*%wew)-2*t(wew)%*%colMeans(rets)
    a<-as.numeric(num/den)
    
    shrinkage<-pmax(0, pmin(1, a))

    return(shrinkage)
    } else {
      
      c3<-as.numeric(((t-n-1)*(t-n-4))/(t*(t-2)))
      covc3 <- ((1/t)*(t(D)%*%D))/c3
      
      
      EW<-matrix(0, n, n)
      diag(EW)<-n*colMeans(rets)
      wew<-rep(1/n, n)
      c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%MASS::ginv(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      num<-(t(wew)%*%cov(rets)%*%wew)-(1+(1/c1))*t(wew)%*%colMeans(rets)+(1/c1)*SR_u
      den<-(1/c1)*(SR_u+(n/t))+(t(wew)%*%cov(rets)%*%wew)-2*t(wew)%*%colMeans(rets)
      a<-as.numeric(num/den)
      
      shrinkage<-pmax(0, pmin(1, a))

      return(shrinkage)
    }
  }
  
  
  .mp_cI <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
    D <- rets-colMeans(rets)
    t<-nrow(rets)
    n<-ncol(rets)
    
    if (t>n+4) {
      
    SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
    c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
    
    covc <- ((1/t)*(t(D)%*%D))/c
    c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
    I<-matrix(0, n, n)
    diag(I)<-1
    lambda<-t(colMeans(rets))%*%colMeans(rets)
    Qhat<-(t/(t-1))*(t(colMeans(rets))%*%cov(rets)%*%colMeans(rets))
    
    num<-(1/c1)*(SR_u^2/(SR_u+(n/t)))-lambda+Qhat-(1/c1)*(SR_u/(SR_u+(n/t)))*lambda
    den<-(1/c1)*((n/t)+SR_u)+Qhat-2*(1/c1)*(SR_u/(SR_u+(n/t)))*lambda
    a<-as.numeric(num/den)

    shrinkage<-pmax(0, pmin(1, a))
    
    return(shrinkage)
    
    } else {
      
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%MASS::ginv(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
      
      covc <- ((1/t)*(t(D)%*%D))/c
      c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
      I<-matrix(0, n, n)
      diag(I)<-1
      lambda<-t(colMeans(rets))%*%colMeans(rets)
      Qhat<-(t/(t-1))*(t(colMeans(rets))%*%cov(rets)%*%colMeans(rets))
      
      num<-(1/c1)*(SR_u^2/(SR_u+(n/t)))-lambda+Qhat-(1/c1)*(SR_u/(SR_u+(n/t)))*lambda
      den<-(1/c1)*((n/t)+SR_u)+Qhat-2*(1/c1)*(SR_u/(SR_u+(n/t)))*lambda
      a<-as.numeric(num/den)
      
      shrinkage<-pmax(0, pmin(1, a))
      
      return(shrinkage)
    }
  }
  
  .mp_cEW <- function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
    D <- rets-colMeans(rets)
    t<-nrow(rets)
    n<-ncol(rets)
    
    if (t>n+4){
      
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%solve(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
      covc <- ((1/t)*(t(D)%*%D))/c
      
      EW<-matrix(0, n, n)
      diag(EW)<-n*colMeans(rets)
      wew<-rep(1/n, n)
      c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
      num<-(1/c1)*(SR_u^2/(SR_u+(n/t)))+(t(wew)%*%cov(rets)%*%wew)-(1/c1)*(SR_u/(SR_u+(n/t)))*t(wew)%*%colMeans(rets)-t(wew)%*%colMeans(rets)
      den<-(1/c1)*(n/t+SR_u)+(t(wew)%*%cov(rets)%*%wew)-2*(1/c1)*(SR_u/(SR_u+(n/t)))*t(wew)%*%colMeans(rets)
      a<-as.numeric(num/den)
      
      shrinkage<-pmax(0, pmin(1, a))
      
      return(shrinkage)
    } else {
      
      SR_u<-as.numeric(((t-n-2)*(t(colMeans(rets))%*%MASS::ginv(((1/t)*(t(D)%*%D)))%*%colMeans(rets))+n)/t)
      c<-as.numeric((((t-n-1)*(t-n-4))/(t*(t-2)))*(SR_u/(SR_u+(n/t))))
      covc <- ((1/t)*(t(D)%*%D))/c
      
      EW<-matrix(0, n, n)
      diag(EW)<-n*colMeans(rets)
      wew<-rep(1/n, n)
      c1<-((t-2)*(t-n-2))/((t-n-1)*(t-n-4))
      num<-(1/c1)*(SR_u^2/(SR_u+(n/t)))+(t(wew)%*%cov(rets)%*%wew)-(1/c1)*(SR_u/(SR_u+(n/t)))*t(wew)%*%colMeans(rets)-t(wew)%*%colMeans(rets)
      den<-(1/c1)*(n/t+SR_u)+(t(wew)%*%cov(rets)%*%wew)-2*(1/c1)*(SR_u/(SR_u+(n/t)))*t(wew)%*%colMeans(rets)
      a<-as.numeric(num/den)
      
      shrinkage<-pmax(0, pmin(1, a))
      
      return(shrinkage)
  }
}
  
  .lwSigmaInv_mkt<-function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
    t<-nrow(rets)-1
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
    
    return(solve(Sigma))
  }
  
  .lwSigmaInv_I<-function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
    t<-nrow(rets)-1
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
    
    return(solve(Sigma))
  }
  
  .nlsSigmaInv<-function(rets) {
    
    if (missing(rets)) {
      stop("rets is missing")
    }
    if (!is.matrix(rets)) {
      stop("rets must be a (T x n) matrix")
    }
    
    t<-nrow(rets)-1
    n<-ncol(rets)
    m<-colMeans(rets)
    D<-rets-m
    sample<-cov(rets)
    
    S <- crossprod(D)/t
    S_eigen <- eigen(sample)
    U <- S_eigen$vectors
    lambda <- S_eigen$values
    lambdasort <- sort(lambda)
    lambdaorder <- order(lambda)
    tau_est <- nlshrink:::tau_estimate(X = D, k = 0, method = "nlminb")
    nlshrink_tau <- nlshrink:::nlshrink_est(nlshrink:::QuEST(tau_est, t))
    
    Sigma <- U %*% (diag(nlshrink_tau[lambdaorder]) %*% t(U))
    
    return(solve(Sigma))
  }
  