# 1D Ising model
library(scales)

#nProps <- R.utils::cmdArg("nProps",42L)

groverProbs <- function(N=2^10,M=1) {
  Iter  <- ceiling(pi*sqrt(N/M)/4)
  marks <- c(rep(1,M),rep(0,N-M))
  sqrtProbs <- matrix(0,Iter,N)
  probSucceed <- rep(0,Iter)
  sqrtProbs[1,] <- 1/sqrt(N)
  probSucceed[1] <- M/N
  
  for(i in 2:Iter) {
    sqrtProbs[i,]  <- -marks*sqrtProbs[i-1,] + (1-marks)*sqrtProbs[i-1,]
    sqrtProbs[i,]  <- -sqrtProbs[i,] + 2*mean(sqrtProbs[i,])
    probSucceed[i] <- sum(marks*sqrtProbs[i,]^2)
  }
  return(probSucceed)
}


chebPoly <- function(x,L){
  if(abs(x)<=1){
    return( cos(L*acos(x)) )
  } else if (x>1) {
    return( cosh(L*acosh(x)) )
  } else {
    stop("Error")
  }
}

RPhiTheta <- function(theta,phi) {
  Z <- matrix(c(1,0,0,-1),2,2)
  X <- matrix(c(0,1,1,0),2,2)
  output <- complexplus::matexp( -0.5*1i*theta*(cos(phi)*Z+sin(phi)*X) )
  return(output)
}

genGroverProbs <- function(N=2^12,M=1,delta=sqrt(0.1),L=NULL,w=NULL) {
  if(is.null(w)) w <- 1/N
  if(is.null(L))  L <- ceiling( log(2/delta)/sqrt(w) )
  if(L %% 2 == 0) L <- L + 1
  l <- (L-1)/2 # iterations
  alphas <- rep(0,l)
  betas  <- rep(0,l)
  gamma  <- 1/chebPoly(1/delta,1/L)
  
  sqrtProbs <- matrix(0,l+1,2)
  probSucceed <- rep(0,l+1)
  sqrtProbs[1,] <- c(sqrt(1-M/N),sqrt(M/N))
  probSucceed[1] <- M/N
  phi    <- 2*asin(sqrt(M/N))
  
  for(j in 1:l) {
    alphas[j]    <- 2*atan( 1/(tan(2*pi*j/L)*sqrt(1-gamma^2)) )
    betas[l-j+1] <- -alphas[j]
  }
  
  for(i in 2:(l+1)) {
    sqrtProbs[i,]  <- matrix(c(1,0,0,exp(1i*betas[i-1])),2,2) %*% sqrtProbs[i-1,]
    sqrtProbs[i,]  <- -matrix(c(1-(1-exp(-1i*alphas[i-1]))*(1-M/N),
                                -(1-exp(-1i*alphas[i-1]))*sqrt((1-M/N)*M/N),
                                -(1-exp(-1i*alphas[i-1]))*sqrt((1-M/N)*M/N),
                                1-(1-exp(-1i*alphas[i-1]))*(M/N)),2,2) %*% sqrtProbs[i,]
    probSucceed[i] <- Mod(sqrtProbs[i,2])^2
  }
  return( list(probSucceed,L-1) )
}

expoSearch <- function(N=2^10,M=1,marks=NULL) {
  success <- FALSE
  m <- 1
  grovIter <- 0
  lmbd <- 6/5
  if(is.null(marks)) {
    cutoff <- sqrt(N)
    marks <- c(rep(1,M),rep(0,N-M))
  } else {
    M <- sum(marks==1)
    N <- length(marks)
    cutoff <- 9*sqrt(N)/4 # in mcmc apps
  }
  while((!success) & (grovIter<cutoff)) {
    j <- sample(size=1,x=0:(m-1))
    sqrtProbs <- rep(1/sqrt(N),N)
    i <- 1
    while(i<=j) {
      sqrtProbs  <- -marks*sqrtProbs + (1-marks)*sqrtProbs
      sqrtProbs  <- -sqrtProbs + 2*mean(sqrtProbs)
      i <- i + 1
    }
    draw <- sample(size=1,x=1:N,prob=sqrtProbs^2)
    grovIter <- grovIter + j
    if(draw<=M) {
      success <- TRUE
    } else {
      m <- min(lmbd*m,sqrt(N))
    }
  }
  return(list(grovIter,draw,success))
}

fixedPointSearch <- function(N=2^10,M=1,marks=NULL,delta=sqrt(0.1)) {
  # 2D representation
  if(is.null(marks)) {
    marks <- c(rep(1,M),rep(0,N-M))
  } else {
    M <- sum(marks==1)
    N <- length(marks)
  }
  w <- 1/N # width
  L <- ceiling( log(2/delta)/sqrt(w) )
  if (L %% 2 == 0) L <- L + 1
  l <- (L-1)/2 # iterations
  alphas <- rep(0,l)
  betas  <- rep(0,l)
  gamma  <- 1/chebPoly(1/delta,1/L)
  
  s         <- rep(sqrt(1/N),N)
  S         <- outer(s,s)
  tMat      <- outer(marks,marks)
  sqrtProbs <- s
  probSucceed <- rep(0,l+1)
  probSucceed[1] <- M/N
  phi    <- 2*asin(sqrt(M/N))
  
  for(j in 1:l) {
    alphas[j]    <- 2*atan( 1/(tan(2*pi*j/L)*sqrt(1-gamma^2)) )
    betas[l-j+1] <- -alphas[j]
  }
  
  for(j in 1:l) {
    mult1 <- exp(1i*betas[j]) #/ Mod(1-exp(1i*betas[j]))
    mult2 <- exp(-1i*alphas[j]) #/ Mod(1-exp(-1i*alphas[j]))
    sqrtProbs <-  (1-(1-mult1)*marks)*sqrtProbs
    sqrtProbs  <- -sqrtProbs + (1-mult2)*mean(sqrtProbs)
    # sqrtProbs  <- #(diag(N) - (1-exp(1i*betas[j]))*tMat/M) %*% sqrtProbs
    #  sqrtProbs  <- #((1-exp(-1i*alphas[j]))*S - diag(N)) %*% sqrtProbs
    probSucceed[j+1] <- sum(marks*Mod(sqrtProbs)^2)
  }
  draw <- sample(size=1,x=1:N,prob=Mod(sqrtProbs)^2)
  success <- draw <= M
  return(list(L-1,draw,success,probSucceed))
}

quantumMin <- function(field, y=NULL) {
  N <- length(field)
  done <- FALSE
  if(is.null(y)) {
    y <- sample(size=1,x=1:N)
  }
  oracleCalls <- 0
  while (!done) {
    marks <- as.numeric(field < field[y])
    res <- expoSearch(marks = marks)
    oracleCalls <- oracleCalls + res[[1]]
    if (res[[3]]) {
      y <- res[[2]]
    } else {
      done <- TRUE
    }
  }
  return(list(y,oracleCalls))
}

squaredError <- function(Y,X,mu0,mu1) {
  out <- 0
  for(i in 1:dim(Y)[1]) {
    for(j in 1:dim(Y)[2]) {
      out <- out + (X[i,j]<0)*(Y[i,j]-mu0)^2 +
        (X[i,j]>0)*(Y[i,j]-mu1)^2
    }
  }
  return(out)
}

multiPropGibbs <- function(Y=NULL,
                           L=1000,
                           beta=0.1,
                           nProps=100,
                           nIts=1000,
                           chainThin=100,
                           logProbThin=10) {
  if( is.null(Y) ) stop("Data required.")
  initial    <- (Y > 100) * 2 - 1
  
  precision0  <- 2/10000
  #  precision0s <- precision0
  #  precision1s <- precision1
  
  mu0        <- 1
  # mu0s       <- mu0
  mu1        <- 255
  #  mu1s       <- mu1
  # logProbs   <- 0
  logProb    <- 0
  state_0   <- matrix(initial,L,L)
  buffered  <- matrix(0,L+2,L+2)
  buffered[2:(L+1),2:(L+1)] <- state_0
  
  buffY                  <- matrix(0,L+2,L+2)
  buffY[2:(L+1),2:(L+1)] <- Y 
  Y                      <- buffY

    # precompute terms for mu updates
  n0 <- sum(buffered<0)
  n1 <- sum(buffered>0)
  sumY1 <- sum((buffered>0)*Y)
  sumY0 <- sum((buffered<0)*Y)
  
  
  chain <- list()
  chain[[1]] <- buffered
  l <- 1
  currentState <- buffered
  currentIndex <- c(2,2)
  currentIndices <- matrix(0,nIts,2)
  currentIndices[1,] <- currentIndex
  accept <- rep(0,nIts-1)
  oracleCalls <- 0
  for (i in 2:nIts) {
    propIndices <- cbind(sample(x=2:(L+1),size=nProps,replace=TRUE),
                         sample(2:(L+1),size=nProps,replace=TRUE))
    
    targets <- rep(0,nProps+1)
    targets[0] <- 0
    for(j in 1:nProps) {
      nbhdVls <- rep(0,4)
      prpsdVl    <- -1 * currentState[propIndices[j,1],propIndices[j,2]]
      nbhdVls[1] <- currentState[propIndices[j,1]+1,propIndices[j,2]] 
      nbhdVls[2] <- currentState[propIndices[j,1]-1,propIndices[j,2]] 
      nbhdVls[3] <- currentState[propIndices[j,1],propIndices[j,2]+1] 
      nbhdVls[4] <- currentState[propIndices[j,1],propIndices[j,2]-1] 
      targets[j+1] <- beta*sum(prpsdVl*nbhdVls)
      if ( prpsdVl==1 ) {
        targets[j+1] <- targets[j+1] -
          precision0/2*(Y[propIndices[j,1],propIndices[j,2]]- mu1)^2 +
          precision0/2*(Y[propIndices[j,1],propIndices[j,2]]- mu0)^2
      } else {
        targets[j+1] <- targets[j+1] +
          precision0/2*(Y[propIndices[j,1],propIndices[j,2]]- mu1)^2 -
          precision0/2*(Y[propIndices[j,1],propIndices[j,2]]- mu0)^2
      }
    }
    
    selectionProbs  <- targets 
    lambdas         <- -selectionProbs + log(rexp(n=nProps+1))
    currentAndProps <- rbind(c(1,1),propIndices)
    currentAndProps <- currentAndProps[order(lambdas),]
    targets         <- targets[order(lambdas)]
    rank1           <- rank(lambdas)[1]
    qmin            <- quantumMin(field=1:(nProps+1),y=rank1)
    propIndex       <- qmin[[1]]
    oracleCalls     <- oracleCalls + qmin[[2]]
    logProb         <- logProb + targets[propIndex]
    
    currentState[currentAndProps[propIndex,1],currentAndProps[propIndex,2]] <-
      - currentState[currentAndProps[propIndex,1],currentAndProps[propIndex,2]]
    
    if(currentState[currentAndProps[propIndex,1],currentAndProps[propIndex,2]]==1) {
      sumY1   <- sumY1 + Y[currentAndProps[propIndex,1],currentAndProps[propIndex,2]]
      sumY0   <- sumY0 - Y[currentAndProps[propIndex,1],currentAndProps[propIndex,2]]
      n1      <- n1 + 1
      n0      <- n0 - 1
    } else {
      sumY1   <- sumY1 - Y[currentAndProps[propIndex,1],currentAndProps[propIndex,2]]
      sumY0   <- sumY0 + Y[currentAndProps[propIndex,1],currentAndProps[propIndex,2]]
      n1      <- n1 - 1
      n0      <- n0 + 1
    }
    
    if(any(currentAndProps[propIndex,] != currentIndex)) {
      accept[i] <- 1
      currentIndex <- currentAndProps[propIndex,]
    }
    currentIndices[i,] <- currentIndex
    
    # # update means
    # # prior mean 255/2, prior precision 1 # 83068.88 (half all divided by 100) 
    # mu0Prop <- 257
    # while(mu0Prop > 255 | mu0Prop < 0) {
    #   mu0Prop <- rnorm(n=1, mean=(10*255/2+precision0*sumY0)/(10+n0*precision0),
    #                    sd=1/sqrt(10+n0*precision0))
    # }
    # mu0 <- mu0Prop
    # 
    # mu1Prop <- 257
    # while(mu1Prop > 255 | mu1Prop < 0) {
    #   mu1Prop <- rnorm(n=1, mean=(10*255/2+precision0*sumY1)/(10+n1*precision0),
    #                    sd=1/sqrt(10+n1*precision0))
    # }
    # mu1 <- mu1Prop
    # browser()
    
    mu0 <- truncnorm::rtruncnorm(n=1, mean=sumY0/n0,sd=1/sqrt(precision0*n0),
                                 a=0,b=255.00001)
    mu1 <- truncnorm::rtruncnorm(n=1, mean=sumY1/n1,sd=1/sqrt(precision0*n1),
                                 a=0,b=255.00001)
    if(sumY1/n1>255) stop("somethings wrong")

    
    #update precisions rarely because expensive
    if(i %% logProbThin == 0) {

      sqrdRrr0 <- sum((currentState<0)*(Y-mu0)^2)
      sqrdRrr1 <- sum((currentState>0)*(Y-mu1)^2)
      
      precision0 <- rgamma(n=1,shape=(n0+n1+1)/2,
                           rate=(sqrdRrr0+sqrdRrr1+1)/2)
    }
    
    
    
    
    
    if(i %% 100 == 0){
      cat("Iteration: ",i,"\n")
    }
    if(i %% chainThin == 0) {
      l <- l + 1
      chain[[l]] <- currentState
    }
    if(i %% logProbThin == 0) {
      cat(i, logProb, mu0,mu1, precision0,
          oracleCalls, "\n", append=TRUE,
          file="~/qpMCMC/blackHoleResultsCenteredPrior.txt")
    }
  }
  return(list(chain,accept,currentIndices))
}

# ################################################################################
# #
# ####
# ####### blackhole experiment
# ####
# #

# library(ggplot2)
# #library(reshape2)
# # import grayscale image flips vertically FSR
# img <- jpeg::readJPEG("~/qpMCMC/figures/eso2208-eht-mwa.jpg")
# img <- img[dim(img)[1]:1,,1] * 255 #> 100
# saveRDS(img, "~/qpMCMC/figures/blackHoleIntensity.rds")
# imgLong <- melt(img)
# ggplot(imgLong, aes(x = Var2, y = Var1)) +
#   geom_raster(aes(fill=value))

# implement for beta=1.2

img <- readRDS("~/qpMCMC/figures/blackHoleIntensity.rds")

set.seed(1)
nIts <- 100000000
beta <- 1.2
thin <- 200000
out <- multiPropGibbs(Y=img,L=4076,nIts=nIts, beta=beta, nProps = 1024,
                      chainThin=thin, logProbThin = thin/10)

saveRDS(out,file="~/qpMCMC/blackHoleCenteredPrior.rds")

# imgLong <- reshape2::melt(out[[1]][[3]])
# #remove(results)
# ggplot(imgLong, aes(x = Var2, y = Var1)) +
#   geom_raster(aes(fill=value))





