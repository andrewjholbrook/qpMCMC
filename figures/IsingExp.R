# 1D Ising model
library(scales)


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

multiProp <- function(L=1000,
                      beta=0.1,
                      nProps=100,
                      nIts=1000,
                      thin=100) {
  worst  <- numeric()
  for (i in 1:L) {
    worst <- c(worst, (-1)^i * rep(c(-1,1),L/2) )
  }
  #logProb   <- rep(0,nIts)
  state_0   <- matrix(worst,L,L)
  buffered  <- matrix(0,L+2,L+2)
  buffered[2:(L+1),2:(L+1)] <- state_0
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
    propIndices <- cbind(sample(x=2:(L+1),size=nProps,replace=TRUE),sample(2:(L+1),size=nProps,replace=TRUE))

    targets <- rep(0,nProps+1)
    for(j in 1:nProps) {
      nbhdVls <- rep(0,4)
      prpsdVl    <- -1 * currentState[propIndices[j,1],propIndices[j,2]]
      nbhdVls[1] <- currentState[propIndices[j,1]+1,propIndices[j,2]] 
      nbhdVls[2] <- currentState[propIndices[j,1]-1,propIndices[j,2]] 
      nbhdVls[3] <- currentState[propIndices[j,1],propIndices[j,2]+1] 
      nbhdVls[4] <- currentState[propIndices[j,1],propIndices[j,2]-1] 
      targets[j+1] <- beta*sum(prpsdVl*nbhdVls)
    }
    
    selectionProbs  <- targets
    lambdas         <- -selectionProbs + log(rexp(n=nProps+1))
    currentAndProps <- rbind(c(1,1),propIndices)
    currentAndProps <- currentAndProps[order(lambdas),]
    rank1           <- rank(lambdas)[1]
    qmin            <- quantumMin(field=1:(nProps+1),y=rank1)
    propIndex       <- qmin[[1]]
    oracleCalls     <- oracleCalls + qmin[[2]]
    
    #chain[[i]] <- chain[[i-1]]
    currentState[currentAndProps[propIndex,1],currentAndProps[propIndex,2]] <-
      - currentState[currentAndProps[propIndex,1],currentAndProps[propIndex,2]]
    if(any(currentAndProps[propIndex,] != currentIndex)) {
      accept[i] <- 1
      currentIndex <- currentAndProps[propIndex,]
    }
    currentIndices[i,] <- currentIndex
    
    if(i %% 100 == 0){
      cat("Iteration: ",i,"\n")
    }
    if(i %% thin == 0) {
      l <- l + 1
      chain[[l]] <- currentState
    }
  }
  return(list(chain,accept,currentIndices,oracleCalls))
}

################################################################################
#
####
####### scaling experiment
####
#
set.seed(1)
nIts <- 1000000
beta <- 1
thin <- 1000
for(nProps in c(4,8,16,32,64,128,256,512,1024,2048)) {
  out <- multiProp(L=500,nIts=nIts, beta=beta, nProps = nProps,thin=thin) 
  distances <- rep(0,(nIts/thin))
  for(i in 1:(nIts/thin)) {
    distances[i] <- sum(out[[1]][[1]] != out[[1]][[i]] )
  }
  esss <- coda::effectiveSize(distances[(nIts/(thin*2)+1):(nIts/thin)])
  cat(nProps,esss,out[[4]],"\n",file="~/qpMCMC/Ising2D.txt",append=TRUE)
  
  for(i in 1:(nIts/thin)) {
    cat(nProps,distances[i],"\n",append=TRUE,file = "~/qpMCMC/Ising2dDistances.txt")
  }
}




