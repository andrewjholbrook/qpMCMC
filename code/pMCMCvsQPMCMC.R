setwd("~/qpMCMC/")

library(ggplot2)
library(wesanderson)
library(forcats)
pal <- wes_palette("Zissou1", 5, type = "discrete")

################################################################################

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

# groverProbsExp <- function(N=2^10,M=1) {
#   Iter  <- ceiling(pi*sqrt(N/M)/4)
#   marks <- c(rep(1,M),rep(0,N-M))
#   sqrtProbs <- matrix(0,Iter,N)
#   probSucceed <- rep(0,Iter) 
#   sqrtProbs[1,] <- 1/sqrt(N)
#   probSucceed[1] <- M/N
#   
#   for(i in 2:Iter) {
#     sqrtProbs[i,]  <- sqrtProbs[i-1,] - 2*marks*mean(sqrtProbs[i-1,marks])  
#     sqrtProbs[i,]  <- -sqrtProbs[i,] + 2*mean(sqrtProbs[i,])
#     probSucceed[i] <- sum(marks*sqrtProbs[i,]^2)
#   }
#   return(probSucceed)
# }

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

proposal <- function(currentState, sigma=1, nProps=100) {
  # centered Gaussian multiproposal
  dimensions <- length(currentState)
  center <- currentState + rnorm(dimensions,sd=sigma)
  props <- matrix(center,nProps,dimensions,byrow=TRUE)
  props <- props + rnorm(n=nProps*dimensions,sd=sigma)
  return(props)
}

target <- function(state,tau=1) {
  # isotropic gaussian pdf
  return( -sum(state^2)/(2*tau^2) )
}

delta <- function(n) {
  return( n^(-0.5) )
}

pMCMC <- function (nIts=1000,nProps=2000,dims=10,sigma=3,targetAccept=0.5,initial=100) {
  
  states <- matrix(0,nIts,dims)
  states[1,] <- initial
  numberLess <- rep(0,nIts)
  batches    <- 1
  
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 20   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  
  for(i in 2:nIts) {
    props           <- proposal(states[i-1,], nProps=nProps, sigma = sigma)
    currentAndProps <- rbind(states[i-1,],props)
    targetVals      <- apply(currentAndProps,MARGIN = 1,FUN=target)
    selectionProbs  <- targetVals
    lambdas         <- selectionProbs - log(rexp(n=nProps+1))
    numberLess[i]   <- sum(lambdas < lambdas[1])
    propIndex       <- which(lambdas==max(lambdas))  #sample(1:length(selectionProbs),size=1,prob=selectionProbs)
    states[i,]      <- currentAndProps[propIndex,]
    SampCount       <- SampCount + 1
    
    if(propIndex != 1) {
      Acceptances <- Acceptances + 1
    }    
    if(i %% 100 == 0) cat("Iteration ", i, " complete.\n")
    if (SampCount == SampBound) { 
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        sigma <- sigma * min(2,(1 + delta(batches)))
      } else {
        sigma <- sigma * max(0.5,(1 - delta(batches)))
      }
      batches <- batches + 1
      cat("Acceptance ratio: ", AcceptRatio, ". sigma: ", sigma, "\n")
      
      SampCount <- 0
      Acceptances <- 0
    }
  }
  return(list(states,numberLess))
}

qpMCMC <- function (nIts=1000,nProps=2000,dims=10,sigma=3,targetAccept=0.5,
                    initial=100,earlyStop=FALSE) {
  
  states <- matrix(0,nIts,dims)
  states[1,] <- initial
  numberLess <- rep(0,nIts)
  batches    <- 1
  oracleCalls <- 0
  
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 20   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  
  for(i in 2:nIts) {
    props           <- proposal(states[i-1,], nProps=nProps, sigma = sigma)
    currentAndProps <- rbind(states[i-1,],props)
    targetVals      <- apply(currentAndProps,MARGIN = 1,FUN=target)
    selectionProbs  <- targetVals
    lambdas         <- -selectionProbs + log(rexp(n=nProps+1))
    currentAndProps <- currentAndProps[order(lambdas),]
    rank1           <- rank(lambdas)[1]
    qmin            <- quantumMin(field=1:(nProps+1),y=rank1)
    propIndex       <- qmin[[1]]
    oracleCalls     <- oracleCalls + qmin[[2]]
    states[i,]      <- currentAndProps[propIndex,]
    SampCount       <- SampCount + 1
    
    if(any(states[i,] != states[i-1,])) {
      Acceptances <- Acceptances + 1
    }    
    if(i %% 1000 == 0) cat("Iteration ", i, " complete.\n")
    if (SampCount == SampBound) { 
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        sigma <- sigma * min(2,(1 + delta(batches)))
      } else {
        sigma <- sigma * max(0.5,(1 - delta(batches)))
      }
      batches <- batches + 1
      cat("Acceptance ratio: ", AcceptRatio, ". sigma: ", sigma, "\n")
      
      SampCount <- 0
      Acceptances <- 0
    }
    
    if (earlyStop & i>101 & i %% 1000 == 0) {
      ess0 <- min(coda::effectiveSize(states[101:i,]))
      cat(ess0,"\n")
      if(ess0 >100) {
        return(list(i,oracleCalls))
      }
    }
  }
  return(list(states,oracleCalls))
}

randomWalk <- function(N, x0, maxIt=10000,earlyStop=FALSE) {
  if(N!=length(x0)) stop("Dimension mismatch.")
  
  chain <- matrix(0,maxIt,N)
  
  sigma <- 2.4/sqrt(N)
  Ct <- sigma^2*diag(N) 
  xbar <- x0
  
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    xStar <- rnorm(N,sd=sigma) + chain[i-1,]
    if(log(runif(1)) < target(xStar) -
       target(chain[i-1,])) {
      accept[i] <- 1
      chain[i,] <- xStar
    } else {
      chain[i,] <- chain[i-1,]
    }
    
    
    if(i %% 1000 == 0) cat(i,"\n")
    if (earlyStop & i>101 & i %% 1000 == 0) {
      ess0 <- min(coda::effectiveSize(chain[101:i,]))
      cat(ess0,"\n")
      if(ess0 > 100) {
        return(list(i,i))
      }
    }
  }
  ratio <- sum(accept)/(maxIt-1)
  cat("Acceptance ratio: ", ratio,"\n")
  return(list(chain,ratio,sigma,diag(N)))
}

################################################################################
library(coda)
seed <- R.utils::cmdArg("seed",42L)
set.seed(seed=seed)

for(props in c(4,8,16,32,64,128) ) {
  pMCMC_Output <- pMCMC(nIts=10000,nProps=props, initial=0)
  qpMCMC_Output <- qpMCMC(nIts = 10000,nProps=props,initial=0)
  
  meanP  <- mean(effectiveSize(pMCMC_Output[[1]]))
  meanQP <- mean(effectiveSize(qpMCMC_Output[[1]]))
  
  relErrorMean <- abs(meanP-meanQP)/meanP
  
  minP  <- min(effectiveSize(pMCMC_Output[[1]]))
  minQP <- min(effectiveSize(qpMCMC_Output[[1]]))
  
  relErrorMin <- abs(minP-minQP)/minP
  
  cat(props,relErrorMean,relErrorMin,"\n",file="data/pMCMCvsQPMCMC.txt",
      append=TRUE)
}

