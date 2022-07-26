library(ggplot2)
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
  return(probSucceed)
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

fixedPointSearch <- function(N=2^10,M=1,marks=NULL,delta=sqrt(0.995)) {
  # 2D representation
  grovIter <- 0
  if(is.null(marks)) {
    marks <- c(rep(1,M),rep(0,N-M))
  } else {
    M <- sum(marks==1)
    N <- length(marks)
  }
  w <- 1/N # width
  L <- ceiling( log(2/delta)/sqrt(w) )
  if (L %% 2 == 0) L <- L + 1
  
  
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

pMCMC <- function (nIts=1000,nProps=2000,dims=10,sigma=3,targetAccept=0.5) {
  
  states <- matrix(0,nIts,dims)
  states[1,] <- 100
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
        sigma <- sigma * (1 + delta(batches))
      } else {
        sigma <- sigma * (1 - delta(batches))
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
    browser()
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
        sigma <- sigma * (1 + delta(batches))
      } else {
        sigma <- sigma * (1 - delta(batches))
      }
      batches <- batches + 1
    #  cat("Acceptance ratio: ", AcceptRatio, ". sigma: ", sigma, "\n")
      
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
set.seed(1)
N <- 2^c(10:14)
df <- data.frame()
for(n in N) {
  for(M in c(20,15,10,5,1)) {
    for (k in 1:500) {
      df <- rbind(df,c(n,M,expoSearch(N=n,M=M)[[1]]))
      if(k %% 100==0)   cat(k,"\n")
    }
    cat(M,"\n")
  }  
}


colnames(df) <- c("N","Solutions","GroverIts")
df$Solutions <- factor(df$Solutions)
df$N <- factor(df$N)

df$N <- fct_rev(df$N)

gg <- ggplot(df,aes(x=Solutions,y=GroverIts,fill=N)) +
  geom_hline(aes(yintercept = 9/4*sqrt(2^(9+6-as.numeric(N))), colour=N),size=1.5,alpha=0.7) +
  geom_boxplot() +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  ylab("Oracle evaluations")+
  xlab("Solutions (M)") +
  ggtitle("Quantum exponential searching algorithm") +
  theme_bw()
gg

ggsave(gg,filename = "expoSearch.pdf",device = "pdf",path = "qpMCMC/figures/",dpi = "retina",
       width=8.14,height=4)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/expoSearch.pdf",
                    "~/qpMCMC/figures/expoSearch.pdf")
)
################################################################################
set.seed(1)
N <- 2^c(10:14)
df <- data.frame()
for(n in N) {
  for (k in 1:500) {
    for(y in 1:5) {
      qmin <- quantumMin(field=1:n,y=y)
      df <- rbind(df,c(n,y,qmin[[1]],qmin[[2]]))
    }
    if(k %% 100==0)   cat(k,"\n")
  }  
}
# mean(df[,3]!=1) = 0.00928
saveRDS(df,file="qpMCMC/quantMinExp.rds")

colnames(df) <- c("N","StartingRank","Success","GroverIts")
df$N <- factor(df$N)
df$N <- fct_rev(df$N)

df$StartingRank <- factor(df$StartingRank)
gg <- ggplot(df,aes(x=StartingRank,y=GroverIts,fill=N)) +
  geom_boxplot() +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  ylab("Oracle evaluations")+
  xlab("Starting rank") +
  ggtitle("Quantum minimization algorithm") +
  theme_bw()
gg

ggsave(gg,filename = "qMinAlg.pdf",device = "pdf",path = "qpMCMC/figures/",dpi = "retina",
       width=8.14,height=4)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/qMinAlg.pdf",
                    "~/qpMCMC/figures/qMinAlg.pdf")
)
################################################################################
set.seed(1)
D <- c(150,300,600,1200,2400)
df  <- data.frame()
for(d in D) {
  out <- pMCMC(nIts=2000,nProps = 2000, dims = D)
  for(i in 2:2000){
    qmin <- quantumMin(field=1:2000,y=2001-out[[2]][i])
    df <- rbind(df,c(i,d,qmin[[1]],qmin[[2]]))
  }
  if(i %% 100==0)   cat(i,"\n")
}
saveRDS(df,"~/qpMCMC/mcmcIsLessThan.rds")
df <- readRDS("~/qpMCMC/mcmcIsLessThan.rds")

colnames(df) <- c("Iteration","Dimension","Rank","GroverIts")
df$`Target dimension` <- factor(df$Dimension)

gg <- ggplot(df,aes(x=Iteration,y=GroverIts,color=`Target dimension`)) +
  geom_hline(yintercept = 2000,linetype=2) +
  annotate("text", x = 400, y = 1950, label = "Conventional implementation") +
  geom_line(alpha=0.1) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 5), se = FALSE)+
  scale_color_manual(values = pal) +
  ggtitle("Evaluations and burn-in for quantum parallel MCMC") +
  ylab("Oracle evaluations") +
  xlab("MCMC iteration") +
  theme_bw()
gg

gg2 <- ggplot(df,aes(x=Iteration,y=GroverIts,color=`Target dimension`)) +
  geom_line(alpha=0.1) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 5), se = FALSE)+
  scale_color_manual(values = pal) +
  ggtitle("") +
  ylab(NULL) +
  xlab("MCMC iteration") +
  theme_bw()
gg2

mean(df$Rank!=1) # 0.0053
sum(df$GroverIts)/(2000*2000*5)#  0.074

library(ggpubr)

ggsave(ggarrange(gg, gg2, ncol = 2, nrow = 1, common.legend = TRUE, 
                 legend = "bottom"),filename = "mcmcIterations.pdf",device = "pdf",path = "~/qpMCMC/figures/",dpi = "retina",
       width=12,height=5)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/mcmcIterations.pdf",
                    "~/qpMCMC/figures/mcmcIterations.pdf")
)

#################################################################################
set.seed(1)
out <- qpMCMC(nIts = 100000,nProps = 2000, dims = 100)
out[[2]]/(100000*2000) # 0.07
saveRDS(out,file="~/qpMCMC/qqplotExample.rds")
out <- readRDS(file="~/qpMCMC/qqplotExample.rds")
out[[1]] <- out[[1]][floor(seq(from=2000,to=100000,by=10)),]

# qqnorm(out[[1]][2000:100000,3])
# qqline(out[[1]][2000:100000,3],col="red")
y <- c()
for (d in 1:100) {
  out[[1]][,d] <- out[[1]][order(out[[1]][,d]),d]
  ytemp <- qnorm(p=((1:9801-0.5)/9801)) + d
  #ytemp <- ytemp[order(ytemp)]
  y <- c(y,ytemp)
}

pal <- wes_palette("Zissou1", 5, type = "discrete")

df1 <- data.frame(Samples=as.vector(out[[1]]), DirectSamples=y,
                  Dimension=rep(1:100,each=9801))
df1$Dimension <- factor(df1$Dimension)

gg <- ggplot(df1,aes(x=Samples,y=DirectSamples,color=Dimension),alpha=0.7) +
  geom_abline(slope=1,intercept = 1:100,color="grey")  + 
  geom_point() +
  scale_color_manual(values=rep(c(pal[1],pal[3],pal[5],pal[2],pal[4]),20)) +
  xlim(c(-5,5)) + ylim(c(-5,105)) +
  ylab("Theoretical quantiles plus target dimension") + xlab("Sample quantiles") + 
  ggtitle("QQ plot for 100D Gaussian target") +
  theme_bw() +
  theme(legend.position="none")
gg

ggsave(gg,filename = "qqPlot.png",device = "png",path = "~/qpMCMC/figures/",dpi = "retina",
       width=5,height=12)


################################################################################
# ess per evaluation 



target <- function(X,B=0.1) { 
  # B is bananicity constant
  N <- length(X)
  output <-  - X[2]^2/200 - 0.5*(X[1]+B*X[2]^2-100*B)^2 - 0.5*sum(X[3:N]^2)
  return(output)
}


set.seed(1)
nProps <- seq(from=100,to=1000,by=100)#,12800,25600)
df  <- data.frame()
nIts <- 10000
#for (k in 1:10) {
for(d in nProps) {
  out <- qpMCMC(nIts=nIts,nProps = d, dims = 20,initial = 0, sigma=0.5)
  df <- rbind(df,c(d,min(coda::effectiveSize(out[[1]][1001:nIts,])),out[[2]]))
}

rwm <- randomWalk(N=20,x0=rep(0,20),maxIt = nIts)
df <- rbind(df,c(1,min(coda::effectiveSize(rwm[[1]][1001:nIts,])),nIts))
#}
df[,3] <- df[,3]*((nIts-1000)/nIts)
df[,4] <- df[,3]/df[,2]
colnames(df) <- c("Proposals", "ESS","Evaluations","EvalsPerESS")

saveRDS(df,"~/qpMCMC/essPerEvalBanana.rds")
df <- readRDS("~/qpMCMC/essPerEvalBanana.rds")

df$Proposals <- factor(df$Proposals)

df2 <- by(data=df$ESS,INDICES=list(df$Proposals),
                        FUN=quantile,probs=c(0.5,0.05,0.95))
df3 <- data.frame()
for(i in 1:length(df2)) df3 <- rbind(df3,df2[[i]])

df2 <- by(data=df$EvalsPerESS,INDICES=list(df$Proposals),
          FUN=quantile,probs=c(0.5,0.05,0.95))
df4 <- data.frame()
for(i in 1:length(df2)) df4 <- rbind(df4,df2[[i]])

df5 <- cbind(df3,df4)

df6 <- data.frame()
for(i in 1:11) df6 <- rbind(df6, df5[1,4:6]/df5[i,4:6])
df <- cbind(df5,df6)

print(xtable::xtable(df),booktabs=TRUE)

##############################


# target <- function(X,B=0.1) { 
#   # B is bananicity constant
#   N <- length(X)
#   output <-  - X[2]^2/200 - 0.5*(X[1]+B*X[2]^2-100*B)^2 #- 0.5*sum(X[3:N]^2)
#   return(output)
# }
# 
# set.seed(1)
# nProps <- seq(from=10,to=100,by=10)#,12800,25600)
# df  <- data.frame()
# nIts <- 1000000
# for (k in 1:20) {
# for(d in nProps) {
#   out <- qpMCMC(nIts=nIts,nProps = d, dims = 2,initial = 0, sigma=0.5,earlyStop=TRUE)
#   df <- rbind(df,c(d,out[[1]],out[[2]]))
# }
# 
# rwm <- randomWalk(N=2,x0=rep(0,2),maxIt = 1000000,earlyStop = TRUE)
# df <- rbind(df,c(1,rwm[[1]],rwm[[2]]))
# }
# 
# 
# colnames(df) <- c("Proposals", "Iterations","Evaluations")
# 
# saveRDS(df,"~/qpMCMC/essPerEvalBananaTimeTo.rds")
# df <- readRDS("~/qpMCMC/essPerEvalBananaTimeTo.rds")
# 
# df$Proposals <- as.factor(df$Proposals)
# by(data=df$Evaluations,INDICES=list(df$Proposals),FUN=quantile,probs=c(0.5,0.05,0.95))

# target <- function(X,B=0.1) {
#   # B is bananicity constant
#     output <- log( exp( - X[2]^2/200 - 0.5*(X[1]+B*X[2]^2-100*B)^2) +
#                      exp( - (X[1]-1000)^2/200 - 0.5*(X[2]+500+B*(X[1]-1000)^2-100*B)^2) +
#                      exp( - (X[2])^2/200 - 0.5*(X[1]-2000-B*(X[2])^2-100*B)^2 ) +
#                      exp( - (X[1]-500)^2/200 - 0.5*(X[2]-500 -B*(X[1]-500)^2-100*B)^2 ) +
#                      exp( - (X[1]-1500)^2/200 - 0.5*(X[2]-500 -B*(X[1]-1500)^2-100*B)^2 ) )
# 
#   return(output)
# }
# 
# set.seed(1)
# nProps <- seq(from=2000,to=10000,by=2000)#,12800,25600)
# df  <- data.frame()
# nIts <- 100000
#   for(d in nProps) {
#     out <- qpMCMC(nIts=nIts,nProps = d, dims = 2,initial = 0, sigma=0.5)
#     df <- rbind(df,c(d,mean(coda::effectiveSize(out[[1]][1001:nIts,])),out[[2]]))
# }
# df[,3] <- df[,3]*((nIts-1000)/nIts)
# df[,4] <- df[,3]/df[,2]
# colnames(df) <- c("Proposals", "ESS","Evaluations","EvalsPerESS")
# 
# saveRDS(df,"~/qpMCMC/essPerEvalMob.rds")
# df <- readRDS("~/qpMCMC/essPerEvalMob.rds")
# 
# df$Proposals <- factor(df$Proposals)
# 
# df2 <- by(data=df$ESS,INDICES=list(df$Proposals),
#           FUN=quantile,probs=c(0.5,0.05,0.95))
# df3 <- data.frame()
# for(i in 1:length(df2)) df3 <- rbind(df3,df2[[i]])
# 
# df2 <- by(data=df$EvalsPerESS,INDICES=list(df$Proposals),
#           FUN=quantile,probs=c(0.5,0.05,0.95))
# df4 <- data.frame()
# for(i in 1:length(df2)) df4 <- rbind(df4,df2[[i]])
# 
# df5 <- cbind(df3,df4)
# 
# df6 <- data.frame()
# for(i in 1:11) df6 <- rbind(df6, df5[1,4:6]/df5[i,4:6])
# df <- cbind(df5,df6)
# 

# get number of evals
set.seed(1)
df <- data.frame("Proposals","Iterations","Evaluations","Speedup")
nProps <- 10000
for(k in 0:4){
  ranks <- read_csv(paste0("qpMCMC/qpMCMC10KRanks",k,".txt"), 
                    col_names = FALSE)
  ranks <- unlist(ranks)
  ranks <- ranks + 1
  ranks <- as.numeric(ranks)
  
  oracleCalls <- 0
  L <- length(ranks) # 5998
  for(i in 1:L) {
    qmin <- quantumMin(1:(nProps+1),ranks[i])
    oracleCalls <- oracleCalls + qmin[[2]]
    if(i%%100==0) cat(i,"\n")
  }
  df <- rbind(df,c(nProps, L,oracleCalls, (L*nProps)/oracleCalls))
}


# get number of evals
set.seed(1)

nProps <- 5000
for(k in 0:4){
  ranks <- read_csv(paste0("qpMCMC/qpMCMC5KRanks",k,".txt"), 
                    col_names = FALSE)
  ranks <- unlist(ranks)
  ranks <- ranks + 1
  ranks <- as.numeric(ranks)
  
  oracleCalls <- 0
  L <- length(ranks) # 5998
  for(i in 1:L) {
    qmin <- quantumMin(1:(nProps+1),ranks[i])
    oracleCalls <- oracleCalls + qmin[[2]]
    if(i%%100==0) cat(i,"\n")
  }
  df <- rbind(df,c(nProps, L,oracleCalls, (L*nProps)/oracleCalls))
}


set.seed(1)

nProps <- 1000
for(k in 0:4){
  ranks <- read_csv(paste0("qpMCMC/qpMCMC1KRanks",k,".txt"), 
                    col_names = FALSE)
  ranks <- unlist(ranks)
  ranks <- ranks + 1
  ranks <- as.numeric(ranks)
  
  oracleCalls <- 0
  L <- length(ranks) # 5998
  for(i in 1:L) {
    qmin <- quantumMin(1:(nProps+1),ranks[i])
    oracleCalls <- oracleCalls + qmin[[2]]
    if(i%%100==0) cat(i,"\n")
  }
  df <- rbind(df,c(nProps, L,oracleCalls, (L*nProps)/oracleCalls))
}


saveRDS(df,"evalsTillResults.rds")

colnames(df) <- c("Proposals","Iterations","Evaluations","Speedup")
df <- df[-1,]
df[,1] <- as.numeric(df[,1])
df[,2] <- as.numeric(df[,2])
df[,3] <- as.numeric(df[,3])
df[,4] <- as.numeric(df[,4])

df$Gain <- df$Evaluations[11:15] / df$Evaluations
#df$MeanSpeedup <- df$Speedup

df2 <- rbind(colMeans(df[df$Proposals==10000,]),colMeans(df[df$Proposals==5000,]),colMeans(df[df$Proposals==1000,]))
df3 <- rbind(apply(df[df$Proposals==10000,],2,min),apply(df[df$Proposals==5000,],2,min),apply(df[df$Proposals==1000,],2,min))

df4 <- rbind(apply(df[df$Proposals==10000,],2,max),apply(df[df$Proposals==5000,],2,max),apply(df[df$Proposals==1000,],2,max))

colnames(df3) <- c("Proposals","MinIts","MinEvals","MinSpeedup","MinGain")
colnames(df4) <- c("Proposals","MaxIts","MaxEvals","MaxSpeedup","MaxGain")
df5 <- merge(df2,df3)
df6 <- merge(df5,df4)
df6 <- df6[,c(1,2,6,10,3,7,11,4,8,12,5,9,13)]
print(xtable::xtable(df6),booktabs=TRUE)
