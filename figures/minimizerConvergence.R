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

expoSearch <- function(N=2^10,M=1,marks=NULL) {
  success <- FALSE
  m <- 1
  grovIter <- 0
  lmbd <- 6/5
  if(is.null(marks)) {
    marks <- c(rep(1,M),rep(0,N-M))
  } else {
    M <- sum(marks==1)
    N <- length(marks)
  }
  while(!success & grovIter<(9*sqrt(N)/4)) {
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
  # multiproposal from normal distribution centered at current state
  dimensions <- length(currentState)
  props <- matrix(currentState,nProps,dimensions,byrow=TRUE)
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

qpMCMC <- function (nIts=1000,nProps=2000,dims=10,sigma=3,targetAccept=0.5) {
  
  states <- matrix(0,nIts,dims)
  states[1,] <- 100
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
  return(list(states,oracleCalls))
}


################################################################################
set.seed(1)
N <- 2^c(10:14)
df <- data.frame()
for(n in N) {
  for(M in 5:1) {
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

library(ggplot2)
library(wesanderson)
library(forcats)
pal <- wes_palette("Zissou1", 5, type = "discrete")

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
saveRDS(df,"qpMCMC/mcmcIsLessThan.rds")
df <- readRDS("qpMCMC/mcmcIsLessThan.rds")

colnames(df) <- c("Iteration","Dimension","Rank","GroverIts")
df$`Target\ndimension` <- factor(df$Dimension)

gg <- ggplot(df,aes(x=Iteration,y=GroverIts,color=`Target\ndimension`)) +
  geom_line(alpha=0.1) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 5), se = FALSE)+
  scale_color_manual(values = pal) +
  ggtitle("Evaluations and burn-in for quantum parallel MCMC") +
  ylab("Oracle evaluations") +
  xlab("MCMC iteration") +
  theme_bw()
gg

mean(df$Rank!=1) # 0.0054
sum(df$GroverIts)/(2000*2000*5)#  0.07212715

ggsave(gg,filename = "mcmcIterations.pdf",device = "pdf",path = "qpMCMC/figures/",dpi = "retina",
       width=8.14,height=4)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/mcmcIterations.pdf",
                    "~/qpMCMC/figures/mcmcIterations.pdf")
)

#################################################################################
set.seed(1)
out <- qpMCMC(nIts = 100000,nProps = 2000, dims = 100)
out[[2]]/(100000*2000) # 0.07
#saveRDS(out,file="qpMCMC/qqplotExample.rds")
out <- readRDS(file="qpMCMC/qqplotExample.rds")
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

ggsave(gg,filename = "qqPlot.png",device = "png",path = "qpMCMC/figures/",dpi = "retina",
       width=5,height=12)
