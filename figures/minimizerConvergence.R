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

bubbleBath <- function (nIts=100,nProps=100,dims=2,sigma=1) {
  
  states <- matrix(0,nIts,dims)
  states[1,] <- 3
  numberGreater <- rep(0,nIts)
  
  for(i in 2:nIts) {
    props           <- proposal(states[i-1,], nProps=nProps, sigma = sigma)
    currentAndProps <- rbind(states[i-1,],props)
    targetVals      <- apply(currentAndProps,MARGIN = 1,FUN=target)
    prpPrbs         <- propProbsProd(currentAndProps, sigma = sigma)
    selectionProbs  <- targetVals + prpPrbs
    lambdas         <- selectionProbs - log(rexp(n=nProps+1))
    numberGreater[i]<- sum(lambdas > lambdas[1])
    propIndex       <- which(lambdas==max(lambdas))  #sample(1:length(selectionProbs),size=1,prob=selectionProbs)
    states[i,]      <- currentAndProps[propIndex,]
    cat(i,"\n")
  }
  return(list(states,numberGreater))
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
out <- bubbleBath(nIts=1000,nProps = 2000)
N <- 2^16
#M <- c(256,128,64,32,16,12,8:1)
probSucceed <- groverProbs(N=N,M=M[1])
# df <- data.frame(Iterations=1:(ceiling(pi*sqrt(N/M[1])/4)),M=M[1],SuccessProb=probSucceed)
# 
# for(i in 2:length(M)) {
#   probSucceed <- selectionProbs(N=N,M=M[i])
#   dfTemp <- data.frame(Iterations=1:(ceiling(pi*sqrt(N/M[i])/4)),M=M[i],SuccessProb=probSucceed)
#   df <- rbind(df,dfTemp)
# }

initialState <- c(10,10)


library(ggplot2)
library(wesanderson)
library(forcats)
pal <- wes_palette("Zissou1", length(M), type = "continuous")

df$Solutions <- fct_rev(factor(df$M,ordered = TRUE))

gg <- ggplot(df,aes(x=Iterations,y=SuccessProb,color=Solutions)) +
  geom_line() +
  scale_color_manual(values = pal) +
  ylab("Probability of success") +
  xlab("Oracle evaluations") +
  ggtitle("Grover search over ~16k items") +
  theme_bw()

gg

ggsave(gg,filename = "groverCurves2.pdf",device = "pdf",path = "qpMCMC/figures/",dpi = "retina",
       width=8.14,height=4)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/groverCurves2.pdf",
                    "~/qpMCMC/figures/groverCurves2.pdf")
)

