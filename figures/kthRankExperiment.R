library(ggplot2)

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

# biased to take kth rank
pMCMC <- function (nIts=1000,nProps=2000,dims=10,sigma=30,acceptQuant=0.5) {
  
  states <- matrix(0,nIts,dims)
  states[1,] <- 0
  
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)

  for(i in 2:nIts) {
    props           <- proposal(states[i-1,], nProps=nProps, sigma = sigma)
    currentAndProps <- rbind(states[i-1,],props)
    targetVals      <- apply(currentAndProps,MARGIN = 1,FUN=target)
    selectionProbs  <- targetVals
    lambdas         <- selectionProbs - log(rexp(n=nProps+1))
    propIndex       <- which.min(abs(lambdas-quantile(lambdas,probs=acceptQuant)))
    states[i,]      <- currentAndProps[propIndex,]

    if(propIndex != 1) {
      Acceptances <- Acceptances + 1
    }    
    if(i %% 100 == 0) cat("Iteration ", i, " complete.\n")
  
  }
  return(states)
}


out <- pMCMC(nIts=10000,nProps=2000,dims=2,acceptQuant = 1,sigma=30)
out <- as.data.frame(out)
colnames(out) <- c("X","Y")

gg <- ggplot(out,aes(x=X,y=Y)) +
  geom_point() + ggtitle("Quantile = 1") + theme_bw()
gg

out <- pMCMC(nIts=10000,nProps=2000,dims=2,acceptQuant = 0.9,sigma=30)
out <- as.data.frame(out)
colnames(out) <- c("X","Y")

gg2 <- ggplot(out,aes(x=X,y=Y)) +
  geom_point() + ggtitle("Quantile = 0.9") + theme_bw()
gg2

out <- pMCMC(nIts=10000,nProps=2000,dims=2,acceptQuant = 0.8,sigma=30)
out <- as.data.frame(out)
colnames(out) <- c("X","Y")

gg3 <- ggplot(out,aes(x=X,y=Y)) +
  geom_point() + ggtitle("Quantile = 0.8") + theme_bw()
gg3

out <- pMCMC(nIts=10000,nProps=2000,dims=2,acceptQuant = 0.7,sigma=30)
out <- as.data.frame(out)
colnames(out) <- c("X","Y")

gg4 <- ggplot(out,aes(x=X,y=Y)) +
  geom_point() + ggtitle("Quantile = 0.7") + theme_bw()
gg4

out <- pMCMC(nIts=10000,nProps=2000,dims=2,acceptQuant = 0.6,sigma=30)
out <- as.data.frame(out)
colnames(out) <- c("X","Y")

gg5 <- ggplot(out,aes(x=X,y=Y)) +
  geom_point() + ggtitle("Quantile = 0.6") + theme_bw()
gg5

out <- pMCMC(nIts=10000,nProps=2000,dims=2,acceptQuant = 0.5,sigma=30)
out <- as.data.frame(out)
colnames(out) <- c("X","Y")

gg6 <- ggplot(out,aes(x=X,y=Y)) +
  geom_point() + ggtitle("Quantile = 0.5") + theme_bw()
gg6

library(ggpubr)

ggsave(ggarrange(gg,gg2,gg3,gg4,gg5,gg6,
                 ncol=2,nrow=3),
       filename = "~/qpMCMC/figures/acceptanceQuantiles.pdf",
       device = "pdf",
       width=12,height=12)
