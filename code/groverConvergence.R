
selectionProbs <- function(N=2^10,M=1,Iter=NULL) {
  if(is.null(Iter)){
    Iter  <- ceiling(pi*sqrt(N/M)/4)
  }
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

################################################################################

N <- 2^14
M <- c(256,128,64,32,16,12,8:1)
probSucceed <- selectionProbs(N=N,M=M[1])
df <- data.frame(Iterations=1:(ceiling(pi*sqrt(N/M[1])/4)),M=M[1],SuccessProb=probSucceed)

for(i in 2:length(M)) {
  probSucceed <- selectionProbs(N=N,M=M[i])
  dfTemp <- data.frame(Iterations=1:(ceiling(pi*sqrt(N/M[i])/4)),M=M[i],SuccessProb=probSucceed)
  df <- rbind(df,dfTemp)
}

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

ggsave(gg,filename = "groverCurves.pdf",device = "pdf",path = "qpMCMC/figures/",dpi = "retina",
       width=8.14,height=4)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/groverCurves.pdf",
                    "~/qpMCMC/figures/groverCurves.pdf")
)

################################################################################

N <- 2^14
M <- c(256,128,64,32,16,12,8:1)
probSucceed <- selectionProbs(N=N,M=M[1],Iter=100)
df <- data.frame(Iterations=1:100,M=M[1],SuccessProb=probSucceed)

for(i in 2:length(M)) {
  probSucceed <- selectionProbs(N=N,M=M[i],Iter=100)
  dfTemp <- data.frame(Iterations=1:100,M=M[i],SuccessProb=probSucceed)
  df <- rbind(df,dfTemp)
}

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

ggsave(gg,filename = "groverCurvesDeath.pdf",device = "pdf",path = "qpMCMC/figures/",dpi = "retina",
       width=8.14,height=4)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/groverCurvesDeath.pdf",
                    "~/qpMCMC/figures/groverCurvesDeath.pdf")
)