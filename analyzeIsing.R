setwd("~/qpMCMC/")
library(ggplot2)

df <- readr::read_table2("~/qpMCMC/Ising2dLogProbs.txt", col_names = FALSE)

df <- df[,1:2]

df$Proposals <- factor(df$X1)
df$LogProbs  <- df$X2
df$Iterations <- rep(seq.int(from=1,to=1000000),10)
df <- df[,3:5]
df <- df[df$Iterations %in% seq.int(from=1,to=1000000,by=1000),]

gg <- ggplot(df,aes(x=Iterations,y=LogProbs,color=Proposals)) +
  geom_line()
gg

gg2 <- ggplot(df,aes(x=Iterations,y=Distance,color=Proposals)) +
  geom_line() +
  ylim(c(124000,126000)) + xlim(c(200000,1000000))
gg2

coda::effectiveSize(df$Distance[df$Proposals==4 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==8 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==16 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==32 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==64 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==128 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==256 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==512 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==1024 & df$Iterations>=200000])
coda::effectiveSize(df$Distance[df$Proposals==2048 & df$Iterations>=200000])


