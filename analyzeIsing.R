setwd("~/qpMCMC/")
library(ggplot2)

df <- read_table2("~/qpMCMC/Ising2dDistances.txt", col_names = FALSE)

df <- df[,1:2]

df$Proposals <- factor(df$X1)
df$Distance  <- df$X2
df$Iterations <- rep(seq.int(from=0,to=999000,by=1000),10)
df <- df[,3:5]

gg <- ggplot(df,aes(x=Iterations,y=Distance,color=Proposals)) +
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


