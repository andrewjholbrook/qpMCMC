setwd("~/qpMCMC/")
library(ggplot2)

df <- readr::read_table2("~/qpMCMC/Ising2dLogProbs.txt", col_names = FALSE)

df <- df[,1:2]

# temporary
df <- df[df$X1 %in% c(4,8,16,32,64,128,256,512,1024),]


df$Proposals <- factor(df$X1)
df$`Parallel\nMCMC\nproposals` <- forcats::fct_rev(df$Proposals)
df$LogProbs  <- df$X2
df$Iterations <- rep(seq.int(from=1,to=10000000,by=1000),9)
df <- df[,3:6]

pal <- wesanderson::wes_palette("Zissou1",9,type="continuous")

gg <- ggplot(df,aes(x=Iterations,y=LogProbs,color=`Parallel\nMCMC\nproposals`)) +
  geom_line() +
  scale_color_manual(values = pal) +
  ylab("Unnormalized log-probabilities") + xlab("MCMC iterations") +
  ggtitle("Convergence for a 2D Ising model on 500-by-500 lattice") +
  theme_bw()
gg

ggsave("~/qpMCMC/figures/Ising2dFig.pdf",device = "pdf",width = 7,height=4)
system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/Ising2dFig.pdf",
                    "~/qpMCMC/figures/Ising2dFig")
)



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


