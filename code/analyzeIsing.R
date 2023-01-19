setwd("~/qpMCMC/")
library(ggplot2)
library(ggpmisc)

df <- readr::read_table2("~/qpMCMC/data/Ising2dLogProbs.txt", col_names = FALSE)

df <- df[,1:2]

# temporary
#df <- df[df$X1 %in% c(4,8,16,32,64,128,256,512,1024),]


df$Proposals <- factor(df$X1)
df$`Parallel\nMCMC\nproposals` <- forcats::fct_rev(df$Proposals)
df$LogProbs  <- df$X2
df$Iterations <- rep(seq.int(from=1,to=10000000,by=1000),10)
df <- df[,3:6]

pal <- wesanderson::wes_palette("Zissou1",10,type="continuous")

df2 <- readr::read_table2("~/qpMCMC/data/Ising2D.txt", col_names = FALSE)
df2 <- df2[,c(1,3)]
colnames(df2) <- c("Proposals","Target evaluations")
df2$Speedup   <- (df2$Proposals*10000000)/df2$`Target evaluations`
df2$Speedup   <- round(df2$Speedup,digits = 2)
df2$`Target evaluations` <- signif(df2$`Target evaluations`,3)
df2$`Target evaluations` <- formatC(df2$`Target evaluations`,format="e",digits=2)


gg <- ggplot(df,aes(x=Iterations,y=LogProbs,color=`Parallel\nMCMC\nproposals`)) +
  geom_line() +
  scale_color_manual(values = pal) +
  ylab("Unnormalized log-probabilities") + xlab("MCMC iterations") +
  ggtitle("Convergence for an Ising model on a 500-by-500 lattice") +
  geom_table(data = df2, 
                 label = list(df2),
                 x = 9000000, y = 400000,
                 table.theme = ttheme_gtbw,size=3) +
  theme_bw()
gg

ggsave("~/qpMCMC/figures/Ising2dFig.pdf",device = "pdf",width = 7,height=4)
system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/Ising2dFig.pdf",
                    "~/qpMCMC/figures/Ising2dFig")
)



