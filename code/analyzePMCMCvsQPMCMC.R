setwd("~/qpMCMC/")

library(forcats)
library(ggplot2)
library(reshape2)
library(wesanderson)

df <- readr::read_table2("data/pMCMCvsQPMCMC.txt", col_names = FALSE)
df <- df[,1:3]
colnames(df) <- c("Proposals","Mean","Minimum")

df <- melt(df,measure.vars = 2:3)
colnames(df) <- c("Proposals","Statistic","Relative error")
df$Statistic <- factor(df$Statistic)
df$Proposals <- factor(df$Proposals)

pal <- wesanderson::wes_palette("Zissou1",10,type="continuous")


gg <- ggplot(df,aes(y=`Relative error`,x=Proposals,fill=Statistic)) +
  geom_boxplot(position=position_dodge()) +
  scale_fill_manual(values=c(pal[1],pal[10])) +
  ggtitle("Comparing ESS for pMCMC and QPMCMC") +
  ylab("Relative differences") +
  theme_bw()
  
gg

ggsave(gg,filename = "figures/pMCMCvsQPMCMC.pdf",device = "pdf",
       width = 6,height=3)
system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/pMCMCvsQPMCMC.pdf",
                    "~/qpMCMC/figures/pMCMCvsQPMCMC.pdf")
)