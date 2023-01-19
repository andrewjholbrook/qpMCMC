setwd("~/qpMCMC")

library(readr)
library(coda)
library(ggplot2)
library(reshape2)
# 
df <- read_table2("data/blackHoleResults.txt",col_names = FALSE)
df <- as.matrix(df[,1:6])

effectiveSize(df[2501:5000,]) # remove first half as burning
# X1          X2          X3          X4 
# 0.000000  257.077281 1578.662332  257.608159 
# X5          X6 
# 2500.000000    1.500075 


#plot(df[2501:5000,5],type="l")


1024*20000000/df[5000,6] # 10.36 speedup






results <- readRDS("~/qpMCMC/data/blackHoleResults2.rds")

# get mode
imgLong <- melt(results[[1]])
for(i in 2:26) {
  imgLong <- imgLong + melt(results[[i]])
}
imgLong <- round(imgLong / 26)
remove(results)
gg2 <- ggplot(imgLong, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value),show.legend = FALSE) +
  scale_fill_gradient(low = "black",high="white") +
  ggtitle("Pixelwise posterior mode") +
  scale_x_continuous(expand = c(0,0) ) +
  scale_y_continuous(expand = c(0,0) ) +
  theme_void()
gg2
remove(imgLong)

img <- readRDS("~/qpMCMC/data/blackHoleIntensity.rds")
imgLong <- melt(img)
remove(img)
gg1 <- ggplot(imgLong, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value),show.legend = FALSE) +
  scale_fill_gradient(low = "black",high="white") +
  ggtitle("Intensities from Sagittarius A* image") +
  scale_x_continuous(expand = c(0,0) ) +
  scale_y_continuous(expand = c(0,0) ) +
  theme_void()
gg1

library(ggpubr)

ggsave(ggarrange(gg1, gg2, ncol = 2, nrow = 1),filename = "blackHole.pdf",device = "pdf",path = "~/qpMCMC/figures/",dpi = "retina",
       width=9,height=3)

system2(command = "pdfcrop",
        args    = c("~/qpMCMC/figures/blackHole.pdf",
                    "~/qpMCMC/figures/blackHole.pdf")
)
