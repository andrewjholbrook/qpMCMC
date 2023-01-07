setwd("~/qpMCMC")

library(readr)
library(coda)

df <- read_table2("blackHoleResults.txt",col_names = FALSE)
df <- as.matrix(df[,1:7])
# for(i in 1:7){
#   df[,i] <- unlist(df[,i])
# }
#effectiveSize(df[4000:8000,4])
plot(df[,6],type="l")

df <- read_table2("blackHoleResultsSinglePrec.txt",col_names = FALSE)
df <- as.matrix(df[,1:6])
# for(i in 1:7){
#   df[,i] <- unlist(df[,i])
# }
#effectiveSize(df[4000:8000,4])
plot(df[,5],type="l")

# results <- readRDS("~/qpMCMC/blackHole.rds")
# 
# 1024*20000000/results[[4]]
# plot(results[[8]],type="l")
# 
# imgLong <- melt(results[[1]][[1000]])
# remove(results)
# ggplot(imgLong, aes(x = Var2, y = Var1)) +
#   geom_raster(aes(fill=value))
