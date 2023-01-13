setwd("~/qpMCMC")

library(readr)
library(coda)
# 
df <- read_table2("blackHoleResultsCenteredPrior.txt",col_names = FALSE)
df <- as.matrix(df[,1:6])
# for(i in 1:7){
#   df[,i] <- unlist(df[,i])
# }
#effectiveSize(df[4000:8000,4])
plot(df[500:856,5],type="l")

effectiveSize(df[500:856,])

df <- read_table2("blackHoleResultsSinglePrec.txt",col_names = FALSE)
df <- as.matrix(df[,1:5])
# for(i in 1:7){
#   df[,i] <- unlist(df[,i])
# }
#effectiveSize(df[4000:8000,4])
plot(df[,4],type="l")





# results <- readRDS("~/qpMCMC/blackHole.rds")
# 
# 1024*20000000/results[[4]]
# plot(results[[8]],type="l")
# 
# imgLong <- melt(results[[1]][[1000]])
# remove(results)
# ggplot(imgLong, aes(x = Var2, y = Var1)) +
#   geom_raster(aes(fill=value))
