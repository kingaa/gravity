file <- commandArgs(trailingOnly=TRUE)

library(doParallel)

cl <- makeCluster(210,type="MPI")

registerDoParallel(cl)

source(file,echo=TRUE,keep.source=TRUE)

stopCluster(cl)
