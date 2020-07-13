
#how to update R
#################

install.packages("installr"); library(installr)
updateR()


#Add ggplot (figure editor)
###########################

install.packages("ggplot2")
library(ggplot2)


#Some important insall manager? for bio tools? like DESeq2?
###########################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("DESeq2")

a
#To view documentation of things
################################

browseVignettes("DESeq2")

or ?DESeq2

###############################
#PCA explorer. Trying new things for replicate correlation

BiocManager::install("pcaExplorer")
