##-----------------------------------------------------------------
## hclust_dendro.R is an Rscript used to 
## rpkm normalized values and generate clustered output in
## "De novo assembly of the green ash transcriptome and 
## identification of genes responding to abiotic and biotic stress"
## Lane et. al. 2015
##
## hclust_dendro.R assumes a directory in home called "hclust_dendro" exists
## hclust_dendro.R also assumes a sub-directory in home called rpkm
## containing "all_final_cut55.tsv" a file of rpkm normalized values
## accross 55 green ash libraries. This file can be produced using the code
## in the rpkm directory
##
##
## hclust_dendro.R will create a single output Tiff "hclust_dendro.tiff"
##

# set up R environment
library("Hmisc")
library("ape")

# load in file as a matrix
setwd("~/hclust_dendro/rpkm")
rpkm <- as.matrix(read.table("all_final_cut55.tsv", sep="\t", header=TRUE))
rpkmone <- rpkm + 1

# log2 transform
logrpkm <- log2(rpkmone)

# run pearson correlation
pearsonrpkm <- cor(logrpkm, method = "pearson")

# create distance object
mydata <- 1 - pearsonrpkm
d <- as.dist(mydata)

# cluster
cluster <- hclust(d, method = "average", members = NULL)

# plot basic tree
setwd("~/hclust_dendro")
tiff("hclust_dendro.tiff", width = 4, height = 4, units = 'in', res = 300)
plot(as.phylo(cluster), cex = 0.25, label.offset = 0.005)
dev.off()
