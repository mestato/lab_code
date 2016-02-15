##-----------------------------------------------------------------
## degs.R (differentially expressed genes) is an Rscript used to 
## calculate and construct some of the figures/tables in the manuscript
## "De novo assembly of the green ash transcriptome and 
## identification of genes responding to abiotic and biotic stress"
## Lane et. al. 2015
##
## degs.R assumes a directory in home called "degs" exists
## degs.R also assumes a sub-directory in home called counts-sorted
## exists and contains HTSeq output from green ash libraries. See
## counts-sorted.zip. 
##
## degs.R will create two directories for outputs: output and images
##

##-----------------------------------------------------------------
## Load required packages
##-----------------------------------------------------------------
library("DESeq2")
library("RColorBrewer")
library("ggplot2")
library("gplots")
library("grid")
library("VennDiagram")
library("pheatmap")

##-----------------------------------------------------------------
## Setup sub-directories for output
##-----------------------------------------------------------------
mainDir <- "~/degs/"
subDirOne  <- "output"
subDirTwo  <- "images"
dir.create(file.path(mainDir, subDirOne), showWarnings = FALSE)
dir.create(file.path(mainDir, subDirTwo), showWarnings = FALSE)

##-----------------------------------------------------------------
## Cold Stress 
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- c("GA-COL",
                 "GA-COP",
                 "GA-COR",
                 #"24_GA-CL",
                 #"24_GA-CP",
                 #"24_GA-CR",
                 "GA-CL",
                 "GA-CP",
                 "GA-CR")
tissue<-rep( c("Leaf",
               "Petiole",
               "Root"),2)
treatment <- c(rep("control",3),
               rep("Cold24",3))
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  tissue=tissue,
  treatment=treatment)
directory<-("~/degs/counts-sorted")

# run differential expression
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ tissue + treatment)
dds$treatment<-relevel(dds$treatment, "control")
dds        <- DESeq(dds)

# gather results
res        <- results(dds)
resOrdered <- res[order(res$padj),]
resSig     <- subset(resOrdered, padj < 0.01)
resSigUp   <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)

# output results
setwd("~/degs/output")
write.csv(as.data.frame(resSigUp),
          file="Cold.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Cold.Dw.DE.csv")

##-----------------------------------------------------------------
## Drought Stress
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- c("GA-COL",
                 "GA-COP",
                 "GA-COR",
                 "GA-DL",
                 "GA-DP",
                 "GA-DR")
tissue<-rep( c("Leaf",
               "Petiole",
               "Root"),2)
treatment <- c(rep("control",3),
               rep("Drought",3))
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  tissue=tissue,
  treatment=treatment)
directory<-("~/degs/counts-sorted/")

# run differential expression
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ tissue + treatment)
dds$treatment<-relevel(dds$treatment, "control")
dds<-DESeq(dds)

# gather results
res        <- results(dds)
resOrdered <- res[order(res$padj),]
resSig     <- subset(resOrdered, padj < 0.01)
resSigUp   <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)

# output results
setwd("~/degs/output")
write.csv(as.data.frame(resSigUp),
          file="Drought.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Drought.Dw.DE.csv")

##-----------------------------------------------------------------
## Heat
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- c("GA-COL",
                 "GA-COP",
                 "GA-COR",
                 "GA-HL",
                 "GA-HP",
                 "GA-HR")
tissue<-rep( c("Leaf",
               "Petiole",
               "Root"),2)
treatment <- c(rep("control"),3),
               rep("Heat",3))
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  tissue=tissue,
  treatment=treatment)
directory<-("~/degs/counts-sorted/")

# run differential expression
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ tissue + treatment)
dds$treatment<-relevel(dds$treatment, "control")
dds<-DESeq(dds)

# gather results
res        <- results(dds)
resOrdered <- res[order(res$padj),]
resSig     <- subset(resOrdered, padj < 0.01)
resSigUp   <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)

# output results
setwd("~/degs/output")
write.csv(as.data.frame(resSigUp),
          file="Heat.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="Heat.Dw.DE.csv")

##-----------------------------------------------------------------
## EAB Treatment vs. Control (TvC)
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- grep("PE-",list.files(),value=TRUE)
genotype<-c(rep("19",2),
            rep("21",2),
            rep("22",2),
            rep("24",2),
            rep("36",2),
            rep("SM",2))
treatment<-rep(c("control", "treatment"),6)

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  genotype=genotype,
  treatment=treatment)

# run differential expression
directory<-("~/degs/counts-sorted")
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ genotype + treatment)
dds$treatment<-relevel(dds$treatment, "control")
# wald test
dds<-DESeq(dds)

# gather results
res        <- results(dds)
resOrdered <- res[order(res$padj),]
resSig     <- subset(resOrdered, padj < 0.01)
resSigUp   <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)

# output results
setwd("~/degs/output")
write.csv(as.data.frame(resSigUp),
          file="EAB_TvC.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="EAB_TvC.Dw.DE.csv")
rld <- rlog(dds)

##-----------------------------------------------------------------
## EAB - PRINCIPAL COMPONENT ANALYSIS
##-----------------------------------------------------------------

# principal component analysis of the samples
#plotPCA(rld, intgroup=c("treatment", "genotype"))
data <- plotPCA(rld, intgroup=c("treatment", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
setwd("~/degs/images")
tiff("PCA_EAB.tiff", width = 10, height = 10, units = 'in', res = 100)
ggplot(data, aes(PC1, PC2, color=genotype, shape=treatment)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

##-----------------------------------------------------------------
## EAB - RESISTANT VS SUSCEPTIBLE - Initial Time Point
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- c("PE-19",
                 "PE-21",
                 "PE-22",
                 "PE-24",
                 "PE-36",
                 "PE-SUM")
susceptibility<-c(rep("Resistant",4), 
                  rep("Susceptible", 2))
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  susceptibility=susceptibility)

# run differential expression
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ susceptibility)
dds$treatment<-relevel(dds$treatment, "Susceptible")
dds<-DESeq(dds)

# gather results
res        <- results(dds)
resOrdered <- res[order(res$padj),]
resSig     <- subset(resOrdered, padj < 0.01)
resSigUp   <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)

# output results
setwd("~/degs/output/")
write.csv(as.data.frame(resSigUp),
          file="EAB_RvS.InitialTime.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="EAB_RvS.InitialTime.Dw.DE.csv")

##-----------------------------------------------------------------
## EAB - RESISTANT VS SUSCEPTIBLE - Later Time Point
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- c("PE-19_egg",
                 "PE-21_egg",
                 "PE-22_egg",
                 "PE-24_egg",
                 "PE-36_egg",
                 "PE-SUM_egg")
susceptibility<-c(rep("Resistant",4), 
                  rep("Susceptible", 2))
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  susceptibility=susceptibility)

# run differential expression
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ susceptibility)
dds$treatment<-relevel(dds$treatment, "Susceptible")
dds<-DESeq(dds)

# gather results
res        <- results(dds)
resOrdered <- res[order(res$padj),]
resSig     <- subset(resOrdered, padj < 0.01)
resSigUp   <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)

# output results
setwd("~/degs/output/")
write.csv(as.data.frame(resSigUp),
          file="EAB_RvS.LaterTime.Up.DE.csv")
write.csv(as.data.frame(resSigDown),
          file="EAB_RvS.LaterTime.Dw.DE.csv")

##-----------------------------------------------------------------
## Mechanical Wounding
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- c("GA-WL-1",
                 "GA-WL",
                 "GA-COL",
                 "GA-WP-1",
                 "GA-WP",
                 "GA-COP",
                 "GA_M-1",
                 "GA_M-2",
                 "GA_M-3",
                 "GA_M-4")
pooled<-c(rep("pooled",6),
          rep("not_pooled",4))
tissue<-c("Leaf",
          "Leaf",
          "Leaf",
          "Petiole",
          "Petiole",
          "Petiole",
          "Leaf",
          "Leaf",
          "Twig",
          "Twig")
treatment <- c("5 hours",
               "24 hours",
               "control",
               "5 hours",
               "24 hours",
               "control",
               "control",
               "24 hours",
               "control",
               "24 hours")
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  tissue=tissue,
  pooled=pooled,
  treatment=treatment)
directory<-("~/degs/counts-sorted")

# run differential expression
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ tissue + pooled + treatment)
dds$treatment <- relevel(dds$treatment, "control")
dds           <- DESeq(dds)
resultsNames(dds)

##----------------------------------------------
## Mechanical Wounding - principal component analysis of the samples
##----------------------------------------------

# setup PCA inputs
rld <- rlog(dds)
#plotPCA(rld, intgroup=c("treatment", "tissue"))
data <- plotPCA(rld, intgroup=c("treatment", "tissue"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

# output results
setwd("~/degs/images/")
tiff("PCA_MW.tiff", width = 10, height = 10, units = 'in', res = 100)
ggplot(data, aes(PC1, PC2, color=tissue, shape=treatment)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

##----------------------------------------------
## Mechanical Wounding - 5 hour vs control 
##----------------------------------------------

# gather results
res5hr        <- results(dds, contrast=c("treatment","control","5 hours"))
res5hrOrdered <- res5hr[order(res5hr$padj),]
res5hrSig     <- subset(res5hrOrdered, padj < 0.01)
res5hrSigUp   <- subset(res5hrSig, log2FoldChange > 0)
res5hrSigDown <- subset(res5hrSig, log2FoldChange < 0)

# output results
setwd("~/degs/output/")
write.csv(as.data.frame(res5hrSigUp),
          file="MW_5hr_vs_Control.Up.DE.csv")
write.csv(as.data.frame(res5hrSigDown),
          file="MW_5hr_vs_Control.Dw.DE.csv")

##----------------------------------------------
## Mechanical Wounding - 24 hour vs control 
##----------------------------------------------

# gather results
res24hr        <- results(dds, contrast=c("treatment","control","24 hours"))
res24hrOrdered <- res24hr[order(res24hr$padj),]
res24hrSig     <- subset(res24hrOrdered, padj < 0.01)
res24hrSigUp   <- subset(res24hrSig, log2FoldChange > 0)
res24hrSigDown <- subset(res24hrSig, log2FoldChange < 0)

# output results
setwd("~/degs/output/")
write.csv(as.data.frame(res24hrSigUp),
          file="MW_24hr_vs_Control.Up.DE.csv")
write.csv(as.data.frame(res24hrSigDown),
          file="MW_24hr_vs_Control.Dw.DE.csv")

##-----------------------------------------------------------------
## Ozone
##-----------------------------------------------------------------

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- grep("ozone",list.files(),value=TRUE)
sampleFiles <- c("GA1-ozone-7hrs-control.txt",
                 "GA2-ozone-7hrs-80ppb.txt",
                 "GA3-ozone-7hrs-125ppb.txt",
                 "GA4-ozone-7hrs-225ppb.txt",
                 "GA5-ozone-14days-control.txt",
                 "GA6-ozone-14days-80ppb.txt",
                 "GA7-ozone-14days-125ppb.txt",
                 "GA8-ozone-14days-225ppb.txt",
                 "GA9-ozone-28days-control.txt",
                 "GA10-ozone-28days-80ppb.txt",
                 "GA11-ozone-28days-125ppb.txt",
                 "GA12-ozone-28days-225ppb.txt")
time=c( rep("07",4),
        rep("336",4),
        rep("672",4))
treatment=rep( c("0","80", "125", "225"), 3)
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  time=time,
  treatment=treatment)
directory<-("~/degs/counts-sorted")

# run differential expression
dds_trmt<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                     directory=directory,
                                     design= ~time + treatment)
dds$treatment <- relevel(dds$treatment, "control")
dds$time      <- relevel(dds$time, "7")
dds_trmt      <- DESeq(dds_trmt, test = c("LRT"), reduced = ~time)

# gather results
res_trmt        <- results(dds_trmt)
resOrdered_trmt <- res_trmt[order(res_trmt$padj),]
resSig_trmt     <- subset(resOrdered_trmt, padj < 0.01)

# output results
setwd("~/degs/output/")
write.csv(as.data.frame(resSig_trmt),
          file="Ozone.DE.csv")

##------------------------------------------
## Ozone - heatmap 
##------------------------------------------

# prepare input plotting data from ozone deseq2 results
betas <- coef(dds_trmt)
colnames(betas)
topGenes <- head(order(res_trmt$padj),350)
mat <- betas[topGenes, -c(1)]
thr <- 3 # threshold for plotting
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
matOrder <- mat[,c(1,2,5,3,4)]

# plotHeatmap
setwd("~/degs/images/")
tiff("OzoneHeatmap.tiff", width = 4, height = 10, units = 'in', res = 100)
pheatmap(matOrder, breaks=seq(from=-thr, to=thr, length=101), show_rownames=FALSE,cluster_col=FALSE, cutree_rows=8, labels_col=c("14 days vs 7 hours","28 days vs 7 hours","80ppb vs normal atmospheric ozone","125ppb vs normal atmospheric ozone","225ppb vs normal atmospheric ozone"))
dev.off()

#
#
#
#----------------------------------------------
# Mechanical Wounding + Ozone
#----------------------------------------------
# The following produces no results
# empty output files are not generated
#

# setup sampleTable
setwd("~/degs/counts-sorted/")
sampleFiles <- c("GA_OW-9",
                 "GA_OW-10",
                 "GA_OW-11",
                 "GA_OW-12")
treatment <- c("control", "80ppb","125ppb","225ppb")
sampleTable<-data.frame(
  sampleName=sampleFiles,
  filename=sampleFiles,
  treatment=treatment)
directory<-("~/degs/counts-sorted")

# run differential expression
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                directory=directory,
                                design= ~ treatment)
dds$treatment<-relevel(dds$treatment, "control")

#----------------------------------------------
# Contrast: 0 vs 80
#----------------------------------------------
dds<-DESeq(dds)
res80<-results(dds, contrast=c("treatment","control","80ppb"))
res80Ordered<-res80[order(res80$padj),]
res80Sig <- subset(res80Ordered, padj < 0.01)

# gather results
res80SigUp <- subset(res80Sig, log2FoldChange > 0)
res80SigDown <- subset(res80Sig, log2FoldChange < 0)

#----------------------------------------------
# Contrast: 0 vs 125
#----------------------------------------------
dds<-DESeq(dds)
res125<-results(dds, contrast=c("treatment","control","125ppb"))
res125Ordered<-res125[order(res125$padj),]
res125Sig <- subset(res125Ordered, padj < 0.01)

# gather results
res125SigUp <- subset(res125Sig, log2FoldChange > 0)
res125SigDown <- subset(res125Sig, log2FoldChange < 0)

#----------------------------------------------
# Contrast: 0 vs 225
#----------------------------------------------
dds<-DESeq(dds)
res225<-results(dds, contrast=c("treatment","control","225ppb"))
res225Ordered<-res225[order(res225$padj),]
res225Sig <- subset(res225Ordered, padj < 0.01)

# gather results
res225SigUp <- subset(res225Sig, log2FoldChange > 0)
res225SigDown <- subset(res225Sig, log2FoldChange < 0)

#----------------------------------------------
# LRT
#----------------------------------------------
dds<-DESeq(dds, test = c("LRT"), reduced = ~1)
res<-results(dds)

# gather results
resOrdered<-res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.01)
