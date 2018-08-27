######################################################
# This is the script for analyzing count-matrix data #
#     from the RNAcluster pipeline on killdevil.     #
######################################################


#===========================================================================================================#
#                                                Initialize                                                 #
#===========================================================================================================#


#------------------------------------------------ Settings ---------------------------------------------------#

################
# RUN SETTINGS #
################
# Make PDFs or run in Rstudio only?
makePDF = T

# Read in from raw data?
readRaw = T

# Recreate dds and results files?
runDESeq = F
runShrink = F

# Make which plots/files?
txiPlot = T
PCAint = F
PCAplot = T
countPlot = T
MAplot = T
GOplot = T




##################
# COLOR SETTINGS #
##################
# Assign sample colors for PCA, plot counts
sample.pal=parula(8)

# Assign GO term colors
BPcol <- "dodgerblue2"
CCcol <- "firebrick2"
MFcol <- "grey"


#################
# PLOT SETTINGS #
#################
par(mar=c(1,1,1,1))
par(mfrow=c(1,1))


#####################
# VARIABLE SETTINGS #
#####################
# Assign project variables
proj = 'LIMA'
projNum = '2-3.1.1'
sort = 'trmt'

transformType = "var"

# Assign cutoffs
p.thr = .01
lfc.thr = 1

# Assign output directory for PDFs, CSVs, etc.
outputName = paste(proj, projNum, Sys.Date(), sep="_")
outputDir = file.path("/Users/phanstiel3/Desktop", outputName)
#outputDir = file.path("/Users/phanstiel3/Research/Data/Projects", proj, "rna/diff")

# Make output directory (in case it's necessary)
dirCmd = paste("mkdir", outputDir, sep=" ")
system(dirCmd)

# Assign GO marts
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host= "www.ensembl.org")


##################
# LOAD LIBRARIES #
##################
# Load any required libraries
library(limma)
library(biomaRt)
library(org.Hs.eg.db)
library(GO.db)

library(DESeq2)
library(Sushi2)
library(colorspace)
library(pheatmap)
library(ggplot2)
library(pcaExplorer)
library(ReportingTools)
library(factoextra)
library(NbClust)
library(pals)


#------------------------------------------------ Functions ------------------------------------------------#

# General function for rounding a number to a variable base
mround <- function(x,base, mode="std")
{ 
  if(mode == "std"){
    out = base*round(x/base) 
  }
  
  if(mode == "floor"){
    out = base*floor(x/base) 
  }
  
  if(mode == "ceiling"){
    out = base*ceiling(x/base) 
  }
  
  return(out)
}


# Function for shortening the decimal ENSEMBL name from dds, res, etc. to the standard ENSEMBL gene ID
nameAbridge <- function(longName)
{
  split <- strsplit(longName, "[.]")
  shortName <- unlist(split)[2*(1:length(longName))-1]
}


# Function for extending the ENSEMBL name to the decimal version in dds, res, etc.
nameExtend <- function(shortName, longKey = longGeneKey, abridgedKey = abridgedGeneKey)
{
  longNames=list()
  for (name in shortName){
    longName = longKey[which(abridgedKey==name)]
    longNames = append(longNames, longName)
  }
  return(longNames)
}


# Function for circling (and optionally labeling) a list of genes on an MA plot
geneGroupPlot <- function(results, geneName, geneID, abridgedKey = abridgedGeneKey)
{
  if (length(geneName) != length(geneID))
  {
    stop("Gene name and gene ID lists are not same length.")
  }
  for (n in seq(1, length(geneName)))
  {
    bMean = results[which(abridgedKey == geneID[n]),1]
    foldChange = results[which(abridgedKey == geneID[n]),2]
    
    points(bMean, foldChange, pch=1)
    text(bMean, foldChange, labels=geneName[n], cex=.7, pos=3)
  }
}


# Function for plotting all MA-plots
MAplotter <- function(results, col1="dodgerblue2", col2="firebrick2", title1="Neg", title2="Pos", geneNameIDList=NULL, pval = .1)
{
  resSig <- results[which(results$padj <= pval),]
  resPos <- resSig[which(resSig$log2FoldChange > 0),]
  resNeg <- resSig[which(resSig$log2FoldChange < 0),]
  
  min.LFC = mround(min(results[,2][!is.na(results[,2])]), 5, "floor")
  max.LFC = mround(max(results[,2][!is.na(results[,2])]),5, "ceiling")
  
  max.counts = mround(max(results[,2][!is.na(results[,2])]), 1000000, "ceiling")
  
  plot(results[,1], results[,2], log="x", xlim=c(.15, max.counts), ylim=c(min.LFC, max.LFC), cex=.25, col="gray38", pch=19, xlab="", ylab="")
  points(resNeg[,1], resNeg[,2], pch=19, cex=.25, col=col1)
  points(resPos[,1], resPos[,2], pch=19, cex=.25, col=col2)
  lines(c(.15, max.counts), c(0,0), col=alpha("indianred1", .5), lwd=3)
  
  main.title = paste("MA-Plot:", title1, "vs", title2, sep=" ")
  title(main=main.title)
  title(xlab="mean of normalized counts", line=2)
  title(ylab="log2 fold change", line=2.5)
  
  if(!is.null(geneNameIDList)){
    geneNames <- geneNameIDList[,1]
    geneIDs <- geneNameIDList[,2]
    geneGroupPlot(results, geneNames, geneIDs)
  }
  
  legtxt = c(paste0(title1," (n=",nrow(resPos), ")"), paste0(title2," (n=",nrow(resNeg), ")"))
  legend("topright", legend=legtxt, pch=19, col=c(col1, col2), cex=.8)
}  


# limma-based GO enrich function
GOenrich <- function(set,background, mart=mart)
{
  bkgd <- c(set,background)
  bkgd <- bkgd[!duplicated(bkgd)]
  
  ensembl2entrez <- getBM(filters="ensembl_gene_id",attributes=c("ensembl_gene_id", "entrezgene"), values=bkgd, mart=mart)
  
  set_entrez <- ensembl2entrez[which(ensembl2entrez[,1] %in% set),2]
  bkgd_entrez <- ensembl2entrez[which(ensembl2entrez[,1] %in% bkgd),2]
  
  set_entrez = set_entrez[complete.cases(set_entrez)]
  bkgd_entrez = bkgd_entrez[complete.cases(bkgd_entrez)]
  
  enriched = goana(de=set_entrez, universe=bkgd_entrez,species="Hs")
  enriched = enriched[order(enriched$P.DE),]
  
  return(enriched)
}


# Wrap text functions for plotting GO Enrichments
wrap.it <- function(x, len)
{
  sapply(x, function(y) paste(strwrap(y, len), collapse="\n"), USE.NAMES=F)
}

wrap.labels <- function(x, len)
{
  lapply(x,wrap.it,len)
}


# Function for running GO enrichments for up/down regulated genes
GOrunner <- function(results, dir)
{
  # Extract genes enriched in either sample (neg: higher in first [het]; pos: higher in second [los])
  pos <- rownames(results[which(results$padj < .05 & results$log2FoldChange > 0),])
  neg <- rownames(results[which(results$padj < .05 & results$log2FoldChange < 0),])
  back <- rownames(results)
  
  
  # Convert gene names to those without the decimal points for compatibility with BioMart
  posSplit <- strsplit(pos, "[.]")
  posAb <- unlist(posSplit)[2*(1:length(pos))-1]
  
  negSplit <- strsplit(neg, "[.]")
  negAb <- unlist(negSplit)[2*(1:length(neg))-1]
  
  backSplit <- strsplit(back, "[.]")
  backAb <- unlist(backSplit)[2*(1:length(back))-1]
  
  
  # Run the GO Enrichment function
  if(dir="pos"){
    enrich <- GOenrich(posAb,backAb)
  }
  
  if(dir="neg"){
    enrich <- GOenrich(negAb,backAb)
  }
  
  return(enrich)
}


# Function for assigning colors to different GO groups
GOcolor <- function(list)
{
  list[list=="BP"] <- BPcol
  list[list=="CC"] <- CCcol
  list[list=="MF"] <- MFcol
  
  return(list)
}


# Function for combining replicates in a 16x(many genes) matrix 
combineReps <- function(matrix)
{
  newMatrix = c()
  for (i in seq(1,ncol(matrix)/2))
  {
    tempMatrix = matrix(data=c(matrix[,i], matrix[,i+ncol(matrix)/2]), ncol=2, byrow=F)
    newCol = apply(tempMatrix,1,median)
    newMatrix  = cbind(newMatrix,newCol)
  }
  rownames(newMatrix) = rownames(matrix)
  
  return(newMatrix)
}


# Function for writing data to a csv
geneWrite <- function(res, geneList, symbolList, count1_1, count1_2, count2_1, count2_2, sampleName1, sampleName2, filePath, pvalCutoff = .01)
{
  geneListLong = nameExtend(geneList)
  
  # Make list of the fold change (straight from res)
  foldChange <- res$log2FoldChange[which(rownames(res) %in% geneListLong)]
  
  # Make list of the adjusted pvalues (straight from res as well)
  pVal <- res$padj[which(rownames(res) %in% geneListLong)]
  
  # Make list based on p-values and foldChange that indicates which sample a gene is higher in, if any
  higherIn <- rep("", times=length(geneListLong))
  for (n in seq(1, length(geneListLong))){
    if (is.na(pVal[n])){
      higherIn[n] = "Neither"
    } else{
      if ((pVal[n] <= pvalCutoff) && (foldChange[n] > 0)){
        higherIn[n] = sampleName2
      }
      if ((pVal[n] <= pvalCutoff) && (foldChange[n] < 0)){
        higherIn[n] = sampleName1
      }
      if (pVal[n] > pvalCutoff){
        higherIn[n] = "Neither"
      }
    }
  }
  
  # Put it all in one data frame
  data <- data.frame(geneList, symbolList, AP1ind, count1_1, count1_2, count2_1, count2_2, foldChange, higherIn, pVal, 
                     stringsAsFactors = FALSE)
  
  # Write to csv
  write.csv(data, file=filePath, row.names = F)
}



#===========================================================================================================#
#                                                   Run                                                     #
#===========================================================================================================#

#----------------------------------------------- Read in ---------------------------------------------------#

################
# READ IN DATA #
################
# Read in the txi file, as defined in the config file -----------------------------------------------------------> *** READRAW ***
if(readRaw == T){
  # Read in config file
  configPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_', projNum, '_samples.csv')
  config <- read.csv(configPath, header=T)
  
  # Read in count matrix
  txiPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_', projNum, '_txi.rds')
  txi <- readRDS(txiPath)
}


# Plot read info based on the txi file --------------------------------------------------------------------------> *** TXIPLOT ***
if(txiPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "readInfo.pdf"), width=8, height=8)
  }
  
  # Plot count info
  par(mfrow=c(2,2))
  
  txi.x.labels = substr(colnames(txi$counts), 16, 25)
  
  # Count boxplots for each sample (gene avg)
  boxplot(txi$counts,outline=F,xaxt='n', main="Counts")
  axis(side=1,las=2,at=1:16,labels=txi.x.labels)
  
  # Total summed counts for each sample
  bp = barplot(colSums(txi$counts),xaxt='n', main="Summed counts per sample")
  axis(side=1,las=2,at=bp,labels=txi.x.labels)
  
  # Histogram of total counts per gene
  hist(rowSums(txi$counts), main="Summed counts per gene")
  rug(rowSums(txi$counts))
  
  # Histogram of total counts per gene for a single sample
  hist(txi$counts[,3], main="Counts per gene for 60m")
  rug(txi$counts[,3])
  
  if(makePDF == T){
    dev.off()
  }
  
}



#############
# RUN DESEQ #
#############
# Create dds and transformed dds  ------------------------------------------------------------------------------> *** RUNDESEQ ***
if(runDESeq == T){
  # Create data frame with sample names and whatever condition you want to distinguish them by (config$Condition, config$Cell.Type, etc.)
  colData <- data.frame(condition = factor(config$Condition), rep = factor(config$BioRep), row.names = config$Name)
  
  
  # Create DESeq Data Set from txImport, using LRT for p-values
  ddsTxI <- DESeqDataSetFromTximport(txi, colData=colData, design =~ rep + condition)
  
  dds <- DESeq(ddsTxI, test="LRT", full=~rep + condition, reduced = ~rep)
  
  # Transform the dds (rlog, vst or ntd) --------------------------------------------------------------------> *** TRANSFORMTYPE ***
  if(transformType == "norm"){
    dds.trans <- normTransform(dds)
  }
  
  if(transformType == "rlog"){
    dds.trans <- rlog(dds, blind=FALSE)
  }
  
  if(transformType == "var"){
    dds.trans <- vst(dds)
  }
}  

# Find shrunken LFC for each comparison vs 0 ------------------------------------------------------------------> *** RUNSHRINK ***
if(runShrink == T){
  res.030 <- lfcShrink(dds, coef="condition_30_vs_0", type="apeglm")
  res.060 <- lfcShrink(dds, coef="condition_60_vs_0", type="apeglm")
  res.090 <- lfcShrink(dds, coef="condition_90_vs_0", type="apeglm")
  res.0120 <- lfcShrink(dds, coef="condition_120_vs_0", type="apeglm")
  res.0240 <- lfcShrink(dds, coef="condition_240_vs_0", type="apeglm")
  res.0360 <- lfcShrink(dds, coef="condition_360_vs_0", type="apeglm")
  res.01440 <- lfcShrink(dds, coef="condition_1440_vs_0", type="apeglm")
}



##############
# OTHER DATA #
##############
# Make long + abridged gene list from dds 
longGeneKey <- rownames(assay(dds))
abridgedGeneKey <- nameAbridge(longKey)


# Read in excel of genes of interest (gene name and ENSG ID)
GoI <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(GoI) <- c("gene", "ID")

# Add long ID based on dds
for (n in 1:nrow(GoI))
  {
  GoI$long.ID[n] = nameExtend(GoI$ID[n])[[1]]
}




#------------------------------------------------ Analyze ----------------------------------------------------#

#################
# GLOBAL TRENDS #
#################

#-------------------#
# INTERACTIVE PLOTS #
#-------------------#
# PCA explorer (featuring pca2go, which requires the names be shortened to regular ENSEMBL IDs) ------------------> *** PCAINT ***
if(PCAint == T){
  dds.rename = dds
  rownames(dds.rename) <- nameAbridge(rownames(dds.rename))
  pcaExplorer(dds=dds.rename)
}


#----------#
# PCA PLOT #
#----------#
# Plot PCA of transformed counts --------------------------------------------------------------------------------> *** PCAPLOT ***
if(PCAplot == T){
  if(makePDF == T){
    pdf(file=file.path(outputDir, "PCAplot.pdf"), width=8, height=8)
  }
  
  # Plot PCA of transformed counts
  par(mfrow=c(1,1))
  plotPCA(dds.trans, intgroup="condition")+labs(title="Transformed Counts PCA")+scale_color_manual(values=sample.pal)
  
  if(makePDF == T){
    dev.off()
  }
}


#-------------#
# PLOT COUNTS #
#-------------#
# Plot counts of genes of interest ----------------------------------------------------------------------------> *** COUNTPLOT ***
if(countPlot == T){
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countPlot.pdf"), width=8, height=8)
  }
  
  # Plot settings
  n=nrow(GoI)
  par(mar=c(3,2,1,2))
  par(mfrow=c(3,mround(n,3)/3))
  
  # Go through genes of interest (GoI) and plot the counts of each gene listed (note: plotCounts output is the same as if you manually plotted the vst dds assay values)
  for(n in 1:nrow(GoI)){
    plotCounts(dds, GoI$long.ID[n], main=GoI$gene[n], pch=16, col=rep(sample.pal,times=2),
              intgroup="condition")
  }
  
  if(makePDF == T){
    dev.off()
  }
}



#####################
# 1 v 1 COMPARISONS #
#####################

#----------#
# MA PLOTS #
#----------#
# Make MA Plots for each comparison vs 0 -------------------------------------------------------------------------> *** MAPLOT ***
if(MAplot == T){
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countPlot.pdf"), width=12, height=6)
  }
  
  # Plot settings
  par(mar=c(3,2,1,2))
  par(mfrow=c(2,4))
  
  # Plot each with shrunken LFC, LRT p-val
  MAplotter(res.030, col1=sample.pal[1], col2=sample.pal[2],title1="0000", title2="0030", geneNameIDList = ap1)
  MAplotter(res.060, col1=sample.pal[1], col2=sample.pal[3],title1="0000", title2="0060", geneNameIDList = ap1)
  MAplotter(res.090, col1=sample.pal[1], col2=sample.pal[4],title1="0000", title2="0090", geneNameIDList = ap1)
  MAplotter(res.0120, col1=sample.pal[1], col2=sample.pal[5],title1="0000", title2="0120", geneNameIDList = ap1)
  MAplotter(res.0240, col1=sample.pal[1], col2=sample.pal[6],title1="0000", title2="0240", geneNameIDList = ap1)
  MAplotter(res.0360, col1=sample.pal[1], col2=sample.pal[7],title1="0000", title2="0360", geneNameIDList = ap1)
  MAplotter(res.01440, col1=sample.pal[1], col2=sample.pal[8],title1="0000", title2="1440", geneNameIDList = ap1)
  
  if(makePDF == T){
    dev.off()
  }
}


#----------------#
# GO ENRICHMENTS #
#----------------#
# Make GO graphs for each comparison vs 0 ------------------------------------------------------------------------> *** GOPLOT ***
if(GOplot == T){
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countPlot.pdf"), width=12, height=6)
  }
  
  # Plot settings
  par(mar=c(3,2,1,2))
  par(mfrow=c(2,4))
  
  # Plot each with shrunken LFC, LRT p-val
  MAplotter(res.030, col1=sample.pal[1], col2=sample.pal[2],title1="0000", title2="0030", geneNameIDList = ap1)
  MAplotter(res.060, col1=sample.pal[1], col2=sample.pal[3],title1="0000", title2="0060", geneNameIDList = ap1)
  MAplotter(res.090, col1=sample.pal[1], col2=sample.pal[4],title1="0000", title2="0090", geneNameIDList = ap1)
  MAplotter(res.0120, col1=sample.pal[1], col2=sample.pal[5],title1="0000", title2="0120", geneNameIDList = ap1)
  MAplotter(res.0240, col1=sample.pal[1], col2=sample.pal[6],title1="0000", title2="0240", geneNameIDList = ap1)
  MAplotter(res.0360, col1=sample.pal[1], col2=sample.pal[7],title1="0000", title2="0360", geneNameIDList = ap1)
  MAplotter(res.01440, col1=sample.pal[1], col2=sample.pal[8],title1="0000", title2="1440", geneNameIDList = ap1)
  
  if(makePDF == T){
    dev.off()
  }
}

##############
# SUBSETTING #
##############

#----------------------#
# SHRUNKEN LFC + P-VAL #
#----------------------#
# subset significant genes according to LRT p-values
res.030.sig <-  res.030[which(res.030$padj < p.thr),]
res.060.sig <-  res.060[which(res.060$padj < p.thr),]
res.090.sig <-  res.090[which(res.090$padj < p.thr),]
res.0120.sig <-  res.0120[which(res.0120$padj < p.thr),]
res.0240.sig <-  res.0240[which(res.0240$padj < p.thr),]
res.0360.sig <-  res.0360[which(res.0360$padj < p.thr),]
res.01440.sig <-  res.01440[which(res.01440$padj < p.thr),]


# subset genes changing over 2-fold according to shrunken LFC
# Note: negative fold-changes indicate higher in second sample (i.e. 0 vs *___*)
res.030.diff <- res.030.sig[which(abs(res.030.sig$log2FoldChange) >= lfc.thr),]
res.060.diff <- res.060.sig[which(abs(res.060.sig$log2FoldChange) >= lfc.thr),]
res.090.diff <- res.090.sig[which(abs(res.090.sig$log2FoldChange) >= lfc.thr),]
res.0120.diff <- res.0120.sig[which(abs(res.0120.sig$log2FoldChange) >= lfc.thr),]
res.0240.diff <- res.0240.sig[which(abs(res.0240.sig$log2FoldChange) >= lfc.thr),]
res.0360.diff <- res.0360.sig[which(abs(res.0360.sig$log2FoldChange) >= lfc.thr),]
res.01440.diff <- res.01440.sig[which(abs(res.01440.sig$log2FoldChange) >= lfc.thr),]

# combine all genes into one, non-redundant list
genes.sig.diff = c(rownames(res.030.diff), rownames(res.060.diff), rownames(res.090.diff), rownames(res.0120.diff), rownames(res.0240.diff), rownames(res.0360.diff), rownames(res.01440.diff))
genes.sig.diff = rownames(res.LRT.full)[rownames(res.LRT.full) %in% genes.sig.diff]



