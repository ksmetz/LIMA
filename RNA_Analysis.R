######################################################
# This is the script for analyzing count-matrix data #
#     from the RNAcluster pipeline on killdevil.     #
######################################################

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

# Function for shortening the decimal ENSEMBL name from dds, res, etc. to the standard ENSEMBL gene ID
nameAbridge <- function(longName){
  split <- strsplit(longName, "[.]")
  shortName <- unlist(split)[2*(1:length(longName))-1]
}


# Function for extending the ENSEMBL name to the decimal version in dds, res, etc.
nameExtend <- function(results, shortName, abridgedKey = abridgedKey)
{
  longNames=list()
  for (name in shortName){
    longName = rownames(results[which(abridgedKey==name),])
    longNames = append(longNames, longName)
  }
  return(longNames)
}


# Function for circling (and optionally labeling) a list of genes on an MA plot
geneGroupPlot <- function(results, geneName, geneID, abridgedKey = abridgedKey)
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


# Functions for spliting DESeq2 results into sig up/downreg genes
MAsplitter <- function(results, dir)
{
  resSig <- results[which(results$padj <= .1),]
  resPos <- resSig[which(resSig$log2FoldChange > 0),]
  resNeg <- resSig[which(resSig$log2FoldChange < 0),]
  
  if(dir = "pos"){
    return(resPos)
  }
  
  if(dir = "neg"){
    return(resNeg)
  }
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



#------------------------------------------------ Initialize ------------------------------------------------#

#####################
# VARIABLE SETTINGS #
#####################
# Assign project variables
proj = 'LIMA'
projNum = '2-3.1.1'
sort = 'trmt'
outputDir = file.path("/Users/phanstiel3/Research/Data/Projects", proj, "rna/diff")

# Assign GO marts
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host= "www.ensembl.org")


################
# RUN SETTINGS #
################


##################
# COLOR SETTINGS #
##################
# Assign GO term colors
BPcol <- "dodgerblue2"
CCcol <- "firebrick2"
MFcol <- "grey"



################
# READ IN DATA #
################
# Read in config file
configPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_', projNum, '_samples.csv')
config <- read.csv(configPath, header=T)

# Read in count matrix
txiPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_', projNum, '_txi.rds')
txi <- readRDS(txiPath)




# Make long + abridged gene list from results file
longKey <- rownames(results)
abridgedKey <- nameAbridge(longKey)