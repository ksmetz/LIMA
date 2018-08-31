######################################################
# This is the script for analyzing count-matrix data #
#     from the RNAcluster pipeline on killdevil.     #
#----------------------------------------------------#
# It will read in the <Proj>_counts.csv files made   #
#   by txImporter.R through the RNAcluster.py        #
#   pipeline run on killdevil, that were then moved  #
#   to a local location in /Data/Projects/<Proj>/    #
#   RNA_Seq/proc/tximport.                           #
# It will then use DESeq2 to explore the data. You   #
#   can change the DESeq2 portions depending on the  #
#   figures or manipulations you want. When adding   #
#   project-specific details, try to save as         #
#   <Proj>Analyze.R.                                 #
######################################################

# Load any required libraries
library(DESeq2)
library(Sushi2)
library(colorspace)
library(biomaRt)
library(org.Hs.eg.db)
library(GO.db)
library(limma)
library(pheatmap)
library(ggplot2)
library(pcaExplorer)
library(ReportingTools)
library(factoextra)
library(NbClust)
library(pals)

#------------------------------------------------ Functions ------------------------------------------------#

# Function for circling (and optionally labeling) a list of genes on an MA plot
geneGroupPlot <- function(results, geneName, geneNumber)
{
  if (length(geneName) != length(geneNumber))
  {
    stop("Gene name and gene ID lists are not same length.")
  }
  for (n in seq(1, length(geneName)))
  {
    resSplit <- strsplit(rownames(results), "[.]")
    resAb <- unlist(resSplit)[2*(1:length(results[,1]))-1]
    bMean = results[which(resAb == geneNumber[n]),1]
    foldChange = results[which(resAb == geneNumber[n]),2]
    
    points(bMean, foldChange, pch=1)
    text(bMean, foldChange, labels=geneName[n], cex=.7, pos=3)
  }
}


# Function for spliting DESeq2 results into sig up/downreg genes
MAsplitter <- function(results)
{
  resSig <- results[which(results$padj <= .1),]
  resPos <<- resSig[which(resSig$log2FoldChange > 0),]
  resNeg <<- resSig[which(resSig$log2FoldChange < 0),]
}


# limma-based GO enrich function
GOenrich <- function(set,background)
{
  bkgd <- c(set,background)
  bkgd <- bkgd[!duplicated(bkgd)]
  
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host= "www.ensembl.org")
  
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
GOrunner <- function(results)
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
  posEnrich <<- GOenrich(posAb,backAb)
  negEnrich <<- GOenrich(negAb,backAb)
}


# Function for assigning colors to different GO groups
GOcolor <- function(list)
{
  list[list=="BP"] <- "dodgerblue2"
  list[list=="CC"] <- "firebrick2"
  list[list=="MF"] <- "grey"
  
  return(list)
}


# Function for extending the ENSEMBL name to the decimal version in dds, res, etc.
nameExtend <- function(shortName)
{
  longNames=list()
  for (name in shortName){
    longName = rownames(res[which(resAb==name),])
    longNames = append(longNames, longName)
  }
  return(longNames)
}


# Function for shortening the decimal ENSEMBL name from dds, res, etc. to the standard ENSEMBL gene ID
nameAbridge <- function(longName){
  split <- strsplit(longName, "[.]")
  shortName <- unlist(split)[2*(1:length(longName))-1]
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

# Settings
cmdLine = F
proj = 'LIMA'
sort = 'trmt'

sample.pal=parula(8)

cat("\n\n")
cat("------------\n")
cat("Initializing\n")
cat("------------\n")
cat("\n\n")


##############
# USER INPUT #
##############

if (cmdLine == T)
{
  
  # Command line inputs for (1) project name and (2) what sort of groups you are comparing [cell type or treatment]
  args <- commandArgs(trailingOnly=T)
  proj <- (args[1])
  sort <- (args[2])
  
}


# Errors for no input given
if (cmdLine == T && length(args)!=2) 
{
  stop("Please provide (1) a project name [see SampleInfo.csv] and (2) a sorting parameter [cell or trmt]")
}


# Restate user input
projStmt <- paste("Project selected:", proj, sep=" ")
cat(projStmt)

cat("\n\n")

sortStmt <- paste("Comparing by:", sort, sep=" ")
cat(sortStmt)

cat ("\n\n\n")




#########################
# MAKE OUTPUT DIRECTORY #
#########################


cat("Making output directory ...")
cat("\n\n")

# Make directories for output
outputDir = file.path("/Users/phanstiel3/Research/Data/Projects", proj, "rna/diff")
dirCmd = paste("mkdir", outputDir, sep=" ")
system(dirCmd)




##########################
# READ IN TXIMPORT FILES #
##########################


# Read in config file
configPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_2-3.1.1_samples.csv')

configStmt <- paste("Accessing config file at: ", configPath, sep=" ")
cat(configStmt)
cat("\n\n")

config <- read.csv(configPath, header=T)


# Read in count matrix
txiPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_2-3.1.1_txi.rds')

txiStmt <- paste("Accessing txi file at: ", txiPath, sep=" ")
cat(txiStmt)
cat("\n\n")

txi <- readRDS(txiPath)


#-------------------------------------------------- Read Info ----------------------------------------------------#

#pdf("~/Desktop/LIMA.pdf", height=8, width=8)
#par(mfrow=c(1,2))

boxplot(txi$counts,outline=F,xaxt='n', main="Counts")
axis(side=1,las=2,at=1:16,labels=colnames(txi$counts))
bp = barplot(colSums(txi$counts),xaxt='n', main="Summed counts per sample")
axis(side=1,las=2,at=bp,labels=colnames(txi$counts))
hist(rowSums(txi$counts), main="Summed counts per gene")
rug(rowSums(txi$counts))
hist(txi$counts[,3], main="Counts per gene for 60m")
rug(txi$counts[,3])


#-------------------------------------------------- Run DESeq2 ----------------------------------------------------#

# Create data frame with sample names and whatever condition you want to distinguish them by (config$Condition, config$Cell.Type, etc.)
colData <- data.frame(condition = factor(config$Condition), rep = factor(config$BioRep), row.names = config$Name)

# Create DESeq Data Set from txImport and save results to object
ddsTxI <- DESeqDataSetFromTximport(txi, colData=colData, design =~ rep + condition)

dds <- DESeq(ddsTxI)
dds$condition <- factor(dds$condition, levels=c("0","30", "60", '90', '120', '240', '360', '1440'))
dds$rep <- factor(dds$rep, levels=c("2", "3"))
res <- results(dds, contrast=c("condition", "0", "1440"))




resSplit <- strsplit(rownames(res), "[.]")
resAb <- unlist(resSplit)[2*(1:length(res[,1]))-1]



#-------------------------------------------------- MA Plots ----------------------------------------------------#

par(mfrow=c(1,1))

# Make MA Plot
# res (contrast) = condition <red> vs <blue>
MAsplitter(res)
plot(res$baseMean, res$log2FoldChange, log = "x", xlim=c(.15, 1000000), ylim=c(-15,15), cex=.25, col="gray38", pch=19, xlab="", ylab="")
points(resPos$baseMean, resPos$log2FoldChange, pch=19, cex=.25, col="firebrick2")
points(resNeg$baseMean, resNeg$log2FoldChange, pch=19, cex=.25, col="dodgerblue2")
lines(c(.08,5000000), c(0, 0), col=alpha("indianred1", .5), lwd=3)
title(main="MA-Plot: 0000 vs 1440")
title(xlab = "mean of normalized counts", line=2)
title(ylab = "log2 fold change", line=2.5)

geneNames <- ap1$gene
geneIDs <- ap1$ID
geneGroupPlot(res, geneNames, geneIDs)

legtxt = c(paste0("Pos (n=",nrow(resPos), ")"), paste0("Neg (n=",nrow(resNeg), ")"))
legend("topright", legend=legtxt, pch=19, col=c("firebrick2", "dodgerblue2"), cex=.8)




# ------------------------------------------------- PCA ---------------------------------------------- #

par(mfrow=c(1,2))

# Transformations
ntd <- normTransform(dds)
rld <- rlog(dds, blind=FALSE)

plotPCA(ntd, intgroup="condition")+labs(title="Normalized Counts Transformation PCA")+scale_color_manual(values=sample.pal)
plotPCA(rld, intgroup="condition")+labs(title="RLD")+scale_color_manual(values=sample.pal)

# pdf("~/Desktop/LIMA_PCA.pdf")
# plotPCA(rld, intgroup="condition")+scale_color_manual(values=sample.pal)
# dev.off()


# ------------------------------------------------- Plot Counts ---------------------------------------------- #

par(mfrow=c(3,1))

#pdf("~/Desktop/LIMA_plotCounts.pdf", height=6,width=8)

# Read in genes of interest
ap1 <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(ap1) <- c("gene", "ID")
for (n in 1:nrow(ap1)){
  ap1$long.ID[n] = nameExtend(ap1$ID[n])[[1]]
}

# Plot genes of interest
mround <- function(x,base){ 
  base*round(x/base) 
}

n=nrow(ap1)
par(mar=c(3,2,1,2))
par(mfrow=c(3,mround(n,3)/3))
for(n in 1:nrow(ap1)){
  plotCounts(dds, ap1$long.ID[n], main=ap1$gene[n], pch=16, col=rep(sample.pal,times=2),
             intgroup="condition")
}

#dev.off()

# # IL1B
# plotCounts(dds, as.character(nameExtend("ENSG00000125538")), main="IL1B", pch=16, col=rep(c("red","orange","yellow",'green','cyan','blue','purple', 'magenta'),times=2),
#            intgroup="condition")
# # FOS
# plotCounts(dds, as.character(nameExtend("ENSG00000170345")), main="FOS", pch=16, col=rep(c("red","orange","yellow",'green','cyan','blue','purple', 'magenta'),times=2))
# 
# #col=c("red","orange","yellow",'green','cyan','blue','purple')
# # JUN
# plotCounts(dds, as.character(nameExtend("ENSG00000177606")), main="JUN", pch=16, col=rep(c("red","orange","yellow",'green','cyan','blue','purple', 'magenta'),times=2))







# ------------------------------------------------- DESeq2 Clustering Heatmaps ---------------------------------------------- #
par(mfrow=c(1,2))

#
# Count matrix heatmaps
#


select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing=TRUE)[1:1000]
allgenes <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
df <- as.data.frame(colData(dds)[,"condition"])
rownames(df) <- colnames(ntd)
colnames(df) <- "condition"

pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames = FALSE, cluster_cols = TRUE, annotation_col=df)



#
# Median-normalized: Top n genes
#

# Normalize by the mean of each row
dds.norm = assay(ntd)
dds.normMean = (dds.norm-rowMeans(dds.norm))

# remove rows with NA
dds.normMean = dds.normMean[complete.cases(dds.normMean),]

# Calculate variance per gene and reorder by it
variancevalues = apply(dds.normMean,1,var)
dds.normorder = order(variancevalues,decreasing=TRUE)[1:1000]
dds.normSelect = dds.normMean[dds.normorder,]

# Trucnate values
maxval = 3
dds.normSelect[which(dds.normSelect>maxval)] = maxval
dds.normSelect[which(dds.normSelect<(-maxval))] = (-maxval)

# Set colors
pos_breaks <- c(seq(0,maxval,length.out=100))
neg_breaks <- c(seq(-maxval,0, length.out=100))
neg_breaks <- neg_breaks[1:(length(neg_breaks)-1)]
mycolors <- c(colorRampPalette(colors=c("red", "white"))(length(neg_breaks)-1),
              colorRampPalette(colors=c("white", "blue"))(length(pos_breaks)))

# Plot heatmap
row.names(dds.normSelect) <- rep("", times=nrow(dds.normSelect))
colnames(dds.normSelect) <- c()
heatmap(dds.normSelect, col=mycolors, breaks=c(neg_breaks,pos_breaks),scale="none",Rowv = T,Colv = NA, 
        ColSideColors = rep(sample.pal,times=2),
        labCol=rep(c("0","30", "60", '90', '120', '240', '360', '1440'), times=2))


#dev.off()



#=======================================================================#
#=======================================================================#
#=======================================================================#


#=============================================#
#                                             #
#                    LRT                      #
#                                             #
#=============================================#

####################
# Make the LRT DDS #
####################

# ---------- Make DDS with LRT test ---------- #

dds.LRT <- dds

#dds.new$condition <- factor(c(0,30,60,90,2*60,4*60,6*60, 24*60))
#dds.new$rep <- factor(rep(c(2,3), each=8))
colData(dds.LRT)

design(dds.LRT) <- ~rep + condition

dds.LRT <- DESeq(dds.LRT, test="LRT", full=~rep + condition, reduced = ~rep)

# Plot genes of interest
mround <- function(x,base){ 
  base*round(x/base) 
}

n=nrow(ap1)
par(mfrow=c(3,mround(n,3)/3))
for(n in 1:nrow(ap1)){
  plotCounts(dds.LRT, ap1$long.ID[n], main=ap1$gene[n], pch=16, col=rep(sample.pal,times=2),
             intgroup="condition")
}

# Alt:
# dds.new <- DESeq(ddsTxI, test="LRT", full=~rep + condition, reduced = ~rep)

res.LRT.full <- results(dds.LRT, contrast=c("condition", "0","30"))
res.LRT.sub <- res.LRT.full[which(res.LRT.full$padj < .01),]

nrow(res.LRT.full)
nrow(res.LRT.sub)

# ---------- Generate + subset shrunken LFC for each comparison vs 0m ---------- #

res.030 <- lfcShrink(dds.LRT, coef="condition_30_vs_0", type="apeglm")
res.060 <- lfcShrink(dds.LRT, coef="condition_60_vs_0", type="apeglm")
res.090 <- lfcShrink(dds.LRT, coef="condition_90_vs_0", type="apeglm")
res.0120 <- lfcShrink(dds.LRT, coef="condition_120_vs_0", type="apeglm")
res.0240 <- lfcShrink(dds.LRT, coef="condition_240_vs_0", type="apeglm")
res.0360 <- lfcShrink(dds.LRT, coef="condition_360_vs_0", type="apeglm")
res.01440 <- lfcShrink(dds.LRT, coef="condition_1440_vs_0", type="apeglm")


# subset significant genes
p.thr = .01
res.030 <-  res.030[which(res.030$padj < p.thr),]
res.060 <-  res.060[which(res.060$padj < p.thr),]
res.090 <-  res.090[which(res.090$padj < p.thr),]
res.0120 <-  res.0120[which(res.0120$padj < p.thr),]
res.0240 <-  res.0240[which(res.0240$padj < p.thr),]
res.0360 <-  res.0360[which(res.0360$padj < p.thr),]
res.01440 <-  res.01440[which(res.01440$padj < p.thr),]


# subset genes changing over 2-fold
# Note: negative fold-changes indicate higher in second sample (i.e. 0 vs *___*)
lfc.thr = 1
res.030.diff <- res.030[which(abs(res.030$log2FoldChange) >= lfc.thr),]
res.060.diff <- res.060[which(abs(res.060$log2FoldChange) >= lfc.thr),]
res.090.diff <- res.090[which(abs(res.090$log2FoldChange) >= lfc.thr),]
res.0120.diff <- res.0120[which(abs(res.0120$log2FoldChange) >= lfc.thr),]
res.0240.diff <- res.0240[which(abs(res.0240$log2FoldChange) >= lfc.thr),]
res.0360.diff <- res.0360[which(abs(res.0360$log2FoldChange) >= lfc.thr),]
res.01440.diff <- res.01440[which(abs(res.01440$log2FoldChange) >= lfc.thr),]


# # Check out "MA plots" for each
# par(mfrow=c(7,1))
# plot(x=res.030.diff$baseMean, y=res.030.diff$log2FoldChange)
# plot(x=res.060.diff$baseMean, y=res.060.diff$log2FoldChange)
# plot(x=res.090.diff$baseMean, y=res.090.diff$log2FoldChange)
# plot(x=res.0120.diff$baseMean, y=res.0120.diff$log2FoldChange)
# plot(x=res.0240.diff$baseMean, y=res.0240.diff$log2FoldChange)
# plot(x=res.0360.diff$baseMean, y=res.0360.diff$log2FoldChange)
# plot(x=res.01440.diff$baseMean, y=res.01440.diff$log2FoldChange)


# ---------- Look into genes meeting LFC, pvalue threshold ---------- #

# count number of "differential" genes
nrow(res.030.diff)
nrow(res.060.diff)
nrow(res.090.diff)
nrow(res.0120.diff)
nrow(res.0240.diff)
nrow(res.0360.diff)
nrow(res.01440.diff)

# combine all genes into one, non-redundant list
geneList = c(rownames(res.030.diff), rownames(res.060.diff), rownames(res.090.diff), rownames(res.0120.diff), rownames(res.0240.diff), rownames(res.0360.diff), rownames(res.01440.diff))
geneList = rownames(res.LRT.full)[rownames(res.LRT.full) %in% geneList]

length(geneList)



# #####################
# # Interactive Views #
# #####################
# 
# # PCA explorer (featuring pca2go, which requires the names be shortened to regular ENSEMBL IDs)
# dds.rename = dds.LRT
# rownames(dds.rename) <- nameAbridge(rownames(dds.rename))
# pcaExplorer(dds=dds.rename)
# 
# # Reporting Tools (creates an html)
# des2Report <- HTMLReport(shortName='RNAseq_analysis_with_DESeq2', title = 'RNA-seq analysis of differential expression using DESeq2', reportDirectory = "./Desktop/reports")
# publish(dds.LRT, des2Report, pvalueCutoff=.01, annotation.db='org.Hm.eg.db', factor=colData(dds.LRT)$condition, reportDir="./Desktop/reports")
# finish(des2Report)



#=============================================#
#                                             #
#                    LFC                      #
#                                             #
#=============================================#


########################
# LFC Ordered Heatmaps #
########################

# ---------- Build matrix with LFC ---------- #

# Build a matrix with one row per gene, and one column per comparison
# For every gene that changes in at least one sample, include all LFC
geneMatrix = matrix(nrow=length(geneList), ncol=7)
geneMatrix[,1] = res.030$log2FoldChange[which(rownames(res.030) %in% geneList)]
geneMatrix[,2] = res.060$log2FoldChange[which(rownames(res.060) %in% geneList)]
geneMatrix[,3] = res.090$log2FoldChange[which(rownames(res.090) %in% geneList)]
geneMatrix[,4] = res.0120$log2FoldChange[which(rownames(res.0120) %in% geneList)]
geneMatrix[,5] = res.0240$log2FoldChange[which(rownames(res.0240) %in% geneList)]
geneMatrix[,6] = res.0360$log2FoldChange[which(rownames(res.0360) %in% geneList)]
geneMatrix[,7] = res.01440$log2FoldChange[which(rownames(res.01440) %in% geneList)]

rownames(geneMatrix) = geneList
colnames(geneMatrix) = c("30", "60", "90", "120", "240", "360", "1440")

head(geneMatrix)


# -------------------------------------------------------- #
# ---------- RUN FROM HERE TO GENERATE HEATMAPS ---------- #
# -------------------------------------------------------- #

# ---------- Order Gene Matrix ---------- #

# Write functions to apply to each row...
# value.function: find the number of time points where LFC > 2
value.function = function(matrix.row){length(which(abs(matrix.row[1:7]) >= lfc.thr))}
# max.function: find the time point where the LFC is greatest
max.function = function(matrix.row){which.max(matrix.row[1:7])}
# min.function: find the earliest time point where LFC > 2
min.function = function(matrix.row){min(which(abs(matrix.row[1:7]) >= lfc.thr))}

# Apply functions to each row, add a new column with it. Also add a column with baseMean counts (in case we want it?)
geneMatrix.order = cbind(geneMatrix, apply(geneMatrix, 1, value.function))
geneMatrix.order = cbind(geneMatrix.order, apply(geneMatrix.order, 1, max.function))
geneMatrix.order = cbind(geneMatrix.order, apply(geneMatrix.order, 1, min.function))
geneMatrix.order = cbind(geneMatrix.order, res.LRT.full$baseMean[which(rownames(res.LRT.full) %in% geneList)])

colnames(geneMatrix.order)[8:11] = c("tp.span", "max.tp", "first.tp", "base.means")

# reorder by these new rows in different ways
geneMatrix.order.8.9 = geneMatrix.order[order(geneMatrix.order[,8], geneMatrix.order[,9]),]
geneMatrix.order.8.10 = geneMatrix.order[order(geneMatrix.order[,8], geneMatrix.order[,10]),]
geneMatrix.order.9.8 = geneMatrix.order[order(geneMatrix.order[,9], geneMatrix.order[,8]),]
geneMatrix.order.9.10 = geneMatrix.order[order(geneMatrix.order[,9], geneMatrix.order[,10]),]
geneMatrix.order.10.8 = geneMatrix.order[order(geneMatrix.order[,10], geneMatrix.order[,8]),]
geneMatrix.order.10.9 = geneMatrix.order[order(geneMatrix.order[,10], geneMatrix.order[,9]),]

# Identify gaps between sections (using formula below)
#   tail(which(geneMatrix.order.8.9[,8] == 1))

# LFC = 4
# gaps8 = c(1132, 1502, 1940, 2104, 2257, 2363)
# gaps9 = c(873, 901, 929, 1069, 1346, 1639)
# gaps10 = c(144,289,473,680,1273,1541)

# LFC = 2
gaps8 = c(2624, 3798, 4945, 5308, 5560, 5772)
gaps9 = c(980, 1563, 2097, 2656, 3201, 3893)
gaps10 = c(89, 370, 725, 1350, 2999, 3723)


# ---------- Truncate Gene Matrix for plotting ---------- #

maxval = 5
# Trucnate values
geneMatrix.plot.8.9 = geneMatrix.order.8.9[,1:7]
geneMatrix.plot.8.9[which(geneMatrix.plot.8.9>maxval)] = maxval
geneMatrix.plot.8.9[which(geneMatrix.plot.8.9<(-maxval))] = (-maxval)
# Trucnate values
geneMatrix.plot.8.10 = geneMatrix.order.8.10[,1:7]
geneMatrix.plot.8.10[which(geneMatrix.plot.8.10>maxval)] = maxval
geneMatrix.plot.8.10[which(geneMatrix.plot.8.10<(-maxval))] = (-maxval)
# Trucnate values
geneMatrix.plot.9.8 = geneMatrix.order.9.8[,1:7]
geneMatrix.plot.9.8[which(geneMatrix.plot.9.8>maxval)] = maxval
geneMatrix.plot.9.8[which(geneMatrix.plot.9.8<(-maxval))] = (-maxval)
# Trucnate values
geneMatrix.plot.9.10 = geneMatrix.order.9.10[,1:7]
geneMatrix.plot.9.10[which(geneMatrix.plot.9.10>maxval)] = maxval
geneMatrix.plot.9.10[which(geneMatrix.plot.9.10<(-maxval))] = (-maxval)
# Trucnate values
geneMatrix.plot.10.8 = geneMatrix.order.10.8[,1:7]
geneMatrix.plot.10.8[which(geneMatrix.plot.10.8>maxval)] = maxval
geneMatrix.plot.10.8[which(geneMatrix.plot.10.8<(-maxval))] = (-maxval)
# Trucnate values
geneMatrix.plot.10.9 = geneMatrix.order.10.9[,1:7]
geneMatrix.plot.10.9[which(geneMatrix.plot.10.9>maxval)] = maxval
geneMatrix.plot.10.9[which(geneMatrix.plot.10.9<(-maxval))] = (-maxval)


# ---------- Plot heatmaps of ordered Gene Matrix ---------- #

# Read in table with genes of interest
ap1 <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(ap1) <- c("gene", "ID")
for (n in 1:nrow(ap1)){
  ap1$long.ID[n] = nameExtend(ap1$ID[n])[[1]]
}

# Plot heatmaps using Sushi2 hotmap function
#pdf("~/Desktop/LIMA_LFC_clustering.pdf")
par(mar=c(5,4,4,2))
par(mfrow=c(1,1))

Sushi2::hotmap(geneMatrix.plot.8.9, labrow=F, labcol=T, gaps=gaps8, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene)
mtext(side=3,line=1.0,font=2,text="Number of TP with LFC above threshold > TP of max LFC",cex=.75)

Sushi2::hotmap(geneMatrix.plot.8.10, labrow=F, labcol=T, gaps=gaps8, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene)
mtext(side=3,line=1.0,font=2,text="Number of TP with LFC above threshold > Earliest LFC above threshold",cex=.75)

Sushi2::hotmap(geneMatrix.plot.9.8, labrow=F, labcol=T, gaps=gaps9, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene)
mtext(side=3,line=1.0,font=2,text="TP of max LFC > Number of TP with LFC above threshold",cex=.75)

Sushi2::hotmap(geneMatrix.plot.9.10, labrow=F, labcol=T, gaps=gaps9, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene)
mtext(side=3,line=1.0,font=2,text="TP of max LFC > Earliest LFC above threshold",cex=.75)

Sushi2::hotmap(geneMatrix.plot.10.8, labrow=F, labcol=T, gaps=gaps10, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene)
mtext(side=3,line=1.0,font=2,text="Earliest LFC above threshold > Number of TP with LFC above threshold",cex=.75)

Sushi2::hotmap(geneMatrix.plot.10.9, labrow=F, labcol=T, gaps=gaps10, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene)
mtext(side=3,line=1.0,font=2,text="Earliest LFC above threshold > TP of max LFC",cex=.75)
#dev.off()




#####################
# LFC K-means plots #
#####################

# ---------- Prep data for clustering ---------- #

# I believe this is centering and scaling the data
geneMatrix.norm <- (geneMatrix - rowMeans(geneMatrix))/rowSds(geneMatrix + .5)
#geneMatrix.norm <- geneMatrix/rowSds(geneMatrix+.5)


# ---------- Determine optimal number of clusters ---------- #

# # Method 1: WSS
# fviz_nbclust(geneMatrix.norm, kmeans, method = "wss")
# 
# # Method 2: Silhouette
# fviz_nbclust(geneMatrix.norm, kmeans, method = "silhouette")
# 
# # Method 3: Nb clustering
# nb <- NbClust(geneMatrix.norm, distance = "euclidean", min.nc = 2,
#               max.nc = 10, method = "kmeans")
# fviz_nbclust(nb)


# ---------- Plot the clusters ---------- #
#pdf("~/Desktop/LIMA_LFC-kmeans.pdf")

# r.num = sample(1:1000, 1)
# r.num
# prev. seed = 5
set.seed(733)

# Set the number of clusters
k <- 8
k.color <- brewer.dark2(8)

# # Perform hierarchical clustering
# hc <- hclust(dist(geneMatrix.norm))
# cut <- cutree(hc, k=k)

# Or perform kmeans clustering
fit = kmeans(geneMatrix.norm,centers=k)
cut = fit$cluster

# make the clusters
par(mfrow=c(2,k/2))
for (i in 1:k) {
  
  # make empty plot
  plot(colMeans(geneMatrix.norm[cut==i,]), type="n",main=paste("n=",table(cut)[i],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative LFC")
  
  # and transparent grey lines
  apply(geneMatrix.norm[cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
  
  # add median line
  lines(colMeans(geneMatrix.norm[cut==i,]),col=k.color[i],lwd=2)
  
  # add axis
  axis(1, at=1:7, labels=c(".5", "1", "1.5", "2", "4", "6", "24"))
}


#######################
# LFC K-means heatmap #
#######################

# ---------- Ordering results by cluster ---------- #

# Sort mat and cut according to highest --> lowest overall change (max value in each row); mixed clusters still at this point
geneMatrix.norm.max.order <- geneMatrix.norm[order(rowMax(geneMatrix.norm), decreasing=T),]
cut.max.order <- cut[rownames(geneMatrix.norm.max.order)]

# Set the desired order of clusters, based on their timing
desired.order <- c(8,1,2,4,3,6,5,7)

# Reorganize mat according to the cluster order desired; now sorted first by cluster, then by max value per row (highest to lowest)
mat.order <- geneMatrix.norm.max.order[order(match(cut.max.order,desired.order)),]


# ---------- K-means heatmap ---------- #

# Trucnate values
mat.order.plot = mat.order
maxval = 3
mat.order.plot[which(mat.order.plot>maxval)] = maxval
mat.order.plot[which(mat.order.plot<(-maxval))] = (-maxval)

colnames(mat.order.plot) <- c("30", "60", '90', '120', '240', '360', '1440')

# Read in genes of interest
ap1 <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(ap1) <- c("gene", "ID")
for (n in 1:nrow(ap1)){
  ap1$long.ID[n] = nameExtend(ap1$ID[n])[[1]]
}

# Assign cluster colors
cluster.col.list = k.color[desired.order]
cluster.col = c(rep(cluster.col.list[1], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[1]))))),
                rep(cluster.col.list[2], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[2]))))),
                rep(cluster.col.list[3], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[3]))))),
                rep(cluster.col.list[4], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[4]))))),
                rep(cluster.col.list[5], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[5]))))),
                rep(cluster.col.list[6], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[6]))))),
                rep(cluster.col.list[7], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[7]))))),
                rep(cluster.col.list[8], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[8]))))))

# Add cluster colors to the plot
mat.sushi.order.plot = cbind(mat.order.plot, cluster.col)

# Assign gaps between clusters
gaps.k=c()
for (i in 1:(length(desired.order)-1)){
  gaps.k = c(gaps.k, (tail(which(mat.sushi.order.plot[,ncol(mat.sushi.order.plot)] == cluster.col.list[i]))[6]))
}

# Plot heatmap
par(mfrow=c(1,1))
par(mar=c(5,4,4,2))
Sushi2::hotmap(mat.order.plot, col=mycolors, labrow=F, labcol=T, gaps=gaps.k, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene, rowcolors=cluster.col)
mtext(side=3,line=1.0,font=2,text="K-means LFC clusters",cex=.75)
#dev.off()

# ##################
# # GO Enrichments #
# ##################
# 
# # Make GO plots
# #pdf("~/Desktop/LIMA_LRT_GOenr.pdf")
# par(mfrow=c(1,1))
# 
# GOrunner(res.030)
# GOval <- c(-log(posEnrich[5:1,5]),log(negEnrich[5:1,5]))
# GOhit <- c(posEnrich[5:1,1], negEnrich[5:1,1])
# GOcol <- c(posEnrich[5:1,2], negEnrich[5:1,2])
# GOcol <- GOcolor(GOcol)
# plot <- barplot(GOval, horiz=T, main="0 v 30 m GO Enrichments", col=c(rep(sample.pal[2], 5), rep(sample.pal[1], 5)), xlim=c(-150,150), xlab="log(P-value)")
# text(x=rep(c(2,-2),each=5), 
#      y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
#      pos=rep(c(4,2),each=5))
# 
# GOrunner(res.060)
# GOval <- c(-log(posEnrich[5:1,5]),log(negEnrich[5:1,5]))
# GOhit <- c(posEnrich[5:1,1], negEnrich[5:1,1])
# GOcol <- c(posEnrich[5:1,2], negEnrich[5:1,2])
# GOcol <- GOcolor(GOcol)
# plot <- barplot(GOval, horiz=T, main="0 v 60 m GO Enrichments", col=c(rep(sample.pal[3], 5), rep(sample.pal[1], 5)), xlim=c(-150,150), xlab="log(P-value)")
# text(x=rep(c(2,-2),each=5), 
#      y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
#      pos=rep(c(4,2),each=5))
# 
# GOrunner(res.090)
# GOval <- c(-log(posEnrich[5:1,5]),log(negEnrich[5:1,5]))
# GOhit <- c(posEnrich[5:1,1], negEnrich[5:1,1])
# GOcol <- c(posEnrich[5:1,2], negEnrich[5:1,2])
# GOcol <- GOcolor(GOcol)
# plot <- barplot(GOval, horiz=T, main="0 v 90 m GO Enrichments", col=c(rep(sample.pal[4], 5), rep(sample.pal[1], 5)), xlim=c(-150,150), xlab="log(P-value)")
# text(x=rep(c(2,-2),each=5), 
#      y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
#      pos=rep(c(4,2),each=5))
# 
# GOrunner(res.0120)
# GOval <- c(-log(posEnrich[5:1,5]),log(negEnrich[5:1,5]))
# GOhit <- c(posEnrich[5:1,1], negEnrich[5:1,1])
# GOcol <- c(posEnrich[5:1,2], negEnrich[5:1,2])
# GOcol <- GOcolor(GOcol)
# plot <- barplot(GOval, horiz=T, main="0 v 120 m GO Enrichments", col=c(rep(sample.pal[5], 5), rep(sample.pal[1], 5)), xlim=c(-150,150), xlab="log(P-value)")
# text(x=rep(c(2,-2),each=5), 
#      y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
#      pos=rep(c(4,2),each=5))
# 
# GOrunner(res.0240)
# GOval <- c(-log(posEnrich[5:1,5]),log(negEnrich[5:1,5]))
# GOhit <- c(posEnrich[5:1,1], negEnrich[5:1,1])
# GOcol <- c(posEnrich[5:1,2], negEnrich[5:1,2])
# GOcol <- GOcolor(GOcol)
# plot <- barplot(GOval, horiz=T, main="0 v 240 m GO Enrichments", col=c(rep(sample.pal[6], 5), rep(sample.pal[1], 5)), xlim=c(-150,150), xlab="log(P-value)")
# text(x=rep(c(2,-2),each=5), 
#      y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
#      pos=rep(c(4,2),each=5))
# 
# GOrunner(res.0360)
# GOval <- c(-log(posEnrich[5:1,5]),log(negEnrich[5:1,5]))
# GOhit <- c(posEnrich[5:1,1], negEnrich[5:1,1])
# GOcol <- c(posEnrich[5:1,2], negEnrich[5:1,2])
# GOcol <- GOcolor(GOcol)
# plot <- barplot(GOval, horiz=T, main="0 v 360 m GO Enrichments", col=c(rep(sample.pal[7], 5), rep(sample.pal[1], 5)), xlim=c(-150,150), xlab="log(P-value)")
# text(x=rep(c(2,-2),each=5), 
#      y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
#      pos=rep(c(4,2),each=5))
# 
# GOrunner(res.01440)
# GOval <- c(-log(posEnrich[5:1,5]),log(negEnrich[5:1,5]))
# GOhit <- c(posEnrich[5:1,1], negEnrich[5:1,1])
# GOcol <- c(posEnrich[5:1,2], negEnrich[5:1,2])
# GOcol <- GOcolor(GOcol)
# plot <- barplot(GOval, horiz=T, main="0 v 1440 m GO Enrichments", col=c(rep(sample.pal[8], 5), rep(sample.pal[1], 5)), xlim=c(-150,150), xlab="log(P-value)")
# text(x=rep(c(2,-2),each=5), 
#      y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
#      pos=rep(c(4,2),each=5))
# #dev.off()

#=============================================#
#                                             #
#                     VST                     #
#                                             #
#=============================================#

######################
# K-means clustering #
######################

# ---------- Prep data for clustering ---------- #

# Perform variance stabalizing transformation. I think this basically mean or median normalizes each row etc.
vsd.LRT <- vst(dds.LRT)

# select significant genes: KATIE CHANGED P-VAL CUTOFF from 0.05 to 0.01
#vsd.LRT.sub <- vsd.LRT[which(res.LRT.full$padj < .01),]

# select significant genes based on p-value and LFC cut-offs from LRT
vsd.LRT.sub <- vsd.LRT[which(rownames(vsd.LRT) %in% geneList),]

# make matrix
mat.vsd.LRT.sub <- assay(vsd.LRT.sub)

# I believe this is centering and scaling the data
# # Method 1
# mat.vsd.LRT.sub.norm <- (mat.vsd.LRT.sub - rowMeans(mat.vsd.LRT.sub))/rowSds(mat.vsd.LRT.sub + .5)

# Maybe try this way?
# # Method 2
# mat.vsd.LRT.sub.norm <- mat.vsd.LRT.sub-rowMin(mat.vsd.LRT.sub)

# Method 3 - No scaling
mat.vsd.LRT.sub.norm <- mat.vsd.LRT.sub

# ---------- Determine optimal number of clusters ---------- #


# # Method 1: WSS
# fviz_nbclust(geneMatrix.norm, kmeans, method = "wss")
# 
# # Method 2: Silhouette
# fviz_nbclust(geneMatrix.norm, kmeans, method = "silhouette")
# 
# # Method 3: Nb clustering
# nb <- NbClust(geneMatrix.norm, distance = "euclidean", min.nc = 2,
#               max.nc = 10, method = "kmeans")
# fviz_nbclust(nb)



# ---------- Plot the clusters ---------- #

# r.num = sample(1:1000, 1)
# r.num
set.seed(275)

# Set the number of clusters
k <- 8
k.color <- brewer.dark2(8)

# # Perform hierarchical clustering
# hc <- hclust(dist(mat.vsd.LRT.sub.norm))
# cut <- cutree(hc, k=k)

# Or perform kmeans clustering
fit = kmeans(mat.vsd.LRT.sub.norm,centers=k)
cut = fit$cluster

# Combine replicates
mat.plot <- combineReps(mat.vsd.LRT.sub.norm)

# make the clusters
#pdf("~/Desktop/LIMA_LRT_counts_clustering.pdf")
par(mfrow=c(2,k/2))
for (i in 1:k) {
  
  # make empty plot
  plot(colMeans(mat.plot[cut==i,]), type="n",main=paste("n=",table(cut)[i],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative expression")
  
  # and transparent grey lines
  apply(mat.plot[cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
  
  # add median line
  lines(colMeans(mat.plot[cut==i,]),col=k.color[i],lwd=2)
  
  # add axis
  axis(1, at=1:8, labels=c("0",".5", "1", "1.5", "2", "4", "6", "24"))
}


# ---------- Ordering results by cluster ---------- #

# Sort mat and cut according to highest --> lowest overall change (max value in each row); mixed clusters still at this point
mat.max.order <- mat.plot[order(rowMax(mat.plot), decreasing=T),]
cut.max.order <- cut[rownames(mat.max.order)]

# Set the desired order of clusters, based on their timing
desired.order <- c(8,1,2,3,5,6,4,7)

# Reorganize mat according to the cluster order desired; now sorted first by cluster, then by max value per row (highest to lowest)
mat.order <- mat.max.order[order(match(cut.max.order,desired.order)),]



# ---------- K-means heatmap ---------- #

# Trucnate values
mat.order.plot = mat.order
# maxval = 3
# mat.order.plot[which(mat.order.plot>maxval)] = maxval
# mat.order.plot[which(mat.order.plot<(-maxval))] = (-maxval)

#colnames(mat.order.plot) <- c("0","30", "60", '90', '120', '240', '360', '1440')
colnames(mat.order.plot) <- c("0","0.5", "1", '1.5', '2', '4', '6', '24')

# Read in genes of interest
ap1 <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(ap1) <- c("gene", "ID")
for (n in 1:nrow(ap1)){
  ap1$long.ID[n] = nameExtend(ap1$ID[n])[[1]]
}
ap1$color = c(rep("tan2", times=7),rep("steelblue", times=3), "black")


# Assign cluster colors
cluster.col.list = k.color[desired.order]
cluster.col = c(rep(cluster.col.list[1], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[1]))))),
                rep(cluster.col.list[2], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[2]))))),
                rep(cluster.col.list[3], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[3]))))),
                rep(cluster.col.list[4], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[4]))))),
                rep(cluster.col.list[5], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[5]))))),
                rep(cluster.col.list[6], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[6]))))),
                rep(cluster.col.list[7], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[7]))))),
                rep(cluster.col.list[8], times=length(which(rownames(mat.order) %in% names(which(cut.max.order == desired.order[8]))))))

# Add cluster colors to the plot
mat.sushi.order.plot = cbind(mat.order.plot, cluster.col)

# Assign gaps between clusters
gaps.k=c()
for (i in 1:(length(desired.order)-1)){
  gaps.k = c(gaps.k, (tail(which(mat.sushi.order.plot[,ncol(mat.sushi.order.plot)] == cluster.col.list[i]))[6]))
}

mycolors <- c(colorRampPalette(colors=c("firebrick2", "black"))(length(neg_breaks)-1),
              colorRampPalette(colors=c("black", "dodgerblue2"))(length(pos_breaks)))

# Plot heatmap
par(mfrow=c(1,1))
par(mar=c(5,4,2,5))
Sushi2::hotmap(mat.order.plot, col=mycolors, labrow=F, labcol=T, gaps=gaps.k, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene, selectylabs.col = ap1$color)
#, rowcolors=cluster.col

Sushi::addlegend(c(-3,3), palette=colorRampPalette(colors=c("firebrick2", "black","dodgerblue2")), title="Normalized Transcript Counts", bottominset=.5, xoffset=.11, title.offset = .07)
#mtext(side=3,line=1.0,font=2,text="K-means counts clusters",cex=.75)
#dev.off()


geneMatrix.norm[which(rownames(geneMatrix.norm) %in% ap1$long.ID),]
mat.order.plot[which(rownames(mat.order.plot) %in% ap1$long.ID),]





#=======================================================================#
#=======================================================================#
#=======================================================================#

######################
# K-Means Clustering #
######################

#---------- Load Libraries ---------- #
library(splines)
library(matrixStats)

# Get ensembl ids for our favorite genes
FOS = "ENSG00000170345"
IL1B = "ENSG00000125538"

dds.k <- dds

#---------- Read in data -----------#
dds.k$time <- c(0,30,60,90,2*60,4*60,6*60, 24*60)
dds.k$rep <- rep(c(2,3), each=8)
colData(dds.k)

#---------- Perform Differential Analysis ---------- #


# spline fit across time with two degrees of freedom.  What is the tilda all about?
design(dds.k) <- ~ns(rep + time, df=4)

# saving degrees of freedom for later.  I tried to use this variable above but got an errror
dgr <- 4

# Estimate fize factors, estimate dispersion, negative binomical GLM fitting. How can we understand reduced?
dds.k <- DESeq(dds.k, test="LRT", reduced=~rep)

# perform differental analysis
res.k <- results(dds.k)

# I think this extracts equations (coefficients) for each gene.
betas <- coef(dds.k)

# print out number of significant genes
sum(res.k$padj < .01, na.rm=TRUE)


# check for our favs
res.k$padj[grep(FOS,rownames(res.k))]
res.k$padj[grep(IL1B,rownames(res.k))]


#---------- Plot 16 significant genes (across range of p values) ---------- #
par(mfrow=c(4,4), mar=c(4,4,1,1))

for (i in (1+100*0:15)) {
 
 # get the index of interest
 idx <- order(res.k$pvalue)[i]
 
 # get normailized counts for all genes
 ncts <- counts(dds.k, normalized=TRUE)
 
 # get time points 
 time <- dds.k$time
 
 # fit the data
 design <- model.matrix(~ns(time, df=dgr))
 fitted <- betas %*% t(design)
 
 # plot the actual normalized data
 plot(log2(ncts[idx,]+.5) ~ time, pch=16, cex=1.5, xlab="time",
      ylab="log2 norm. counts", main=paste0(i,"-th sig gene"))
 
 # add the fitted points
 points(fitted[idx,] ~ time, col="blue")
 
 # draw a line throug the fitted points
 fit <- lm(fitted[idx,] ~ ns(time, df=dgr))
 df <- data.frame(time=seq(0,360))
 curve <- predict(fit, newdata=df)
 lines(df$time, curve, col="blue")
}

# Perform variance stabalizing transformation. I think this basically mean or median normalizes each row etc.
vsd <- vst(dds.k)

# select significant genes: KATIE CHANGED P-VAL CUTOFF from 0.05 to 0.01
vsd.sub <- vsd[which(res.k$padj < .01),]

# make matrix
mat.vsd.sub <- assay(vsd.sub)

# plot the k different clusters
par(mfrow=c(1,1))
plot(rowMeans(mat.vsd.sub), rowSds(mat.vsd.sub))
mat <- (mat.vsd.sub - rowMeans(mat.vsd.sub))/rowSds(mat.vsd.sub + .5)

# # determine k
# fviz_nbclust(mat, kmeans, method = "wss")
# fviz_nbclust(mat, kmeans, method = "silhouette")
# 
# nb <- NbClust(mat, distance = "euclidean", min.nc = 2,
#               max.nc = 10, method = "kmeans")
# fviz_nbclust(nb)

# make the clusters
hc <- hclust(dist(mat))
k <- 8
cut <- cutree(hc, k=k)
par(mfrow=c(2,k/2))
for (i in 1:k) {
 plot(colMeans(mat[cut==i,]), type="b",main=paste("n=",table(cut)[i],sep=""), xaxt="n", xlab="Hours", ylab="Relative expression")
 axis(1, at=1:7, labels=c("0",".5", "1", "1.5", "2", "4", "6"))
}

# which clusters are our favorite genes in?
cut[grep(FOS,rownames(mat))]
cut[grep(IL1B,rownames(mat))]

# interestingly the normalized counts and variance stablizing transformed counts tell a different story regarding when upregulation really occurs for IL1B
geneofinterest = IL1B
idx = grep(geneofinterest,rownames(ncts))
plot(ncts[idx,],type="b",col="dodgerblue2",pch=19)
plot(mat[grep(geneofinterest,rownames(mat)),],type="b",col="firebrick2",pch=19)

ap1 <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(ap1) <- c("gene", "ID")
for (n in 1:nrow(ap1)){
 ap1$cluster[n] = cut[grep(ap1[n,2],rownames(mat))]
}



###############################
# Ordering results by cluster #
###############################

# KATIE: Sort mat AND cut by maximum change value; reorganize based on group (order of choice based on timing)

# Sort mat and cut according to highest --> lowest overall change (max value in each row); mixed clusters still at this point
mat.max.order <- mat[order(rowMax(mat), decreasing=T),]
cut.max.order <- cut[rownames(mat.max.order)]

# Set the desired order of clusters, based on their timing
desired.order <- c(7,4,2,8,6,3,1,5)

# Reorganize mat according to the cluster order desired; now sorted first by cluster, then by max value per row (highest to lowest)
mat.order <- mat.max.order[order(match(cut.max.order,desired.order)),]

#saveRDS(mat.order, file="~/Desktop/ordered.matrix.RDS")
#saveRDS(cut.max.order, file="~/Desktop/cluster.id.RDS")

#test <-readRDS("~/Desktop/ordered.matrix.RDS")


###################
# K-means heatmap #
###################

# Trucnate values
mat.order.plot = mat.order
maxval = 3
mat.order.plot[which(mat.order.plot>maxval)] = maxval
mat.order.plot[which(mat.order.plot<(-maxval))] = (-maxval)

colnames(mat.order.plot) <- rep(c("0","30", "60", '90', '120', '240', '360', '1440'), times=2)

# # Set colors (heatmap function only)
# pos_breaks <- c(seq(0,maxval,length.out=100))
# neg_breaks <- c(seq(-maxval,0, length.out=100))
# neg_breaks <- neg_breaks[1:(length(neg_breaks)-1)]
# 
# mycolors <- c(colorRampPalette(colors=c("red", "black"))(length(neg_breaks)-1),
#               colorRampPalette(colors=c("black", "green"))(length(pos_breaks)))
# 
# # Remove names
# row.names(mat.order.plot) <- rep("", times=nrow(mat.order.plot))
# colnames(mat.order.plot) <- c()
# 
# # Plot heatmap
# heatmap(mat.order.plot, col=mycolors, breaks=c(neg_breaks,pos_breaks),scale="none",Rowv = NA,Colv = NA, 
#         ColSideColors = rep(c("red","orange","yellow",'green','cyan','blue','purple', 'magenta'),times=2),
#         labCol=rep(c("0","30", "60", '90', '120', '240', '360', '1440'), times=2))


ap1 <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(ap1) <- c("gene", "ID")
for (n in 1:nrow(ap1)){
  ap1$long.ID[n] = nameExtend(ap1$ID[n])[[1]]
}

# Assign cluster colors
cluster.col = c(rep("red", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 7))))),
                 rep("orange", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 4))))),
                 rep("yellow", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 2))))),
                 rep("green", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 8))))),
                 rep("cyan", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 6))))),
                 rep("blue", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 3))))),
                 rep("purple", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 1))))),
                 rep("magenta", times=length(which(rownames(mat.order) %in% names(which(cut.max.order == 5))))))

# Add cluster colors to the plot
mat.sushi.order.plot = cbind(mat.order.plot, cluster.col)

# Assign gaps between clusters
gaps.k = c(701, 2764, 4625, 5187, 6138, 7416, 9861, 11676)

# # Temporary fix for gene labeling: Swap ENSEMBL ID rownames for gene names for genes of interest
# for (n in rownames(mat.order.plot)[which(rownames(mat.order.plot) %in% ap1$long.ID)]){
#   rownames(mat.order.plot)[rownames(mat.order.plot) == n] <- toString(ap1$gene[which(n == ap1$long.ID)])
# }

# Plot heatmap
par(mfrow=c(1,1))
par(mar=c(5,4,4,2))
Sushi2::hotmap(mat.order.plot, col=mycolors, labrow=F, labcol=F, gaps=gaps.k, selectylabs=ap1$long.ID, selectylabs.label = ap1$gene, rowcolors=cluster.col)





#------------------------------------------------ Excel Overview ------------------------------------------------#

# Generic Info
# ============

# Make list of all genes based on VSD Results (should be sorted by ENSEMBL ID #)
ENSEMBLnames <- nameAbridge(rownames(mat.vsd.LRT.sub))

# Make matching list of all gene names, using BioMart to retrieve them... 
NamesDb <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', values=ENSEMBLnames, mart=mart)
# which(duplicated(MGInamesDb$ensembl_gene_id)==TRUE) --> 2657, 12757, 42071 are duplicates and need to be removed
# MGInamesDb <- MGInamesDb[-c(2657, 12757, 42071),]

# TEMPORARY FIX: Remove rows that didn't have ensembl_gene_id's fetched with BioMart, so for loop can work
ENSEMBLnames = ENSEMBLnames[(ENSEMBLnames %in% NamesDb$ensembl_gene_id)]

# ... and then a for loop to order them according to the ENSEMBL name list (BioMart comes out a different order than results)
HGNC_symbols <- rep("", times=length(ENSEMBLnames))
for (n in seq(1, length(ENSEMBLnames))){
  HGNC_symbols[n] = NamesDb$hgnc_symbol[which(NamesDb$ensembl_gene_id == ENSEMBLnames[n])]
}

# Make a list indicating which genes are of interest (called ap1 above)
AP1ind <- rep("", times=length(ENSEMBLnames))
for (n in seq(1, length(ENSEMBLnames))){
  if (ENSEMBLnames[n] %in% ap1$ID){
    AP1ind[n] = "Yes"
  } else{
    AP1ind[n] = "No"
  }
}

# Make lists of all the raw count data (straight from dds because dds and results should be in same order)
ENSEMBLnames.long = nameExtend(ENSEMBLnames)

count0000_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),1]
count0000_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),9]

count0030_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),2]
count0030_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),10]

count0060_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),3]
count0060_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),11]

count0090_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),4]
count0090_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),12]

count0120_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),5]
count0120_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),13]

count0240_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),6]
count0240_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),14]

count0360_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),7]
count0360_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),15]

count01440_1 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),8]
count01440_2 <- assay(dds)[which(rownames(assay(dds)) %in% ENSEMBLnames.long),16]



# Results-specific Info
# =====================

geneWrite(res=res.030, geneList=ENSEMBLnames, symbolList=HGNC_symbols, count1_1 = count0000_1, count1_2 = count0000_2, 
          count2_1 = count0030_1, count2_2 = count0030_2, sampleName1 = "0", sampleName2 = "30", filePath = "~/Desktop/030data.csv")

geneWrite(res=res.060, geneList=ENSEMBLnames, symbolList=HGNC_symbols, count1_1 = count0000_1, count1_2 = count0000_2, 
          count2_1 = count0060_1, count2_2 = count0060_2, sampleName1 = "0", sampleName2 = "60", filePath = "~/Desktop/060data.csv")

geneWrite(res=res.090, geneList=ENSEMBLnames, symbolList=HGNC_symbols, count1_1 = count0000_1, count1_2 = count0000_2, 
          count2_1 = count0090_1, count2_2 = count0090_2, sampleName1 = "0", sampleName2 = "90", filePath = "~/Desktop/090data.csv")

geneWrite(res=res.0120, geneList=ENSEMBLnames, symbolList=HGNC_symbols, count1_1 = count0000_1, count1_2 = count0000_2, 
          count2_1 = count0120_1, count2_2 = count0120_2, sampleName1 = "0", sampleName2 = "120", filePath = "~/Desktop/0120data.csv")

geneWrite(res=res.0240, geneList=ENSEMBLnames, symbolList=HGNC_symbols, count1_1 = count0000_1, count1_2 = count0000_2, 
          count2_1 = count0240_1, count2_2 = count0240_2, sampleName1 = "0", sampleName2 = "240", filePath = "~/Desktop/0240data.csv")

geneWrite(res=res.0360, geneList=ENSEMBLnames, symbolList=HGNC_symbols, count1_1 = count0000_1, count1_2 = count0000_2, 
          count2_1 = count0360_1, count2_2 = count0360_2, sampleName1 = "0", sampleName2 = "360", filePath = "~/Desktop/0360data.csv")

geneWrite(res=res.01440, geneList=ENSEMBLnames, symbolList=HGNC_symbols, count1_1 = count0000_1, count1_2 = count0000_2, 
          count2_1 = count01440_1, count2_2 = count01440_2, sampleName1 = "0", sampleName2 = "1440", filePath = "~/Desktop/01440data.csv")

# Save LFC matrix to csv as well
geneMatrix.order.sub = geneMatrix.order[which(rownames(geneMatrix.order) %in% ENSEMBLnames.long),]
LFCdata <- data.frame(HGNC_symbols, geneMatrix.order.sub, AP1ind, 
                   stringsAsFactors = FALSE)
write.csv(LFCdata, file="~/Desktop/LFCdata.csv", row.names = T)


