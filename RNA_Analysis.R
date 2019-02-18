######################################################
# This is the script for analyzing count-matrix data #
#     from the RNAcluster pipeline on longleaf.      #
######################################################


#############################################################################################################
#                                                                                                           #
#                                            I N I T I A L I Z E                                            #
#                                                                                                           #
#############################################################################################################

#-------------------------------------------------------------------------------------#
#                                       Libraries                                     #
#-------------------------------------------------------------------------------------#

# Load any required libraries
library(limma)
library(biomaRt)
library(org.Hs.eg.db)
library(GO.db)
library(DESeq2)
library(Sushi2)
library(ggplot2)
library(pcaExplorer)
library(ReportingTools)
library(NbClust)
library(pals)
library(readr)


#------------------------------------------------------------------------------------#
#                                       Settings                                     #
#------------------------------------------------------------------------------------#

#                    ################                    #
#--------------------# RUN SETTINGS #--------------------#
#                    ################                    #

# Make PDFs or run in Rstudio only?
makePDF = F

# Use directory specifically for today's figures?
newDir = T

# Read in from raw data?
readRaw = T

# Recreate dds and results files?
runDESeq = T
runShrink = T

# Make which plots/files?
txiPlot = F
PCAint = F
PCAplot = F
countPlot = F
MAplot = F
GOplot = F
LFCclusterOpt = F
LFCkclustPlot = F
LFCkclustHeatPlot = F
LFCorderPlot = F
countClusterOpt = F
countKclustPlot = F
countKclustHeatPlot = F
tsvFull = F
tsvSub = F


#                   ##################                   #
#-------------------# COLOR SETTINGS #-------------------#
#                   ##################                   #

# Assign sample colors for PCA, plot counts
sample.pal=parula(8)

# Assign GO term colors
BPcol <- "dodgerblue2"
CCcol <- "firebrick2"
MFcol <- "grey"

# Assign colors for clusters
k.colors <- brewer.dark2(8)

# Assign colors for heatmaps
heatmap.min = "steelblue1"
heatmap.mid = "black"
heatmap.max = "gold"



#                  #####################                  #
#------------------# VARIABLE SETTINGS #------------------# 
#                  #####################                  #

# Assign project variables
proj = 'LIMA'
projNum = '2-3.1.1'
sort = 'trmt'

transformType = "var"

# Assign cutoffs for subsetting
p.thr = .01
lfc.thr = 1

# Set default output directory
outputDir = file.path("/Users/phanstiel3/Desktop")

# Assign + make output directory for PDFs, CSVs, etc.
if(newDir == T & makePDF == T){
  
  # Assign new output dir
  outputName = paste(proj, projNum, Sys.Date(), sep="_")
  outputDir = file.path("/Users/phanstiel3/Desktop", outputName)
  #outputDir = file.path("/Users/phanstiel3/Research/Data/Projects", proj, "rna/diff")
  
  # Make dir if necessary
  dirCmd = paste("mkdir", outputDir, sep=" ")
  system(dirCmd)
  
}




#--------------------------------------------------------------------------------------#
#                                      Functions                                       #
#--------------------------------------------------------------------------------------#

#                       #########                       #
#-----------------------# BASIC #-----------------------#
#                       #########                       #

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


# Function for combining replicates in a 16x(many genes) matrix; r1/r1/r1/r2/r2/r2 structure
combineReps <- function(matrix, new.colnames)
{
  newMatrix = c()
  for (i in seq(1,ncol(matrix)/2))
  {
    tempMatrix = matrix(data=c(matrix[,i], matrix[,i+ncol(matrix)/2]), ncol=2, byrow=F)
    newCol = apply(tempMatrix,1,median)
    newMatrix  = cbind(newMatrix,newCol)
  }
  rownames(newMatrix) = rownames(matrix)
  
  if(!is.null(new.colnames)){
    colnames(newMatrix) = new.colnames
  }
  
  return(newMatrix)
}


# Function for coloring heatmap based on breaks
makeBreaks <- function(maxval, num, min.col, mid.col, max.col)
{
  
  # Make breaks according to absolute max value (negative needs one less)
  pos.breaks = c(seq(0, maxval, length.out=num))
  neg.breaks = c(seq(-maxval, 0, length.out=num))
  neg.breaks = neg.breaks[1:(length(neg.breaks)-1)]
  
  # Make color ramp (length = 2*num)
  colors = c(colorRampPalette(colors=c(min.col, mid.col))(length(neg.breaks)-1),
             colorRampPalette(colors=c(mid.col, max.col))(length(pos.breaks)))
  
}



#                      ############                      #
#----------------------# GOI PLOT #----------------------#
#                      ############                      #

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


#                      ############                      #
#----------------------# MA PLOTS #----------------------#
#                      ############                      #
# Function for plotting all MA-plots
MAplotter <- function(results, col1="dodgerblue2", col2="firebrick2", title1="Neg", title2="Pos", geneNameIDList=NULL, pval = .1)
{
  # Find significantly pos/neg genes
  resSig <- results[which(results$padj <= pval),]
  resPos <- resSig[which(resSig$log2FoldChange > 0),]
  resNeg <- resSig[which(resSig$log2FoldChange < 0),]
  
  # Set graph bounds
  min.LFC = mround(min(results[,2][!is.na(results[,2])]), 5, "floor")
  max.LFC = mround(max(results[,2][!is.na(results[,2])]),5, "ceiling")
  
  max.counts = mround(max(results[,2][!is.na(results[,2])]), 1000000, "ceiling")
  
  # Plot all points + highlight significant points
  plot(results[,1], results[,2], log="x", xlim=c(.15, max.counts), ylim=c(min.LFC, max.LFC), cex=.25, col="gray38", pch=19, xlab="", ylab="")
  points(resNeg[,1], resNeg[,2], pch=19, cex=.25, col=col1)
  points(resPos[,1], resPos[,2], pch=19, cex=.25, col=col2)
  lines(c(.15, max.counts), c(0,0), col=alpha("indianred1", .5), lwd=3)
  
  # Add titles
  main.title = paste("MA-Plot:", title1, "vs", title2, sep=" ")
  title(main=main.title)
  title(xlab="mean of normalized counts", line=2)
  title(ylab="log2 fold change", line=2.5)
  
  # Circle genes of interest
  if(!is.null(geneNameIDList)){
    geneNames <- geneNameIDList[,1]
    geneIDs <- geneNameIDList[,2]
    geneGroupPlot(results, geneNames, geneIDs)
  }
  
  # Add legend
  legtxt = c(paste0(title1," (n=",nrow(resPos), ")"), paste0(title2," (n=",nrow(resNeg), ")"))
  legend("topright", legend=legtxt, pch=19, col=c(col1, col2), cex=.8)
}  


#                   ##################                   #
#-------------------# GO ENRICHMENTS #-------------------#
#                   ##################                   #
# limma-based GO enrich function
GOenrich <- function(set,bkgd, ensembl.entrez, speciesName = "Hs")
{
  set_entrez <- ensembl.entrez[which(ensembl.entrez[,1] %in% set),2]
  bkgd_entrez <- ensembl.entrez[which(ensembl.entrez[,1] %in% bkgd),2]
  
  set_entrez = set_entrez[complete.cases(set_entrez)]
  bkgd_entrez = bkgd_entrez[complete.cases(bkgd_entrez)]
  
  enriched = goana(de=set_entrez, universe=bkgd_entrez,species=speciesName)
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

# Function for assigning colors to different GO terms
GOcolor <- function(list)
{
  list[list=="BP"] <- BPcol
  list[list=="CC"] <- CCcol
  list[list=="MF"] <- MFcol
  
  return(list)
}

# Function for running GO enrichments for up/down regulated genes and plotting the top n
GOplotter <- function(results, pval = .05, LFC = 0, n=5, color1=NULL, color2=NULL, title1="Pos", title2="Neg")
{
  # Extract genes enriched in either sample (neg: higher in first [het]; pos: higher in second [los])
  pos <- rownames(results[which(results$padj < pval & results$log2FoldChange > LFC),])
  neg <- rownames(results[which(results$padj < pval & results$log2FoldChange < LFC),])
  back <- rownames(results)
  
  
  # Convert gene names to those without the decimal points for compatibility with BioMart
  pos <- nameAbridge(pos)
  neg <- nameAbridge(neg)
  back <- nameAbridge(back)
  
  
  # Run the GO Enrichment function
  posEnrich <- GOenrich(pos,back,ensembl2entrez)
  negEnrich <- GOenrich(neg,back,ensembl2entrez)
  
  # Make list of top 5 p-values, names, colors
  GOval <- c(-log(posEnrich[n:1,5]), log(negEnrich[n:1,5]))
  GOhit <- c(posEnrich[n:1,1], negEnrich[n:1,1])
  
  if(any(is.null(color1), is.null(color2))){
    GOcol <- c(posEnrich[n:1,2], negEnrich[n:1,2])
    GOcol <- GOcolor(GOcol)
  }
  else{
    GOcol <- c(rep(c(color2, color1), each=n))
  }
  
  # Set graph bounds
  min.val = mround(min(GOval), 5, "floor")
  max.val = mround(max(GOval), 5, "ceiling")
  
  bound = max(abs(min.val), abs(max.val))
  
  # Plot the top n enrichments
  main.title = paste(title1, "vs", title2, "GO enrichments", sep=" ")
  plot <- barplot(GOval, horiz=T, main=main.title, col=GOcol, xlim=c(-bound, bound), xlab="log(P-value)")
  text(x=rep(c(2,-2),each=n),
       y=plot, labels=wrap.labels(GOhit,30), cex=.8, font=2,
       pos=rep(c(4,2),each=n))
}



#                      ##############                    #
#----------------------# SUBSETTING #--------------------#
#                      ##############                    #
# Function for subsetting all genes according to p-value and LFC thresholds for 1 v 1 comparisons
subsetter <- function(results1, results2, results3, results4, results5, results6, results7, p.val=p.thr, lfc=lfc.thr){
  
  # Subset significant genes according to LRT p-values
  res1.sig = results1[which(results1[,5] < p.val),]
  res2.sig = results2[which(results2[,5] < p.val),]
  res3.sig = results3[which(results3[,5] < p.val),]
  res4.sig = results4[which(results4[,5] < p.val),]
  res5.sig = results5[which(results5[,5] < p.val),]
  res6.sig = results6[which(results6[,5] < p.val),]
  res7.sig = results7[which(results7[,5] < p.val),]
  
  # subset genes according to shrunken LFC
  # Note: negative fold-changes indicate higher in second sample (i.e. 0 vs *___*)
  res1.diff <- res1.sig[which(abs(res1.sig[,2]) >= lfc),]
  res2.diff <- res2.sig[which(abs(res2.sig[,2]) >= lfc),]
  res3.diff <- res3.sig[which(abs(res3.sig[,2]) >= lfc),]
  res4.diff <- res4.sig[which(abs(res4.sig[,2]) >= lfc),]
  res5.diff <- res5.sig[which(abs(res5.sig[,2]) >= lfc),]
  res6.diff <- res6.sig[which(abs(res6.sig[,2]) >= lfc),]
  res7.diff <- res7.sig[which(abs(res7.sig[,2]) >= lfc),]
  
  # combine all genes into one, non-redundant list
  geneList.sig.diff = c(rownames(res1.diff), rownames(res2.diff), rownames(res3.diff), rownames(res4.diff), rownames(res5.diff), rownames(res6.diff), rownames(res7.diff))
  geneList.sig.diff = rownames(results1)[rownames(results1) %in% geneList.sig.diff]
  
  return(geneList.sig.diff) 
}



#                     ###############                    #
#---------------------# LFC SORTING #--------------------#
#                     ###############                    #
# Function to find the number of time points where LFC > 2
span.function <- function(matrix.row){
  length(which(abs(matrix.row[1:7]) >= lfc.thr))
}

# Function to find the time point where the LFC is greatest
max.function <- function(matrix.row){
  which.max(matrix.row[1:7])
}

# Function to find the earliest time point where LFC > 2
first.function <- function(matrix.row){
  min(which(abs(matrix.row[1:7]) >= lfc.thr))
}

# Function to make the LFC matrix from shrunken results objects
LFCmatrixMaker <- function(geneList, results1, results2, results3, results4, results5, results6, results7){
  
  # Make empty matrix
  LFCmatrix = matrix(nrow=length(geneList), ncol=7)
  
  # Fill with LFC from 1 v 1 results
  LFCmatrix[,1] = results1[,2][which(rownames(results1) %in% geneList)]
  LFCmatrix[,2] = results2[,2][which(rownames(results2) %in% geneList)]
  LFCmatrix[,3] = results3[,2][which(rownames(results3) %in% geneList)]
  LFCmatrix[,4] = results4[,2][which(rownames(results4) %in% geneList)]
  LFCmatrix[,5] = results5[,2][which(rownames(results5) %in% geneList)]
  LFCmatrix[,6] = results6[,2][which(rownames(results6) %in% geneList)]
  LFCmatrix[,7] = results7[,2][which(rownames(results7) %in% geneList)]
  
  # Set NA's to 0's
  LFCmatrix[is.na(LFCmatrix)]=0
  
  # Apply functions for sorting
  LFCmatrix = cbind(LFCmatrix, apply(LFCmatrix, 1, span.function))
  LFCmatrix = cbind(LFCmatrix, apply(LFCmatrix, 1, max.function))
  LFCmatrix = cbind(LFCmatrix, apply(LFCmatrix, 1, first.function))
  LFCmatrix = cbind(LFCmatrix, results1$baseMean[which(rownames(results1) %in% geneList)])
  
  # Set rownames to geneID, colnames to sample (compared to 0) + sorting parameter
  rownames(LFCmatrix) = geneList
  colnames(LFCmatrix) = c("30", "60", "90", "120", "240", "360", "1440", "tp.span", "tp.max", "tp.first", "base.means")
  
  return(LFCmatrix)
}

# Function for reordering and plotting LFC matrix
LFCplot <- function(LFC.matrix, cluster, order, truncate=5,geneNameIDList=NULL, k=7, maxval=NULL, main.title=NULL){
  
  # Convert name inputs for cluster/order into column numbers (based on build of LFC.matrix)
  conversion = c("tp.span"=8, "tp.max"=9, "tp.first"=10, "base.means"=11)
  
  if(is.character(cluster)){
    cluster = conversion[cluster]
  }
  
  if(is.character(order)){
    order = conversion[order]
  }
  
  # Make titles
  if(is.null(main.title)){
    if(cluster =="tp.span" | cluster == 8){
      title1 = "Number of TP with LFC above threshold"
    }
    if(cluster =="tp.max" | cluster == 9){
      title1 = "TP of max LFC"
    }
    if(cluster =="tp.first" | cluster == 10){
      title1 = "First TP with LFC above threshold"
    }
    if(cluster =="base.means" | cluster == 11){
      title1 = "Mean counts"
    }
    
    
    if(order =="tp.span" | order == 8){
      title2 = "Number of TP with LFC above threshold"
    }
    if(order =="tp.max" | order == 9){
      title2 = "TP of max LFC"
    }
    if(order =="tp.first" | order == 10){
      title2 = "First TP with LFC above threshold"
    }
    if(order =="base.means" | order == 11){
      title2 = "Mean counts"
    }
    main.title = paste0(title1, " > ", title2)
  }
  
  # Reorder LFC matrix according to (1) cluster and (2) order parameters
  LFC.matrix.order = LFC.matrix[order(LFC.matrix[,cluster], LFC.matrix[,order]),]
  
  # Identify gaps between sections
  gapList = c()
  for(i in 1:(k-1)){
    gapList = c(gapList, which(LFC.matrix.order[,cluster] == i)[length(which(LFC.matrix.order[,cluster] == i))])
  }
  
  # Truncate for plotting
  if(!is.null(maxval)){
    LFC.matrix.order[which(LFC.matrix.order>maxval)] = maxval
    LFC.matrix.order[which(LFC.matrix.order<(-maxval))] = (-maxval)
  }
  
  # Read in genes of interest to highlight
  geneNames = c()
  geneIDs = c()
  
  if(!is.null(geneNameIDList)){
    geneNames <- geneNameIDList[,1]
    geneIDs <- nameExtend(geneNameIDList[,2])
  }
  
  # Plot using Sushi2
  Sushi2::hotmap(LFC.matrix.order[,1:k], labrow=F, labcol=T, gaps=gapList, selectylabs=geneIDs, selectylabs.label = geneNames)
  mtext(side=3,line=1.0,font=2,text=main.title,cex=.5)
}



#                      ###########                        #
#----------------------# K-MEANS #------------------------#
#                      ###########                        #
# Function for k clustering
kclust <- function(matrix, k, type){
  if(type=="hclust" | type=="hierarchical" | type=="h"){
    clust = hclust(dist(matrix))
    cut = cutree(clust, k=k)
  }
  
  if(type=="kmeans" | type=="k"){
    clust = kmeans(matrix, centers=k)
    cut = clust$cluster
  }
  
  return(cut)
}

# Function for plotting kmeans heatmap
kHotmap <- function(kmatrix.norm, cut, cluster.order, k.colors=NULL, maxval=NULL, geneNameIDList=NULL, title=""){
  
  # Identify number of samples
  n = ncol(kmatrix.norm)
  
  # Identify number of clusters
  k = max(cut)
  
  # Order matrix and cluster assignment list in order of most --> least changed
  kmatrix.norm.order = kmatrix.norm[order(rowMax(kmatrix.norm), decreasing=T),]
  cut.order = cut[rownames(kmatrix.norm.order)]
  
  # Separate clusters in matrix according to desired order
  kmatrix.norm.order = kmatrix.norm.order[order(match(cut.order, cluster.order)),]
  cut.order = cut[rownames(kmatrix.norm.order)]
  
  # Truncate values
  if(!is.null(maxval)){
    kmatrix.norm.order[which(kmatrix.norm.order>maxval)] = maxval
    kmatrix.norm.order[which(kmatrix.norm.order<(-maxval))] = (-maxval)
  }
  
  # Assign cluster colors, if selected
  cluster.cols=c()
  if(!is.null(k.colors)){
    cluster.cols = k.colors[cut.order]
  }
  
  # Combine LFC and cluster assignments into one matrix to find gaps
  kmatrix.norm.order.gaps = cbind(kmatrix.norm.order, cut.order)
  
  # Identify gaps between clusters
  gaps.k = c()
  for(i in 1:(k-1)){
    gaps.k = c(gaps.k, which(kmatrix.norm.order.gaps[,(n+1)] == cluster.order[i])[length(which(kmatrix.norm.order.gaps[,(n+1)] == cluster.order[i]))])
  }
  
  # Set heatmap colors
  map.ramp = makeBreaks(maxval=maxval, num=100, min.col=heatmap.min, mid.col=heatmap.mid, max.col=heatmap.max)
  
  # Circle genes of interest
  geneNames=c()
  geneIDs=c()
  if(!is.null(geneNameIDList)){
    geneNames <- geneNameIDList[,1]
    geneIDs <- nameExtend(geneNameIDList[,2])
  }
  
  # Plot heatmap
  Sushi2::hotmap(kmatrix.norm.order[,1:n], col=map.ramp, labrow=F, labcol=T, gaps=gaps.k, selectylabs=geneIDs, selectylabs.label=geneNames, rowcolors=cluster.cols)
  Sushi::addlegend(c(-3,3), palette=colorRampPalette(colors=c(heatmap.min, heatmap.mid, heatmap.max)), title="Normalized Transcript Counts", bottominset=.5, xoffset=.11, title.offset = .07)
  mtext(side=3,line=1.0,font=2,text=title,cex=2)
}



#                        #######                        #
#------------------------# CSV #------------------------#
#                        #######                        #
# Function for determining from LFC and p-values which sample is higher
higherIn <- function(results, sampleName1, sampleName2, pval.cutoff){
  
  foldChange = results$log2FoldChange
  pVal = results$padj
  
  indicator = rep("", times=nrow(results))
  for(n in seq(1, length(geneList))){
    if(is.na(pVal[n])){
      indicator[n] = "Neither"
    } else{
      if((pVal[n] <= pval.cutoff) && (foldChange[n] > 0)){
        indicator[n] = sampleName2
      }
      if((pVal[n] <= pval.cutoff) && (foldChange[n] < 0)){
        indicator[n] = sampleName1
      }
      if(pVal[n] > pval.cutoff){
        indicator[n] = "Neither"
      }
    }
  }
  
  return(indicator)
}

# # Function for extracting count info for given samples
# countExtract <- function(dds.obj, sampleName1, sampleName2, sample1, sample2){
#   
#   count1.1 = assay(dds.obj)[,sample1]
#   count1.2 = assay(dds.obj)[,(sample1+8)]
#   
#   count2.1 = assay(dds.obj)[,sample2]
#   count2.2 = assay(dds.obj)[,(sample2+8)]
#   
#   name.count1.1 = paste("count", sampleName1, ".1")
#   name.count1.2 = paste("count", sampleName1, ".2")
#   
#   name.count2.1 = paste("count", sampleName2, ".1")
#   name.count2.2 = paste("count", sampleName2, ".2")
#   
#   count.df = data.frame(name.count1.1 = count1.1, name.count1.2 = count1.2, name.count2.1=count2.1, name.count2.2=count2.2)
#   
#   return(count.df)
# }

# Function for determining which genes are part of a set
goiCheck <- function(matrix, GoI){
  
  geneIDlist = nameAbridge(rownames(matrix))
  
  indicator = rep("", times=length(geneIDlist))
  for (n in seq(1, length(geneIDlist))){
    if (geneIDlist[n] %in% GoI){
      indicator[n] = "Y"
    } else{
      indicator[n] = "N"
    }
  }
  
  return(indicator)
  
}

# Function for obtaining a list of matched ENSEMBL + HGNC IDs
ENSEMBL.HGNC.extract <- function(geneList){
  
  # convert gene IDs to abridged, BioMaRT-compatible ENSEMBL IDs
  ENSEMBL = nameAbridge(rownames(geneList))
  
  # make ensembl-hgnc db with BioMaRT
  nameDb = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', values=ENSEMBL, mart=mart)
  
  # remove duplicated rows
  if(length(which(duplicated(nameDb$ensembl_gene_id))) > 0){
    nameDb = nameDb[-which(duplicated(nameDb$ensembl_gene_id)),]
  }
  
  # TEMPORARY FIX: Remove rows that didn't have ensembl_gene_id's fetched with BioMart, so for loop can work
  ENSEMBL = ENSEMBL[(ENSEMBL %in% nameDb$ensembl_gene_id)]
  
  # order according to ENSEMBL name list (BioMaRt comes out a different order than results)
  HGNC = rep("", times=length(ENSEMBL))
  for(n in seq(1, length(ENSEMBL))){
    HGNC[n] = nameDb$hgnc_symbol[which(nameDb$ensembl_gene_id == ENSEMBL[n])]
  }
  
  ENSEMBL.HGNC = data.frame(ENSEMBL, HGNC)
  return(ENSEMBL.HGNC)
}

# # Function for writing data to a csv
# geneWrite <- function(res, geneList, symbolList, count1_1, count1_2, count2_1, count2_2, sampleName1, sampleName2, filePath, pvalCutoff = .01)
# {
#   geneListLong = nameExtend(geneList)
#   
#   # Make list of the fold change (straight from res)
#   foldChange <- res$log2FoldChange[which(rownames(res) %in% geneListLong)]
#   
#   # Make list of the adjusted pvalues (straight from res as well)
#   pVal <- res$padj[which(rownames(res) %in% geneListLong)]
#   
#   # Make list based on p-values and foldChange that indicates which sample a gene is higher in, if any
#   higherIn <- rep("", times=length(geneListLong))
#   for (n in seq(1, length(geneListLong))){
#     if (is.na(pVal[n])){
#       higherIn[n] = "Neither"
#     } else{
#       if ((pVal[n] <= pvalCutoff) && (foldChange[n] > 0)){
#         higherIn[n] = sampleName2
#       }
#       if ((pVal[n] <= pvalCutoff) && (foldChange[n] < 0)){
#         higherIn[n] = sampleName1
#       }
#       if (pVal[n] > pvalCutoff){
#         higherIn[n] = "Neither"
#       }
#     }
#   }
#   
#   # Put it all in one data frame
#   data <- data.frame(geneList, symbolList, AP1ind, count1_1, count1_2, count2_1, count2_2, foldChange, higherIn, pVal, 
#                      stringsAsFactors = FALSE)
#   
#   # Write to csv
#   write.csv(data, file=filePath, row.names = F)
# }



#############################################################################################################
#                                                                                                           #
#                                                 R U N                                                     #
#                                                                                                           #
#############################################################################################################

#------------------------------------------------------------------------------------#
#                                       Read in                                      #
#------------------------------------------------------------------------------------#

#                    ################                    #
#--------------------# READ IN DATA #--------------------#
#                    ################                    #

# Read in the txi file, as defined in the config file - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** READRAW ***
if(readRaw == T){
  # Read in config file
  configPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_', projNum, '_samples.csv')
  config <- read.csv(configPath, header=T)
  
  # Read in count matrix
  txiPath <- paste0('/Users/phanstiel3/Research/Data/Projects/', proj, '/rna/proc/tximport/', proj, '_', projNum, '_txi.rds')
  txi <- readRDS(txiPath)
}


# Plot read info based on the txi file - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** TXIPLOT ***
if(txiPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "readInfo.pdf"), width=8, height=8)
  }
  
  # Plot count info
  par(mfrow=c(2,2))
  par(mar=c(5,4,4,2))
  
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



#                     #############                     #
#---------------------# RUN DESEQ #---------------------#
#                     #############                     #

# Create dds and transformed dds - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** RUNDESEQ ***
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

# Find shrunken LFC for each comparison vs 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** RUNSHRINK ***
if(runShrink == T){
  res.030 <- lfcShrink(dds, coef="condition_30_vs_0", type="apeglm")
  res.060 <- lfcShrink(dds, coef="condition_60_vs_0", type="apeglm")
  res.090 <- lfcShrink(dds, coef="condition_90_vs_0", type="apeglm")
  res.0120 <- lfcShrink(dds, coef="condition_120_vs_0", type="apeglm")
  res.0240 <- lfcShrink(dds, coef="condition_240_vs_0", type="apeglm")
  res.0360 <- lfcShrink(dds, coef="condition_360_vs_0", type="apeglm")
  res.01440 <- lfcShrink(dds, coef="condition_1440_vs_0", type="apeglm")
}



#                    ##############                    #
#--------------------# OTHER DATA #--------------------#
#                    ##############                    #

# Make long + abridged gene list from dds 
longGeneKey <- rownames(assay(dds))
abridgedGeneKey <- nameAbridge(longGeneKey)


# Read in excel of genes of interest (gene name and ENSG ID)
GoI <- read.csv("~/Research/Documents/Projects/LIMA/RNAseq/AP1.csv", header=F)
colnames(GoI) <- c("gene", "ID")

# Add long ID based on dds
for (n in 1:nrow(GoI))
{
  GoI$long.ID[n] = nameExtend(GoI$ID[n])[[1]]
}




#------------------------------------------------------------------------------------#
#                                       Analyze                                      #
#------------------------------------------------------------------------------------#

#                    #################                    #
#--------------------# GLOBAL TRENDS #--------------------#
#                    #################                    #

#-------------------#
# INTERACTIVE PLOTS #
#-------------------#
# PCA explorer (featuring pca2go, which requires the names be shortened to regular ENSEMBL IDs) - - - - - - - - - > *** PCAINT ***
if(PCAint == T){
  dds.rename = dds
  rownames(dds.rename) <- nameAbridge(rownames(dds.rename))
  pcaExplorer(dds=dds.rename)
}


#----------#
# PCA PLOT #
#----------#
# Plot PCA of transformed counts - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** PCAPLOT ***
if(PCAplot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "PCAplot.pdf"), width=8, height=8)
  }
  
  # Plot PCA of transformed counts
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2))
  
  print(plotPCA(dds.trans, intgroup="condition")+labs(title="Transformed Counts PCA")+scale_color_manual(values=sample.pal))
  
  if(makePDF == T){
    dev.off()
  }
}


#-------------#
# PLOT COUNTS #
#-------------#
# Plot counts of genes of interest - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** COUNTPLOT ***
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



#                  #####################                  #
#------------------# 1 v 1 COMPARISONS #------------------# 
#                  #####################                  #

#----------#
# MA PLOTS #
#----------#
# Make MA Plots for each comparison vs 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** MAPLOT ***
if(MAplot == T){
  if(makePDF == T){
    pdf(file=file.path(outputDir, "MAPlot.pdf"), width=12, height=6)
  }
  
  # Plot settings
  par(mar=c(5,4,4,2))
  par(mfrow=c(2,4))
  
  # Plot each with shrunken LFC, LRT p-val
  MAplotter(res.030, col1=sample.pal[1], col2=sample.pal[2],title1="0000", title2="0030", geneNameIDList = GoI)
  MAplotter(res.060, col1=sample.pal[1], col2=sample.pal[3],title1="0000", title2="0060", geneNameIDList = GoI)
  MAplotter(res.090, col1=sample.pal[1], col2=sample.pal[4],title1="0000", title2="0090", geneNameIDList = GoI)
  MAplotter(res.0120, col1=sample.pal[1], col2=sample.pal[5],title1="0000", title2="0120", geneNameIDList = GoI)
  MAplotter(res.0240, col1=sample.pal[1], col2=sample.pal[6],title1="0000", title2="0240", geneNameIDList = GoI)
  MAplotter(res.0360, col1=sample.pal[1], col2=sample.pal[7],title1="0000", title2="0360", geneNameIDList = GoI)
  MAplotter(res.01440, col1=sample.pal[1], col2=sample.pal[8],title1="0000", title2="1440", geneNameIDList = GoI)
  
  if(makePDF == T){
    dev.off()
  }
}


#----------------#
# GO ENRICHMENTS #
#----------------#
# Make GO graphs for each comparison vs 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** GOPLOT ***
if(GOplot == T){
  
  # Assign GO marts
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host= "www.ensembl.org")
  
  # Set background
  background = nameAbridge(rownames(res.030))
  ensembl2entrez <- getBM(filters="ensembl_gene_id",attributes=c("ensembl_gene_id", "entrezgene"), values=background, mart=mart)
  
  # Make PDF(?) and plot
  if(makePDF == T){
    pdf(file=file.path(outputDir, "GOplot.pdf"), width=12, height=6)
  }
  
  # Plot settings
  par(mar=c(5,4,4,2))
  par(mfrow=c(2,4))
  
  # Plot graphs of top 5 enrichments
  GOplotter(res.030, color1=sample.pal[1], color2=sample.pal[2], title1="0", title2="30")
  GOplotter(res.060, color1=sample.pal[1], color2=sample.pal[3], title1="0", title2="60")
  GOplotter(res.090, color1=sample.pal[1], color2=sample.pal[4], title1="0", title2="90")
  GOplotter(res.0120, color1=sample.pal[1], color2=sample.pal[5], title1="0", title2="120")
  GOplotter(res.0240, color1=sample.pal[1], color2=sample.pal[6], title1="0", title2="240")
  GOplotter(res.0360, color1=sample.pal[1], color2=sample.pal[7], title1="0", title2="360")
  GOplotter(res.01440, color1=sample.pal[1], color2=sample.pal[8], title1="0", title2="1440")
  
  if(makePDF == T){
    dev.off()
  }
}



#                    ##############                    #
#--------------------# SUBSETTING #--------------------# 
#                    ##############                    #

#----------------------#
# SHRUNKEN LFC + P-VAL #
#----------------------#
genes.sig.diff = subsetter(res.030, res.060, res.090, res.0120, res.0240, res.0360, res.01440)



#                   ###################                   #
#-------------------# CLUSTERING: LFC #-------------------#
#                   ###################                   #

# Build a matrix of subsetted genes with one row per gene, one column per comparison, and extra columns with sorting parameters
LFCmatrix.cat = LFCmatrixMaker(genes.sig.diff, res.030, res.060, res.090, res.0120, res.0240, res.0360, res.01440)

# Make another matrix without the sorting parameters
LFCmatrix.k = LFCmatrix.cat[,1:7]

#---------#
# K-MEANS #
#---------#
# Center and scale data
LFCmatrix.k.norm <- (LFCmatrix.k - rowMeans(LFCmatrix.k))/rowSds(LFCmatrix.k + .5)

# Find optimal number of clusters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** LFCCLUSTEROPT ***
if(LFCclusterOpt == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "LFCclusterOpt.pdf"), width=8, height=8)
  }
  
  # Method 1: WSS
  fviz_nbclust(LFCmatrix.k.norm, kmeans, method = "wss")
  
  # Method 2: Silhouette
  fviz_nbclust(LFCmatrix.k.norm, kmeans, method = "silhouette")
  
  # Method 3: Nb clustering
  nb <- NbClust(LFCmatrix.k.norm, distance = "euclidean", min.nc = 2,
                max.nc = 10, method = "kmeans")
  fviz_nbclust(nb)
  
  if(makePDF == T){
    dev.off()
  }
  
}

# Set the number of clusters
LFC.k <- 8

# Set seed to preserve manual ordering
set.seed(733)

# Perform clustering
LFC.cut = kclust(LFCmatrix.k.norm, k=LFC.k, type="k")

# Plot the individual k clusters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** LFCKCLUSTPLOT ***
if(LFCkclustPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "LFCkclustersPlot.pdf"), width=12, height=6)
  }
  
  # Plot each cluster
  par(mar=c(3,2,1,2))
  par(mfrow=c(2,LFC.k/2))
  
  for (i in 1:LFC.k) {
    
    # make empty plot
    plot(colMeans(LFCmatrix.k.norm[LFC.cut==i,]), type="n",main=paste("n=",table(LFC.cut)[i],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative LFC")
    
    # and transparent grey lines
    apply(LFCmatrix.k.norm[LFC.cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
    
    # add median line
    lines(colMeans(LFCmatrix.k.norm[LFC.cut==i,]),col=k.colors[i],lwd=2)
    
    # add axis
    axis(1, at=1:(LFC.k-1), labels=c(".5", "1", "1.5", "2", "4", "6", "24"))
  }
  
  if(makePDF == T){
    dev.off()
  }
}

# Plot heatmap of k clusters  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** LFCKCLUSTHEATPLOT ***
if(LFCkclustHeatPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "LFCkclustersHeatPlot.pdf"), width=8, height=8)
  }
  
  # Set the desired order of clusters, based on their timing (manual; changes with seed)
  LFC.cluster.order <- c(8,1,2,4,3,6,5,7)
  
  # Plot heatmap
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,5))
  
  kHotmap(LFCmatrix.k.norm, LFC.cut, LFC.cluster.order, maxval=3, geneNameIDList=GoI, title="LFC K-Means Clustering")
  
  if(makePDF == T){
    dev.off()
  }
}



#-------------#
# BY CATEGORY #
#-------------#
# Plot LFC heatmaps clustered by one parameter and sorted by the other - - - - - - - - - - - - - - - - - - > *** LFCORDERPLOT ***
if(LFCorderPlot == T){
  if(makePDF == T){
    pdf(file=file.path(outputDir, "LFCorderPlot.pdf"), width=12, height=8)
  }
  
  par(mfrow=c(2,3))
  par(mar=c(5,4,4,2))
  maxval=5
  LFCplot(LFCmatrix.cat, cluster="tp.span", order="tp.max", geneNameIDList = GoI, maxval=maxval)
  LFCplot(LFCmatrix.cat, cluster="tp.span", order="tp.first", geneNameIDList = GoI, maxval=maxval)
  LFCplot(LFCmatrix.cat, cluster="tp.max", order="tp.span", geneNameIDList = GoI, maxval=maxval)
  LFCplot(LFCmatrix.cat, cluster="tp.max", order="tp.first", geneNameIDList = GoI, maxval=maxval)
  LFCplot(LFCmatrix.cat, cluster="tp.first", order="tp.span", geneNameIDList = GoI, maxval=maxval)
  LFCplot(LFCmatrix.cat, cluster="tp.first", order="tp.max", geneNameIDList = GoI, maxval=maxval)
  
  if(makePDF == T){
    dev.off()
  }
}



#            ##################################           #
#------------# CLUSTERING: TRANSFORMED COUNTS #-----------#
#            ##################################           #

# Build a matrix of subsetted genes with one row per gene, one column per comparison
countMatrix = dds.trans[which(rownames(dds.trans) %in% genes.sig.diff),]
countMatrix = assay(countMatrix)

#---------#
# K-MEANS #
#---------#

# Center and scale data
countMatrix.norm <- (countMatrix - rowMeans(countMatrix))/rowSds(countMatrix + .5)

# Find optimal number of clusters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** CLUSTEROPT ***
if(countClusterOpt == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countClusterOpt.pdf"), width=8, height=8)
  }
  
  # Method 1: WSS
  fviz_nbclust(countMatrix.norm, kmeans, method = "wss")
  
  # Method 2: Silhouette
  fviz_nbclust(countMatrix.norm, kmeans, method = "silhouette")
  
  # Method 3: Nb clustering
  nb <- NbClust(countMatrix.norm, distance = "euclidean", min.nc = 2,
                max.nc = 10, method = "kmeans")
  fviz_nbclust(nb)
  
  if(makePDF == T){
    dev.off()
  }
  
}

# Set the number of clusters
count.k <- 8

# Set seed to preserve manual ordering
set.seed(733)

# Perform clustering
count.cut = kclust(countMatrix.norm, k=count.k, type="k")

# Combine replicates
countMatrix.norm.combo <- combineReps(countMatrix.norm, new.colnames=c("0", "30", "60", "90", "120", "240", "360", "1440"))

# Plot the individual k clusters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** COUNTKCLUSTPLOT ***
if(countKclustPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countKclustersPlot.pdf"), width=12, height=6)
  }
  
  # Plot each cluster
  par(mar=c(3,2,1,2))
  par(mfrow=c(2,count.k/2))
  
  for (i in 1:count.k) {
    
    # make empty plot
    plot(colMeans(countMatrix.norm.combo[count.cut==i,]), type="n",main=paste("n=",table(count.cut)[i],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative expression")
    
    # and transparent grey lines
    apply(countMatrix.norm.combo[count.cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
    
    # add median line
    lines(colMeans(countMatrix.norm.combo[count.cut==i,]),col=k.colors[i],lwd=2)
    
    # add axis
    axis(1, at=1:count.k, labels=c("0",".5", "1", "1.5", "2", "4", "6", "24"))
  }
  
  if(makePDF == T){
    dev.off()
  }
}

# Plot heatmap of k clusters  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** COUNTKCLUSTHEATPLOT ***
if(countKclustHeatPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countKclustersHeatPlotCommittee.pdf"), width=8, height=8)
  }
  
  # Set the desired order of clusters, based on their timing (manual; changes with seed)
  count.cluster.order <- c(8,3,2,4,6,5,1,7)
  
  # Plot heatmap
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,5))
  
  kHotmap(countMatrix.norm.combo, count.cut, count.cluster.order, maxval=3, geneNameIDList=GoI, title="Normalized Counts Clustering")
  
  if(makePDF == T){
    dev.off()
  }
}


#          #######################################         #
#----------# TSV OVERVIEWS: LFC, COUNTS, CLUSTER #---------#
#          #######################################         #

# ENSG ID
# hgnc symbol
# GoI?
# Count1
# Count2
# Count3
# Count4
# Count5
# Count6
# Count7
# Count8
# LFC1
# LFC2
# LFC3
# LFC4
# LFC5
# LFC6
# LFC7
# tp.span
# tp.max
# tp.first
# base.mean
# cluster(counts)

# Make a TSV of significant/differential genes (incl. transformed counts, cluster, LFC, LFC order) - - - - - - - > *** TSVSUB ***
if(tsvSub==T){
  # Prep transformed count matrix (from countMatrix, selected for pval/LFC)
  tsv.countMatrix.combo = combineReps(countMatrix, new.colnames=c("0", "30", "60", "90", "120", "240", "360", "1440"))
  # Find ENSG + HGNC pairs
  tsv.sub.ENSG.HGNC = ENSEMBL.HGNC.extract(tsv.countMatrix.combo)
  # Combine counts and clusters, only for genes with ENSG-HGNC pairs
  tsv.countMatrix.combo = data.frame(tsv.countMatrix.combo, count.cut)
  tsv.countMatrix.combo = tsv.countMatrix.combo[which(nameAbridge(rownames(tsv.countMatrix.combo)) %in% tsv.sub.ENSG.HGNC$ENSEMBL),]
  
  # Prep LFC matrix
  tsv.LFCmatrix = LFCmatrix.cat[which(nameAbridge(rownames(LFCmatrix.cat)) %in% tsv.sub.ENSG.HGNC$ENSEMBL),]
  
  # Identify GoI
  tsv.sub.GOI = goiCheck(tsv.countMatrix.combo, GoI$ID)
  
  # Put all the data in one dataframe
  tsv.sub <- data.frame(tsv.sub.ENSG.HGNC, tsv.sub.GOI, tsv.countMatrix.combo, tsv.LFCmatrix, 
                        stringsAsFactors = FALSE)
  rownames(tsv.sub) = c()
  colnames(tsv.sub) = c("ENSEMBL", "HGNC", "GoI?", "0.norm.count", "30.norm.count", "60.norm.count", "90.norm.count",
                        "120.norm.count", "240.norm.count", "360.norm.count", "1440.norm.count", "cluster", "30.LFC", 
                        "60.LFC", "90.LFC", "120.LFC", "240.LFC", "360.LFC", "1440.LFC", "tp.span", "tp.max", "tp.first",
                        "base.mean")
  
  # Write to tsv
  write_tsv(tsv.sub, path=file.path(outputDir, "LIMA_RNA_Subset-genes.tsv"))
}

# Make a TSV of all genes (incl. transformed counts, cluster, LFC, LFC order) - - - - - - - - - - - - - - - - - > *** TSVFULL ***
if(tsvFull==T){ 

  # Prep transformed count matrix (from dds, no selection for pval/LFC)
  tsv.full.combo = combineReps(assay(dds.trans), new.colnames=c("0", "30", "60", "90", "120", "240", "360", "1440"))
  # Find ENSG + HGNC pairs
  tsv.full.ENSG.HGNC = ENSEMBL.HGNC.extract(tsv.full.combo)
  # Make full list of clusters (including cluster 0 for non-significant)
  full.cut = rep(0, times=nrow(tsv.full.combo))
  names(full.cut) = rownames(tsv.full.combo)
  for(gene in names(full.cut)){
    if(gene %in% names(count.cut)){
      full.cut[gene] = count.cut[gene]
    }
  }
  
  # Combine counts and clusters, only for genes with ENSG-HGNC pairs
  tsv.full.combo = data.frame(tsv.full.combo, full.cut)
  tsv.full.combo = tsv.full.combo[which(nameAbridge(rownames(tsv.full.combo)) %in% tsv.full.ENSG.HGNC$ENSEMBL),]
  
  # Prep LFC matrix
  full.list = rownames(tsv.full.combo)
  tsv.LFCmatrix.full = LFCmatrixMaker(full.list, res.030, res.060, res.090, res.0120, res.0240, res.0360, res.01440)
  
  # Identify GoI
  tsv.full.GOI = goiCheck(tsv.full.combo, GoI$ID)
  
  # Put all the data in one dataframe
  tsv.full <- data.frame(tsv.full.ENSG.HGNC, tsv.full.GOI, tsv.full.combo, tsv.LFCmatrix.full, 
                        stringsAsFactors = FALSE)
  rownames(tsv.full) = c()
  colnames(tsv.full) = c("ENSEMBL", "HGNC", "GoI?", "0.norm.count", "30.norm.count", "60.norm.count", "90.norm.count",
                        "120.norm.count", "240.norm.count", "360.norm.count", "1440.norm.count", "cluster", "30.LFC", 
                        "60.LFC", "90.LFC", "120.LFC", "240.LFC", "360.LFC", "1440.LFC", "tp.span", "tp.max", "tp.first",
                        "base.mean")
  
  # Write to tsv
  write_tsv(tsv.full, path=file.path(outputDir, "LIMA_RNA.tsv"))
}




# #######################################################################################################################
# 
# ## For committee meeting figure:
# 
# #pdf(file=file.path(outputDir, "countKclustersHeatPlotCommittee.pdf"), width=8, height=8)
# png(filename=file.path(outputDir, "countKclustersHeatPlotCommittee.png"), width=8, height=8, units="in", res=300)
# 
# # Set the desired order of clusters, based on their timing (manual; changes with seed)
# count.cluster.order <- c(8,3,2,4,6,5,1,7)
# 
# # Plot heatmap
# par(mfrow=c(1,1))
# par(mar=c(5,4,4,5))
# 
# kmatrix.norm = countMatrix.norm.combo
# cut = count.cut
# cluster.order = count.cluster.order
# maxval = 3
# geneNameIDList = GoI
# title=""
# k.colors=NULL
# 
# # Identify number of samples
# n = ncol(kmatrix.norm)
# 
# # Identify number of clusters
# k = max(cut)
# 
# # Order matrix and cluster assignment list in order of most --> least changed
# kmatrix.norm.order = kmatrix.norm[order(rowMax(kmatrix.norm), decreasing=T),]
# cut.order = cut[rownames(kmatrix.norm.order)]
# 
# # Separate clusters in matrix according to desired order
# kmatrix.norm.order = kmatrix.norm.order[order(match(cut.order, cluster.order)),]
# cut.order = cut[rownames(kmatrix.norm.order)]
# 
# # Truncate values
# if(!is.null(maxval)){
#   kmatrix.norm.order[which(kmatrix.norm.order>maxval)] = maxval
#   kmatrix.norm.order[which(kmatrix.norm.order<(-maxval))] = (-maxval)
# }
# 
# # Assign cluster colors, if selected
# cluster.cols=c()
# if(!is.null(k.colors)){
#   cluster.cols = k.colors[cut.order]
# }
# 
# # Combine LFC and cluster assignments into one matrix to find gaps
# kmatrix.norm.order.gaps = cbind(kmatrix.norm.order, cut.order)
# 
# # Identify gaps between clusters
# gaps.k = c()
# for(i in 1:(k-1)){
#   gaps.k = c(gaps.k, which(kmatrix.norm.order.gaps[,(n+1)] == cluster.order[i])[length(which(kmatrix.norm.order.gaps[,(n+1)] == cluster.order[i]))])
# }
# 
# # Set heatmap colors
# map.ramp = makeBreaks(maxval=maxval, num=100, min.col=heatmap.min, mid.col=heatmap.mid, max.col=heatmap.max)
# 
# # Circle genes of interest
# geneNames=c()
# geneIDs=c()
# if(!is.null(geneNameIDList)){
#   geneNames <- geneNameIDList[,1]
#   geneIDs <- nameExtend(geneNameIDList[,2])
# }
# 
# colnames(kmatrix.norm.order) = c("0", "0.5", "1", "1.5", "2", "4", "6", "24")
# 
# # Plot heatmap
# Sushi2::hotmap(kmatrix.norm.order[,1:n], col=map.ramp, labrow=F, labcol=T, gaps=gaps.k, selectylabs=geneIDs, selectylabs.label=geneNames, selectylabs.col = c(rep("darkorange3", times=6), rep("royalblue3", times=3), "darkorange3"), rowcolors=cluster.cols)
# Sushi::addlegend(c(-3,3), palette=colorRampPalette(colors=c(heatmap.min, heatmap.mid, heatmap.max)), title="Normalized Transcript Counts", bottominset=.5, xoffset=.11, title.offset = .07)
# mtext(side=3,line=1.0,font=2,text=title,cex=2)
# 
# dev.off()





# #######################################################################################################################
# 
# # Playing with clusters:

# Set time points, DoF
time <- c(0,30,60,90,2*60,4*60,6*60, 24*60)
time = 1:8
dgr <- 3

# # Re-design dds with lm + polynomial formula???
# design <- model.matrix(~poly(time, degree=dgr))
# dds.new = dds
# dds.new$time = time
# design(dds.new) <- ~poly(time, degree=4)
# dds.new <- DESeq(dds.new, test="LRT", reduced=~1)
# 
# betas <- coef(dds.new)
# fitted <- betas %*% t(design)

identical(names(count.cut), rownames(countMatrix.norm.combo))

# Cut the normalized + combined count Matrix into 8 clusters, find mean of each
clusterExtract <- function(k){
  clust.means = matrix(nrow=0, ncol=8)
  for (i in 1:k){
    clust = countMatrix.norm.combo[count.cut == i,]
    clust.mean = colMeans(clust)
    clust.means = rbind(clust.means, clust.mean)
  }
  colnames(clust.means) = time
  return(clust.means)
}

means = clusterExtract(8)

# Plot the cluster means, and the predicted poly curve based on only the means for each cluster
clusterPlot <- function(k){
  df = data.frame(time=1:8)
  for (i in 1:k){
    fit = lm(means[i,] ~ poly(time, degree=dgr))
    curve = predict(fit, newdat=df)
    plot(time, means[i,])
    lines(df$time, curve, col="blue")
  }
}

par(mfrow=c(2,4))
clusterPlot(8)

# Extract the betas + residuals for each curve, to compare with betas for other clusters from other data sets (???)
clusterBetas <- function(k){
  betas = matrix(nrow=0, ncol=(dgr+1))
  for (i in 1:k){
    fit = lm(means[i,] ~ poly(time, degree=dgr))
    betas=rbind(betas, as.vector(coef(fit)))
  }
  colnames(betas) = names(coef(fit))
  return(betas)
}

clusterResid <- function(k){
  resid = matrix(nrow=0, ncol=8)
  for (i in 1:k){
    fit = lm(means[i,] ~ poly(time, degree=dgr))
    resid=rbind(resid, as.vector(resid(fit)))
  }
  return(resid)
}

betas = clusterBetas(8)
resid = clusterResid(8)













