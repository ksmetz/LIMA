#####################################################
# This is the script for analyzing peak count  data #
#    from the ChIPcluster pipeline on longleaf.     #
#####################################################


#############################################################################################################
#                                                                                                           #
#                                            I N I T I A L I Z E                                            #
#                                                                                                           #
#############################################################################################################

#-------------------------------------------------------------------------------------#
#                                       Libraries                                     #
#-------------------------------------------------------------------------------------#

library(readr)
library(DESeq2)
library(pals)


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
FRiPplot = F
PCAplot = F
MAplot = F
LFCclusterOpt = F
KclustPlot = F
KclustHeatPlot = F
tsvSub = F
tsvFull = F



#                   ##################                   #
#-------------------# COLOR SETTINGS #-------------------#
#                   ##################                   #

# Assign sample colors for PCA, plot counts
sample.pal=parula(8)

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
projNum = '1-2.1.1'
target = 'H3K27ac'

tp = c("0", "30", "60", "90", "120", "240", "360", "1440")

transformType = "var"

# Assign cutoffs for subsetting
p.thr = .01
lfc.thr = 1

# Assign + make output directory for PDFs, CSVs, etc.
outputDir = file.path("/Users/phanstiel3/Research/Data/Projects", proj, "chip/h3k27ac/diff")
if(newDir == T){
  
  # Assign new output dir
  outputName = paste(proj, "ChIP", target, projNum, Sys.Date(), sep="_")
  outputDir = file.path("/Users/phanstiel3/Desktop", outputName)
  
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

# Function for combining replicates in a 16x(many genes) matrix; r1/r2/r1/r2/r1/r2 structure
combineReps <- function(matrix, new.colnames)
{
  newMatrix = c()
  for (i in seq(1,ncol(matrix), by=2))
  {
    tempMatrix = matrix(data=c(matrix[,i], matrix[,i+1]), ncol=2, byrow=F)
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

# Read in peak read counts, FRiP table and adjust for analysis - - - - - - - - - - - - - - - - - - - - - - - - - -> *** READRAW ***
if(readRaw == T){
  # Read in the counts.tsv file containing the reads (values) at each peak (rows) for each sample (column)
  countTable <- data.frame(read_tsv('/Users/phanstiel3/Research/Data/Projects/LIMA/chip/h3k27ac/proc/peaks/counts.tsv', col_names=F))
  
  # Set columns to shortened sample names, bed columns; set rows to combination of coordinates and subset to remove coordinates
  colnames(countTable) <- c("chr", "start", "end",
                            "0000_r1", "0000_r2",
                            "0030_r1", "0030_r2",
                            "0060_r1", "0060_r2",
                            "0090_r1", "0090_r2",
                            "0120_r1", "0120_r2",
                            "0240_r1", "0240_r2",
                            "0360_r1", "0360_r2",
                            "1440_r1", "1440_r2")
  rownames(countTable) = paste(countTable$chr, countTable$start, countTable$end, sep="_")
  
  countMatrix = countTable[,4:19]
  
  # Read in the FRiP table
  FRiPtable <- data.frame(read_delim('/Users/phanstiel3/Research/Data/Projects/LIMA/chip/h3k27ac/proc/frip/FRiP.txt', delim=" ", col_names=F))
  
  # Convert to a 2-column, where each row is sample (column 1 = reads in peaks, column 2 = total reads)
  FRiPmatrix = data.frame(matrix(FRiPtable$X2, ncol=2, byrow = T))
  rownames(FRiPmatrix) = colnames(countMatrix)
  
  # Calculate FRiP from columns 1 and 2 
  FRiPmatrix$FRiP = FRiPmatrix$X1/FRiPmatrix$X2
}

#                     #############                     #
#---------------------# RUN DESEQ #---------------------#
#                     #############                     #

# Create dds and transformed dds - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** RUNDESEQ ***

if(runDESeq == T){
  # Create data frame with sample names and whatever condition you want to distinguish them by (config$Condition, config$Cell.Type, etc.)
  colData <- data.frame(condition = factor(rep(tp, each=2)), rep = factor(rep(c(1,2),times=8)), row.names = colnames(countMatrix))
  
  # Create DESeq Data Set from matrix, using LRT for p-values
  dds <- DESeqDataSetFromMatrix(countMatrix, colData=colData, design =~ rep + condition)
  
  dds <- DESeq(dds, test="LRT", full=~rep + condition, reduced = ~rep)

  # Transform the dds (rlog, vst or ntd) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** TRANSFORMTYPE ***
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

# Find shrunken LFC for each comparison vs 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** RUNSHRINK ***
if(runShrink == T){
  res.030 <- lfcShrink(dds, coef="condition_30_vs_0", type="apeglm")
  res.060 <- lfcShrink(dds, coef="condition_60_vs_0", type="apeglm")
  res.090 <- lfcShrink(dds, coef="condition_90_vs_0", type="apeglm")
  res.0120 <- lfcShrink(dds, coef="condition_120_vs_0", type="apeglm")
  res.0240 <- lfcShrink(dds, coef="condition_240_vs_0", type="apeglm")
  res.0360 <- lfcShrink(dds, coef="condition_360_vs_0", type="apeglm")
  res.01440 <- lfcShrink(dds, coef="condition_1440_vs_0", type="apeglm")
}


#------------------------------------------------------------------------------------#
#                                       Analyze                                      #
#------------------------------------------------------------------------------------#

#                    #################                    #
#--------------------# GLOBAL TRENDS #--------------------#
#                    #################                    #

#-----------#
# FRiP PLOT #
#-----------#
# Plot Fraction of Reads in Peaks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** FRIPPLOT ***
if(FRiPplot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "FRiPplot.pdf"), width=8, height=8)
  }

  # Plot FRiP for each sample + rep
  bp = barplot(FRiPmatrix$FRiP, main="FRiP", col=rep(sample.pal, each=2))
  axis(side=1,las=2,at=bp,labels=rownames(FRiPmatrix))
  
  if(makePDF == T){
    dev.off()
  }

}


#----------#
# PCA PLOT #
#----------#
# Plot PCA of transformed counts - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** PCAPLOT ***
if(PCAplot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "PCAplot.pdf"), width=8, height=8)
  }

  # Plot PCA of transformed counts
  print(plotPCA(dds.trans, intgroup="condition")+labs(title="Transformed Counts PCA")+scale_color_manual(values=sample.pal))
  
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
# Make MA Plots for each comparison vs 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** MAPLOT ***
if(MAplot == T){
  if(makePDF == T){
    pdf(file=file.path(outputDir, "MAPlot.pdf"), width=12, height=6)
  }
  
  # Plot settings
  par(mar=c(5,4,4,2))
  par(mfrow=c(2,4))
  
  # Plot each with shrunken LFC, LRT p-val
  MAplotter(res.030, col1=sample.pal[1], col2=sample.pal[2],title1="0000", title2="0030")
  MAplotter(res.060, col1=sample.pal[1], col2=sample.pal[3],title1="0000", title2="0060")
  MAplotter(res.090, col1=sample.pal[1], col2=sample.pal[4],title1="0000", title2="0090")
  MAplotter(res.0120, col1=sample.pal[1], col2=sample.pal[5],title1="0000", title2="0120")
  MAplotter(res.0240, col1=sample.pal[1], col2=sample.pal[6],title1="0000", title2="0240")
  MAplotter(res.0360, col1=sample.pal[1], col2=sample.pal[7],title1="0000", title2="0360")
  MAplotter(res.01440, col1=sample.pal[1], col2=sample.pal[8],title1="0000", title2="1440")
  
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
peaks.sig.diff = subsetter(res.030, res.060, res.090, res.0120, res.0240, res.0360, res.01440, p.val=.01, lfc=1)



#                   ###################                   #
#-------------------# CLUSTERING: LFC #-------------------#
#                   ###################                   #

# Build a matrix of subsetted genes with one row per gene, one column per comparison
subsetMatrix = dds.trans[which(rownames(dds.trans) %in% peaks.sig.diff),]
subsetMatrix = assay(subsetMatrix)

#---------#
# K-MEANS #
#---------#

# Center and scale data
subsetMatrix.norm <- (subsetMatrix - rowMeans(subsetMatrix))/rowSds(subsetMatrix + .5)

# Find optimal number of clusters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** LFCCLUSTEROPT ***
if(LFCclusterOpt == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "LFCclusterOpt.pdf"), width=8, height=8)
  }
  
  # Method 1: WSS
  fviz_nbclust(subsetMatrix.norm, kmeans, method = "wss")
  
  # Method 2: Silhouette
  fviz_nbclust(subsetMatrix.norm, kmeans, method = "silhouette")
  
  # Method 3: Nb clustering
  nb <- NbClust(subsetMatrix.norm, distance = "euclidean", min.nc = 2,
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
count.cut = kclust(subsetMatrix.norm, k=count.k, type="k")

# Combine replicates
subsetMatrix.norm.combo <- combineReps(subsetMatrix.norm, new.colnames=c("0", "30", "60", "90", "120", "240", "360", "1440"))


# Plot the individual k clusters - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - > *** KCLUSTPLOT ***
if(KclustPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countKclustersPlot.pdf"), width=12, height=6)
  }

  # Plot each cluster
  par(mar=c(3,2,1,2))
  par(mfrow=c(2,count.k/2))
  
  for (i in 1:count.k) {
    
    # make empty plot
    plot(colMeans(subsetMatrix.norm.combo[count.cut==i,]), type="n",main=paste("n=",table(count.cut)[i],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative expression")
    
    # and transparent grey lines
    apply(subsetMatrix.norm.combo[count.cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
    
    # add median line
    lines(colMeans(subsetMatrix.norm.combo[count.cut==i,]),col=k.colors[i],lwd=2)
    
    # add axis
    axis(1, at=1:count.k, labels=c("0",".5", "1", "1.5", "2", "4", "6", "24"))
  }

  if(makePDF == T){
    dev.off()
  }
  
}


# Plot heatmap of k clusters  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -> *** KCLUSTHEATPLOT ***
if(KclustHeatPlot == T){
  
  if(makePDF == T){
    pdf(file=file.path(outputDir, "countKclustersHeatPlot.pdf"), width=8, height=8)
  }
  
  # Set the desired order of clusters, based on their timing (manual; changes with seed)
  count.cluster.order <- c(1,2,3,4,5,6,7,8)
  
  # Plot heatmap
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,5))
  
  kHotmap(subsetMatrix.norm.combo, count.cut, count.cluster.order, maxval=3, title="Normalized Counts Clustering", k.colors=k.colors)
  
  if(makePDF == T){
    dev.off()
  }

}



#          #######################################         #
#----------# TSV OVERVIEWS: LFC, COUNTS, CLUSTER #---------#
#          #######################################         #

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
# cluster(counts)

# Make a TSV of significant/differential genes (incl. transformed counts, cluster, LFC, LFC order) - - - - - - - - > *** TSVSUB ***
if(tsvSub==T){
  # Prep transformed count matrix (from subsetMatrix, selected for pval/LFC)
  tsv.subsetMatrix.combo = combineReps(subsetMatrix, new.colnames=c("0", "30", "60", "90", "120", "240", "360", "1440"))

  # Combine counts and clusters
  tsv.subsetMatrix.combo = data.frame(tsv.subsetMatrix.combo, count.cut)

  # Split rownames into 3 columns (bed format coords)
  coords = strsplit(rownames(tsv.subsetMatrix.combo), "_")
  coords = unlist(coords)
  
  chr = coords[seq(1, length(coords), by=3)]
  srt = coords[seq(2, length(coords), by=3)]
  end = coords[seq(3, length(coords), by=3)]

  # Put all the data in one dataframe
  tsv.sub <- data.frame(chr, srt, end, tsv.subsetMatrix.combo,
                        stringsAsFactors = FALSE)
  rownames(tsv.sub) = c()
  colnames(tsv.sub) = c("chr", "start", "end", "0.norm.count", "30.norm.count", "60.norm.count", "90.norm.count",
                        "120.norm.count", "240.norm.count", "360.norm.count", "1440.norm.count", "cluster")

  # Write to tsv
  write_tsv(tsv.sub, path=file.path(outputDir, "LIMA_ChIP_H3k27ac_Subset-peaks.tsv"))
}

# Make a TSV of all genes (incl. transformed counts, cluster, LFC, LFC order) - - - - - - - - - - - - - - - - - - > *** TSVFULL ***
if(tsvFull==T){

  # Prep transformed count matrix (from dds, no selection for pval/LFC)
  tsv.full.combo = combineReps(assay(dds.trans), new.colnames=c("0", "30", "60", "90", "120", "240", "360", "1440"))

  # Make full list of clusters (including cluster 0 for non-significant)
  full.cut = rep(0, times=nrow(tsv.full.combo))
  names(full.cut) = rownames(tsv.full.combo)
  for(peak in names(full.cut)[(names(full.cut) %in% names(count.cut))]){
    full.cut[peak] = count.cut[peak]
  }
  # Combine counts and clusters
  tsv.full.combo = data.frame(tsv.full.combo, full.cut)
  
  # Split rownames into 3 columns (bed format coords)
  coords = strsplit(rownames(tsv.full.combo), "_")
  coords = unlist(coords)
  
  chr = coords[seq(1, length(coords), by=3)]
  srt = coords[seq(2, length(coords), by=3)]
  end = coords[seq(3, length(coords), by=3)]

  # Put all the data in one dataframe
  tsv.full <- data.frame(chr, srt, end, tsv.full.combo,
                         stringsAsFactors = FALSE)
  rownames(tsv.full) = c()
  colnames(tsv.full) = c("chr", "start", "end", "0.norm.count", "30.norm.count", "60.norm.count", "90.norm.count",
                         "120.norm.count", "240.norm.count", "360.norm.count", "1440.norm.count", "cluster")

  # Write to tsv
  write_tsv(tsv.full, path=file.path(outputDir, "LIMA_ChIP_H3K27ac.tsv"))
}






# ######################################################################################################################
# 
# ## For committee meeting figure:
# 
# #pdf(file=file.path(outputDir, "countKclustersHeatPlotCommittee.pdf"), width=8, height=8)
# png(filename=file.path(outputDir, "countKclustersHeatPlotCommittee.png"), width=8, height=8, units="in", res=300)
# 
# # Set the desired order of clusters, based on their timing (manual; changes with seed)
# count.cluster.order <- c(7,6,8,3,4,1,2,5)
# 
# # Plot heatmap
# par(mfrow=c(1,1))
# par(mar=c(5,4,4,5))
# 
# kmatrix.norm = subsetMatrix.norm.combo
# cut = count.cut
# cluster.order = count.cluster.order
# maxval = 3
# geneNameIDList = NULL
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
# Sushi::addlegend(c(-3,3), palette=colorRampPalette(colors=c(heatmap.min, heatmap.mid, heatmap.max)), title="Normalized Peak Read Counts", bottominset=.5, xoffset=.11, title.offset = .07)
# mtext(side=3,line=1.0,font=2,text=title,cex=2)
# 
# dev.off()
# 
