# This is the very first look into the ChIP data.
# I copied the counts.tsv file from longleaf and will just run DESeq2 on it to see if the format works.
# Using copied code snippets from RNA_Analysis.R

library(readr)
library(DESeq2)
library(pals)

k.colors <- brewer.dark2(8)
sample.pal=parula(8)
tp = c("0000", "0030", "0060", "0090", "0120", "0240", "0360", "1440")
heatmap.min = "red"
heatmap.mid = "white"
heatmap.max = "blue"


# Functions
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

# Function for combining replicates in a 16x(many genes) matrix 
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
  Sushi2::hotmap(kmatrix.norm.order[,1:7], col=map.ramp, labrow=F, labcol=T, gaps=gaps.k, selectylabs=geneIDs, selectylabs.label=geneNames, rowcolors=cluster.cols)
  Sushi::addlegend(c(-3,3), palette=colorRampPalette(colors=c(heatmap.min, heatmap.mid, heatmap.max)), title="Normalized Transcript Counts", bottominset=.5, xoffset=.11, title.offset = .07)
  mtext(side=3,line=1.0,font=2,text=title,cex=2)
}


countTable <- data.frame(read_tsv('/Users/phanstiel3/Desktop/counts.tsv', col_names=F))
colnames(countTable) <- c("chr", "start", "end",
                          "0000_r1", "0000_r2",
                          "0030_r1", "0030_r2",
                          "0060_r1", "0060_r2",
                          "0090_r1", "0090_r2",
                          "0120_r1", "0120_r2",
                          "0240_r1", "0240_r2",
                          "0360_r1", "0360_r2",
                          "1440_r1", "1440_r2")
head(countTable)

rownames(countTable) = paste(countTable$chr, countTable$start, countTable$end, sep="_")

countMatrix = countTable[,4:19]

head(countMatrix)



# Create data frame with sample names and whatever condition you want to distinguish them by (config$Condition, config$Cell.Type, etc.)
colData <- data.frame(condition = factor(rep(tp, each=2)), rep = factor(rep(c(1,2),times=8)), row.names = colnames(countMatrix))

# Create DESeq Data Set from matrix, using LRT for p-values
dds <- DESeqDataSetFromMatrix(countMatrix, colData=colData, design =~ rep + condition)

dds <- DESeq(dds, test="LRT", full=~rep + condition, reduced = ~rep)

dds.trans <- vst(dds)

pdf(file=file.path("/Users/phanstiel3/Desktop/ChIPplots.pdf"), width=8, height=8)

print(plotPCA(dds.trans, intgroup="condition")+labs(title="Transformed Counts PCA")+scale_color_manual(values=sample.pal))


res.030 <- lfcShrink(dds, coef="condition_0030_vs_0000", type="apeglm")
res.060 <- lfcShrink(dds, coef="condition_0060_vs_0000", type="apeglm")
res.090 <- lfcShrink(dds, coef="condition_0090_vs_0000", type="apeglm")
res.0120 <- lfcShrink(dds, coef="condition_0120_vs_0000", type="apeglm")
res.0240 <- lfcShrink(dds, coef="condition_0240_vs_0000", type="apeglm")
res.0360 <- lfcShrink(dds, coef="condition_0360_vs_0000", type="apeglm")
res.01440 <- lfcShrink(dds, coef="condition_1440_vs_0000", type="apeglm")


MAplotter(res.030, col1=sample.pal[1], col2=sample.pal[2],title1="0000", title2="0030")
MAplotter(res.060, col1=sample.pal[1], col2=sample.pal[3],title1="0000", title2="0060")
MAplotter(res.090, col1=sample.pal[1], col2=sample.pal[4],title1="0000", title2="0090")
MAplotter(res.0120, col1=sample.pal[1], col2=sample.pal[5],title1="0000", title2="0120")
MAplotter(res.0240, col1=sample.pal[1], col2=sample.pal[6],title1="0000", title2="0240")
MAplotter(res.0360, col1=sample.pal[1], col2=sample.pal[7],title1="0000", title2="0360")
MAplotter(res.01440, col1=sample.pal[1], col2=sample.pal[8],title1="0000", title2="1440")


peaks.sig.diff = subsetter(res.030, res.060, res.090, res.0120, res.0240, res.0360, res.01440, p.val=.01, lfc=1)

subsetMatrix = dds.trans[which(rownames(dds.trans) %in% peaks.sig.diff),]
subsetMatrix = assay(subsetMatrix)

# Center and scale data
subsetMatrix.norm <- (subsetMatrix - rowMeans(subsetMatrix))/rowSds(subsetMatrix + .5)


# Set the number of clusters
count.k <- 8

# Perform clustering
count.cut = kclust(subsetMatrix.norm, k=count.k, type="k")

# Combine replicates
subsetMatrix.norm.combo <- combineReps(subsetMatrix.norm, new.colnames=c("0", "30", "60", "90", "120", "240", "360", "1440"))


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

# Set the desired order of clusters, based on their timing (manual; changes with seed)
count.cluster.order <- c(1,2,3,4,5,6,7,8)

# Plot heatmap
par(mfrow=c(1,1))
par(mar=c(5,4,4,5))

kHotmap(subsetMatrix.norm.combo, count.cut, count.cluster.order, maxval=3, title="Normalized Counts Clustering", k.colors=k.colors)

dev.off()
