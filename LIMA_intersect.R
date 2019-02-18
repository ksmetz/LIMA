# This is the very start of the code intersecting RNA-Seq, H3K27ac ChIP-Seq, and HiC (starting with prelim from Macrophages) data.

library(readr)
library(DESeq2)
library(pals)
library(Sushi2)
library(bedtoolsr)

# Function for shortening the decimal ENSEMBL name from dds, res, etc. to the standard ENSEMBL gene ID
nameAbridge <- function(longName)
{
  split <- strsplit(longName, "[.]")
  shortName <- unlist(split)[2*(1:length(longName))-1]
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

# Assign colors for heatmaps
heatmap.min = "steelblue1"
heatmap.mid = "black"
heatmap.max = "gold"

# Assign colors for clusters
k.colors <- brewer.dark2(8)

chip <- data.frame(read_tsv("/Users/phanstiel3/Dropbox/Work/Research/Data/Projects/LIMA/chip/h3k27ac/diff/LIMA_ChIP_H3K27ac.tsv"))
chipMerge <- data.frame(chr = chip$chr, start = chip$start, end = chip$end, count = paste(chip$X0.norm.count, chip$X30.norm.count, 
                                                                                          chip$X60.norm.count, chip$X90.norm.count, 
                                                                                          chip$X120.norm.count, chip$X240.norm.count,
                                                                                          chip$X360.norm.count, chip$X1440.norm.count,
                                                                                          sep=","), 
                        stringsAsFactors = F)
write_tsv(chipMerge, path="/Users/phanstiel3/Desktop/chip.bed", col_names=F)


rna <- data.frame(read_tsv("/Users/phanstiel3/Dropbox/Work/Research/Data/Projects/LIMA/rna/diff/LIMA_RNA.tsv"))
cluster = data.frame(gene = rna$ENSEMBL, cluster = rna$cluster)

hic <- data.frame(read_tsv("/Users/phanstiel3/Dropbox/Work/Research/Data/Projects/CI/hic/loops/CI_THP1_A_0.0.0.loops.10Kb.bedpe", col_names=F))

anchor1 = data.frame(chr=hic$X1, start=hic$X2, end=hic$X3)
write_tsv(anchor1, path="/Users/phanstiel3/Desktop/anchor1.bed", col_names=F)
anchor2 = data.frame(chr=hic$X4, start=hic$X5, end=hic$X6)
write_tsv(anchor2, path="/Users/phanstiel3/Desktop/anchor2.bed", col_names=F)

prom <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promList.bed", col_names=F))
prom$X5 = nameAbridge(prom$X5)
prom = prom[order(prom$X5),]

promSub = data.frame(chr = prom$X1, start = prom$X2, stop = prom$X3, name = prom$X5, strand = prom$X4, stringsAsFactors = F)

promSize = 5000

prom0 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 0]),]
prom1 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 1]),]
prom2 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 2]),]
prom3 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 3]),]
prom4 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 4]),]
prom5 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 5]),]
prom6 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 6]),]
prom7 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 7]),]
prom8 = promSub[which(promSub$name %in% cluster$gene[cluster$cluster == 8]),]

prom0 = prom0[prom0$start > 0,]

write_tsv(prom0, path="/Users/phanstiel3/Desktop/promoter_cluster0.bed", col_names=F)
write_tsv(prom1, path="/Users/phanstiel3/Desktop/promoter_cluster1.bed", col_names=F)
write_tsv(prom2, path="/Users/phanstiel3/Desktop/promoter_cluster2.bed", col_names=F)
write_tsv(prom3, path="/Users/phanstiel3/Desktop/promoter_cluster3.bed", col_names=F)
write_tsv(prom4, path="/Users/phanstiel3/Desktop/promoter_cluster4.bed", col_names=F)
write_tsv(prom5, path="/Users/phanstiel3/Desktop/promoter_cluster5.bed", col_names=F)
write_tsv(prom6, path="/Users/phanstiel3/Desktop/promoter_cluster6.bed", col_names=F)
write_tsv(prom7, path="/Users/phanstiel3/Desktop/promoter_cluster7.bed", col_names=F)
write_tsv(prom8, path="/Users/phanstiel3/Desktop/promoter_cluster8.bed", col_names=F)


# # in python:
# chip = "/Users/phanstiel3/Desktop/chip.bed"
# for cluster in ["0", "1", "2", "3", "4", "5", "6", "7", "8"]:
#   promBed = "/Users/phanstiel3/Desktop/promoter_cluster" + cluster + ".bed"
# intersect = "/Users/phanstiel3/Desktop/promoter_chip" + cluster + ".bed"
# intersectCmd = "bedtools intersect -wa -a " + chip + " -b " + promBed + " > " + intersect
# os.system(intersectCmd)


chip_prom0 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip0.bed", col_names = F))
chip_prom1 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip1.bed", col_names = F))
chip_prom2 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip2.bed", col_names = F))
chip_prom3 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip3.bed", col_names = F))
chip_prom4 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip4.bed", col_names = F))
chip_prom5 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip5.bed", col_names = F))
chip_prom6 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip6.bed", col_names = F))
chip_prom7 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip7.bed", col_names = F))
chip_prom8 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chip8.bed", col_names = F))


out = matrix(as.numeric(unlist(strsplit(chip_prom0$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom0 = cbind(chip_prom0, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom1$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom1 = cbind(chip_prom1, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom2$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom2 = cbind(chip_prom2, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom3$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom3 = cbind(chip_prom3, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom4$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom4 = cbind(chip_prom4, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom5$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom5 = cbind(chip_prom5, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom6$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom6 = cbind(chip_prom6, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom7$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom7 = cbind(chip_prom7, out)

out = matrix(as.numeric(unlist(strsplit(chip_prom8$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_prom8 = cbind(chip_prom8, out)





# # in python:
# # Intersect with HiC data
# anchor1 = "/Users/phanstiel3/Desktop/anchor1.bed"
# anchor2 = "/Users/phanstiel3/Desktop/anchor2.bed"
# for cluster in ["0", "1", "2", "3", "4", "5", "6", "7", "8"]:
#   promBed = "/Users/phanstiel3/Desktop/promoter_cluster" + cluster + ".bed"
# intersect1 = "/Users/phanstiel3/Desktop/promoter_loop1_" + cluster + ".bed"	
# intersect2 = "/Users/phanstiel3/Desktop/promoter_loop2_" + cluster + ".bed"
# intersectCmd1 = "bedtools intersect -c -a " + anchor1 + " -b " + promBed + " > " + intersect1
# intersectCmd2 = "bedtools intersect -c -a " + anchor2 + " -b " + promBed + " > " + intersect2
# os.system(intersectCmd1)
# os.system(intersectCmd2)

loop1_prom0 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_0.bed", col_names = F))
loop1_prom1 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_1.bed", col_names = F))
loop1_prom2 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_2.bed", col_names = F))
loop1_prom3 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_3.bed", col_names = F))
loop1_prom4 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_4.bed", col_names = F))
loop1_prom5 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_5.bed", col_names = F))
loop1_prom6 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_6.bed", col_names = F))
loop1_prom7 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_7.bed", col_names = F))
loop1_prom8 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_8.bed", col_names = F))

loop2_prom0 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_0.bed", col_names = F))
loop2_prom1 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_1.bed", col_names = F))
loop2_prom2 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_2.bed", col_names = F))
loop2_prom3 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_3.bed", col_names = F))
loop2_prom4 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_4.bed", col_names = F))
loop2_prom5 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_5.bed", col_names = F))
loop2_prom6 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_6.bed", col_names = F))
loop2_prom7 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_7.bed", col_names = F))
loop2_prom8 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_8.bed", col_names = F))

loop_prom0 <- cbind(loop1_prom0, loop2_prom0)
loop_prom1 <- cbind(loop1_prom1, loop2_prom1)
loop_prom2 <- cbind(loop1_prom2, loop2_prom2)
loop_prom3 <- cbind(loop1_prom3, loop2_prom3)
loop_prom4 <- cbind(loop1_prom4, loop2_prom4)
loop_prom5 <- cbind(loop1_prom5, loop2_prom5)
loop_prom6 <- cbind(loop1_prom6, loop2_prom6)
loop_prom7 <- cbind(loop1_prom7, loop2_prom7)
loop_prom8 <- cbind(loop1_prom8, loop2_prom8)

loop_prom0 <- loop_prom0[which(loop_prom0[,4] > 0 | loop_prom0[,8] > 0),]
loop_prom1 <- loop_prom1[which(loop_prom1[,4] > 0 | loop_prom1[,8] > 0),]
loop_prom2 <- loop_prom2[which(loop_prom2[,4] > 0 | loop_prom2[,8] > 0),]
loop_prom3 <- loop_prom3[which(loop_prom3[,4] > 0 | loop_prom3[,8] > 0),]
loop_prom4 <- loop_prom4[which(loop_prom4[,4] > 0 | loop_prom4[,8] > 0),]
loop_prom5 <- loop_prom5[which(loop_prom5[,4] > 0 | loop_prom5[,8] > 0),]
loop_prom6 <- loop_prom6[which(loop_prom6[,4] > 0 | loop_prom6[,8] > 0),]
loop_prom7 <- loop_prom7[which(loop_prom7[,4] > 0 | loop_prom7[,8] > 0),]
loop_prom8 <- loop_prom8[which(loop_prom8[,4] > 0 | loop_prom8[,8] > 0),]

loop1_prom0 = loop_prom0[,1:4]
loop1_prom1 = loop_prom1[,1:4]
loop1_prom2 = loop_prom2[,1:4]
loop1_prom3 = loop_prom3[,1:4]
loop1_prom4 = loop_prom4[,1:4]
loop1_prom5 = loop_prom5[,1:4]
loop1_prom6 = loop_prom6[,1:4]
loop1_prom7 = loop_prom7[,1:4]
loop1_prom8 = loop_prom8[,1:4]

loop2_prom0 = loop_prom0[,5:8]
loop2_prom1 = loop_prom1[,5:8]
loop2_prom2 = loop_prom2[,5:8]
loop2_prom3 = loop_prom3[,5:8]
loop2_prom4 = loop_prom4[,5:8]
loop2_prom5 = loop_prom5[,5:8]
loop2_prom6 = loop_prom6[,5:8]
loop2_prom7 = loop_prom7[,5:8]
loop2_prom8 = loop_prom8[,5:8]

enh0 <- rbind(loop1_prom0[which(loop1_prom0[,4] == 0),], loop2_prom0[which(loop2_prom0[,4] == 0),])
enh1 <- rbind(loop1_prom1[which(loop1_prom1[,4] == 0),], loop2_prom1[which(loop2_prom1[,4] == 0),])
enh2 <- rbind(loop1_prom2[which(loop1_prom2[,4] == 0),], loop2_prom2[which(loop2_prom2[,4] == 0),])
enh3 <- rbind(loop1_prom3[which(loop1_prom3[,4] == 0),], loop2_prom3[which(loop2_prom3[,4] == 0),])
enh4 <- rbind(loop1_prom4[which(loop1_prom4[,4] == 0),], loop2_prom4[which(loop2_prom4[,4] == 0),])
enh5 <- rbind(loop1_prom5[which(loop1_prom5[,4] == 0),], loop2_prom5[which(loop2_prom5[,4] == 0),])
enh6 <- rbind(loop1_prom6[which(loop1_prom6[,4] == 0),], loop2_prom6[which(loop2_prom6[,4] == 0),])
enh7 <- rbind(loop1_prom7[which(loop1_prom7[,4] == 0),], loop2_prom7[which(loop2_prom7[,4] == 0),])
enh8 <- rbind(loop1_prom8[which(loop1_prom8[,4] == 0),], loop2_prom8[which(loop2_prom8[,4] == 0),])

write_tsv(enh0, path="/Users/phanstiel3/Desktop/enhancer_cluster0.bed", col_names=F)
write_tsv(enh1, path="/Users/phanstiel3/Desktop/enhancer_cluster1.bed", col_names=F)
write_tsv(enh2, path="/Users/phanstiel3/Desktop/enhancer_cluster2.bed", col_names=F)
write_tsv(enh3, path="/Users/phanstiel3/Desktop/enhancer_cluster3.bed", col_names=F)
write_tsv(enh4, path="/Users/phanstiel3/Desktop/enhancer_cluster4.bed", col_names=F)
write_tsv(enh5, path="/Users/phanstiel3/Desktop/enhancer_cluster5.bed", col_names=F)
write_tsv(enh6, path="/Users/phanstiel3/Desktop/enhancer_cluster6.bed", col_names=F)
write_tsv(enh7, path="/Users/phanstiel3/Desktop/enhancer_cluster7.bed", col_names=F)
write_tsv(enh8, path="/Users/phanstiel3/Desktop/enhancer_cluster8.bed", col_names=F)

# # in python:
# # Intersect with ChIP data
# chip = "/Users/phanstiel3/Desktop/chip.bed"
# for cluster in ["0", "1", "2", "3", "4", "5", "6", "7", "8"]:
#   enhBed = "/Users/phanstiel3/Desktop/enhancer_cluster" + cluster + ".bed"
# intersect = "/Users/phanstiel3/Desktop/enhancer_chip" + cluster + ".bed"
# intersectCmd = "bedtools intersect -wa -a " + chip + " -b " + enhBed + " > " + intersect
# os.system(intersectCmd)


chip_enh0 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip0.bed", col_names = F))
chip_enh1 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip1.bed", col_names = F))
chip_enh2 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip2.bed", col_names = F))
chip_enh3 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip3.bed", col_names = F))
chip_enh4 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip4.bed", col_names = F))
chip_enh5 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip5.bed", col_names = F))
chip_enh6 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip6.bed", col_names = F))
chip_enh7 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip7.bed", col_names = F))
chip_enh8 <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chip8.bed", col_names = F))


out = matrix(as.numeric(unlist(strsplit(chip_enh0$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh0 = cbind(chip_enh0, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh1$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh1 = cbind(chip_enh1, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh2$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh2 = cbind(chip_enh2, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh3$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh3 = cbind(chip_enh3, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh4$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh4 = cbind(chip_enh4, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh5$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh5 = cbind(chip_enh5, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh6$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh6 = cbind(chip_enh6, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh7$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh7 = cbind(chip_enh7, out)

out = matrix(as.numeric(unlist(strsplit(chip_enh8$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
chip_enh8 = cbind(chip_enh8, out)




pdf("/Users/phanstiel3/Desktop/Promoter_ChIP.pdf", width=12, height=6)
par(mfrow=c(2,4))

# # make empty plot
# plot(colMeans(chip_prom0[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# # and transparent grey lines
# apply(chip_prom0[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# # add median line
# lines(colMeans(chip_prom0[,5:12]),col="black",lwd=2)
# # add axis
# axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom1[,5:12]), type="n", main=paste("n=",nrow(chip_prom1),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom1[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom1[,5:12]),col=k.colors[1],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom2[,5:12]), type="n", main=paste("n=",nrow(chip_prom2),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom2[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom2[,5:12]),col=k.colors[2],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom3[,5:12]), type="n", main=paste("n=",nrow(chip_prom3),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom3[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom3[,5:12]),col=k.colors[3],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom4[,5:12]), type="n", main=paste("n=",nrow(chip_prom4),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom4[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom4[,5:12]),col=k.colors[4],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom5[,5:12]), type="n", main=paste("n=",nrow(chip_prom5),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom5[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom5[,5:12]),col=k.colors[5],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom6[,5:12]), type="n", main=paste("n=",nrow(chip_prom6),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom6[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom6[,5:12]),col=k.colors[6],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom7[,5:12]), type="n", main=paste("n=",nrow(chip_prom7),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom7[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom7[,5:12]),col=k.colors[7],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_prom8[,5:12]), type="n", main=paste("n=",nrow(chip_prom8),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_prom8[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_prom8[,5:12]),col=k.colors[8],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

dev.off()







pdf("/Users/phanstiel3/Desktop/Enhancer_ChIP.pdf", width=12, height=6)
par(mfrow=c(2,4))

# # make empty plot
# plot(colMeans(chip_enh0[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# # and transparent grey lines
# apply(chip_enh0[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# # add median line
# lines(colMeans(chip_enh0[,5:12]),col="black",lwd=2)
# # add axis
# axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh1[,5:12]), type="n", main=paste("n=",nrow(chip_enh1),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh1[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh1[,5:12]),col=k.colors[1],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh2[,5:12]), type="n", main=paste("n=",nrow(chip_enh2),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh2[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh2[,5:12]),col=k.colors[2],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh3[,5:12]), type="n", main=paste("n=",nrow(chip_enh3),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh3[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh3[,5:12]),col=k.colors[3],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh4[,5:12]), type="n", main=paste("n=",nrow(chip_enh4),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh4[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh4[,5:12]),col=k.colors[4],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh5[,5:12]), type="n", main=paste("n=",nrow(chip_enh5),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh5[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh5[,5:12]),col=k.colors[5],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh6[,5:12]), type="n", main=paste("n=",nrow(chip_enh6),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh6[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh6[,5:12]),col=k.colors[6],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh7[,5:12]), type="n", main=paste("n=",nrow(chip_enh7),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh7[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh7[,5:12]),col=k.colors[7],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

# make empty plot
plot(colMeans(chip_enh8[,5:12]), type="n", main=paste("n=",nrow(chip_enh8),sep=""), ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# and transparent grey lines
apply(chip_enh8[,5:12],1,lines,col=adjustcolor("grey",alpha.f=0.1))
# add median line
lines(colMeans(chip_enh8[,5:12]),col=k.colors[8],lwd=2)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))

dev.off()




rna.norm = as.matrix(rna[,4:11])
rna.norm = (rna.norm - rowMeans(rna.norm))/rowSds(rna.norm + .5)

# par(mfrow=c(2,4))
# for (i in 1:8){
#   # make empty plot
#   plot(colMeans(rna.norm[rna$cluster==i,]), type="n", main=paste("n=",table(rna$cluster)[i+1],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative expression")
#   # and transparent grey lines
#   apply(rna.norm[rna$cluster==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
#   # add median line
#   lines(colMeans(rna.norm[rna$cluster==i,]),col=k.colors[i],lwd=2)
#   # add axis
#   axis(1, at=1:8, labels=c("0",".5", "1", "1.5", "2", "4", "6", "24"))
# }

pdf("/Users/phanstiel3/Desktop/RNA_ChIP_sig.pdf", width=12, height=6)
par(mfrow=c(2,4))

# make empty plot
plot(colMeans(chip_enh1[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==1,]),col=k.colors[1],lwd=2)
# add prom median line
lines(colMeans(chip_prom1[,5:12]),col=k.colors[1],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh1[,5:12]),col=k.colors[1],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[2],")",sep=""), 
                paste("prom (",nrow(chip_prom1),")",sep=""), 
                paste("enh (",nrow(chip_enh1),")",sep="")), 
       col=k.colors[1], lty=c(1, 2, 3), lwd=2)

# make empty plot
plot(colMeans(chip_enh2[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==2,]),col=k.colors[2],lwd=2)
# add prom median line
lines(colMeans(chip_prom2[,5:12]),col=k.colors[2],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh2[,5:12]),col=k.colors[2],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[3],")",sep=""), 
                paste("prom (",nrow(chip_prom2),")",sep=""), 
                paste("enh (",nrow(chip_enh2),")",sep="")), 
       col=k.colors[2], lty=c(1, 2, 3), lwd=2)

# make empty plot
plot(colMeans(chip_enh3[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==3,]),col=k.colors[3],lwd=2)
# add prom median line
lines(colMeans(chip_prom3[,5:12]),col=k.colors[3],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh3[,5:12]),col=k.colors[3],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[4],")",sep=""), 
                paste("prom (",nrow(chip_prom3),")",sep=""), 
                paste("enh (",nrow(chip_enh3),")",sep="")), 
       col=k.colors[3], lty=c(1, 2, 3), lwd=2)

# make empty plot
plot(colMeans(chip_enh4[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==4,]),col=k.colors[4],lwd=2)
# add prom median line
lines(colMeans(chip_prom4[,5:12]),col=k.colors[4],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh4[,5:12]),col=k.colors[4],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[5],")",sep=""), 
                paste("prom (",nrow(chip_prom4),")",sep=""), 
                paste("enh (",nrow(chip_enh4),")",sep="")), 
       col=k.colors[4], lty=c(1, 2, 3), lwd=2)

# make empty plot
plot(colMeans(chip_enh5[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==5,]),col=k.colors[5],lwd=2)
# add prom median line
lines(colMeans(chip_prom5[,5:12]),col=k.colors[5],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh5[,5:12]),col=k.colors[5],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[6],")",sep=""), 
                paste("prom (",nrow(chip_prom5),")",sep=""), 
                paste("enh (",nrow(chip_enh5),")",sep="")), 
       col=k.colors[5], lty=c(1, 2, 3), lwd=2)

# make empty plot
plot(colMeans(chip_enh6[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==6,]),col=k.colors[6],lwd=2)
# add prom median line
lines(colMeans(chip_prom6[,5:12]),col=k.colors[6],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh6[,5:12]),col=k.colors[6],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[7],")",sep=""), 
                paste("prom (",nrow(chip_prom6),")",sep=""), 
                paste("enh (",nrow(chip_enh6),")",sep="")), 
       col=k.colors[6], lty=c(1, 2, 3), lwd=2)

# make empty plot
plot(colMeans(chip_enh7[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==7,]),col=k.colors[7],lwd=2)
# add prom median line
lines(colMeans(chip_prom7[,5:12]),col=k.colors[7],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh7[,5:12]),col=k.colors[7],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[8],")",sep=""), 
                paste("prom (",nrow(chip_prom7),")",sep=""), 
                paste("enh (",nrow(chip_enh7),")",sep="")), 
       col=k.colors[7], lty=c(1, 2, 3), lwd=2)

# make empty plot
plot(colMeans(chip_enh8[,5:12]), type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna median line
lines(colMeans(rna.norm[rna$cluster==8,]),col=k.colors[8],lwd=2)
# add prom median line
lines(colMeans(chip_prom8[,5:12]),col=k.colors[8],lwd=2, lty=2)
# add enh median line
lines(colMeans(chip_enh8[,5:12]),col=k.colors[8],lwd=2, lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c(paste("RNA (",table(rna$cluster)[9],")",sep=""), 
                paste("prom (",nrow(chip_prom8),")",sep=""), 
                paste("enh (",nrow(chip_enh8),")",sep="")), 
       col=k.colors[8], lty=c(1, 2, 3), lwd=2)

dev.off()


# Cross correlation test

clust1 = data.frame(
  enh = colMeans(chip_enh1[,5:12]),
  pro = colMeans(chip_prom1[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==1,])
)
clust2 = data.frame(
  enh = colMeans(chip_enh2[,5:12]),
  pro = colMeans(chip_prom2[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==2,])
)
clust3 = data.frame(
  enh = colMeans(chip_enh3[,5:12]),
  pro = colMeans(chip_prom3[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==3,])
)
clust4 = data.frame(
  enh = colMeans(chip_enh4[,5:12]),
  pro = colMeans(chip_prom4[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==4,])
)
clust5 = data.frame(
  enh = colMeans(chip_enh5[,5:12]),
  pro = colMeans(chip_prom5[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==5,])
)
clust6 = data.frame(
  enh = colMeans(chip_enh6[,5:12]),
  pro = colMeans(chip_prom6[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==6,])
)
clust7 = data.frame(
  enh = colMeans(chip_enh7[,5:12]),
  pro = colMeans(chip_prom7[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==7,])
)
clust8 = data.frame(
  enh = colMeans(chip_enh8[,5:12]),
  pro = colMeans(chip_prom8[,5:12]),
  rna = colMeans(rna.norm[rna$cluster==8,])
)

pdf("/Users/phanstiel3/Desktop/RNA_ChIP_sig_ccf.pdf", width=12, height=6)
par(mfrow=c(2,4))
ccf(x=clust1$enh, y=clust1$rna)
abline(v=0, col="red", lty=3)
ccf(x=clust2$enh, y=clust2$rna)
abline(v=0, col="red", lty=3)
ccf(x=clust3$enh, y=clust3$rna)
abline(v=0, col="red", lty=3)
ccf(x=clust4$enh, y=clust4$rna)
abline(v=0, col="red", lty=3)
ccf(x=clust5$enh, y=clust5$rna)
abline(v=0, col="red", lty=3)
ccf(x=clust6$enh, y=clust6$rna)
abline(v=0, col="red", lty=3)
ccf(x=clust7$enh, y=clust7$rna)
abline(v=0, col="red", lty=3)
ccf(x=clust8$enh, y=clust8$rna)
abline(v=0, col="red", lty=3)
dev.off()

# # rna-enh
# plot(clust4$rna, type="n")
# lines(x=1:8, y=clust4$rna)
# lines(x=1:8, y=clust4$enh)
# lines(x=1:8, y=clust4$pro)
# plot(x=clust4$rna, y=clust4$enh)
# #plot(x=clust4$rna, y=clust4$pro)
# model <- lm(clust4$rna ~ clust4$enh)
# coef(model)
# abline(a = coef(model)[1], b = coef(model)[2], lty = 2, lwd = 2, col = "red")
# cor(clust4[,1:2])[1,2]
# 
# # shift +1
# model <- lm(clust4$rna[1:7] ~ clust4$enh[2:8])
# plot(x=clust4$rna[1:7], y=clust4$enh[2:8])
# plot(x=clust4$rna[1:6], y=clust4$enh[3:8])
# plot(x=clust4$rna[1:5], y=clust4$enh[4:8])




# BCL6; gene in massive hub from Phanstiel 2017
# MAFB doesn't show any contacts? Expand the range for looking for overlapping anchors?
# HES1; gene in loop, see HiC2 regions from rotation
# MALT1; differential region (see "FinalLoops" pres from rotation)
# DUSP10; differential region
# AFF3; decreasing gene from diff region

# IL1B
IL1B_rna = rna[which(rna$HGNC == "IL1B"),]
IL1B.rna.norm = as.matrix(IL1B_rna[,4:11])
IL1B.rna.norm = (IL1B.rna.norm - rowMeans(IL1B.rna.norm))/rowSds(IL1B.rna.norm + .5)

IL1B_pro = promSub[which(promSub$name == IL1B_rna$ENSEMBL),]
write_tsv(IL1B_pro, path="/Users/phanstiel3/Desktop/promoter_IL1B.bed", col_names=F)

# # in python:
# chip = "/Users/phanstiel3/Desktop/chip.bed"
# promBed = "/Users/phanstiel3/Desktop/promoter_IL1B.bed"
# intersect = "/Users/phanstiel3/Desktop/promoter_chipIL1B.bed"
# intersectCmd = "bedtools intersect -wa -a " + chip + " -b " + promBed + " > " + intersect
# os.system(intersectCmd)

IL1B_pro_chip <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_chipIL1B.bed", col_names = F))

out = matrix(as.numeric(unlist(strsplit(IL1B_pro_chip$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
IL1B_pro_chip = cbind(IL1B_pro_chip, out)

# # in python:
# # Intersect with HiC data
# anchor1 = "/Users/phanstiel3/Desktop/anchor1.bed"
# anchor2 = "/Users/phanstiel3/Desktop/anchor2.bed"
# promBed = "/Users/phanstiel3/Desktop/promoter_IL1B.bed"
# intersect1 = "/Users/phanstiel3/Desktop/promoter_loop1_IL1B.bed"
# intersect2 = "/Users/phanstiel3/Desktop/promoter_loop2_IL1B.bed"
# intersectCmd1 = "bedtools intersect -c -a " + anchor1 + " -b " + promBed + " > " + intersect1
# intersectCmd2 = "bedtools intersect -c -a " + anchor2 + " -b " + promBed + " > " + intersect2
# os.system(intersectCmd1)
# os.system(intersectCmd2)

loop1_proIL1B <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop1_IL1B.bed", col_names = F))
loop2_proIL1B <- data.frame(read_tsv("/Users/phanstiel3/Desktop/promoter_loop2_IL1B.bed", col_names = F))

loop_proIL1B <- cbind(loop1_proIL1B, loop2_proIL1B)

loop_proIL1B <- loop_proIL1B[which(loop_proIL1B[,4] > 0 | loop_proIL1B[,8] > 0),]

loop1_proIL1B = loop_proIL1B[,1:4]
loop2_proIL1B = loop_proIL1B[,5:8]

enh_IL1B <- rbind(loop1_proIL1B[which(loop1_proIL1B[,4] == 0),], loop2_proIL1B[which(loop2_proIL1B[,4] == 0),])

write_tsv(enh_IL1B, path="/Users/phanstiel3/Desktop/enhancer_IL1B.bed", col_names=F)
enh_IL1B

# # in python:
# # Intersect with ChIP data
# chip = "/Users/phanstiel3/Desktop/chip.bed"
# enhBed = "/Users/phanstiel3/Desktop/enhancer_IL1B.bed"
# intersect = "/Users/phanstiel3/Desktop/enhancer_chipIL1B.bed"
# intersectCmd = "bedtools intersect -wa -a " + chip + " -b " + enhBed + " > " + intersect
# os.system(intersectCmd)

IL1B_enh_chip <- data.frame(read_tsv("/Users/phanstiel3/Desktop/enhancer_chipIL1B.bed", col_names = F))

out = matrix(as.numeric(unlist(strsplit(IL1B_enh_chip$X4, split=","))), ncol=8, byrow=T)
out <- (out - rowMeans(out))/rowSds(out + .5)
IL1B_enh_chip = cbind(IL1B_enh_chip, out)

pdf("~/Desktop/RNA_ChIP_sig_IL1B.pdf", width=12, height=4)
par(mfrow=c(1,3))
# make empty plot
plot(1:8, type="n", ylim=c(-2, 2), xaxt="n", xlab="Hours after LPS treatment", ylab="Normalized counts")
# add rna line
lines(x=1:8, y=IL1B.rna.norm,col="darkgrey",lwd=2)
# add prom lines
apply(IL1B_pro_chip[,5:12],1,lines,col=adjustcolor("black",alpha.f=0.1), lty=2)
lines(x=1:8, y=colMeans(IL1B_pro_chip[,5:12]),col="black",lty=2)
# add enh lines
apply(IL1B_enh_chip[,5:12],1,lines,col=adjustcolor("red",alpha.f=0.1), lty=3)
lines(x=1:8, y=colMeans(IL1B_enh_chip[,5:12]),col="red",lty=3)
# add axis
axis(1, at=1:8, labels=c("0", ".5", "1", "1.5", "2", "4", "6", "24"))
legend("topleft", 
       legend=c("RNA", "prom", "enh"), 
       col=c("darkgrey", "black", "red"), lty=c(1,2,3), lwd=2)


IL1B_data = data.frame(
  enh = colMeans(IL1B_enh_chip[,5:12]),
  pro = colMeans(IL1B_pro_chip[,5:12]),
  rna = as.vector(IL1B.rna.norm))

ccf(x=IL1B_data$enh, y=IL1B_data$rna)
abline(v=0, col="red", lty=3)
ccf(x=IL1B_data$pro, y=IL1B_data$rna)
abline(v=0, col="red", lty=3)
dev.off()




# Clustering promoters/enhancers

chip_prom = rbind(chip_prom0, chip_prom1, chip_prom2, chip_prom3, chip_prom4, chip_prom5, chip_prom6, chip_prom7, chip_prom8)
chip_prom[,1:3] = trimws(chip_prom[,1:3])
# names = apply(chip_prom[,1:3], 1, paste, collapse="_")
names = c()
for (row in 1:nrow(chip_prom)){
  name = paste(chip_prom[row,1:3], "_")
  names = c(names,name)
  if ((row %% 1000 ) == 0){
    print(paste0("On row ", row))
  }
}
rownames(chip_prom) = names
chip_prom = as.matrix(chip_prom[,-(1:4)])

chip_enh = rbind(chip_enh0, chip_enh1, chip_enh2, chip_enh3, chip_enh4, chip_enh5, chip_enh6, chip_enh7, chip_enh8)
names = c()
for (row in 1:nrow(chip_enh)){
  name = paste(chip_enh[row,1:3], "_")
  names = c(names,name)
}
rownames(chip_enh) = names
chip_enh = as.matrix(chip_enh[,-(1:4)])


# Center and scale data
chip_prom.norm <- (chip_prom - rowMeans(chip_prom))/rowSds(chip_prom + .5)
rownames(chip_prom.norm) = 1:nrow(chip_prom.norm)
chip_enh.norm <- (chip_enh - rowMeans(chip_enh))/rowSds(chip_enh + .5)
rownames(chip_enh.norm) = 1:nrow(chip_enh.norm)

# Set the number of clusters
count.k <- 8

# Set seed to preserve manual ordering
set.seed(733)

# Perform clustering
prom.count.cut = kclust(chip_prom.norm, k=count.k, type="k")
enh.count.cut = kclust(chip_enh.norm, k=count.k, type="k")

# Plot the individual k clusters

#png('/Users/phanstiel3/Desktop/promoter_clusters.png', width=8, height=8, units="in", res=300)
# Plot each cluster
par(mar=c(3,2,1,2))
par(mfrow=c(2,count.k/2))

for (i in 1:count.k) {
  
  # make empty plot
  plot(colMeans(chip_prom.norm[prom.count.cut==i,]), type="n",main=paste("n=",table(prom.count.cut)[i],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative expression")
  
  # and transparent grey lines
  apply(chip_prom.norm[prom.count.cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
  
  # add median line
  lines(colMeans(chip_prom.norm[prom.count.cut==i,]),col=k.colors[i],lwd=2)
  
  # add axis
  axis(1, at=1:count.k, labels=c("0",".5", "1", "1.5", "2", "4", "6", "24"))
}

# Set the desired order of clusters, based on their timing (manual; changes with seed)
count.cluster.order <- c(1,2,3,4,5,6,7,8)

# Plot heatmap
par(mfrow=c(1,1))
par(mar=c(5,4,4,5))

kHotmap(chip_prom.norm, prom.count.cut, count.cluster.order, k.colors = k.colors, maxval=3, title="Normalized Counts Clustering")

#dev.off()


#png('/Users/phanstiel3/Desktop/enhancer_clusters.png', width=8, height=8, units="in", res=300)
pdf(file="/Users/phanstiel3/Desktop/enhancer_clusters.pdf", width=12, height=6)

par(mar=c(3,2,1,2))
par(mfrow=c(2,count.k/2))

for (i in 1:count.k) {
  
  # make empty plot
  plot(colMeans(chip_enh.norm[enh.count.cut==i,]), type="n",main=paste("n=",table(enh.count.cut)[i],sep=""),ylim=c(-2,2), xaxt="n", xlab="Hours after LPS treatment", ylab="Relative expression")
  
  # and transparent grey lines
  apply(chip_enh.norm[enh.count.cut==i,],1,lines,col=adjustcolor("grey",alpha.f=0.1))
  
  # add median line
  lines(colMeans(chip_enh.norm[enh.count.cut==i,]),col=k.colors[i],lwd=2)
  
  # add axis
  axis(1, at=1:count.k, labels=c("0",".5", "1", "1.5", "2", "4", "6", "24"))
}

# Set the desired order of clusters, based on their timing (manual; changes with seed)
count.cluster.order <- c(1,2,3,4,5,6,7,8)

# Plot heatmap
par(mfrow=c(1,1))
par(mar=c(5,4,4,5))

kHotmap(chip_enh.norm, enh.count.cut, count.cluster.order, k.colors = k.colors, maxval=3, title="Normalized Counts Clustering")

dev.off()

ccf(x=colMeans(chip_enh.norm[enh.count.cut==1,]), y=colMeans(chip_prom.norm[prom.count.cut==5,]))
abline(v=0, col="red", lty=3)
ccf(x=IL1B_data$pro, y=IL1B_data$rna)
abline(v=0, col="red", lty=3)


# #######################################################################################################################
# 
# # Playing with HiC hubs:


hist(anchor1$start)
hist(anchor2$start)

length(anchor1$start[anchor1$start %in% anchor2$start])
head(hic)
head(anchor1, n=50)
head(anchor2)

length(unique(anchor1$start))
length(anchor1$start[!duplicated.default(anchor1$start)])

anchor1.unique = anchor1$start[!duplicated.default(anchor1$start)]
anchor1.unique.index = !duplicated.default(anchor1$start)

count=c()
maxval=c()
for (anchor in anchor1.unique){
  count = c(count, length(anchor2$start[anchor1$start == anchor]))
  maxval = c(maxval, max(anchor2$start[anchor1$start == anchor]))
}

anchor1.unique.multi = anchor1.unique[count > 1]
anchor1.unique.multi.2 = maxval[count > 1]

hist(count)
hist(maxval)

anchor1$end[anchor1.unique.index]

count = c()
for (anchor in unique(anchor1$start)){
  count = c(count, length(anchor1$start[anchor1$start == anchor]))
}




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

# # store for later; from RNA_Analysis.R
# rna.count.cut = count.cut
# rna.countMatrix = countMatrix.norm.combo

# use data from RNA_Analysis.R
identical(names(rna.count.cut), rownames(rna.countMatrix))

# Cut the normalized + combined count Matrix into 8 clusters, find mean of each
clusterExtract <- function(k, countMatrix, cut){
  clust.means = matrix(nrow=0, ncol=8)
  for (i in 1:k){
    clust = countMatrix[cut == i,]
    clust.mean = colMeans(clust)
    clust.means = rbind(clust.means, clust.mean)
  }
  colnames(clust.means) = time
  return(clust.means)
}

rna.means = clusterExtract(8, countMatrix.norm.combo, rna.count.cut)

# Plot the cluster means, and the predicted poly curve based on only the means for each cluster
clusterPlot <- function(k, means, color){
  df = data.frame(time=1:8)
  for (i in 1:k){
    fit = lm(means[i,] ~ poly(time, degree=dgr))
    curve = predict(fit, newdat=df)
    plot(time, means[i,])
    lines(df$time, curve, col=color)
  }
}

par(mfrow=c(2,4))
clusterPlot(8, rna.means, "blue")

# Extract the betas + residuals for each curve, to compare with betas for other clusters from other data sets (???)
clusterBetas <- function(k, means){
  betas = matrix(nrow=0, ncol=(dgr+1))
  for (i in 1:k){
    fit = lm(means[i,] ~ poly(time, degree=dgr))
    betas=rbind(betas, as.vector(coef(fit)))
  }
  colnames(betas) = names(coef(fit))
  return(betas)
}

clusterResid <- function(k, means){
  resid = matrix(nrow=0, ncol=8)
  for (i in 1:k){
    fit = lm(means[i,] ~ poly(time, degree=dgr))
    resid=rbind(resid, as.vector(resid(fit)))
  }
  return(resid)
}

betas = clusterBetas(8, rna.means)
resid = clusterResid(8, rna.means)


# CHIP PROMOTERS

# # store for later; from ChIP_Analysis.R
# chip.count.cut = count.cut
# chip.countMatrix = subsetMatrix.norm.combo

# use data from ChIP_Analysis.R
identical(names(chip.count.cut), rownames(chip.countMatrix))

chip.prom.means = clusterExtract(8,chip_prom.norm, prom.count.cut)
chip.enh.means = clusterExtract(8,chip_enh.norm, enh.count.cut)

par(mfrow=c(2,4))
clusterPlot(8, chip.prom.means, "red")
clusterPlot(8, chip.enh.means, "green")

chip.betas = clusterBetas(8, chip.means)
chip.resid = clusterResid(8, chip.means)

match = c(4, 5)


chip_prom.norm
prom.count.cut









