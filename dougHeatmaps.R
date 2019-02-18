library('readr')
library('Sushi')
library('Sushi2')
library('devtools')
library('matrixStats')


rna = data.frame(read_tsv('~/Phanstiel Lab Dropbox/Kathleen Metz/Work/Research/Data/Projects/LIMA/rna/diff/LIMA_RNA_Subset-genes.tsv'))
rna.order = rna[order(rna$cluster),]
rna.order =rbind(rna[which(rna$cluster == 8),],
                 rna[which(rna$cluster == 3),],
                 rna[which(rna$cluster == 2),],
                 rna[which(rna$cluster == 7),],
                 rna[which(rna$cluster == 4),],
                 rna[which(rna$cluster == 6),],
                 rna[which(rna$cluster == 5),],
                 rna[which(rna$cluster == 1),]
)

rna.counts = as.matrix(rna.order[,4:11])
rna.4plot = (rna.counts  - rowMeans(rna.counts))/rowSds(rna.counts + 1)

# set the colors
rna.colors = colorRampPalette(c("grey","black"))(9 + 1)
names(rna.colors) = c(8,3,2,7,4,6,5,1)
colnames=as.character(rna.order$cluster)
rna.color = as.character(rna.colors[colnames])
colnames(rna.4plot) =c(0,.5,1,1.5,2,4,6,24)

chp = data.frame(read_tsv('~/Phanstiel Lab Dropbox/Kathleen Metz/Work/Research/Data/Projects/LIMA/chip/h3k27ac/diff/LIMA_ChIP_H3k27ac_Subset-peaks.tsv'))
#chp.order = chp[order(chp$cluster),]
chp.order =rbind(
  chp[which(chp$cluster == 7),],              
  chp[which(chp$cluster == 6),],
  chp[which(chp$cluster == 3),],
  chp[which(chp$cluster == 2),],
  chp[which(chp$cluster == 8),],
  chp[which(chp$cluster == 4),],
  chp[which(chp$cluster == 5),],
  chp[which(chp$cluster == 1),]
)
chp.counts = as.matrix(chp.order[,4:11])
chp.4plot = (chp.counts  - rowMeans(chp.counts))/rowSds(chp.counts + 1)
colnames(chp.4plot) =c(0,.5,1,1.5,2,4,6,24)

# set the colors
chp.colors = colorRampPalette(c("grey","black"))(9 + 1)
names(chp.colors) = c(7,6,3,2,8,4,5,1,0)
colnames=as.character(chp.order$cluster)
chp.color = as.character(chp.colors[colnames])


rna.4plot[which(rna.4plot > 3)] = 3
rna.4plot[which(rna.4plot < -3)] = -3
chp.4plot[which(chp.4plot > 3)] = 3
chp.4plot[which(chp.4plot < -3)] = -3

scalefactor = 1.0
pdf("~/Desktop/test.pdf",height=5,width=10)
#png("~/Desktop/test.png",width = 10, height = 5,res = 300,units = 'in')
par(mfrow=c(1,2),mgp=c(3,.2,0),mar=c(3,3,2,4))
hotmap(rna.4plot,labrow = FALSE,xlab.cex = 0.75,rowcolors = rna.color)
addlegend(range=c(-3,3),
          palette=colorRampPalette(c("deepskyblue2","black", "gold")),
          title.cex = .5,txt.cex = .5,
          title="Log2 (fold-change compared to mean)",side="right",title.offset = 0.075, 
          bottominset=0.6,topinset=0.00,xoffset=0.025,labelside="right",width=0.025)
mtext(side=3,line=.5,text="Differential Genes (n = 5,848)")
mtext(side=1,line=1.5,text=expression(paste("Hours after LPS/IFN",gamma," Treatment",sep="")))
mtext(side=1,adj=1.03,text="cluster",line=.2,cex=0.75)
labelplot("A")

hotmap(chp.4plot,labrow = FALSE,xlab.cex = 0.75,rowcolors = chp.color)
addlegend(range=c(-3,3),
          palette=colorRampPalette(c("deepskyblue2","black", "gold")),
          title.cex = .5,txt.cex = .5,
          title="Log2 (fold-change compared to mean)",side="right",title.offset = 0.075, 
          bottominset=0.6,topinset=0.00,xoffset=0.025,labelside="right",width=0.025)
mtext(side=3,line=.5,text="Differential H3K27ac Peaks (n = 59,830)")
mtext(side=1,line=1.5,text=expression(paste("Hours after LPS/IFN",gamma," Treatment",sep="")))
mtext(side=1,adj=1.025,text="cluster",line=.2,cex=0.75)
labelplot("B")
dev.off()

