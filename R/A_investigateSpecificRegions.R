###### Find an interesting example of a flipped block #######
# Find an interesting Mb that's flipped in a cell type, and #
# shows some enhancer enrichment, as a specific example of  #
# the general overall trend shown in Figure 4.              #
#############################################################
library("dplyr")
library("scales")
## nb. custom R package,
## devtools::install_github('blmoore', "blmR")
library("blmR")
data(g.dat)
data(h.dat)
data(k.dat)

source("~/hvl/R/hvl.R")
options(scipen=99)

read.feats <- function(fn){
  feats <- read.table(fn)
  colnames(feats) <- c("chr", "start", "end", "block", "chr.e",
                         "start.e", "end.e", "type", "length.e")
  feats
}

h1.open <- read.feats("~/hvl/ice/bedfiles/h.open.bed_cts.out")
h1.open$coords <- with(h1.open, paste0(chr,":", start,"-", end))
h.e <- group_by(h1.open, coords, type) %.% summarise(count=n())

ecounts <- subset(h.e, type =="E")

ecounts[order(ecounts$count, decreasing=T),]
# 1 PAX1

eigs <- data.frame(chr=gsub("-.*", "", rownames(h.dat)),
                   start=gsub(".*?-", "", rownames(h.dat)),
                   h=h.dat$eigen, g=g.dat$eigen, k=k.dat$eigen,
                   hstates=callStates(h.dat$eigen)$state,
                   gstates=callStates(g.dat$eigen)$state,
                   kstates=callStates(k.dat$eigen)$state)
head(eigs)
# chr20:21000000-22000000    E    43
ind <- which(eigs$chr == "chr20" & eigs$start == 21000000)
plot <- eigs[(ind-5):(ind+5),]
plot(as.numeric(as.character(plot$start)), plot$h, type="l", col=muted("orange"))
lines(as.numeric(as.character(plot$start)),plot$k, col=muted("red"))
abline(h=0, lty=2, col="red")
lines(as.numeric(as.character(plot$start)), plot$g, col=muted("blue"))
rect(21e6, -1, 22e6, 1)

# find consecutive open H1 regions:
diff(as.numeric(as.character(eigs$start)))


## Old code:
library("caret")
options(scipen=999)
source("~/hvl/R/hvl.R")
library("Gviz")
library("IRanges")
library("reshape2")
library("RColorBrewer")
devtools::load_all("blmR")

## 1) Get OrderEigs ( from fig4aStratifyByVariability:
options(scipen=999)
justEigs <- data.frame(gm=g.dat$eigen, h1=h.dat$eigen, k5=k.dat$eigen)
justEigs <- cbind(justEigs, mad=apply(justEigs,1,sd))
rownames(justEigs) <- rownames(h.dat)
justEigs[,4] <- apply(justEigs[,1:3],1,function(x) mad(x))
orderedEigs <- justEigs[order(justEigs[,4], decreasing=T),]

## 3) blockHighlight() form hvl.R:
rankToBlock <- function(rank, s=justEigs, oEigs=orderedEigs){
  which(rownames(s) == rownames(oEigs)[rank])
}
rankToBlock(2)

rownameToBlock <- function(rowname, oEigs=orderedEigs){
  which(rownames(oEigs) == rowname)
}
rownameToBlock(rownames(h.dat)[rankToBlock(2)])

# ## 4) Gviz full visualisation
# source("http://bioconductor.org/biocLite.R")
blockPlot <- function(rank, eigTable=justEigs, 
                      ordered=orderedEigs, offset=3){
#     eigTable <- justEigs
#     ordered <- orderedEigs
#     rank <- 2
#     offset <- 3
  gen <- "hg19"
  e.s <- 50000 #amount to add either side of block
  #block <- rownames(eigTable)[rankToBlock(rank, oEigs=ordered)]
  block <- rownames(ordered)[rank]
  csome <- as.character(gsub("chr","",gsub("-.*","",block)))
  start.c <- as.numeric(gsub(".*-","",block))
  start.pos <- if( (start.c - offset*1000000) > 0) start.c - (offset*1000000) else 0 
  end.c <- start.c+1e6
  end.pos <- end.c + offset*1000000
  #pall <- brewer.pal(5, "Dark2")
  pall <- c("#0000ff98", "#FFA50098", "#ff000098")
  ##Ideogram:
  ideo <- IdeogramTrack(genome=gen, chromosome=csome)
  ##GenomeAxis:
  gtrack <- GenomeAxisTrack(range=IRanges(start=start.c, end=end.c, names=" "), 
                            fill.range="yellow", add53=T, add35=T, littleTicks=T)
  ##DataTracks:
  c.only <- eigTable[gsub("chr","",gsub("-.*","",rownames(eigTable))) == csome,]
  #m.c <- melt(c.only)
  dt.s <- as.numeric(gsub(".*-", "", rownames(c.only)))
  dt.e <- dt.s+1e6
  #colnames(c.only)
  H1 <- DataTrack(chromosome=csome, start=dt.s, end=dt.e, showAxis=T, baseline=0,
                  data=c.only[,-ncol(c.only)], name="Compartment eigenvectors", genome=gen, ylim=c(-.4,.4), 
                  col.baseline="grey40", background.title="grey40",
                  col=pall, type=c("g","a","p"), groups=colnames(c.only)[-ncol(c.only)], legend=T,
                  stackedBars=F, col.grid="grey90", cex.legend=.7)
  
  biom <- BiomartGeneRegionTrack(genome=gen, chromosome=csome, start=start.c-e.s, 
                                 end=end.c+e.s, name="ENSEMBL genes", geneSymbols=T,
                                 background.title="grey60", shape="box", stacking="squish",
                                 filters = "entrezgene", mergeGroups=T, lex=.5, cex=.5,
                                 min.distance=1)
  # RNA-seq for each cell line from ENCODE
  hexp <- DataTrack(range=h1.bw, from=start.pos, to=end.pos, genome=gen, 
                    chromosome=csome, type="h", name="H1 expr",
                    transformation=function(x) log2(x+1), col.histogram=pall[2],
                    fill.histogram=pall[2], window=-1, col=pall[2],
                    background.title="grey80", ylim=c(0,12), windowSize=2500)
  gexp <- DataTrack(range=gm.bw, from=start.pos, to=end.pos, genome=gen, 
                    chromosome=csome, type="h", name="GM expr",
                    transformation=function(x) log2(x+1), col.histogram=pall[1],
                    window=-1, fill.histogram=pall[1], col=pall[1],
                    background.title="grey80", ylim=c(0,12), windowSize=2500)
  kexp <- DataTrack(range=k5.bw, from=start.pos, to=end.pos, genome=gen, 
                    chromosome=csome, type="h", name="K5 expr",
                    col.histogram=pall[3], window=-1, col=pall[3],
                    transformation=function(x) log2(x+1), fill.histogram=pall[3],
                    background.title="grey80", ylim=c(0,12), windowSize=2500)
  # Pretty heatmap of mappability
 mappa <- DataTrack(range=map, from=start.pos, to=end.pos, genome=gen, 
                      chromosome=csome, type="heatmap", name="Mapability", col="grey40",
                      window=-1000, windowSize=5000, showAxis=F, strand="+",
                      aggregation="median",
                      background.title="white")  
  glist <- as.data.frame(attr(biom, "range"))
  scl <- GenomeAxisTrack(scale=1000000, labelPos="below")
  plotTracks(list(ideo, gtrack, H1, biom, gexp, hexp, kexp, mappa, scl), 
             from=start.pos, to=end.pos, showId=T,  collapseTranscripts=T,
             stackHeight=0.2, collapse=F, sizes=c(1.5,2,6,7,2,2,2,1,.1),
             transcriptAnnotation = "symbol")
}

source("http://bioconductor.org/biocLite.R")
biocLite()
library("Gviz")
library("IRanges")
library("GenomicRanges")
library("reshape2")
library("rtracklayer")
h1.bw <- "/Volumes/Data/HvL/rnaSeq/wgEncodeRegTxnCaltechRnaSeqH1hescR2x75Il200SigPooled.bw"
gm.bw <- "/Volumes/Data/HvL/rnaSeq/wgEncodeRegTxnCaltechRnaSeqGm12878R2x75Il200SigPooled.bw"
k5.bw <- "/Volumes/Data/HvL/rnaSeq/wgEncodeRegTxnCaltechRnaSeqK562R2x75Il200SigPooled.bw"
gc <- "/Volumes/Data/HvL/rnaSeq/gc5base.wig"
map <- "/Volumes/Data/HvL/rnaSeq/wgEncodeCrgMapabilityAlign50mer.bigWig"

pdf("~/hvl/ice/plots/supplementary/block2_variabiltyeg_v2.pdf", 5.5, 8)
blockPlot(2, justEigs, offset=3, ordered=orderedEigs)
dev.off()

for(i in 1:10){
  png(paste0("~/hvl/ice/plots/supplementary/stratVar_eg", i, ".png"), 5.5, 8, units="in", res=300)
  blockPlot(i, justEigs, offset=3, ordered=orderedEigs) 
  dev.off()
}

transcript(biom)

rownames(justEigs[rankToBlock(2)-1,])
## block chr22 : 22 Mb missing -- why is this...
rownames(h.dat)[grepl("chr22", rownames(h.dat))]

h1o <- read.table("~/hvl/ice/bedfiles/h.open.bed")
h1o <- sapply(1:nrow(h1o), function(i) paste0(h1o[i,1], "-", h1o[i,2]))
ranks <- sapply(h1o, rownameToBlock)
length(ranks)

for(i in 11:35){
  png(paste0("~/hvl/ice/plots/supplementary/hopen/ho_eg", i, ".png"), 5.5, 8, units="in", res=150)
  blockPlot(ranks[i], justEigs, offset=3, ordered=orderedEigs) 
  dev.off()
} 


rownameToBlock(1)
bedFromRownames(table=all.negative,outfile="~/hvl/allSignifNeg.bed")

conservation <- UcscTrack(genome = "hg19", chromosome = "chr1",
                          track = "Conservation", table = "phyloP46wayPlacental",
                          from=1000000, to=2000000, trackType = "DataTrack",
                          start = "start", end = "end", data = "score",
                          name = "Cons")


pdfplot <- function(rank){
  pdf(file=paste("~/hvl/plots/DivPlots/d2Plot_", rank, ".pdf", sep=""), 5.5,7.5)
  blockPlot(rank, s.in, offset=5, ordered=o5)
  dev.off()
}


## UCSC bedGraph files for custom tracks::
l <- rownames(h.dat)
options(scipen=999)

## 1) H1::
head <- gsub("\n", " ", 'track type=bedGraph name=\"H1 compartments\" description=center_label
                visibility=display_mode color=000,000,222 altColor=222,000,000
                priority=priority autoScale=on alwaysZero=on
                gridDefault=on maxHeightPixels=100:48:32
                graphType=bar viewLimits=-.4:.4
                yLineMark=0 yLineOnOff=off
                windowingFunction=mean smoothingWindow=off')
df <- data.frame(chrom=gsub("-.*", "", l),
           chromStart=gsub(".*-", "", l),
           chromEnd=as.numeric(gsub(".*-", "", l)) + 1e6,
           dataValue=h.dat$eigen)

sink("~/hvl/ice/bedfiles/customTrack_H1eigs.bedGraph")
cat(head, "\n")
sink()

## nb: chrom chr5 size (181000000 > 180915260
df[df$chrom == "chr5" & df$chromStart == 180000000,3]  <- 180915260
## chrom chr10 size (136000000 > 135534747)
df[df$chrom == "chr10" & df$chromStart == 135000000,3] <- 135534747
## chrom chr12 size (134000000 > 133851895)
df[df$chrom == "chr12" & df$chromStart == 133000000,3] <- 133851895
## chrom chr15 size (103000000 > 102531392)
df[df$chrom == "chr15" & df$chromStart == 102000000,3] <- 102531392

write.table(df, file="~/hvl/ice/bedfiles/customTrack_H1eigs.bedGraph",
            append=T, quote=F, col.names=F, row.names=F, sep="\t")


head <- gsub("\n", " ", 'track type=bedGraph name=\"GM12878 compartments\" description=center_label
                visibility=display_mode color=000,000,222 altColor=222,000,000
                priority=priority autoScale=on alwaysZero=on
                gridDefault=on maxHeightPixels=100:48:32
                graphType=bar viewLimits=-.4:.4
                yLineMark=0 yLineOnOff=off
                windowingFunction=mean smoothingWindow=off')
df$dataValue <- g.dat$eigen

sink("~/hvl/ice/bedfiles/customTrack_GMeigs.bedGraph")
cat(head, "\n")
sink()
write.table(df, file="~/hvl/ice/bedfiles/customTrack_GMeigs.bedGraph",
            append=T, quote=F, col.names=F, row.names=F, sep="\t")

head <- gsub("\n", " ", 'track type=bedGraph name=\"K562 compartments\" description=center_label
                visibility=display_mode color=000,000,222 altColor=222,000,000
                priority=priority autoScale=on alwaysZero=on
                gridDefault=on maxHeightPixels=100:48:32
                graphType=bar viewLimits=-.4:.4
                yLineMark=0 yLineOnOff=off
                windowingFunction=mean smoothingWindow=off')
df$dataValue <- k.dat$eigen

sink("~/hvl/ice/bedfiles/customTrack_K5eigs.bedGraph")
cat(head, "\n")
sink()
write.table(df, file="~/hvl/ice/bedfiles/customTrack_K5eigs.bedGraph",
            append=T, quote=F, col.names=F, row.names=F, sep="\t")

## ROI (from 6_.*?.R), cross-ref with g.open.bed etc.
head(enhancers[!grepl("none", enhancers$id),], 10)

#g.49:
chr6:11000000-12000000 (chr6:8800000-15500000)
# no: 44, 13, k.25, 
#g.39:
chr5:142000000-143000000 (chr5:140700000-145000000)



