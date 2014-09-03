############### draw TADs diagram ####################
# Show TADs called for a section of a chromosome,    #
# assess how well-conserved they are, draw a pretty  #
# picture relating TADs across the three cell lines. #
######################################################
library("plotrix")
library("reshape2")
library("scales")
library("blmR")

# Call TADs using Bing Ren's software, published with
# Dixon et al. (2012). Currently available from his 
# lab website at: http://yuelab.org/hi-c/download.html
# (under Hi-C Domain Caller).
# TADs were called using recommended settings using 
# 40 kb resolution ICE-corrected contact matrices.

# Results are included:
#   * .di       directionality index (BR script)
#   * .tads     placement of TADs / boundaries (BR script)
#   * x100.bed  100 kb resolution eigenvectors (from ICE)

# Bing Ren's HMM_calls.m MATLAB script can run under
# GNU Octave wiht minor adjustments and v3 of pmtk
  
plot.sizes <- function(tad.file){
  # Look at size distribution of called tads
  tsizes <- (tad.file$V3-tad.file$V2)/1e6
  plot(density((tad.file$V3-tad.file$V2)/1e6), xlab="Size of TAD (MB)", 
       main="Size distribution of H1 TADs")
  abline(v=median(tsizes), col="blue")
  text(median(tsizes)*1.03, y=.7, pos=4, 
       paste0("Median = ", signif(median(tsizes), 3)), col="blue")
  #abline(v=mean(tsizes))
}

tads <- read.table("data/text/h1.tads")
plot.sizes(tads)

drawTads <- function(di.file, tad.file, dat.file, 
                     state.file=NULL, ytop=200, chr="chr2", 
                     range=0:1000, ctcf=NA, ...){
  # build composite plot of DI, tads, A/B
  all.di <- read.table(di.file)
  #all.di <- all.di[order(all.di[,1]),]
  c1.di <- all.di[all.di[,1] == gsub("chr","",chr),]
  colnames(c1.di) <- c("chr", "start", "end", "di")
  tads <- read.table(tad.file)
  t1 <- tads[tads$V1 == chr,]
  xcords <- apply(c1.di[,2:3], 1, mean) / 1e6

  d <- c1.di[,4]
  cols <- ifelse(d < -max(d)/50, muted("red"),
                 ifelse(d > max(d)/50, muted("green"), "grey40"))
  
  plot(xcords[range], d[range], type="n", xlab="Position along chr1 (Mb)",
       ylab="Directionality index", frame=F, ylim=c(-ytop, ytop), axes=F, ...)
  axis(2, at=c(-(ytop*.9), 0, ytop*.9), cex.axis=.8,
       labels=c("Upstream", "0", "Downstream"))
  axis(1, at=seq(round(min(xcords)), round(max(xcords)), by=5), 
       labels=seq(round(min(xcords)), round(max(xcords)), by=5))
  abline(h=0, col="darkgrey")
  gradient.rect(xleft=t1$V2/1e6, ybottom=-ytop*1.2, ytop=ytop*1.2, xright=t1$V3/1e6,
                col=smoothColors(muted("green"), 10,
                                 "white", 10,
                                 muted("red"), alpha=90), border=NA)
  s <- seq(length(xcords)-1)
  segments(x0=xcords[s], x1=xcords[s+1], y0=d[s], y1=d[s+1], col=cols, lwd=1)
  if(!is.null(state.file)){
    #debugging
    sf <- read.table(state.file)
    points(t1[,2]/1e6, rep(30, nrow(t1)), pch="x", col="red")
    points(xcords[range], (sf[range,4]-2)*20)
  }  
  # now get A/B compartments for overlay:
  heigs <- read.table(dat.file)
  heigs <- heigs[heigs$V2 == chr,]
  heigs$V6 <- scale(heigs$V6) * (ytop/4)
  points(((heigs$V3 + heigs$V4)/2)/1e6, smooth(heigs$V6), type="l", col="blue")
  
  if(!is.na(ctcf)){
    allc <- read.table(ctcf)
    allc <- allc[allc[,1] == chr,]
    ycord <- replicate(nrow(allc), -(ytop*.9)+rnorm(1, ytop/40, ytop/40))
    points(((allc[,2] + allc[,3])/2)/1e6,  ycord, 
           pch=15, col=rgb(0,0,0,.3), cex=.4)
    rm(allc)
  }
}

## Zoom to 22Mb highlighted section:
cat("Drawing: figures/f1_ZoomedRegion.pdf\n")
pdf("figures/f1_ZoomedRegion.pdf", 6, 4.5)
csome = "chr2"
region = 1150:1600
par(mfrow=c(3,1), mar=c(2,5,2,2), xaxs="i")
drawTads("data/text/h1.di",
         "data/text/h1.tads",
         "data/bedfiles/h100.bed",
         ytop=400, main="H1 hESC", chr=csome, range=region)
drawTads("data/text/gm.di",
         "data/text/gm.tads",
         "data/bedfiles/g100.bed",
         ytop=30, main="GM12878", chr=csome, range=region)
drawTads("data/text/k5.di",
         "data/text/k5.tads",
         "data/bedfiles/k100.bed",
         ytop=30, main="K5 (new)", chr=csome, range=region)
dev.off()

## END figure 1
