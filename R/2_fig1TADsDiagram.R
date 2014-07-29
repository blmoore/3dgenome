############### draw TADs diagram ####################
# Show TADs called for a section of a chromosome,    #
# assess how well-conserved they are, draw a pretty  #
# picture relating TADs across the three cell lines. #
######################################################
library("scales")
library("plotrix")
library("R.matlab")
library("reshape2")
source("~/hvl/R/hvl.R")

# Call TADs using Bing Ren's software, published with
# Dixon et al. (2012). Currently available from his 
# lab website at: http://yuelab.org/hi-c/download.html
# (under Hi-C Domain Caller).
# TADs were called using recommended settings using 
# 40 kb resolution ICE-corrected contact matrices.

cmax <- read.delim("~/hvl/hg19.chrom.sizes.txt", header=F)

# 1) Build input matrices for calculating D.I.
for(ct in c("h1", "gm", "k5")){
  for (c in 0:22){
    cchr <- ifelse(c == 22, paste0("chr", "X"), paste0("chr", c + 1))
    if(!file.exists(paste0("~/hvl/ice/tads/inMats/", ct, "_", cchr, ".inmat"))){ 
      mf <- paste0("(", c, ", ", c, ")")
      c10 <- readMat(paste0("~/hvl/hiclib/hm/", ct, "40kb/", mf, ".mat"), verbose=F)
      mat <- c10[[1]]
      starts <- seq(0, (nrow(mat)-1)*4e4, by=4e4)
      ends <- seq(39999, nrow(mat)*4e4-1, by=4e4)
      #tail(cbind(starts, ends))
      
      max.chr <- cmax[cmax[,1] == cchr,2]
      if(ends[length(ends)] > max.chr){
        # If the last bin overruns, decrease it. 
        # First check nothing weird is going on:
        if(ends[length(ends)] > (max.chr + 4e4))
          stop("Error: co-ordinates overrun by >1 bin length")
        else
          ends[length(ends)] <- max.chr
      }
      
      options(scipen=99)
      colnames(mat) <- rownames(mat) <- paste0(cchr, ":", starts, "-", ends)
      
      mat <- as.data.frame(mat)
      m <- data.frame(chr   = as.character(cchr), 
                      start = starts,
                      end   = ends,
                      mat)
      cat("Writing", cchr, "of cell type", ct, "...\n")
      write.table(m, file=paste0("~/hvl/ice/tads/inMats/", ct, "_", cchr, ".inmat"),
                  row.names=F, col.names=F, quote=F, sep="\t")
    } else
      cat(paste0("File \"", ct, "_", cchr, ".inmat\" already exists, skipping ... \n"))
  }
}
  
# 2) Run Bing Ren's DI script on each file (takes ~5 min single-threaded)
call.di <- function(mat)
  system(paste0("perl hvl/ice/tads/brenSoft/perl_scripts/DI_from_matrix.pl ",
                mat, " 40000 2000000 ",
                "hvl/ice/tads/brenSoft/ucsc.hg19.fasta.fai > ",
                "hvl/ice/tads/di/", basename(mat)))

call.di <- Vectorize(call.di)
fl <- list.files("hvl/ice/tads/inMats/", full.names=T)
call.di(fl)

# 3) Restructure files and build HMM to call D.I. states
for(ct in c("h1", "gm", "k5")){
  cfl <- list.files("~/hvl/ice/tads/di/", pattern=paste0("^", ct, ".*inmat"), full.names=T)
  options(warn=-1)
  cfl <- cfl[order(as.numeric(gsub(".*chr(.+)\\..*", "\\1", cfl)))]
  options(warn=1)
  # chrX is now required as chr23:
  if(!file.exists(paste0("~/hvl/ice/tads/di/", ct, "_chr23.inmat"))){
    system(paste0("sed \'s/X/23/g\' ", cfl[length(cfl)], " > ", 
                  dirname(cfl[length(cfl)]), "/", ct, "_chr23.inmat"))
  }
  cfl <- cfl[-grep("chrX", cfl)]
  
  # concatenate all files together
  system(paste0("cat ", paste(cfl, collapse=" "), " > ~/hvl/ice/tads/di/", ct, "_all.di"))
}
# --- To use these .di files with Bing Ren HMM_calls.m MATLAB script,
# you need to hardcode input/output filenames for each. Can be run 
# using Octave with a bit of tweaking but I used MATLAB site install
# on cluster, i.e.: matlab -nojvm -nodesktop -r "run HMM_calls.m" ---
# Also runs under GNU Octave with minor adjustments and v3 of pmtk
  
tads <- read.table("~/hvl/tads/h1_tads.final")
tads.2 <- read.table("~/hvl/ice/tads/h1_allChr_tads.final")

plot.sizes <- function(tad.file){
  tsizes <- (tad.file$V3-tad.file$V2)/1e6
  plot(density((tad.file$V3-tad.file$V2)/1e6), xlab="Size of TAD (MB)", 
       main="Size distribution of H1 TADs")
  abline(v=median(tsizes), col="blue")
  text(median(tsizes)*1.03, y=.7, pos=4, 
       paste0("Median = ", signif(median(tsizes), 3)), col="blue")
  #abline(v=mean(tsizes))
}

plot.sizes(tads)
plot.sizes(tads.2)

drawTads <- function(di.file, tad.file, dat.file, state.file=NULL, ytop=200, chr="chr2", 
                     range=0:1000, ctcf=NA, ...){
  library("scales")
  # composite plot of DI, tads, A/B
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
                col=smoothColors(muted("green"),10,"white",10,muted("red"), alpha=90), border=NA)
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

## New version:
dev.off()
pdf("~/hvl/ice/plots/f1_22mbZoom.pdf", 6, 4.5)
csome = "chr2"
region = 1150:1600
par(mfrow=c(3,1), mar=c(2,5,2,2), xaxs="i")
drawTads("~/hvl/ice/tads/hmmout/h1_all.di",
         "~/hvl/ice/tads/h1_allChr_tads.final",
         "~/hvl/ice/h100.bed",
         ytop=400, main="H1 hESC", chr=csome, range=region)
drawTads("~/hvl/ice/tads/hmmout/gm_all.di",
         "~/hvl/ice/tads/gm_allChr_tads.final",
         "~/hvl/ice/g100.bed",
         ytop=30, main="GM12878", chr=csome, range=region)
drawTads("~/hvl/ice/tads/hmmout/k5_all.di",
         "~/hvl/ice/tads/k5_allChr_tads.final",
         "~/hvl/ice/k100.bed",
         ytop=30, main="K5 (new)", chr=csome, range=region)
dev.off()


## END figure 1

## Super-boundaries: is there an enhanced CTCF enrichment at
#  TAD + Compartment boundaries, relative to the average tad
#  boundary. What else is going on at these regions?

h.b <- read.table("~/hvl/tads/h1_allChr_tads.final")
h.b <- read.table("~/hvl/tads/k5_allChr_tads.final")

library(naturalsort)
h.b <- h.b[naturalorder(h.b[,1]),]
all.bounds <- h.b
bounds <- sapply(paste0("chr", 1:22), function(x) {
  h <- h.b[h.b[,1] == x,]
  #between <- h[2:nrow(h),2] - h[1:(nrow(h)-1),3]
  new.row <- cbind(chr=x, start=h[1:(nrow(h)-1),3],
                   end=h[2:nrow(h),2])
  return(new.row)
})

bounds <- as.data.frame(do.call(rbind, bounds))
bounds$start <- as.numeric(as.character(bounds$start))
bounds$end <- as.numeric(as.character(bounds$end))
bounds$diff <- bounds[,3] - bounds[,2]

#nrow(bounds[bounds$diff < 400000,]) / nrow(bounds) # 90% < 400kb
#nrow(bounds[bounds$diff < 50000,]) / nrow(bounds) # 57% < 50kb

tad.b <- bounds[bounds$diff < 400000,]
tad.b$mids <- apply(tad.b[,2:3], 1, mean)
new.end <- tad.b$mids + 20000
new.start <- tad.b$mids - 20000
options(scipen=999)
bed.df <- data.frame(chr=tad.b$chr, start=new.start,
                     end=new.end, id=paste0(tad.b$chr, "-", new.start))
write.table(bed.df, file="~/hvl/tads/k5_tbounds.bed", quote=F, 
            row.names=F, col.names=F, sep="\t")

########################################################################

# compare with bren's tads, remember: different genome
# assemblies, mappers, HMM calls/EM build variance
brentads <- read.csv("~/hvl/tads/brenH1.tads.csv", header=F)
brentads <- brentads[brentads[,1] == "chr2",]
h1 <- h.b[h.b[,1] =="chr2",]
brendiff <- sapply(1:nrow(brentads), function(j){
  return(c(min(abs(brentads[j,2] - h1[,2])),
           min(abs(brentads[j,3] - h1[,3]))))
})
  
## ctcf bounds:
ctcf <- read.table("~/hvl/tads/h1_counds_ctcf.isect")
allc <- read.table("~/hvl/tads/ctcf//wgEncodeOpenChromChipH1hescCtcfPk.narrowPeak")

k100 <- read.table("~/hvl/hicmap/h/h.100k.PC1.txt")

g <- as.data.frame(matrix(0, ncol=4, nrow=nrow(k100)*10))
for( i in 1:nrow(k100) ){
  chr <- gsub("-.*", "", k100[i,1])
  perbins <- seq(from=k100[i,3], to=k100[i,4], length.out=11)[1:10]
  perbine <- perbins + 10000
  id <- paste0(chr, "-", perbins)
  g[i*10-9:(i+9),] <- cbind(chr, perbins, perbine, id)
}


ctcf.dat <- read.table("~/hvl/tads/h1.10k.ctcf.bed")
c2d <- ctcf.dat[ctcf.dat[,1] == "chr2",]
points(((c2d[,2] + c2d[,3]) /2) /1e6, c2d[,5], type="h")

