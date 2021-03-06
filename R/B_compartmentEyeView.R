library("magrittr")
library("R.matlab")
library("RColorBrewer")
library("reshape2")
library("blmR")
options(scipen=999)

## Code from 7_fig5boundaryEnrichments.R ::
# call 1MB boundaries (100kb were used for boundary profiles)
h.1mb <- read.table("data/bedfiles/h1_eigs.bed")
k.1mb <- read.table("data/bedfiles/k5_eigs.bed")
g.1mb <- read.table("data/bedfiles/gm_eigs.bed")

h.1mb$states <- callStates(h.1mb[,6])$state
k.1mb$states <- callStates(k.1mb[,6])$state
g.1mb$states <- callStates(g.1mb[,6])$state

colnames(k.1mb) <- colnames(g.1mb) <- colnames(h.1mb) <- 
  c("id", "chr", "start", "end", "strand", "eig", "states")

callCBounds <- function(sfile){
  compartments <- sapply(paste0("chr", c(1:22, "X")), function(x) {
    cc <- sfile[sfile$chr == x,]
    rle <- rle(cc$states)
    # boundary blocks (i.e., last block of given state)
    bounds <- cc[cumsum(rle$lengths),]   
    bounds[,3:4] <- bounds[,3:4]
    return(bounds)
  }, simplify=F)
  compartments <- do.call(rbind, compartments)
  outbed <- compartments[,c("chr", "start", "end")]
  outbed$names <- paste0(compartments$chr, "-", compartments$start)
  return(outbed)
}

kcb <- callCBounds(k.1mb)
gcb <- callCBounds(g.1mb)
hcb <- callCBounds(h.1mb)

readMbMat <- function(ct=c("gm", "h1", "k5")){
  library("R.matlab")
  hmb <- readMat(paste0("~/hvl/hiclib/hm/", ct, "c_1mb/heatmap.mat"))$heatmap
  ind <- readMat(paste0("~/hvl/hiclib/hm/", ct, "c_1mb/positionIndex.mat"))$positionIndex
  chr <- readMat(paste0("~/hvl/hiclib/hm/", ct, "c_1mb/chromosomeIndex.mat"))$chromosomeIndex
  
  colnames(hmb) <- rownames(hmb) <- paste0("chr", chr+1, ".", ind)
  
  rm <- which(colSums(hmb) == 0, arr.ind=T) 
  hmb <- hmb[-rm, -rm]
  hmb
}

gm.gw <- readMbMat("gm")
h1.gw <- readMbMat("h1")
k5.gw <- readMbMat("k5")

plotLRI <- function(imat, compartments, fn){
  # Remove non-autosomal::
  rm <- substr(colnames(imat), 1, 4) %in% paste0("chr", c("X", "Y", "M"))
  imat <- imat[!rm,!rm]
    
  # format IDs per interaction matrix
  compartments <- paste0(compartments[,1], ".", compartments[,3]) 
  colnames(imat) <- gsub("chr23", "chrX", colnames(imat))
  breaks <- unique(c(1, which(colnames(imat) %in% compartments), nrow(imat)))
  
  # don't allow compartments to run over chromosomes
  csomes <- cumsum(rle(gsub("\\..*", "", colnames(imat)))$lengths)
  breaks <- sort(unique(c(breaks, csomes)))
  
  b <- melt(as.matrix(imat))
  # factor to numeric, bad idea
  b$Var1 <- as.numeric(factor(b$Var1))
  b$Var2 <- as.numeric(factor(b$Var2))
  
  bb <- with(b, tapply(value, list(
    y=cut(Var1, breaks=breaks, include.lowest=T),
    x=cut(Var2, breaks=breaks, include.lowest=T)),
    sum)
  )
  bsq <- melt(bb)
  
  getNum <- . %>%
    # rm brackets
    gsub("\\[|\\(|\\]|\\)", "", .) %>%
    # split digits and convert
    strsplit(",") %>%
    #xleft xright (?)
    unlist %>% as.numeric
  
  y <- t(sapply(bsq[,1], getNum))
  x <- t(sapply(bsq[,2], getNum))
  
  bsq$size <- (y[,2] - y[,1]) * (x[,2] - x[,1])
  bsq$norm <- bsq$value / bsq$size
  bsq$norm[is.na(bsq$norm)] <- 0

  pal <- colorRamp(brewer.pal(8, "RdBu")[2:7], interpolate="spline")
  u <- bsq$norm / max(bsq$norm)
  uni <- ecdf(u)(u)
  
  pdf(paste0("figures/suppl/", fn, "_lri.pdf"), 9, 9)
  par(mar=rep(2,4), mgp=rep(0,3), xaxs="i", yaxs="i")
  plot(1:nrow(imat), 1:nrow(imat), type="n", frame=F, axes=F, xlab="", ylab="")
  rect(ybottom=y[,1], ytop=y[,2],
       xleft=x[,1], xright=x[,2], 
       col=rgb(pal(uni), max=255),
       border=NA)
  
  abline(h=csomes, col="grey60", lwd=2)
  abline(v=csomes, col="grey60", lwd=2)
  axis(1, at=c(mean(c(1, csomes[1])), zoo::rollmean(csomes, 2)), 
       labels=c(1:22, "X"), tick=F, col="grey60", cex.axis=.7)
  axis(2, at=c(mean(c(1, csomes[1])), zoo::rollmean(csomes, 2)), 
       labels=c(1:22, "X"), tick=F, col="grey60", cex.axis=.7)
  axis(1, at=c(1, csomes), labels=NA, col="grey60", lwd=2, tck=-.005)
  axis(2, at=c(1, csomes), labels=NA, col="grey60", lwd=2, tck=-.005)
  dev.off()
  
  invisible(NULL)
}  
  
plotLRI(h1.gw, hcb, fn="h1gw")
plotLRI(k5.gw, kcb, fn="k5gw")
plotLRI(gm.gw, gcb, fn="gmgw")

# Remove non-autosomal::
rm <- substr(colnames(gm.gw), 1, 4) %in% paste0("chr", c("X", "Y", "M"))
gm.gw <- gm.gw[!rm,!rm]

# load compartment calls
#gm.c <- read.table("data/bedfiles/gm_compartmentbounds.bed")
#gm.c <- read.table("data/bedfiles/h1_compartmentbounds.bed")
gm.c <- gcb

# format IDs per interaction matrix
gm.c <- paste0(gm.c[,1], ".", gm.c[,3]) #rowMeans(gm.c[,2:3]))

breaks <- unique(c(1, which(colnames(gm.gw) %in% gm.c), nrow(gm.gw)))

b <- melt(as.matrix(gm.gw))
# factor to numeric, bad idea
b$Var2 <- as.numeric(b$Var2)

bb <- with(b, tapply(value, list(
  y=cut(Var1, breaks=breaks, include.lowest=T),
  x=cut(Var2, breaks=breaks, include.lowest=T)),
  sum)
)

bsq <- melt(bb)

getNum <- . %>%
  # rm brackets
  gsub("\\[|\\(|\\]|\\)", "", .) %>%
  # split digits and convert
  strsplit(",") %>%
  #xleft xright (?)
  unlist %>% as.numeric


y <- t(sapply(bsq[,1], getNum))
x <- t(sapply(bsq[,2], getNum))

bsq$size <- (y[,2] - y[,1]) * (x[,2] - x[,1])
bsq$norm <- bsq$value / bsq$size
bsq$norm[is.na(bsq$norm)] <- 0

## negative!
min(bsq$norm)
stopifnot(min(bsq$norm) < 5)
bsq$norm <- bsq$norm + 5

#pdf("figures/suppl/h1_longRange_c.pdf", 6, 6)
par(mar=rep(1,4))
plot(1:nrow(gm), 1:nrow(gm), type="n", frame=F, axes=F, xlab="", ylab="")
rect(ybottom=y[,1], ytop=y[,2],
     xleft=x[,1], xright=x[,2], 
     col=rgb(colorRamp(c("darkred", "white","steelblue4"), 
                       interpolate="spline")(bsq$norm / max(bsq$norm)), 
             alpha=255*(bsq$norm / max(bsq$norm)), max=255),
     border=NA)
     #border="white")

#dev.off()
library("RColorBrewer")
pal <- colorRamp(brewer.pal(9, "Blues")[c(1:3, 5:7)], interpolate="spline")

pdf("figures/suppl/gm_lri.pdf", 10, 10)
pal <- colorRamp(brewer.pal(8, "RdBu")[2:7], interpolate="spline")
u <- bsq$norm / max(bsq$norm)
uni <- ecdf(u)(u)

par(mar=rep(2,4), mgp=rep(0,3), xaxs="i", yaxs="i")
plot(1:nrow(gm.gw), 1:nrow(gm.gw), type="n", frame=F, axes=F, xlab="", ylab="")
rect(ybottom=y[,1], ytop=y[,2],
     xleft=x[,1], xright=x[,2], 
     col=rgb(pal(uni), max=255),
     border=NA)

csomes <- cumsum(rle(gsub("\\..*", "", colnames(gm.gw)))$lengths)

abline(h=csomes, col="grey50", lwd=2)
abline(v=csomes, col="grey50", lwd=2)


axis(1, at=c(mean(c(1, csomes[1])), zoo::rollmean(csomes, 2)), 
     labels=1:22, tick=F, col="grey60")
axis(2, at=c(mean(c(1, csomes[1])), zoo::rollmean(csomes, 2)), 
     labels=1:22, tick=F, col="grey60")
axis(1, at=c(1, csomes), labels=NA, col="grey60", lwd=2)
axis(2, at=c(1, csomes), labels=NA, col="grey60", lwd=2)
dev.off()


## map normal-ish data to uniform::
hist(bsq$norm)
n1 <- (bsq$norm - mean(bsq$norm)) / sd(bsq$norm)
n1 <- ecdf(n1)(n1)

dev.off()
par(mar=rep(2,4), mgp=rep(0,3), xaxs="i", yaxs="i")
plot(1:nrow(gm.gw), 1:nrow(gm.gw), type="n", frame=F, axes=F, xlab="", ylab="")
rect(ybottom=y[,1], ytop=y[,2],
     xleft=x[,1], xright=x[,2], 
     col=rgb(colorRamp(c("white", "steelblue4"))(n1), max=255),
     border=NA)
