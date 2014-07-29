
library("R.matlab")
source("~/hvl/R/hvl.R")
c18 <- readMat("hvl/hiclib/hm/h140kb/(17, 17).mat")
image(cor(log10(as.matrix(c18$`(17, 17)`)+1)))
dev.off()

chr <- readMat("~/hvl/hiclib/hm/k5_1mb_corr/chromosomeIndex.mat")$chromosomeIndex
bin <- readMat("~/hvl/hiclib/hm/k5_1mb_corr/positionIndex.mat")$positionIndex

chr <- chr + 1
chr[chr == 23]  <- "X"

#paste0("chr", chr, ":", bin, "-", bin+1e6-1)
matchnames <- paste0("chr", chr, "-", bin)

gm.pc <- as.vector(t(read.delim("~/hvl/hiclib/hm/gm_corrected_hm_1Mb.hdf5_PCA.out", sep=" ", header=F))[,1])
h1.pc <- t(read.delim("~/hvl/hiclib/hm/h1_corrected_hm_1Mb.hdf5_PCA.out", sep=" ", header=F))[,1]
k5.pc <- t(read.delim("~/hvl/hiclib/hm/k5_corrected_hm_1Mb.hdf5_PCA.out", sep=" ", header=F))[,1]
gen <- readMat("~/hvl/hiclib/hm/h11mb/heatmap.mat")$heatmap

# get 0 bins (excluded from PC calculation)
rownames(gen) <- colnames(gen) <- matchnames
cut <- rownames(gen)[which(colSums(gen) == 0)]

names(gm.pc) <- matchnames[!matchnames %in% cut]
names(k5.pc) <- matchnames[!matchnames %in% cut]
names(h1.pc) <- matchnames[!matchnames %in% cut]

both <- intersect(names(gm.pc), rownames(h.dat))

plot(gm.pc[names(gm.pc) %in% both], type="l")
par(new=T)
plot(g.dat[rownames(h.dat) %in% both,]$eigen, type="l", col="blue")

testhd <- h.dat
testhd <- h.dat[rownames(h.dat) %in% both,]
testhd$eigen <- h1.pc[names(h1.pc) %in% both]

tmod <- modelEigens(testhd)
plotPredRes.homer(tmod)
barsWError.3(tmod$imp)


testgd <- g.dat
testgd <- g.dat[rownames(g.dat) %in% both,]
testgd$eigen <- gm.pc[names(gm.pc) %in% both]

tmod <- modelEigens(testgd)
plotPredRes.homer(tmod)
barsWError.3(tmod$imp)


testkd <- k.dat
testkd <- k.dat[rownames(k.dat) %in% both,]
testkd$eigen <- k5.pc[names(k5.pc) %in% both]

tmod <- modelEigens(testkd)
plotPredRes.homer(tmod)
barsWError.3(tmod$imp)

# write bed file of new bins:
options(scipen=99)
bed <- data.frame(chr   = gsub("-.*", "", names(k5.pc)),
                  start = as.numeric(gsub(".*-", "", names(k5.pc))),
                  end   = as.numeric(gsub(".*-", "", names(k5.pc))) + 1e6-1,
                  id    = names(k5.pc))
write.table(bed, "~/hvl/ice/bins.bed", sep="\t", quote=F,
            row.names=F, col.names=F)


### 100 kb bins:
g100.pc <- read.pc("~/hvl/hiclib/hm/gm_corrected_hm_100kb.hdf5_PCA.out")
h100.pc <- read.pc("~/hvl/hiclib/hm/h1_corrected_hm_100kb.hdf5_PCA.out")
k100.pc <- read.pc("~/hvl/hiclib/hm/k5_corrected_hm_100kb.hdf5_PCA.out")

pos <- readMat("~/hvl/hiclib/hm/k5100_new/positionIndex.mat")$positionIndex
cstart <- readMat("~/hvl/hiclib/hm/k5100_new/chromosomeStarts.mat")$chromosomeStarts
frags <- readMat("~/hvl/hiclib/hm/k5100_new/frags.mat")$frags
# Rows that have been removes for eigenvector calculation, i.e. those
# rows in the interaction matrix with 0 total mapped contacts. These
# are picked up by which(rowSums(mat) == 0) on the huge 100 kb res
# genome-wide interataction matrix (~30,000**2 cells)
rem <- read.csv("~/hvl/hiclib/hm/nToRemove.csv", row.names=1)
pos[-rem[,1]]

times <- diff(c(cstart, cstart[length(cstart)] + (length(pos) - cstart[length(cstart)])))

bins <- data.frame(chr = rep(paste0("chr", c(1:22, "X")), times),
                   bp  = pos)
bins$id <- paste(bins$chr, bins$bp, sep="-")

pc <- cbind(bins[-rem[,1],], gm=g100.pc, k5=k100.pc, h1=h100.pc)

## Trim end bins is overrun:
chr.sizes <- read.table("~/hvl/hg19.chrom.sizes.txt")[1:24,]
chr.sizes <- chr.sizes[!chr.sizes[,1] %in% "chrY",]
for(i in 1:nrow(chr.sizes)){
  c <- pc[as.character(pc$chr) == chr.sizes[i,1],]
  if((c[nrow(c),]$bp + 1e5) > chr.sizes[i,2]){
    print(paste(c[nrow(c),]$id, chr.sizes[i,2]))
  }
}

pc$end <- pc$bp + 1e5
pc[pc$id == "chr17-81100000",]$end <- 81195210

pc <- pc[,c(3,1,2,7,4,6,5)]
pc$strand <- "+"
kp <- pc[,c(1:4, 8, 7)]
hp <- pc[,c(1:4, 8, 6)]
gp <- pc[,c(1:4, 8, 5)]

options(scipen=99)
write.bed <- function(bed, file)
  write.table(bed, file, quote=F, row.names=F, 
              col.names=F, sep="\t")

write.bed(pc, "~/hvl/ice/100kb_bins.bed")

# For use in 2_figure1TadsDiagram.R
write.bed(kp, "~/hvl/ice/k100.bed")
write.bed(gp, "~/hvl/ice/g100.bed")
write.bed(hp, "~/hvl/ice/h100.bed")


