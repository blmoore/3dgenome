

readMbMat <- function(ct=c("gm", "h1", "k5")){
  A <- read.csv2(paste0("~/hvl/hicmap/", ct, "_gw_1mb.mat"), 
            sep="\t", stringsAsFactors=F, header=T)
  A <- apply(A[,-c(1,2)], 2, as.numeric)
}

gm.mb <- readMbMat("gm")
h1.mb <- readMbMat("h1")
k5.mb <- readMbMat("k5")

str(gm.mb)

# flipped regions of interest

#gm: chr5:158,000,000 (1043)
boi <- which(colnames(gm.mb) == "chr5.158000000")

#h1: chr20:21,000,000 (2752)
boi <- which(colnames(gm.mb) == "chr20.21000000")

#k5: chr8: 106,000,000 (1504)
boi <- which(colnames(k5.mb) == "chr8.106000000")

boi <- 2752

r50 <- function(bin)
  (bin-100):(bin+100)

image(cbind(h1.mb[r50(boi),boi], k5.mb[r50(boi),boi], gm.mb[r50(boi),boi]))

image(cbind(h1.mb[,boi], k5.mb[,boi], gm.mb[,boi]))

cor(k5.mb[r50(boi),boi], gm.mb[r50(boi), boi])
cor(k5.mb[r50(boi),boi], h1.mb[r50(boi), boi])
cor(gm.mb[r50(boi),boi], h1.mb[r50(boi), boi])

plot(k5.mb[r50(boi),boi], type="b", ylim=c(-3,3), col="darkred", pch=20)
lines(h1.mb[r50(boi),boi], type="b", col="orange", pch=20)
lines(gm.mb[r50(boi),boi], type="b", col="darkblue", pch=20)

plot(k5.mb[r50(boi),boi], h1.mb[r50(boi),boi])

# data.frame conversion for rle class
cor(colSums(h1.mb), colSums(gm.mb))

h.dat <- readRDS("hvl/production/3dgenome/data/rds/H1hesc_35Vars.rds")
library("blmR")
hstates <- callStates(h.dat$eigen)
nrow(h1.mb)

offchr <- which( gsub("\\..*", "", colnames(h1.mb)) %in% paste0("chr", c("X", "Y", "M")))

empty <- which(rowSums(h1.mb) == 0)

h1 <- h1.mb[,which(colnames(h1.mb) %in% gsub("-", ".", rownames(h.dat)))]
dim(h1)

sum(h1.mb[,boi])
sum(gm.mb[,boi])
sum(k5.mb[,boi])

library("fields")
image(apply(cbind(h1.mb[,boi], k5.mb[,boi], gm.mb[,boi]), 2, function(x) cut(x, seq(1:length(x), by=10)))

as.data.frame.rle <- function(x, ...) do.call(data.frame, x)

# length of chrs in bins:
h.bins <- as.data.frame(rle(gsub("-.*", "", h1.mb[,1])))
h.bins$total <- cumsum(h.bins$lengths)

k.bins <- as.data.frame(rle(gsub("-.*", "", k5.mb[,1])))
k.bins$total <- cumsum(k.bins$lengths)

cbind(h.bins,
k.bins)


hist(h1.mb)
length(h1.mb[h1.mb == 0]) / length(h1.mb[])
