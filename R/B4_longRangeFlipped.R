library("dplyr")
library("ggplot2")
library("gridExtra")
library("reshape2")
library("blmR")

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


gm.mb <- readMbMat("gm")
h1.mb <- readMbMat("h1")
k5.mb <- readMbMat("k5")

## pF from 6_fig4bFlippedBoxplots.R
pF <- readRDS("data/rds/chromFeaturesFlipped.rds")

enhancers <- pF[pF$feature == "E",]
enhancers <- enhancers[order(enhancers$olap, decreasing=T),]
i <- as.character(head(enhancers[!grepl("none|closed", enhancers$id),], 60)$id)

bins <- data.frame()
for( d in paste0(c("g", "h", "k"), rep(c(".open.bed", ".closed.bed"),3))){
  n <- read.table(paste0("~/hvl/ice/bedfiles/", d))
  bins <- rbind(bins, n) 
}

bins <- bins[bins$V4 %in% i,]
bins <- bins[match(i,bins$V4),]


## get macthing eigs, filter those that aren't -, - vs. +
bins[,c("h", "g", "k")] <- 0
# for x in rows
for (i in 1:nrow(bins)) {
  bref <- paste0(bins[i,1] , "-", bins[i,2])
  eig.index <- which(rownames(h.dat) == bref)
  bins[i,5:7] <- c(h.dat[eig.index,1], g.dat[eig.index,1], k.dat[eig.index,1])
}

b <- cbind(bins, ucsc=paste0(bins[,1], ":", bins[,2] - 1.2e6, "-", bins[,3] + 1.2e6),
           s=sign(bins$h * bins$g * bins$k))
b <- b[b$s == 1,-ncol(b)]
b$ct <- substr(b$V4, 1, 1)

flip <- as.data.frame(group_by(b, ct) %>% mutate(c=1:n()), 105)
# flipped open regions ordered per cell type by number of enhancers

#paste0(flip[flip$ct == "g",1:2], collapse=".")
g.ids <- gsub(" ", "", apply(flip[flip$ct == "h",1:2], 1, paste0, collapse="."))

grabFlipped <- function(mat, ids, dat){
  flip <- mat[,colnames(mat) %in% ids]
  flip <- as.data.frame(flip[,match(ids, colnames(flip))])
  flip$states <- callStates(dat$eigen)$state
  f <- flip %>% group_by(states) %>% summarise_each(funs(mean))
  melt(f, id.var="states")
}
  
flippedCt <- function(celltype=c("h", "g", "k")){
  ids <- gsub(" ", "", apply(flip[flip$ct == celltype,1:2], 1, paste0, collapse="."))
  gf <- grabFlipped(gm.mb, ids, g.dat)
  hf <- grabFlipped(h1.mb, ids, h.dat)
  kf <- grabFlipped(k5.mb, ids, k.dat)
  
  # combine all
  gf$ct <- "Gm12878"
  hf$ct <- "H1 hESC"
  kf$ct <- "K562"
  
  all <- rbind(gf, hf, kf)
  all
}
  
  
g.flip <- flippedCt("g")
h.flip <- flippedCt("h")
k.flip <- flippedCt("k")

## Visualise actual contacts::
# pal <- colorRampPalette(c("white", "steelblue4"))(56)
# par(mfrow=c(14,2), mar=rep(1, 4))
# for(i in 1:ncol(g.flip)){
# #for(i in 1:2){
#   compare <- as.matrix(data.frame(g.flip[,i], 
#                                   h.flip[,i], 
#                                   k.flip[,i]))
#   image(log10(compare + 1),
#         frame=F, axes=F, col=pal)
# }

props <- . %>% group_by(ct, variable, states) %>% 
  summarise(tot=sum(value)) %>% mutate(t=tot/sum(tot))

g.flip <- props(g.flip)
h.flip <- props(h.flip)
k.flip <- props(k.flip)

g.flip$flip <- "Gm12878"
h.flip$flip <- "H1 hESC"
k.flip$flip <- "K562"

flipped <- rbind(g.flip, h.flip, k.flip)

# chrX.longnum -> chrX:10-11
chr <- gsub("\\..*", "", as.character(flipped$variable))
start <- as.numeric(gsub(".*\\.", "", as.character(flipped$variable)))/1e6
end <- start+1
flipped$variable <- paste0(chr, ":", start, "-", end)

# ordering for variable factor
var <- flipped %>% subset(ct == flip & states == 2) %>% 
  group_by(ct) %>% arrange(t)
active <- subset(flipped, states==2)
active$variable <- factor(active$variable, levels=var$variable, ordered=T)

pdf("figures/suppl/flipped_Longrange.pdf", 6, 8)
ggplot(active, aes(x=t, col=ct, y=variable)) + 
  geom_point(size=I(0)) + #theme_bw() +
  #geom_vline(xintercept=.5, col=I("grey40"), linetype=2) +
  # GM average
  geom_vline(xintercept=.466, col=I("#0000ff98"), linetype=2) +
  # H1 average
  geom_vline(xintercept=.513, col=I("#FFA50098"), linetype=2) +
  # K5 average
  geom_vline(xintercept=.434, col=I("#ff000098"), linetype=2) +
  scale_x_continuous(limits=c(.3, .7)) +
  scale_colour_manual(values=c("#0000ff", "#FFA500", "#ff0000")) +
  facet_grid(flip~., scales="free", space="free") +
  xlab("Proportion of interactions with active compartments") +
  ylab("Megabase regions of variable structure") +
  labs(colour="Cell type") + geom_point(size=I(2.5)) + 
  theme(legend.position=c(.85,.8))
dev.off()

## lines on the plot:
# 1) Average B compartment
# 2) Average A compartment
# 3) 50%

## none flipped::
gstates <- callStates(g.dat$eigen)$state
stopifnot(length(gstates) == nrow(gm.mb))
hstates <- callStates(h.dat$eigen)$state
kstates <- callStates(k.dat$eigen)$state

# all contacts with A
plot(density(rowSums(gm.mb[,gstates == 2]) / rowSums(gm.mb)), 
     xlim=c(0,1), ylim=c(0, 8))
lines(density(rowSums(h1.mb[,hstates == 2]) / rowSums(h1.mb)))
lines(density(rowSums(k5.mb[,kstates == 2]) / rowSums(k5.mb)))

# .466
mean(rowSums(gm.mb[,gstates == 2]) / rowSums(gm.mb))
# .513
mean(rowSums(h1.mb[,hstates == 2]) / rowSums(h1.mb))
# .434
mean(rowSums(k5.mb[,kstates == 2]) / rowSums(k5.mb))

library("plotrix")
1.96*std.error(rowSums(gm.mb[,gstates == 2]) / rowSums(gm.mb))
1.96*std.error(rowSums(h1.mb[,hstates == 2]) / rowSums(h1.mb))
1.96*std.error(rowSums(k5.mb[,kstates == 2]) / rowSums(k5.mb))
