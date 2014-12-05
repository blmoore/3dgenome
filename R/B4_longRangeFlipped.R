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
  geom_vline(xintercept=.5, col=I("grey40"), linetype=2) +
  scale_x_continuous(limits=c(.3, .7)) +
  scale_colour_manual(values=c("#0000ff", "#FFA500", "#ff0000")) +
  facet_grid(flip~., scales="free", space="free") +
  xlab("Proportion of interactions with active compartments") +
  ylab("Megabase regions of variable structure") +
  labs(colour="Cell type") + geom_point(size=I(2.5)) + 
  theme(legend.position=c(.85,.8))
dev.off()

## rowSums across all flipped regions

res <- data.frame(rbind(cbind("A", rowSums(gflip[gflip$state == 1,-ncol(gflip)][,1])),
           cbind("B", rowSums(gflip[gflip$state == 2,-ncol(gflip)][,1]))))

res$X2 <- as.numeric(as.character(res$X2))
res$chr <- factor(gsub("\\..*", "", rownames(res)), levels=paste0("chr", 1:22))
res$pos <- as.numeric(gsub(".*\\.", "", rownames(res)))


ggplot(res, aes(x=pos/1e6, y=X2, fill=X1, col=X1)) + 
  geom_bar(stat="identity") +
  facet_wrap(~chr, scales="free_x") +
  theme_bw() + theme(legend.position=c(.9,.1)) 

group_by(res, X1) %>% summarise(m=sum(X2) / n())

## average A / B contacts for an A / B compartment genome-wide
## compare flipped regions to there, do they change their long-
## range contact profiles?
# remember A compartments make lots more contacts

a <- which(callStates(g.dat)$state == 1)
b <- which(callStates(g.dat)$state == 2)
g.a <- gm.mb[,a]
par(mfrow=c(2,1))
# A compartments interaction frequencies with other A compartments::
plot(density(rowMeans(gm.mb[a,a])), col="darkgreen",
     main="Contact profile of active compartment")
# B compartments 
lines(density(rowMeans(gm.mb[a,b])), col="darkred")

g.b <- gm.mb[,-a]
plot(density(rowSums(g.b[a,])), col="darkgreen", xlim=c(1000,5000),
     main="Contact profile of inactive compartment")
lines(density(rowSums(g.b[-a,])), col="darkred", xlim=c(1000,5000))
# B contacts


sum(gm.mb[a,a]) / length(a)**2
sum(gm.mb[-a,-a]) / (nrow(gm.mb) - length(a))**2

