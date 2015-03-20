library("dplyr")
library("ggplot2")
library("gridExtra")
library("reshape2")
library("blmR")
options(scipen=99)

## Generates supplementary long range contacts figure part A and B
## and (as yet unused) "compartments-eye view" of the genome plot.
g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

readMbMat <- function(ct=c("gm", "h1", "k5")){
  library("R.matlab")
  # 1Mb normalised contact matrices as output by hiclib / ICE
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
# 
# ## pF from 6_fig4bFlippedBoxplots.R
pF <- readRDS("data/rds/chromFeaturesFlipped.rds")

enhancers <- pF[pF$feature == "E",]
enhancers <- enhancers[order(enhancers$olap, decreasing=T),]
i <- as.character(enhancers[!grepl("none|closed", enhancers$id),]$id)

bins <- data.frame()
for( d in paste0(c("g", "h", "k"), rep(c(".open.bed", ".closed.bed"),3)) ){
  message(paste0("data/bedfiles/", d))
  n <- read.table(paste0("data/bedfiles/", d))
  bins <- rbind(bins, n) 
}

bins <- bins[bins$V4 %in% i,]
bins <- bins[match(i,bins$V4),]

# ## get macthing eigs, filter those that aren't -, - vs. +
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

flip <- as.data.frame(group_by(b, ct) %>% mutate(c=1:n()))

#------------------------------------------#
# split into long and short range contacts #
#------------------------------------------#
long_range_lab <- "Long range"
near_cis_lab <- "Regional (<2 Mb)"
cutoff <- 2e6

# v2 of grabFlipped() splits contacts into near-cis or long range (>10 Mb)
grabFlipped_longRange <- function(mat, ids, dat, means=T){
  flip <- mat[,colnames(mat) %in% ids]
  flip <- as.data.frame(flip[,match(ids, colnames(flip))])
  flip_states <- callStates(dat$eigen)$state

  chr <- gsub("\\..*", "", rownames(flip))
  bins <- as.numeric(gsub(".*?\\.", "", rownames(flip)))

  df <- data.frame(flipped=character(),
                   contacts=numeric(),
                   type=character(),
                   states=character(),
                   bins=character())
  for(i in 1:ncol(flip)){
    tchr <- gsub("\\..*", "", colnames(flip)[i])
    bin_dist <- abs(bins - as.numeric(gsub(".*?\\.", "", colnames(flip)[i])))
    # if distance > 20 Mb OR different chromosome, "long range", else "near cis"
    type <- ifelse(tchr != chr | bin_dist > cutoff, long_range_lab, near_cis_lab)
    #type[bin_dist == 0] <- NA # 0 at self-interaction doesn't change proportion
    odf <- data.frame(flipped=colnames(flip)[i],
                      contacts=flip[,i],
                      type=type,
                      states=flip_states,
                      bins=rownames(flip))
    df <- rbind(df, odf)
  }
  
  if(means){
    df %>% group_by(flipped, type, states) %>% summarise(m=mean(contacts))
  } 
  
  df
}


flippedCt_longRange <- function(celltype=c("h", "g", "k"), ...){
  ids <- gsub(" ", "", apply(flip[flip$ct == celltype,1:2], 1, paste0, collapse="."))
  gf <- grabFlipped_longRange(gm.mb, ids, g.dat, ...)
  hf <- grabFlipped_longRange(h1.mb, ids, h.dat, ...)
  kf <- grabFlipped_longRange(k5.mb, ids, k.dat, ...)
  
  # combine all
  gf$ct <- "Gm12878"
  hf$ct <- "H1 hESC"
  kf$ct <- "K562"
  
  all <- rbind(gf, hf, kf)
  all
}

# convert counts to proportions
props2 <- . %>% group_by(ct, flipped, type, states) %>% 
  summarise(tot=sum(contacts)) %>% mutate(t=tot/sum(tot))

gflr <- props2(flippedCt_longRange("g"))
hflr <- props2(flippedCt_longRange("h"))
kflr <- props2(flippedCt_longRange("k"))

gflr$flip <- "Gm12878"
hflr$flip <- "H1 hESC"
kflr$flip <- "K562"

flipped_2 <- rbind(gflr, hflr, kflr)

# chrX.longnum -> chrX:10-11
chr <- gsub("\\..*", "", as.character(flipped_2$flipped))
start <- as.numeric(gsub(".*\\.", "", as.character(flipped_2$flipped)))/1e6
end <- start+1
#flipped_2$flipped <- paste0(chr, ":", start, "-", end)
flipped_2$flipped <- paste0(gsub("chr", "", chr), ":", start)
  
# ordering for variable factor
var <- flipped_2 %>% subset(ct == flip & states == 2 & type == "Long range") %>% 
  group_by(ct) %>% arrange(t)
active_2 <- subset(flipped_2, states==2 & type == "Long range")
active_2$flipped <- factor(active_2$flipped, levels=var$flipped, ordered=T)

# means
means <- active_2 %>% group_by(ct, type) %>% summarise(mean(t))
active_2$mean_t <- with(active_2, ifelse(ct == "Gm12878", 0.4366260,
                                         ifelse(ct == "H1 hESC", 0.5137763, 0.4183666)))

#pdf("figures/suppl/flipped_long_excl2_side.pdf", 4, 9)
ggplot(active_2, aes(x=t - mean_t, col=ct, y=flipped)) + 
  geom_point(size=I(0)) + 
  geom_vline(xintercept=0, col="grey70") +
  scale_colour_manual(values=c("#0000ff", "#FFA500", "#ff0000")) +
  facet_grid(flip~., scales="free", space="free") +
  ylab("All individual megabase regions of variable structure") + xlab("") +
  labs(colour="Cell type") + geom_point(size=I(1.2)) +
  theme_bw() +
  theme(legend.position="none", 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major=element_blank()) +
  scale_x_continuous(limits=c(-.22, .22), 
                     breaks=seq(-.2, .2, by=.1),
                     labels=c("-0.2", "-0.1", "0", "+0.1", "+0.2"))
#dev.off()

pdf("figures/suppl/flipped_long_excl2_horiz.pdf", 7, 3.1)
ggplot(active_2, aes(y=t - mean_t, col=ct, x=ct, fill=ct, group=ct)) + 
  geom_vline(yintercept=0, col="black") +
  scale_colour_manual(values=c("#0000ff", "#FFA500", "#ff0000")) +
  scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) +
  facet_grid(~flip, scales="free", space="free") +
  geom_hline(xintercept=0, col="grey80") +
  ggtitle("Cell type in which variable region is active") +
  ylab("Proportion of long-range interactions (>2 Mb)\nwith A compartments, relative to expected") +
  labs(colour="Cell type") + xlab("") +
  geom_boxplot(notch=T, outlier.colour="grey70") +
  theme_minimal() + #coord_flip() +
  theme(legend.position="none",
        axis.title.y=element_text(size=10),
        axis.ticks.x=element_line(size=.2),
        axis.title.x=element_text(size=0),
        axis.text.x=element_text(size=8),
        title=element_text(size=8)) +
  stat_summary(fun.y=median, geom="point", col="white") +
  scale_y_continuous(limits=c(-.22, .22), 
                     breaks=seq(-.2, .2, by=.1),
                     labels=c("-0.2", "-0.1", "0", "+0.1", "+0.2"))
dev.off()




# dashed lines should really be for long/short range, could
# be done with diagonals of interaction matrix:
#
# Near cis::
# x x x . . .
#   x x x . .
#     x x x .
#       x x x
#
# Long range::
# . . . x x x
#   . . . x x
#     . . . x
#       . . .

# contact proportions::
gf <- flippedCt_longRange("g", means=F)
hf <- flippedCt_longRange("h", means=F)
kf <- flippedCt_longRange("k", means=F)

gflr$flip <- "Gm12878"
hflr$flip <- "H1 hESC"
kflr$flip <- "K562"

v5c <- function(ff){
  # Raw numbers of near / far contacts:
  ff = kf
  ff %>% group_by(flipped, type) %>% summarise(m=sum(contacts), mean=mean(contacts))  
  
  ff$chr <- factor(gsub("\\..*", "", ff$bins), levels=paste0("chr", c(1:23)))
  ff$pos <- as.numeric(gsub(".*?\\.", "", ff$bins))
  
  for( j in 1:length(levels(ff$flipped)) ){
    #j = 7
    plot <- levels(ff$flipped)[j]
    hline_data <- data.frame(chr=gsub("\\..*", "", plot), 
                             x=as.numeric(gsub(".*?\\.", "", plot)))
    message(paste(j, plot))
    #pdf("figures/suppl/virt5c_flipped6_gm.pdf", 17, 4)
    print(ggplot(subset(ff, flipped==plot & contacts !=0), #& type == "Near-cis (<5 Mb)"), 
                 aes(y=(log10(contacts+1) / log10(sum(contacts))), x=pos/1e6)) + 
            geom_point(aes(col=states), alpha=I(1), size=I(2)) + 
            #scale_y_log10() + 
            scale_x_continuous(breaks=c()) +
            facet_grid(ct~chr, scales="free", space="free") +
            labs(y="Normalised interaction frequency") +  theme_bw() +
            theme(legend.position="none") + xlab("") +
            geom_vline(data=hline_data, aes(xintercept=x/1e6)))
    #dev.off()
    Sys.sleep(15)
  }
}

v5c(gf)
# g 10

levels(kf$flipped)[9]
h.dat[rownames(h.dat) == sub("\\.", "-", levels(kf$flipped)[7]),]
g.dat[rownames(g.dat) == sub("\\.", "-", levels(kf$flipped)[7]),]
k.dat[rownames(k.dat) == sub("\\.", "-", levels(kf$flipped)[7]),]


head(kf)
kf[kf$flipped == levels(kf$flipped)[9] & kf$bins == levels(kf$flipped)[9],]

kf[kf$contacts == 0 & kf$flipped == levels(kf$flipped)[9] & kf$chr == "chr6",]
