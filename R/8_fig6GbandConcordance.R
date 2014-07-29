###### Compare metaphase G-bands with i-phase structure ######
# Compare giemsa bands by overlap of A/B compartments and by #
# boundary distances between TAD and compartment edges.      #
##############################################################
library("dplyr")
library("ggplot2")
library("gridExtra")
library("RHmm")
library("scales")
source("~/hvl/R/hvl.R")
options(scipen=99)
set.seed(42)
gbands <- read.table("~/hvl/ice/giemsa/ucsc_g.dat", header=T, sep="\t")
chr.sizes <- read.table("~/hvl/hg19.chrom.sizes.txt")[1:24,]
dr <- "~/hvl/ice/giemsa/"

## given g-band what's the nearest TAD / compartment
circDiffs<- function(vec, n, max, bounds, toEdge=0){
  # take a vector of positions and add a constant,
  # if the updated value is greater than the chromosome 
  # size, wrap it.
  vec <- vec + n
  vec <- vec %% max
  if(!toEdge == 0){
    diffs1 <- sapply(bounds+toEdge, function(x) min(abs(x - vec)))
    diffs2 <- sapply(bounds-toEdge, function(x) min(abs(x - vec)))
    diffs3 <- sapply(bounds, function(x) min(abs(x - vec)))
    diffs <- apply(cbind(diffs1, diffs2, diffs3), 1, min)
  } else {
    diffs <- sapply(bounds, function(x) min(abs(x - vec)))
  }
  return(diffs)
}

gToBoundDist <- function(cb, reps=10, edge=0){
  ## pre-allocate data frame
  # 931 g-bands, 942 compartment bounds (*100, +100)
  rn <- nrow(gbands)*reps
  #odf <- data.frame(chr=character(rn), type=character(rn), dist=numeric(rn))
  pb <- txtProgressBar(min = 1, max = rn, style = 3)
  i <- 1
  odf <- data.frame(chr=character(), type=character(), dist=numeric())
  for(c in paste0("chr", 1:22)) {
    curr.bands <- gbands[gbands$chrom == c,]
    g.bands <- c(curr.bands$chromStart[1], curr.bands$chromEnd)
    curr.bounds <- rowMeans(cb[cb[,1] == c, 2:3])
    max <- chr.sizes[chr.sizes[,1] == c, 2]
    
    # calculate permuted distances, from supplied 
    perms <- c(replicate(reps, circDiffs(curr.bounds,         # compartment bounds
                                         sample(1:max, 1),    # degree of shift for this chr
                                         max,                 # max chromosome position 
                                         g.bands,             # giemsa bands (end?)
                                         toEdge=edge)))       # amount of boundary wobble allowed
    odf <- rbind(odf, data.frame(chr  = rep(c, length(perms)), 
                                 type = rep("permute", length(perms)), 
                                 dist = perms))
    
    # calculate differences using actual boundaries
    actual <- circDiffs(curr.bounds,  0, max, g.bands, toEdge=edge)
    odf <- rbind(odf, data.frame(chr  = rep(c, length(actual)), 
                                 type = rep("empirical", length(actual)),
                                 dist = actual))
    i <- i + length(perms) 
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(odf)
}

## 1) Compartments:
gcb <- read.table("~/hvl/ice/boundaries/g_cb.bed")
hcb <- read.table("~/hvl/ice/boundaries/h_cb.bed")
kcb <- read.table("~/hvl/ice/boundaries/k_cb.bed")

gcg <- gToBoundDist(gcb, reps=20)
hcg <- gToBoundDist(hcb, reps=20)
kcg <- gToBoundDist(kcb, reps=20)

gcg$ct <- "GM12878"
hcg$ct <- "H1 hESC"
kcg$ct <- "K562"
comps <- rbind(gcg, hcg, kcg)

## 2) TADs
gtb <- read.table("~/hvl/ice/boundaries/gm_tbounds.bed")
htb <- read.table("~/hvl/ice/boundaries/h1_tbounds.bed")
ktb <- read.table("~/hvl/ice/boundaries/k5_tbounds.bed")

gtg <- gToBoundDist(gtb, reps=20)
htg <- gToBoundDist(htb, reps=20)
ktg <- gToBoundDist(ktb, reps=20)

gtg$ct <- "GM12878"
htg$ct <- "H1 hESC"
ktg$ct <- "K562"
tads <- rbind(gtg, htg, ktg)

comps$test <- "Compartments"
tads$test <- "TADs"
bv2 <- rbind(comps, tads)
bv2$type <- ifelse(as.character(bv2$type) == "permute", "Null model", "Empirical")


pdf("~/hvl/ice/plots/supplementary/s10_GbandDistECDF.pdf", 7, 5)
ggplot(bv2, aes(x=dist+1, col=type)) + #, linetype=ct)) +
  facet_grid(ct~test) + stat_ecdf(geom="line") +
  scale_x_log10(breaks=c(1e3, 1e4, 1e5, 1e6, 1e7),
                labels=c("1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb")) +
  coord_cartesian(xlim=c(10000,3e7)) + theme_bw() +
  labs(list(x="Distance from G-band edge to nearest boundary", y="ECDF",
            linetype="Cell type", col="Distribution")) +
  scale_color_brewer(type="qual", palette=2)
dev.off()

## significance:
for(ct in unique(bv2$ct)) {
  res <- ks.test(bv2[bv2$ct == ct & 
                       bv2$type == "Null model" & 
                       bv2$test == "Compartments", "dist"],
                 bv2[bv2$ct == ct & 
                       bv2$type == "Empirical" & 
                       bv2$test == "Compartments", "dist"])
  cat(ct, "\t", res$statistic, "\t", res$p.value, "\n")
}

res <- with(bv2, ks.test(bv2[ type == "Null model" & 
                                test == "Compartments", "dist"],
                         bv2[type == "Empirical" & 
                               test == "Compartments", "dist"]))

res.2 <- with(bv2, ks.test(bv2[ type == "Null model" & 
                                  test == "TADs", "dist"],
                           bv2[type == "Empirical" & 
                                 test == "TADs", "dist"]))

## Compartments:
res$statistic  # 0.097
res$p.value    # 0
## TADs:
res.2$statistic # 0.04
res.2$p.value   # 0.001

options(scipen=1)
for(ct in unique(bv2$ct)) {
  cat(median(bv2[bv2$ct == ct & 
                   bv2$type == "Null model" & 
                   bv2$test == "Compartments", "dist"]), "\t",
      median(bv2[bv2$ct == ct & 
                   bv2$type == "Empirical" & 
                   bv2$test == "Compartments", "dist"]), "\n")
}

for(ct in unique(bv2$ct)) {
  res <- ks.test(bv2[bv2$ct == ct & 
                       bv2$type == "Null model" & 
                       bv2$test == "TADs", "dist"],
                 bv2[bv2$ct == ct & 
                       bv2$type == "Empirical" & 
                       bv2$test == "TADs", "dist"])
  cat(ct, "\t", res$p.value, "\n")
}

test.dist = 500000
with(bv2, nrow(bv2[type == "Empirical" & 
                     dist <= test.dist & 
                     test == "Compartments" ,]) / 
       nrow(bv2[type == "Empirical",])) 
with(bv2, nrow(bv2[type == "Null model" & 
                     dist <= test.dist & 
                     test == "Compartments",]) /
       nrow(bv2[type == "Null model",]))

## How many TADs are called per CT and chromosome?
ktads <- group_by(ktb, V1) %.% summarise(n()) %.% arrange(as.numeric(gsub("chr","",V1)))
gtads <- group_by(gtb, V1) %.% summarise(n()) %.% arrange(as.numeric(gsub("chr","",V1)))
htads <- group_by(htb, V1) %.% summarise(n()) %.% arrange(as.numeric(gsub("chr","",V1)))

## Note lots more TADs called in H1
barplot(htads[,2], add=F, names.arg=htads[,1], las=3, ylab="Number of TADs called",
        col=rgb(t(col2rgb("orange")), max=255, alpha=150), ylim=c(0,300))
barplot(gtads[,2], add=T, axes=F, 
        col=rgb(t(col2rgb("blue")), max=255, alpha=150))
barplot(ktads[,2], add=T, axes=F, 
        col=rgb(t(col2rgb("red")), max=255, alpha=150))
legend("topright", legend=c("H1 hESC", "GM12878", "K562"),
       fill=c(rgb(t(col2rgb("orange")), max=255, alpha=150),
              rgb(t(col2rgb("blue")), max=255, alpha=150),
              rgb(t(col2rgb("red")), max=255, alpha=150)), bty="n")

######## Segment overlap between compartments and G-bands

# 1 - calculate A/B co-ordinates per chr
# states from (e.g.) hstates, bins from h.dat

getComps <- function(chromosome, states=kstates){
  bin.chr <- gsub("-.*", "", rownames(h.dat))
  cur.chr <- states[which(bin.chr == chromosome),]
  bin.bp <- as.numeric(gsub(".*-", "", rownames(h.dat[which(bin.chr == chromosome),])))
  # nb. need to add 
  ends.ind <- cumsum(rle(cur.chr[,2])$lengths)
  ends <- bin.bp[ends.ind]#+1e6
  starts <- c(bin.bp[1], ends[1:(length(ends)-1)])
  # outfmt:
  # chr   start   end     compartment
  # chr1  100000  5000000 A
  df.rows <- data.frame(chr=chromosome, start=starts, end=ends, 
                        compartment=ifelse(rle(cur.chr[,2])$values == 1, "B", "A"))
}

cwrapper <- function(s, fname){
  options(scipen=99)
  dr <- "~/hvl/ice/giemsa/"
  f <- do.call(rbind, sapply(paste0("chr", c(1:22, "X")), 
                             getComps, s, simplify=F))
  write <- transform(f, id=paste0(rownames(f), ".", f$compartment))
  write.table(write[,c(1:3,5)], paste0(dr, fname), quote=F,
              row.names=F, col.names=F, sep="\t")
  
}

## Call HMM states: hvl.R::callStates(dat$eigen)[,2]
h.comps <- cwrapper(callStates(h.dat$eigen), "hcomps.bed")
k.comps <- cwrapper(callStates(k.dat$eigen), "kcomps.bed")
g.comps <- cwrapper(callStates(g.dat$eigen), "gcomps.bed")

# Giemsa band BED:
write.table(gbands[,c(1:3,5)], paste0(dr, "gbands.bed"), quote=F,
            row.names=F, col.names=F, sep="\t")

getIntersections <- function(comps.fn){
  cmdres <- system(paste0("bedtools intersect -b ", dr, "gbands.bed -a ", dr, comps.fn, " -wo"), wait=T, intern=T)
  isect <- as.data.frame(do.call(rbind, strsplit(cmdres, "\t")))
  colnames(isect) <- c("chr", "start", "end", "id", "chr.B", "start.B", "end.B", "gband", "overlap")
  isect$overlap <- as.numeric(as.character(isect$overlap))
  isect$start <- as.numeric(as.character(isect$start))
  isect$end <- as.numeric(as.character(isect$end))
  isect <- transform(isect, comp=gsub(".*\\.", "", id))
  return(isect)
}

h.isect <- getIntersections("hcomps.bed")
g.isect <- getIntersections("gcomps.bed")
k.isect <- getIntersections("kcomps.bed")

h.isect <- transform(h.isect, ct="H1 hESC")
g.isect <- transform(g.isect, ct="GM12878")
k.isect <- transform(k.isect, ct="K562")

all.isect <- rbind(h.isect, g.isect, k.isect)

pergband <- group_by(h.isect, id, gband) %.% summarise(total = sum(overlap), .drop=F)
percomp <- group_by(all.isect, ct, comp, gband) %.% summarise(total = sum(overlap))
percomp <- subset(percomp, gband != "acen" & gband != "gvar")
percomp$gband <- factor(percomp$gband, levels=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100"))

## Build null:
set.seed(42)
emitStates <- function(indat){
  hmm <- HMMFit(indat, dis="NORMAL", nStates=2, control=list(init="RANDOM", iter=5000))
  out <- data.frame(eigen=indat, state=HMMSim(length(indat), hmm$HMM)$states)
  return(out)
}

overlapPermute <- function(dat, reps=100){
  # dats should be thought of as global, immutable
  ct <- if(identical(dat, h.dat)) "h" else if(identical(dat, k.dat)) "k" else "g"
  rdf <- data.frame(chr=character(), start=integer(),
                    end=integer(), id=character(), chr.B=character(),
                    start.B=character(), gband=character(), overlap=integer(),
                    comp=character(), ct=character(), run=integer())
  
  options(scipen=99)
  for(i in 1:reps){
    rand <- emitStates(dat$eigen)
    r.comps <- cwrapper(rand, paste0(ct, "randcomps.bed"))
    r.isect <- getIntersections(paste0(ct, "randcomps.bed"))
    r.isect <- transform(r.isect, ct=paste0(ct, ".null"), run=i)
    rdf <- rbind(rdf, r.isect)
  }
  return(rdf)
}

h.rdf <- overlapPermute(h.dat)
g.rdf <- overlapPermute(g.dat)
k.rdf <- overlapPermute(k.dat)

rdf <- rbind(h.rdf, g.rdf, k.rdf)

pcr <- group_by(rdf, ct, comp, gband, run) %.% summarise(total = sum(overlap))
pcr <- subset(pcr, gband != "acen" & gband != "gvar")
pcr$gband <- factor(pcr$gband, levels=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100"))

pos <- position_dodge(.9)

pdf("~/hvl/ice/plots/supplementary/s7_gBandIntersectionVnull.pdf", 10, 8)
grid.arrange(
  ggplot(percomp, aes(x=gband, y=total/1e6, fill=comp)) +
    facet_grid(ct~.) + theme_bw() + 
    theme(legend.position=c(.85,.92), legend.background=element_blank()) +
    geom_bar(stat="identity", position=position_dodge()) + 
    labs(list(y="Total intersection across genome (Mb)", 
              x="Giemsa stain", fill="Compartment")) +
    scale_fill_brewer(palette=3, type="qual") +
    ggtitle("Observed overlaps"),
  ggplot(pcr, aes(x=gband, y=total/1e6, fill=comp, group=comp)) +
    facet_grid(ct~.) +  theme_bw() + 
    theme(legend.position=c(.85,.92), legend.background=element_blank()) +
    #geom_bar(aes(x=gband, y=total/1e6), stat="identity", position=pos) +
    stat_summary(fun.y=mean, geom="bar", position=pos, alpha=I(.75)) + 
    stat_summary(fun.data=mean_cl_normal, geom="errorbar", 
                 col="black", position=pos, width=.3) +
    labs(list(y="Total intersection across genome (Mb)", 
              x="Giemsa stain", fill="Compartment")) +
    scale_fill_brewer(palette=3, type="qual") +
    ggtitle("Expected overlaps"),
  ncol=2
)
dev.off()


pc.c <- percomp
pcr$ct <- plyr::revalue(pcr$ct, c("h.null"="H1 hESC", "k.null"="K562", "g.null"="GM12878"))

pos <- position_dodge(.4)

## Better dataviz IMO, place expected behind (alpha'd) the observed
ggplot(pcr, aes(x=gband, y=total/1e6, fill=comp, group=comp)) +
  facet_grid(ct~.) +  theme_bw() +
  stat_summary(fun.y=mean, geom="bar", position="stack", 
               aes(width=.8), alpha=I(.5))+ 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", 
               col="black", position=pos, width=.15) +
  labs(list(y="Total intersection across genome (Mb)", 
            x="Giemsa stain", fill="Compartment")) +
  scale_fill_manual(values=brewer.pal(4, "Paired")[1:2]) +
  geom_bar(aes(x=gband, y=total/1e6, group=comp), stat="identity",
           position="stack", data=pc.c, width=.4, alpha=I(1)) +
  coord_flip()

# sum of A in gneg and gpos 25:
totals <- group_by(pcr, comp, gband) %.% summarise(olap=mean(total) / 1e6)
sum(with(totals, totals[comp == "A" & gband %in% c("gneg", "gpos25"),"olap"]))
665.719 / 2811 # 23.7 % of (genome) in A + gneg/gpos25

sum(with(totals, totals[comp == "B" & gband %in% c("gpos75", "gpos100"),"olap"]))
492.3537 / 2811 # 17.5% of (genome) in B + gpos75/100

100 * (665.719 + 492.3537) / 2811

cell.types <- unique(as.character(pcr$ct))
bands <- unique(as.character(pcr$gband))
cs <- c("A", "B")

options(scipen=0)
for(cell in cell.types){
  for(c in cs){
    for(b in bands){
      #print(subset(percomp, ct == cell & comp == c &
      #               gband == b)$total)
      r <- wilcox.test(
        subset(pcr, ct == cell & comp == c &
                 gband == b)$total,
        mu=subset(percomp, ct == cell & comp == c &
                    gband == b)$total)
      cat(paste(cell, c, b, r$p.value, collapse="\t"), "\n")
    }
  }
}

