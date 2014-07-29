######### Cell-type specific / shared enrichments ###########
# Look at chromHMM / SegWay combined annotations in flipped #
# open, flipped closed compartments and test for enrichment #
# or depletion. For some states compare regions that are    #
# shared between cell types and those that are specific to  #
# one (cell type-specfic vs. shared; e.g. enhancers).       #
#############################################################
# notes = {
#  runtime:       5-10 mins,
#  external deps: [bedtools, chromhmm+segway files] 
# };
library("data.table")
library("dplyr")
library("GenomicRanges")
library("ggplot2")
require("gridExtra")
source("~/hvl/R/hvl.R")

## 1) Calculate intersections between flipped regions
##    and chromHMM / SegWay chromatin states
eigs <- data.frame(g=g.dat$eigen, h=h.dat$eigen, k=k.dat$eigen)
rownames(eigs) <- rownames(g.dat)

eigs <- transform(eigs, g.state = callStates(g)[,2],
                  h.state = callStates(h)[,2],
                  k.state = callStates(k)[,2])

eigs$sum <- rowSums(eigs[,4:6])
eigs$flip <- NA

eigs$flip <- with(eigs, {
  # if all 2 or all 1
  ifelse(sum == 6 | sum == 3, "none",
         ifelse(h.state == 1 & sum == 5, "h.closed",
                ifelse(g.state == 1 & sum == 5, "g.closed",
                       ifelse(k.state == 1 & sum ==5, "k.closed",
                              ifelse(h.state == 2 & sum == 4, "h.open",
                                     ifelse(g.state == 2 & sum == 4, "g.open",
                                            ifelse(k.state == 2 & sum == 4, "k.open", "NA")))))))
})

# Optional: look at distribution of open, closed:
# pie(table(eigs$flip), col=brewer.pal(7, "Reds"))

writeBed <- function(flip=c("h.open", "h.closed",
                            "g.open", "g.closed",
                            "k.open", "k.closed",
                            "none")){
  options(scipen=99)
  rns <- rownames(eigs[eigs$flip == flip,])
  df <- data.frame(chr   = gsub("-.*", "", rns),
                   start = gsub(".*-", "", rns),
                   end   = as.numeric(gsub(".*-", "", rns)) + 1e6,
                   id    = paste0(flip, ".", seq_along(rns)))
  write.table(df, file=paste0("~/hvl/ice/bedfiles/", flip, ".bed"),
              quote=F, row.names=F, col.names=F, sep="\t")
}

### As-is, none.bed:
#   h  g  k  
#   1  2  1  <- g.open
#   2  2  2  <- none 
#   1  2  2  <- h.closed
#   2  1  2

# write all types of flip (+none) to bed files for intersection
writeBed <- Vectorize(writeBed)
writeBed(unique(eigs$flip))

fnames <- c(paste0(c("h", "g", "k"), rep(c(".open.bed", ".closed.bed"), 3)), "none.bed")

for(f in fnames){
  state.file <- ifelse(substr(f, 1, 1) == "h",
                       "chromhmm.segway.h1hesc.comb11.concord4.bed.bed",
                       ifelse(substr(f, 1, 1) == "g", 
                              "chromhmm.segway.gm12878.comb11.concord4.bed.bed",
                              "chromhmm.segway.k562.comb11.concord4.bed.bed"))
  cmd <- paste0("bedtools intersect -a ~/hvl/ice/bedfiles/", f, " -b ~/hvl/ice/bedfiles/", 
                state.file, " -wao > ~/hvl/ice/", f, "_isect.out")
  cat("Running: ", cmd, "\n")
  res <- system(cmd, intern=F)
}

### END 1 ###

## 2) Aggregate and analyse intersection results (no cell type specific / shared distinction)
ctFeats <- function(ct=c("h", "g", "k")){
  hoi <- read.table(paste0("~/hvl/ice/", ct, ".open.bed_isect.out"), stringsAsFactors=F)
  hci <- read.table(paste0("~/hvl/ice/", ct, ".closed.bed_isect.out"), stringsAsFactors=F)
  none <- read.table("~/hvl/ice/none.bed_isect.out", stringsAsFactors=F)
  cnames <- c("chr", "start", "end", "id",
              "chr.b", "start.b", "end.b", "feature",
              "1k", "dot", "start.o", "end.o", "col", "overlap")
  colnames(hoi) <- cnames
  colnames(hci) <- cnames
  colnames(none) <- cnames
  hoi$flip <- "open"
  hci$flip <- "closed"
  none$flip <- "none"
  hi <- rbind(hoi, hci, none)
  
  perFeat <- group_by(hi, flip, id, feature) %.% summarise(olap = sum(overlap))
  perFeat <- subset(perFeat, feature != ".")
  perFeat$ct <- ct
  return(perFeat)
}

hf1 <- ctFeats("h")
gf1 <- ctFeats("g")
kf1 <- ctFeats("k")
f1 <- rbind(hf1, gf1, kf1)

ggplot(f1, aes(x=flip, y=olap/1e3, fill=flip)) +
  facet_wrap(ct~feature, scales="free_y", nrow=3) +
  geom_boxplot(notch=T, position=position_dodge(.8),
               outlier.size=.5) +
  labs(list(y="Annotation coverage per Mb (kb)", x="",
            fill="Flipped state")) +
  scale_fill_brewer(palette="Blues") + theme_bw()


## 3) Build subset of annotations that are preserved in all cell types (?)

class.vec <- c("character", "numeric", "numeric", "character", "numeric", "character",
               "numeric", "numeric", "character")
chrom.gm <- read.table("~/hvl/ice/bedfiles/chromhmm.segway.gm12878.comb11.concord4.bed.bed",
                       stringsAsFactors=F, colClasses=class.vec)
chrom.h1 <- read.table("~/hvl/ice/bedfiles/chromhmm.segway.h1hesc.comb11.concord4.bed.bed",
                       stringsAsFactor=F, colClasses=class.vec)
chrom.k5 <- read.table("~/hvl/ice/bedfiles/chromhmm.segway.k562.comb11.concord4.bed.bed",
                       stringsAsFactors=F, colClasses=class.vec)

c.gm <- with(chrom.gm, GRanges(V1, IRanges(start=V2, end=V3), "+", feat=V4))
c.h1 <- with(chrom.h1, GRanges(V1, IRanges(start=V2, end=V3), "+", feat=V4))
c.k5 <- with(chrom.k5, GRanges(V1, IRanges(start=V2, end=V3), "+", feat=V4))

featSubset <- function(f=unique(c.k5$feat), c1=c("g", "h", "k")){
  #   f = "E"
  #   c1 = "h"
  f1 <- if(c1 == "g") c.gm else if(c1 == "h") c.h1 else c.k5
  f2 <- if(c1 == "g") c.h1 else if(c1 == "h") c.k5 else c.gm
  f3 <- if(c1 == "g") c.k5 else if(c1 == "h") c.gm else c.h1
  f1 <- subset(f1, feat == f)
  f2 <- subset(f2, feat == f)
  f3 <- subset(f3, feat == f)
  
  ## Find overlapping regions:
  f4 <- subsetByOverlaps(f1, f2)
  f4.2 <- subsetByOverlaps(f1, f3)
  f5 <- subsetByOverlaps(f4, f3)
  cat("%", f, "shared in all cell types", 100* length(f5) / length(f1), "\n")
  
  # unique ID for each element
  id.assign <- function(df) return(as.vector(with(df, paste(seqnames, start, feat, sep="."))))
  
  # Get elements in all cell types
  shared <- as.data.frame(f5)
  shared$id <- id.assign(shared)
  
  # Get those in any two (single shared)
  single <- rbind(as.data.frame(f4), as.data.frame(f4.2))
  single$id <- id.assign(single)
  
  #cat(length(shared$id))
  orig <- as.data.frame(f1)
  orig$id <- id.assign(orig)
  orig$status <- ifelse(!orig$id %in% single$id, "in:one", 
                        ifelse(orig$id %in% shared$id, "in:all", "in:some"))
  
  return(orig)
}

gm.feats <- sapply(unique(c.k5$feat), featSubset, c1="g", simplify=F)
gm.feats <- do.call(rbind, gm.feats)

par(mfrow=c(4,4), mar=c(2,4,3,1.5))
options(scipen=99)
pie.dens <- function(feat){
  pie(table(gm.feats[gm.feats$feat == feat,]$status), main=feat)
  plot(ecdf(gm.feats[gm.feats$feat == feat,]$width/1000), main=feat,
       verticals=T)
}  
sapply(unique(gm.feats$feat), pie.dens)

dev.off()
plot(density(gm.feats[gm.feats$feat == "R",]$width), log="x")

h1.feats <- sapply(unique(c.k5$feat), featSubset, c1="h", simplify=F)
h1.feats <- do.call(rbind, h1.feats)

k5.feats <- sapply(unique(c.k5$feat), featSubset, c1="k", simplify=F)
k5.feats <- do.call(rbind, k5.feats)

## then write each of these three types to seperate files for 
## bedtools intersect (i.e. step 1)

splitWrite <- function(feats, ct=c("g", "h", "k")){
  options(scipen=99)
  write.table(subset(feats, status == "in:all")[,c(1:3, 6)],
              file = paste0("~/hvl/ice/bedfiles/", ct, "_shared_chromstates.bed"),
              row.names=F, col.names=F, quote=F, sep="\t")
  write.table(subset(feats, status == "in:one")[,c(1:3, 6)],
              file = paste0("~/hvl/ice/bedfiles/", ct, "_celltypespecific_chromstates.bed"),
              row.names=F, col.names=F, quote=F, sep="\t")
  write.table(subset(feats, status == "in:some")[,c(1:3, 6)],
              file = paste0("~/hvl/ice/bedfiles/", ct, "_partshared_chromstates.bed"),
              row.names=F, col.names=F, quote=F, sep="\t")
  write.table(subset(feats, status %in% c("in:all", "in:some"))[,c(1:3, 6)],
              file = paste0("~/hvl/ice/bedfiles/", ct, "_shared2plus_chromstates.bed"),
              row.names=F, col.names=F, quote=F, sep="\t")
}

splitWrite(gm.feats, ct="g")
splitWrite(h1.feats, ct="h")
splitWrite(k5.feats, ct="k")

for(f in fnames[-length(fnames)]){
  ct <- substr(f, 1, 1)
  flip <- ifelse(grepl("closed", f), "closed", 
                 ifelse(grepl("open", f), "open", "none"))
  cmd1 <- paste0("bedtools intersect -a ~/hvl/ice/bedfiles/", f, " -b ~/hvl/ice/bedfiles/", 
                 ct, "_celltypespecific_chromstates.bed -wao > ~/hvl/ice/bedfiles/", 
                 f, "_cts.out")
  ## This command treats "shared" as present in ALL THREE only:
  #   cmd2 <- paste0("bedtools intersect -a ~/hvl/ice/bedfiles/", f, " -b ~/hvl/ice/bedfiles/", 
  #                  ct, "_shared_chromstates.bed -wao > ~/hvl/ice/bedfiles/", 
  #                  f, "_shared.out")
  ## This one considers "shared" as anything not cell type specific:
  cmd2 <- paste0("bedtools intersect -a ~/hvl/ice/bedfiles/", f, " -b ~/hvl/ice/bedfiles/", 
                 ct, "_shared2plus_chromstates.bed -wao > ~/hvl/ice/bedfiles/", 
                 f, "_shared.out")
  cat("Running: ", cmd1, "\n\n")
  system(cmd1, intern=F)
  cat("Running: ", cmd2, "\n\n")
  system(cmd2, intern=F)
}

## process the blocks with no flip for comparison:
for(ct in c("h", "g", "k")){
  ## Same as above, rm "2plus" for shared == ALL 3 cell type overlap
  cmd3 <- paste0("bedtools intersect -a ~/hvl/ice/bedfiles/none.bed  -b ~/hvl/ice/bedfiles/", 
                 ct, "_shared2plus_chromstates.bed -wao > ~/hvl/ice/bedfiles/", 
                 ct, ".none_shared.out")
  cmd4 <- paste0("bedtools intersect -a ~/hvl/ice/bedfiles/none.bed -b ~/hvl/ice/bedfiles/", 
                 ct, "_celltypespecific_chromstates.bed -wao > ~/hvl/ice/bedfiles/", 
                 ct, ".none_cts.out")
  cat("Running: ", cmd3, "\n\n")
  system(cmd3, intern=F)
  cat("Running: ", cmd4, "\n\n")
  system(cmd4, intern=F)
}

## 4) Read results and analyse, as per #2
ctFeats.p2 <- function(ct=c("h", "g", "k"), overlap=TRUE){
  ## Debugging:
  #   ct = "h"
  #   overlap = FALSE
  
  cnames <- c("chr", "start", "end", "id",
              "chr.b", "start.b", "end.b", "feature", "overlap")
  
  read.feat <- function(file){
    f <- read.table(paste0("~/hvl/ice/bedfiles/", file), stringsAsFactors=F, col.names=cnames)
    f$flip <- ifelse(grepl("open", file), "open",
                     ifelse(grepl("closed", file), "closed", "none"))
    f$type <- ifelse(grepl("shared", file), "shared",
                     ifelse(grepl("cts", file), "cts", "ERROR"))
    return(f)
  }
  
  oc <- read.feat(paste0(ct, ".open.bed_cts.out"))
  cc <- read.feat(paste0(ct, ".closed.bed_cts.out"))
  os <- read.feat(paste0(ct, ".open.bed_shared.out"))
  cs <- read.feat(paste0(ct, ".closed.bed_shared.out"))
  ns <- read.feat(paste0(ct, ".none_shared.out"))
  nc <- read.feat(paste0(ct, ".none_cts.out"))
  i <- rbind(oc, cc, os, cs, ns, nc)
  
  perFeat <- if(overlap == T){
    ## sum base pairs
    group_by(i, type, flip, id, feature) %.% dplyr::summarise(olap = sum(overlap))
  } else {
    ## count number
    group_by(i, type, flip, id, feature) %.% dplyr::summarise(olap = n())
  }
  
  ## nb do this after summarising
  perFeat <- subset(perFeat, feature != ".")
  perFeat$feature <- factor(perFeat$feature)
  perFeat$id <- factor(perFeat$id)
  all.combs <- expand.grid(levels(perFeat$id), unique(perFeat$type), levels(perFeat$feature))
  
  pfdt <- data.table(perFeat)  
  setkey(pfdt, id)
  all.combs <- data.table(all.combs) 
  # nb id is tied to flip
  setnames(all.combs, 1:3, c("id", "type", "feature"))
  setkey(all.combs, id)
  all.combs$flip <- ifelse(grepl("open", all.combs$id), "open",
                           ifelse(grepl("closed", all.combs$id), "closed", "none"))
  m <- merge(all.combs, pfdt, by=c("id", "type", "flip", "feature"), all.x=T, allow.cartesian=T)
  m$olap[is.na(m$olap)] <- 0
  
  perFeat <- as.data.frame(m)
  # "." == NULL feature, i.e. no overlaps, instead add a zero for each feature?
  perFeat$ct <- ct
  return(perFeat)
}

pF.g <- ctFeats.p2("g", overlap=F)
pF.k <- ctFeats.p2("k", overlap=F)
pF.h <- ctFeats.p2("h", overlap=F)

pF <- rbind(pF.g, pF.h, pF.k)
pF$ct <- ifelse(pF$ct == "h", "H1 hESC", 
                ifelse(pF$ct == "g", "GM12878", "K562"))
stopifnot(length(unique(pF$type)) == 2)
#pF$type
pF$type <- ifelse(pF$type == "cts", "Cell type specific", "Shared")

options(scipen=99)
pd <- position_dodge(width=.9)
rt <- subset(pF, feature %in% c("R", "T") & type == "Cell type specific")
rt.counts <- group_by(rt, ct, flip, feature) %.% summarise(top=max(olap)*1.05, count=paste0("n = ", n()))

pf.i <- subset(pF, feature %in% c("R", "T") & type == "Cell type specific")
group_by(pf.i, feature, flip, feature) %.% summarise(p=wilcox.test(olap)$p.value)

f4.plot <- function(chrom.state){
  epf <- subset(pF, feature == chrom.state & type == "Cell type specific")
  epf.shared <- subset(pF, feature == chrom.state & type == "Shared")
  e.counts <- group_by(epf, ct, flip, feature) %.% summarise(top=max(olap)*.7, count=paste0("n = ", n()))
  
  grid.arrange(
    ggplot(epf, aes(x=flip, y=olap, fill=flip,
                    group=interaction(flip,type,feature))) +
      scale_fill_brewer(palette="Blues") + theme_bw() +
      coord_cartesian(ylim = quantile(epf$olap, c(0, 0.95))) +
      geom_boxplot(width=.8, notch=T, outlier.size=0, position=pd,
                   colour="black") +
      stat_summary(fun.y=mean, geom="point", position=pd,
                   colour="grey40", shape=18, size=3) +
      facet_grid(type~ct,shrink=T) +
      labs(list(y="Number of annotated features per Mb", x="",
                fill="Flipped state")) #+ theme(legend.position="none") 
    ,
    ggplot(epf.shared, aes(x=flip, y=olap, fill=flip,
                           group=interaction(flip,type,feature))) +
      scale_fill_brewer(palette="Blues") + theme_bw() +
      coord_cartesian(ylim = quantile(epf.shared$olap, c(0, 0.95))) +
      geom_boxplot(width=.8, notch=T, outlier.size=0, position=pd, colour="black") +
      stat_summary(fun.y=mean, geom="point", position=pd,
                   colour="grey40", shape=18, size=3) +
      facet_grid(type~ct, shrink=T) +
      labs(list(y="Number of annotated features per Mb", x="",
                fill="Flipped state")) #+ theme(legend.position="none")
    , ncol=1)
}

# Figure 4b, flipped open compartments enriched for enhancers
#pdf("~/hvl/ice/plots/f4b_enhancerEnrichFlipped.pdf", 7, 5.5)
f4.plot("E")
#dev.off()

# interesting:
f4.plot("WE")
f4.plot("CTCF")
f4.plot("TSS")

#pdf("~/hvl/ice/plots/supplementary/s5_transcribedFlipped.pdf", 7, 5.5)
f4.plot("T")
#dev.off()

# boring or weird:
f4.plot("R")
f4.plot("PF")

# Open/closed vs. None wilcox tests
cat("Cell type\tFlip\tType\tSignif\tMedian(flipped)\tMedian(none)\n")
for(c in c("GM12878", "H1 hESC", "K562")){
  for(f in c("open", "closed")){
    for(t in c("Cell type specific", "Shared")){
      res <- with(pF, wilcox.test(pF[feature == "E" & flip == f & ct == c & type == t, "olap"],
                                  pF[feature == "E" & flip == "none" & ct == c & type == t, "olap"]))
      cat(c, "\t", f, "\t", t, "\t", res$p.value, "\t", 
          with(pF, median(pF[feature == "E" & flip == f & ct == c & type == t, "olap"])),
          with(pF, median(pF[feature == "E" & flip == "none" & ct == c & type == t, "olap"])),
          "\n")
    }
  }
}

## Open vs. Closed wilcox tests
cat("Cell type\tFlip\tType\tSignif\tMedian(flipped)\tMedian(none)\n")
for(c in c("GM12878", "H1 hESC", "K562")){
  for(t in c("Cell type specific", "Shared")){
    res <- with(pF, wilcox.test(pF[feature == "E" & flip == "open" & ct == c & type == t, "olap"],
                                pF[feature == "E" & flip == "closed" & ct == c & type == t, "olap"]))
    cat(c, "\t", f, "\t", t, "\t", res$p.value, "\t", 
        with(pF, median(pF[feature == "E" & flip == f & ct == c & type == t, "olap"])),
        with(pF, median(pF[feature == "E" & flip == "none" & ct == c & type == t, "olap"])),
        "\n")
  }
}

## Supplementary figure: all tests:
#pdf("~/hvl/ice/plots/s11_allBeans.pdf", 9, 12)
ggplot(pF, aes(x=flip, y=olap, fill=flip,
               group=interaction(flip,type,feature,ct), ymax=0)) +
  scale_fill_brewer(palette="Blues") +
  geom_violin(scale="width", position=pd) + 
  geom_boxplot(width=.1, outlier.size=0, position=pd, fill="black") +
  stat_summary(fun.y=median, geom="point", position=pd,
               fill="white", shape=21, size=3) +
  facet_grid(feature~ct+type, scales="free_y") + theme_bw() +
  labs(list(y="Annotation coverage per Mb (kb)", x="",
            fill="Flipped state")) + theme(legend.position="none")
#dev.off()


## Investigate specific regions

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

options(scipen=99)
cat(paste0(bins[,1], ":", bins[,2] - 2e6, "-", bins[,3] + 2e6, "\n"))
# Yes - interesting, use; N - HMM calls, don't use; U - uninteresting but could use;
# YNNUU Y
bins

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

require("dplyr")
group_by(b, ct) %>% mutate(c=1:n())

