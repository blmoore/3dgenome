## split ctcf and yy1, i.e. is ctcf * yy1 more
## significantly enriched at boundaries than
## ctcf or yy1
library("ggplot2")
library("dplyr")

all.t <- readRDS("data/rds/tad_boundary_features.rds")
# avoid melt
mall <- melt(all.t, id.vars=c("name", "bound", "feat", "ct", "type"))

mall[mall$feat == "Atf3",]$value <- subset(mall, feat == "Ctcf")$value * subset(mall, feat == "Yy1")$value

all.t[all.t$feat=="Atf3", c(4:28)] <- all.t[all.t$feat=="Ctcf", c(4:28)] * all.t[all.t$feat=="Yy1", c(4:28)]

comp.t <- group_by(all.t, ct, type, feat) %>% 
  summarise(bound.mean = mean(X12),
            edge.mean  = mean(c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21)),
            p = wilcox.test(X12, c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21))$p.value)


ggplot(comp.t, aes(x=feat, y=-log10(p), col=type,
                 size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~., scales="free_y") + 
  geom_point(position=position_dodge(.15)) + theme_bw() +
  labs(list(size="Absolute difference\nat boundary",
            y=expression(-log[10](italic(p))),
            x="", col="Boundary type")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(range=c(2,7)) +
  guides(col=guide_legend(override.aes=list(size=4))) +
  geom_hline(yintercept=-log10(.01 / 297), linetype="dashed") +
  scale_colour_brewer(type="qual", palette=6) 


read_peaks <- function(fname){
 # fname = "data/peaks//wgEncodeBroadHistoneGm12878CtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bam.regionPeak"
  tmp <- read.table(fname)
  with(tmp, GRanges(seqnames=V1, ranges=IRanges(start=V2, end=V3)))
}
  
write_overlaps <- function(ct=1:3){
  library("GenomicRanges")
  library("rtracklayer")
  # ct 1:: Gm12878; 2:: H1 hESC; 3:: K562
  cname <- ifelse(ct == 1, "Gm12878", ifelse(ct == 2, "H1hesc", "K562"))
  short <- tolower(substr(cname, 1, 2))
  
  ## read peak calls, intersect, then boundary profile each
  infiles <- data.frame(
    ctcf=c("wgEncodeBroadHistoneGm12878CtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bam.regionPeak",
           "wgEncodeBroadHistoneH1hescCtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneH1hescControlStdAlnRep0.bam.regionPeak",
           "wgEncodeBroadHistoneK562CtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneK562ControlStdAlnRep1.bam.regionPeak"),
    yy1=c("wgEncodeHaibTfbsGm12891Yy1sc281V0416101AlnRep0.bam_VS_wgEncodeHaibTfbsGm12891RxlchPcr1xAlnRep0.bam.regionPeak",
          "wgEncodeHaibTfbsH1hescYy1sc281V0416102AlnRep0.bam_VS_wgEncodeHaibTfbsH1hescRxlchPcr1xAlnRep0.bam.regionPeak",
          "wgEncodeHaibTfbsK562Yy1V0416102AlnRep0.bam_VS_wgEncodeHaibTfbsK562RxlchPcr1xAlnRep0.bam.regionPeak"))
  
  peak_dir <- "data/peaks/"
  ctcf <- read_peaks(paste0(peak_dir, as.character(infiles$ctcf[[ct]])))
  yy1 <- read_peaks(paste0(peak_dir, as.character(infiles$yy1[[ct]])))
  
  # intersections: ctcf v tads, yy1 v tads && ctcf + yy1 v tads
  both <- subsetByOverlaps(ctcf, yy1)
  
  # now intersect with TAD boundaries:
  tads <- read_peaks(paste0("data/bedfiles/", short, "_tadbounds.bed"))
  yy1_tads <- subsetByOverlaps(tads, yy1)
  message(paste0(signif(100*length(yy1_tads) / length(tads), digits=3), 
                 "% of tads have yy1"))
  
  ctcf_tads <- subsetByOverlaps(tads, ctcf)
  # ~61% of TAD boundaries (40kb) have signif CTCF peak::
  message(paste0(signif(100 * length(ctcf_tads) / length(tads), digits=3),
                 "% of tads have ctcf"))
  
  
  marked_tads <- subsetByOverlaps(tads, both)
  # ~12% of TAD boundaries contain both::
  message(paste0(signif(100*length(marked_tads) / length(tads), digits=3),
                 "% of tads have ctcf + yy1"))
  
  which <- ifelse(tads %in% marked_tads, "both", 
                  ifelse(tads %in% yy1_tads, "yy1", 
                         ifelse(tads %in% ctcf_tads, "ctcf", "none")))
  
  message("Writing bed files...") 
  export.bed(tads[which == "ctcf",], paste0("data/bedfiles/", short, "_ctcf_tads.bed"))
  export.bed(tads[which == "both",], paste0("data/bedfiles/", short, "_both_tads.bed"))
  export.bed(tads[which == "yy1",], paste0("data/bedfiles/", short, "_yy1_tads.bed"))
  export.bed(tads[which == "none",], paste0("data/bedfiles/", short, "_none_tads.bed"))
  
  message(length(tads))

  return(as.factor(which))
}
  
# gm12878
gt_status <- write_overlaps(1)

# h1
ht_status <- write_overlaps(2)

# k5
kt_status <- write_overlaps(3)

par(mfrow=c(3,1), mar=rep(1,4))
pie(table(gt_status))
pie(table(ht_status))
pie(table(kt_status))


## now: e.g., running from py/: 
#     python binAroundBed.py ../data/bedfiles/h1_both_tads.bed 500000 > ../data/nonrepo/inbed/hbb.bins


### --- Adapted from 7_fig5boundaryEnrichments.R --- ###

buildBoundariesFull <- function(ct=c("h", "g", "k"), 
                                type=c("b", "c", "y", "n")){

  dir <- "data/nonrepo/split/"
  steps=25
  
  # filenames:: ct (h|g|k) + boundary type ([b]oth|[c]tcf|[y]y1|[n]one) + bins (b)
  # e.g. hyb.bins_HaibTfbs...
  fl <- paste0(dir, list.files(dir, pattern=paste0(ct, type, "b.bins_", ".*")))
  looper <- fl  
  # initiate all.dat df
  a <- read.delim(fl[1], header=F)
  all.dat <- data.frame(row.names=a$V1)
  vars <- colnames(h.dat)[-1]
  
  matches <- lapply(vars, function(x) grep(x, fl)[1])
  m <- fl[unlist(matches)]
  m <- na.omit(m)
  m <- unique(m)
  
  # Inspect matches:
  cbind(gsub(".*\\/", "", m), vars)
  looper <- m
  
  gdf <- data.frame(name=character(), bound=numeric(),
                    value=numeric())
  bins <- matrix(nrow=0, ncol=steps)
  for(i in 1:length(looper)){
    f <- looper[i]
    message(f)
    var <- vars[i]
    feat <- read.table(f)  
    ac <- matrix(feat$V5, ncol=steps, byrow=T)
    ac[is.nan(ac)] <- 0
    outr <- data.frame(name=f, bound=1:nrow(ac), 
                       feat=rep(var, nrow(ac)), ac)
    gdf <- rbind(gdf, outr)
  }
  gdf$ct <- ct
  gdf$type <- type
  return(gdf)
}

types <- c("b", "n", "y", "c")
h.all <- sapply(types, buildBoundariesFull, ct="h", simplify=F)
h.all <- do.call(rbind, h.all)

# g.all <- sapply(types, buildBoundariesFull, ct="g", simplify=F)
# g.all <- do.call(rbind, g.all)

k.all <- sapply(types, buildBoundariesFull, ct="k", simplify=F)
k.all <- do.call(rbind, k.all)

all.t <- rbind(h.all, k.all)

saveRDS(all.t, "data/rds/tad_split_boundaries.rds")

library("dplyr")
library("ggplot2")

comp <- group_by(all.t, ct, type, feat) %.% 
  summarise(bound.mean = mean(X12),
            edge.mean  = mean(c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21)),
            p = wilcox.test(X12, c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21))$p.value)

# Proper nomenclature, caps:
caps2 <- factor(comp$feat)
levels(caps2) <- c("ATF3", "CEBP", "CHD1", "CHD2", "MYC", "CTCF",
                   "DNase", "EGR1", "EZH2", "GABP", "JUND", "MAX",
                   "H2A.Z", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1",
                   "H3K4me2", "H3K4me3", "H3K79me2", "H3K9ac", "H3K9me3",
                   "H4K20me1", "MXI1", "NRSF", "P300",
                   "POL2", "RAD21", "SIX5", "SP1", "TAF1", "TBP", "YY1", 
                   "ZNF143", "GC", "GERP")
comp$feat <- as.character(caps2)

#pdf("figures/f5b_boundaryEnrichmentBubble.pdf", 9, 5)
ggplot(comp, aes(x=feat, y=-log10(p), col=type,
                 size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~.)+#, scales="free_y") + 
  geom_point(position=position_dodge(.15)) + theme_bw() +
  labs(list(size="Absolute difference\nat boundary",
            y=expression(-log[10](italic(p))),
            x="", col="Boundary type")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(range=c(2,7)) +
  guides(col=guide_legend(override.aes=list(size=4))) +
  geom_hline(yintercept=-log10(.01 / 297), linetype="dashed") +
  scale_colour_brewer(type="qual", palette=6) 
#dev.off()


