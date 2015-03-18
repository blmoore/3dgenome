## split ctcf and yy1, i.e. is ctcf * yy1 more
## significantly enriched at boundaries than
## ctcf or yy1
library("ggplot2")
library("dplyr")

g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

#all.t <- readRDS("data/rds/tad_boundary_features.rds")

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
  both <- subsetByOverlaps(yy1, ctcf)
  
  # now intersect with TAD boundaries:
  tads <- read_peaks(paste0("data/bedfiles/", short, "_tadbounds.bed"))
  yy1_tads <- subsetByOverlaps(tads, yy1)
  message(paste0(signif(100*length(yy1_tads) / length(tads), digits=3), 
                 "% of tads have yy1"))
  
  ctcf_tads <- subsetByOverlaps(tads, ctcf)
  # ~61% of TAD boundaries (40kb) have signif CTCF peak::
  message(paste0(signif(100 * length(ctcf_tads) / length(tads), digits=3),
                 "% of tads have ctcf"))
  
  #message(paste0(length(both))
  print(length(yy1))
  marked_tads <- subsetByOverlaps(tads, both)
  # ~12% of TAD boundaries contain both::
  message(paste0(signif(100*length(marked_tads) / length(tads), digits=3),
                 "% of tads have ctcf + yy1"))
  
  # problem is: some have both ctcf and yy1 but the two aren't overlapping
#   which <- ifelse(tads %in% marked_tads, "both", 
#                   ifelse(tads %in% yy1_tads, "yy1",
#                          ifelse(tads %in% ctcf_tads, "ctcf", 
#                                 "none")))
  both_n <- which(tads %in% marked_tads)
  ctcf_n <- which(tads %in% ctcf_tads & !tads %in% marked_tads)
  
  # make sure no ctcf only tads are both
  stopifnot(length(which(ctcf_n %in% both_n)) == 0)
  
  yy1_n <- which(tads %in% yy1_tads & !tads %in% marked_tads)
  
  # should be exclusively yy1
  stopifnot(length(which(yy1_n %in% both_n)) == 0)

  none_n <- which(!tads %in% yy1_tads & !tads %in% ctcf_tads)
  stopifnot(length(which(none_n %in% c(yy1_n, ctcf_n, both_n))) == 0)


  message("Writing bed files...") 
  export.bed(tads[ctcf_n,], paste0("data/bedfiles/", short, "_ctcf_tads.bed"))
  export.bed(tads[both_n,], paste0("data/bedfiles/", short, "_both_tads.bed"))
  export.bed(tads[yy1_n,], paste0("data/bedfiles/", short, "_yy1_tads.bed"))
  export.bed(tads[none_n,], paste0("data/bedfiles/", short, "_none_tads.bed"))
  
  message(length(tads))
  
  return(invisible())
  #return(as.factor(which))
}
  
# gm12878
gt_status <- write_overlaps(1)

# h1
ht_status <- write_overlaps(2)

# k5
kt_status <- write_overlaps(3)

# Run bash script from sh/:
#   sh splitTadBounds.sh
# This generates BED files for input to bigWigAverageOverBed,
# see exampleClusterScript.sh for rough idea.

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
  # tiry up var names
  gdf$ct <- if(ct == "h") "H1 hESC" else 
    if(ct == "g") "GM12878" else "K562"
  gdf$type <- if(type == "b") "CTCF + YY1" else 
    if(type == "y") "YY1" else if(type == "c") "CTCF" else "None"
  return(gdf)
}

types <- c("b", "n", "y", "c")
h.all <- sapply(types, buildBoundariesFull, ct="h", simplify=F)
h.all <- do.call(rbind, h.all)

g.all <- sapply(types, buildBoundariesFull, ct="g", simplify=F)
g.all <- do.call(rbind, g.all)

k.all <- sapply(types, buildBoundariesFull, ct="k", simplify=F)
k.all <- do.call(rbind, k.all)

all.t <- rbind(h.all, g.all, k.all)

#saveRDS(all.t, "data/rds/tad_split_boundaries.rds")
#all.t <- readRDS("data/rds/tad_split_boundaries.rds")

library("dplyr")
library("ggplot2")

comp <- group_by(all.t, ct, type, feat) %>% 
  summarise(bound.mean = mean(X12),
            edge.mean  = mean(c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21)),
            p = wilcox.test(X12, c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21))$p.value)

# comp_all <- group_by(all.t, ct, feat) %>% 
#   summarise(bound.mean = mean(X12),
#             edge.mean  = mean(c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21)),
#             p = wilcox.test(X12, c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21))$p.value)

tidy_comp <- function(comp, typed=T){
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
  
  if(typed == T)
    comp$type <- factor(comp$type, levels=c("CTCF + YY1", "YY1", "CTCF", "None"))
  
  comp
}

comp <- tidy_comp(comp)
comp_all <- tidy_comp(comp_all, F)

## sample sizes per boundary type
sizes <- all.t %>% group_by(ct, type) %>% 
  #filter(type != "CTCF + YY1") %>% 
  summarise(sample=n()/34)
  #summarise(sample=max(bound))

# should be: {Gm: 1662, H1: 2810, K5: 1449};
group_by(sizes, ct) %>% summarise(tot=sum(sample))

sizes$x <- "CHD1"
sizes$y <- c(120, 140, 110, 130) # hack for ordering
sizes$type <- factor(sizes$type, levels=c("CTCF + YY1", "YY1", "CTCF", "None"))

pdf("figures/suppl/splitBoundaryEnrichmentBubble_v3.pdf", 9, 7)
ggplot(comp, aes(x=feat, y=-log10(p), col=type,
                 size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~., scales="free_y") +
  #geom_point(data=comp_all, inherit.aes=F, aes(x=feat, y=-log10(p)),
  #           size=I(3), col=I("grey70")) +
  geom_point(position=position_dodge(.15)) + theme_bw() +
  labs(list(size="Absolute difference\nat boundary",
            y=expression(-log[10](italic(p))),
            x="", col="Boundary type")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(range=c(2,7)) +
  guides(col=guide_legend(override.aes=list(size=4))) +
  geom_hline(yintercept=-log10(.01 / 297), linetype="dashed") +
  scale_colour_brewer(type="qual", palette=6) +
  geom_text(data=sizes, aes(label=paste0("n = ", sample), 
                            x=x, y=y, col=type), inherit.aes=F,
            size=I(3.5))
dev.off()


## individual profiles plot:


buildBoundariesDat <- function(ct=c("h", "g", "k"), 
                               type=c("b", "c", "y", "n")){
  
  dir <- "data/nonrepo/split/"
  steps=25
  
  fl <- paste0(dir, list.files(dir, pattern=paste0(ct, type, "b.bins_", ".*")))
  looper <- fl  
  # initiate all.dat df
  a <- read.delim(fl[1], header=F)
  all.dat <- data.frame(row.names=a$V1)
  vars <- colnames(h.dat)[-c(1, ncol(h.dat))]
  
  matches <- lapply(vars, function(x) grep(x, fl)[1])
  m <- fl[unlist(matches)]
  m <- na.omit(m)
  m <- unique(m)
  
  # Inspect matches:
  # print(cbind(gsub(".*\\/", "", m), vars))
  looper <- m
  
  gdf <- data.frame(name=character(), pos=numeric(),
                    mean=numeric(), mean.l=numeric(), 
                    mean.u=numeric())
  bins <- matrix(nrow=0, ncol=steps)
  for(f in looper){
    
    message(f)
    feat <- read.table(f)  
    ac <- matrix(feat$V5, ncol=steps, byrow=T)
    ac[is.nan(ac)] <- 0
    
    means <- colMeans(ac)
    ## NB Confidence intervals ? *1.96, Std error? as-is
    se <- apply(ac, 2, plotrix::std.error)
    outr <- data.frame(name=f, pos=1:steps, 
                       mean=means, means.l=means-se, means.u=means+se)
    gdf <- rbind(gdf, outr)
  }
  gdf$feat <-  rep(vars, each=steps)
  gdf$ct <- ct
  gdf$type <- type
  return(gdf)
}

types <- c("b", "n", "y", "c")
h_bounds <- suppressMessages(sapply(types, buildBoundariesDat, ct="h", simplify=F))
h_bounds <- do.call(rbind, h_bounds)

# g_bounds <- sapply(types, buildBoundariesDat, ct="g", simplify=F)
# g_bounds <- do.call(rbind, g.bounds)

k_bounds <- suppressMessages(sapply(types, buildBoundariesDat, ct="k", simplify=F))
k_bounds <- do.call(rbind, k_bounds)

all_bounds <- rbind(h_bounds, k_bounds)

# Colours for Fig. 5: TADs: blue; Compartments: red;
col <- brewer.pal(5, "Blues")[-c(1,2)]
colf <- paste0(col, "55")
col2 <- brewer.pal(5, "Reds")[-c(1,2)]
col2f <- paste0(col2, "33")

# Proper nomenclature, caps:
caps <- factor(all_bounds$feat)
levels(caps) <- c("ATF3", "CEBP", "CHD1", "CHD2", "MYC", "CTCF",
                  "DNase", "EGR1", "EZH2", "GABP", "GC", "GERP",
                  "H2A.Z", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1",
                  "H3K4me2", "H3K4me3", "H3K79me2", "H3K9ac", "H3K9me3",
                  "H4K20me1", "JUND", "MAX", "MXI1", "NRSF", "P300",
                  "POL2", "RAD21", "SIX5", "SP1", "TAF1", "TBP", "YY1", 
                  "ZNF143")
all_bounds$feat <- as.character(caps)

ggplot(subset(all_bounds, feat %in% c("CTCF", "YY1")),
       aes(x=pos, y=mean, ymin=means.l, ymax=means.u, 
           col=ct, fill=ct, shape=type, group=ct)) +
  geom_ribbon(alpha=I(.2), aes(linetype=NA)) +
  facet_wrap(type~feat, scale="free_y", ncol=5) + 
  geom_point() + geom_line() + theme_bw() +
  #scale_y_continuous(breaks=seq(0, 1.5, by=.1)) +
  #scale_x_continuous(breaks=c(1, 10.5, 25),
  #                   labels=c("-1.5 Mb", "Boundary", "+1.5 Mb")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5)) +
  labs(list(x="", y="",
            fill="Compartment\ncell type", col="Compartment\ncell type", 
            shape="Compartment\ncell type")) +
  scale_color_manual(values=col2) + scale_fill_manual(values=col2f) 


