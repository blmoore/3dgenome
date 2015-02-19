######## Call + analyse compartment/TAD boundaries ###########
# First pull out those regions comparments switch from open  #
# to closed, and the edge of TADs. Then calculate feature    #
# enrichment or depletion over the different boundary types. #
##############################################################
library("dplyr")
library("ggplot2")
require("gridExtra")
library("naturalsort")
library("plotrix")
require("RColorBrewer")
library("RHmm")
library("scales")
library("blmR")

g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

## External files required: boundaries/
## Once boundaries are generated, build a bins file around each:
##   python binAroundBed.py h_cb.bed 1500000 > h_cb_bins.bed
## Then use these bins with bigWigAverageOverBed to bin all of the
## bigWig files used in this work, i.e. 187 K562 vars, 115 Gm12878 etc.

h.100k <- read.table("data/bedfiles/h100.bed")
k.100k <- read.table("data/bedfiles/k100.bed")
g.100k <- read.table("data/bedfiles/g100.bed")

h.100k$states <- callStates(h.100k[,6])$state
k.100k$states <- callStates(k.100k[,6])$state
g.100k$states <- callStates(g.100k[,6])$state

colnames(k.100k) <- colnames(g.100k) <- colnames(h.100k) <- 
  c("id", "chr", "start", "end", "strand", "eig", "states")

callCBounds <- function(sfile){
  compartments <- sapply(paste0("chr", c(1:22, "X")), function(x) {
    cc <- sfile[sfile$chr == x,]
    rle <- rle(cc$states)
    # boundary blocks (i.e., last block of given state)
    bounds <- cc[cumsum(rle$lengths),]   
    bounds[,3:4] <- bounds[,3:4]+50000
    return(bounds)
  }, simplify=F)
  compartments <- do.call(rbind, compartments)
  outbed <- compartments[,c("chr", "start", "end")]
  outbed$names <- paste0(compartments$chr, "-", compartments$start)
  return(outbed)
}

write.bed <- function(bed, file)
  write.table(bed, file, quote=F, row.names=F, 
              col.names=F, sep="\t")

kb <- callCBounds(k.100k)
gb <- callCBounds(g.100k)
hb <- callCBounds(h.100k)

write.bed(kb, "data/bedfiles/k5_compartmentbounds.bed")
write.bed(gb, "data/bedfiles/gm_compartmentbounds.bed")
write.bed(hb, "data/bedfiles/h1_compartmentbounds.bed")

# Now need to generate bins around these boundaries. This is done 
# via binAroundBed.py, i.e.:
#   python binAroundBed.py h_compartmentbounds.bed 1500000 > h_cb_bins.bed
#   python binAroundBed.py g_compartmentbounds.bed 1500000 > g_cb_bins.bed
#   python binAroundBed.py k_compartmentbounds.bed 1500000 > k_cb_bins.bed

# Once this is done, use bigWigAverageOverBed with these bed files
# and all the 34ish (-gc) features from ENCODE (.bigwig format).
# Don't forget GC content! calc with (e.g.) :
#   bedtools nuc -fi <hg19.fa> -bed gm_40kb_cb_bins.bed | cut -f-4,6-6 | sed '1d' > Gm12878_gc_40kb.gc


## END compartments

## BEGIN TADs

procBounds <- function(b, ct=c("h1", "gm", "k5")){
  # Old / bad code, ideally clean up before release
  # order properly
  b <- b[naturalorder(b[,1]),]
  all.bounds <- b
  bounds <- sapply(paste0("chr", 1:22), function(x) {
    h <- b[b[,1] == x,]
    new.row <- cbind(chr=x, start=h[1:(nrow(h)-1),3],
                     end=h[2:nrow(h),2])
  })
  bounds <- as.data.frame(do.call(rbind, bounds))
  bounds$start <- as.numeric(as.character(bounds$start)) + 1
  bounds$end <- as.numeric(as.character(bounds$end))
  # dist == interspace distance between TADs
  bounds$diff <- bounds[,3] - bounds[,2]

  # Check proportions:
  #nrow(bounds[bounds$diff < 400000,]) / nrow(bounds) # 97% < 400kb
  #nrow(bounds[bounds$diff < 50000,]) / nrow(bounds) # 70% < 50kb
  tad.b <- bounds[bounds$diff < 400000,]
  tad.b$mids <- apply(tad.b[,2:3], 1, mean)
  new.end <- tad.b$mids + 20000
  new.start <- tad.b$mids - 20000
  options(scipen=999)
  bed.df <- data.frame(chr=tad.b$chr, start=new.start,
                       end=new.end, id=paste0(tad.b$chr, "-", new.start))
  write.table(bed.df, file=paste0("data/bedfiles/", ct, "_tadbounds.bed"), quote=F, 
              row.names=F, col.names=F, sep="\t")
  invisible()
}

# Load TAD calls
h.b <- read.table("data/text/h1.tads")
k.b <- read.table("data/text/k5.tads")
g.b <- read.table("data/text/gm.tads")

procBounds(h.b, "h1")
procBounds(k.b, "k5")
procBounds(g.b, "gm")

# User these bed files as above with binAroundBed.py, then
# run bigWigAverageOverBed to bin ENCODE bigwig data.

## End TADs

buildBoundariesDat <- function(ct=c("H1hesc", "Gm12878", "K562"), 
                               type=c("Compartments", "TADs")){
  ## Combine TADs and boundaries into big dataframe
  if(type == "Compartments"){
    dir <- "data/nonrepo/cluster/"
    steps <- 30
  } else {
    dir <- "data/nonrepo/cluster_tad/"
    steps <- 25
  }
  
  fl <- paste0(dir, list.files(dir, pattern=paste0(".*", ct,".*")))
  looper <- fl  
  # initiate all.dat df
  a <- read.delim(fl[1], header=F)
  all.dat <- data.frame(row.names=a$V1)
  vars <- c(colnames(h.dat)[-1], "gerp")
  
  matches <- lapply(vars, function(x) grep(x, fl)[1])
  m <- fl[unlist(matches)]
  m <- na.omit(m)
  m <- unique(m)
  
  # c(Inspect matches:
  cbind(gsub(".*\\/", "", m), vars)
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
    se <- apply(ac, 2, std.error)
    outr <- data.frame(name=f, pos=1:steps, 
                       mean=means, means.l=means-se, means.u=means+se)
    gdf <- rbind(gdf, outr)
  }
  gdf$feat <-  rep(vars, each=steps)
  gdf$ct <- ct
  gdf$type <- type
  return(gdf)
}

## error re: wrong number of replacement rows propbably means you 
## forgot to run GC seperately, see above bedtools cmd.

# Compartment boundaries:
h.c <- buildBoundariesDat("H1hesc", "Compartments")
g.c <- buildBoundariesDat("Gm12878", "Compartments")
k.c <- buildBoundariesDat("K562", "Compartments")

# TADs:
h.t <- buildBoundariesDat("H1hesc", "TADs")
g.t <- buildBoundariesDat("Gm12878", "TADs")
k.t <- buildBoundariesDat("K562", "TADs")

gdf <- rbind(h.c, g.c, k.c, # comps
             h.t, g.t, k.t) # TADs

# Colours for Fig. 5: TADs: blue; Compartments: red;
col <- brewer.pal(5, "Blues")[-c(1,2)]
colf <- paste0(col, "55")
col2 <- brewer.pal(5, "Reds")[-c(1,2)]
col2f <- paste0(col2, "33")

# Proper nomenclature, caps:
caps <- factor(gdf$feat)
levels(caps) <- c("ATF3", "CEBP", "CHD1", "CHD2", "MYC", "CTCF",
                  "DNase", "EGR1", "EZH2", "GABP", "GC", "GERP",
                  "H2A.Z", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1",
                  "H3K4me2", "H3K4me3", "H3K79me2", "H3K9ac", "H3K9me3",
                  "H4K20me1", "JUND", "MAX", "MXI1", "NRSF", "P300",
                  "POL2", "RAD21", "SIX5", "SP1", "TAF1", "TBP", "YY1", 
                  "ZNF143")
gdf$feat <- as.character(caps)

## Figure 5a: Selected average-o-grams for interesting features
pdf("figures/f5a_boundaryEnrichmentProfiles.pdf", 7.36, 4.27)
grid.arrange(
  ggplot(subset(gdf, type == "Compartments" & 
                  feat %in% c("CTCF", "H2A.Z", "YY1")),
         aes(x=pos, y=mean, ymin=means.l, ymax=means.u, 
             col=ct, fill=ct, shape=ct)) +
    geom_ribbon(alpha=I(.2), aes(linetype=NA)) +
    facet_wrap(~feat, scale="free_y", ncol=5) + 
    geom_point() + geom_line() + theme_bw() +
    scale_y_continuous(breaks=seq(0, 1.5, by=.1)) +
    scale_x_continuous(breaks=c(1, 14, 30),
                       labels=c("-1.5 Mb", "Boundary", "+1.5 Mb")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5)) +
    labs(list(x="", y="",
              fill="Compartment\ncell type", col="Compartment\ncell type", 
              shape="Compartment\ncell type")) +
    scale_color_manual(values=col2) + scale_fill_manual(values=col2f) 
  ,
  ## TADs::
  ggplot(subset(gdf, type == "TADs" & 
                  feat %in% c("CTCF", "H3K9me3", "POL2")), 
         aes(x=pos, y=mean, ymin=means.l, ymax=means.u, 
             col=ct, fill=ct, shape=ct)) +
    geom_ribbon(alpha=I(.2), aes(linetype=NA)) +
    facet_wrap(~feat, scale="free_y", ncol=5) + 
    geom_point() + geom_line() + theme_bw() +
    scale_x_continuous(breaks=c(1, 12, 25),
                       labels=c("-500 kb", "Boundary", "+500 kb")) +
    scale_y_continuous(breaks=seq(0, 1.5, by=.1)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5)) +
    labs(list(x="", y="Mean normalised signal per bin",
              fill="TAD cell type", col="TAD cell type", 
              shape="TAD cell type")) +
    scale_color_manual(values=col) + scale_fill_manual(values=colf)
  , nrow=2)
dev.off()


## All compartment boundaries, supplementary figure:
pdf("figures/suppl/s10_compartmentBoundaries.pdf", 12, 12)
ggplot(subset(gdf, type == "Compartments"),
       aes(x=pos, y=mean, ymin=means.l, ymax=means.u, 
           col=ct, fill=ct, shape=ct)) +
  geom_ribbon(alpha=I(.3), aes(linetype=NA)) +
  facet_wrap(~feat, scale="free_y", ncol=5) + 
  geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks=c(1, 14, 30),
                     labels=c("-1.5 Mb", "Boundary", "+1.5 Mb")) +
  #scale_y_continuous(breaks=seq(0,2, by=.1)) +
  labs(list(x="", y="Normalised ChIP-seq signal relative to input per 100 kb bin",
            fill="Cell type", col="Cell type", shape="Cell type")) +
  ggtitle("Compartment boundary feature enrichments") +
  scale_fill_manual(values=c("#0000ff48", "#FFA50048", "#ff000048")) +
  scale_colour_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) 
dev.off()

## All TAD boundaries, supplementary figure:
pdf("figures/suppl/s9_tadBoundaries.pdf", 12, 12)
ggplot(subset(gdf, type == "TADs"), 
       aes(x=pos, y=mean, ymin=means.l, ymax=means.u, 
           col=ct, fill=ct, shape=ct, show_guide=F)) +
  geom_ribbon(alpha=I(.3), aes(linetype=NA)) +
  facet_wrap(~feat, scale="free_y", ncol=5) + 
  geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks=c(1, 12, 25),
                     labels=c("-500 kb", "Boundary", "+500 kb")) +
  #scale_y_continuous(breaks=seq(0,2, by=.1)) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  labs(list(x="", y="Normalised ChIP-seq signal relative to input per 40 kb bin",
            fill="Cell type", col="Cell type", shape="Cell type")) +
  ggtitle("TAD boundary feature enrichments") +
  scale_fill_manual(values=c("#0000ff48", "#FFA50048", "#ff000048")) +
  scale_colour_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) 
dev.off()

## To test, need (row-wise) counts per position per feature!
buildBoundariesFull <- function(ct=c("H1hesc", "Gm12878", "K562"), 
                               type=c("Compartments", "TADs")){
  if(type == "Compartments"){
    dir <- "data/nonrepo/cluster/"
    steps <- 30
  } else {
    dir <- "data/nonrepo/cluster_tad/"
    steps <- 25
  }
  
  fl <- paste0(dir, list.files(dir, pattern=paste0(".*", ct,".*")))
  looper <- fl  
  # initiate all.dat df
  a <- read.delim(fl[1], header=F)
  all.dat <- data.frame(row.names=a$V1)
  vars <- c(colnames(h.dat)[-1], "gerp")
  
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

# Compartments:
hc.full <- buildBoundariesFull("H1hesc", "Compartments")
gc.full <- buildBoundariesFull("Gm12878", "Compartments")
kc.full <- buildBoundariesFull("K562", "Compartments")

full.c <- rbind(hc.full, gc.full, kc.full)

comp.p <- group_by(full.c, ct, type, feat) %.% 
          summarise(bound.mean = mean(X15),
            edge.mean  = mean(c(X1, X2, X3, X4, X5, X30, X29, X28, X27, X26)),
            p = wilcox.test(X15, c(X1, X2, X3, X4, X5, X30, X29, X28, X27, X26))$p.value)

# TADs:
ht.full <- buildBoundariesFull("H1hesc", "TADs")
gt.full <- buildBoundariesFull("Gm12878", "TADs")
kt.full <- buildBoundariesFull("K562", "TADs")

full.t <- rbind(ht.full, gt.full, kt.full)
saveRDS(full.t, "data/rds/tad_boundary_features.rds")

comp.t <- group_by(full.t, ct, type, feat) %.% 
  summarise(bound.mean = mean(X12),
            edge.mean  = mean(c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21)),
            p = wilcox.test(X12, c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21))$p.value)

comp <- rbind(comp.p, comp.t)

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

## Figure 5b; Bubble plot of p-value vs. features, scaled by effect size
pdf("figures/f5b_boundaryEnrichmentBubble.pdf", 9, 5)
ggplot(comp, aes(x=feat, y=-log10(p), col=type,
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
dev.off()
