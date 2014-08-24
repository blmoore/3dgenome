############### Fig 1. Compartment profiles ################
# Generate compartment profiles for both the chromosome    #
# section in figure 1, and the entire genome supplementary #
# figure.                                                  #
############################################################
library("ggplot2")
library("gridExtra")
library("blmR")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
g.dat <- readRDS("data/rds/")
options(scipen=99)

chrs <- gsub("-.*", "", rownames(h.dat))
poss <- as.numeric(as.character(gsub(".*-", "", rownames(h.dat))))
fill <- data.frame(chrs, poss)
missing <- data.frame()

# Checking: is every possible bin covered? Need to fill in
# missing with NAs for plotting
for(k in unique(fill$chrs)){
  ps <- fill[fill$chrs == k, 2]
  o <- seq(min(ps), max(ps), 1e6)
  if(length(o) == length(ps)){
    cat(k, " OK\n")
  } else {
    cat(k, " Missing\n")
    missing <- rbind(missing, cbind(as.character(k), o[!o %in% ps]))
  }
}

missing

# Messy way of adding NAs for geom_line to not connext over 
# centromeres, unmappable regions, poorly covered by Hi-C etc.
eigens <- data.frame(rbind(
  cbind(h.dat$eigen, chrs, poss, "H1 hESC"),
  cbind(g.dat$eigen, chrs, poss, "GM12878"),
  cbind(k.dat$eigen, chrs, poss, "K562")))

missing <- cbind(NA, missing)
colnames(missing) <- NULL
addrows <- rbind(cbind(missing, ct="H1 hESC"),
      cbind(missing, ct="K562"),
      cbind(missing, ct="GM12878"))
colnames(addrows) <- colnames(eigens)
eigens <- rbind(eigens, addrows)

colnames(eigens) <- c("eig", "chr", "pos", "ct")
eigens$pos <- as.numeric(as.character(eigens$pos))
eigens$eig <- as.numeric(as.character(eigens$eig))
eigens$chr <- factor(eigens$chr, levels=paste0("chr", c(1:22, "X")))

# Supplementary figure 1. Compartment eigenvectors for each
# chromosome in three cell types. Stretched laterally to 
# edge of axis to make undulations more clear. Pick out one
# of these chromosomes to make up Figure 1 (e.g. chr2).
bigcs <- eigens[eigens$chr %in% paste0("chr", 1:11),]
smallc <- eigens[!eigens$chr %in% paste0("chr", 1:11),]

pdf("figures/suppl/s1_GenomewideWigglePlots.pdf", 8, 10)
grid.arrange(
  ggplot(bigcs, aes(x=pos/1000000, y=eig, group=ct, col=ct)) + 
    facet_wrap(~chr, nrow=1) + geom_line() + coord_flip() + 
    xlab("Chromosome position (Mb)") + theme_bw() +
    labs(colour = "Cell type") + ylab("") +
    scale_y_continuous(limits=c(-.22,.22), 
                       breaks=c(-.19,0,.19), labels=c("B", "", "A"))
  ,
  ggplot(smallc, aes(x=pos/1000000, y=eig, group=ct, col=ct)) + 
    facet_wrap(~chr, nrow=1) + geom_line() + coord_flip() + 
    xlab("Chromosome position (Mb)") + theme_bw() +
    labs(colour = "Cell type") + ylab("") +
    scale_y_continuous(limits=c(-.3,.3), 
                       breaks=c(-.28,0,.28), labels=c("B", "", "A"))
)
dev.off()

pdf("~/hvl/ice/plots/f1_chr2.pdf", 3, 7)
ggplot(subset(eigens, chr == "chr2"), 
       aes(x=pos/1000000, y=eig*-1, group=ct, col=ct)) + 
  facet_wrap(~chr, nrow=1) + geom_line() + coord_flip() + 
  xlab("Chromosome position (Mb)") + theme_bw() +
  labs(colour = "Cell type") + ylab("") +
  scale_y_continuous(limits=c(-.18,.18), breaks=c(-.16,0,.16), 
                     labels=c("A", "", "B"))
dev.off()
