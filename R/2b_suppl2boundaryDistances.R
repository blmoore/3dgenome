###### Compare metaphase G-bands with i-phase structure ######
# Compare giemsa bands by overlap of A/B compartments and by #
# boundary distances between TAD and compartment edges.      #
##############################################################
library("corrgram")
library("fGarch")
library("ggplot2")
library("gridExtra")
library("reshape2")
library("stats")
library("blmR")
set.seed(42)

csizes <- read.table("data/text/hg19.chrom.sizes.txt")

rand.tads <- function(intb){
  # null TADS, fit to this dist?:
  sdist <- htb[,3] - htb[,2]
  pars <- snormFit(sdist)$par
  cs <- csizes
  ns <- table(intb[,1]) 
  rdf <- data.frame(V1=character(),
                    V2=numeric(),
                    V3=numeric())
  for ( x in 1:length(ns) ){
    max <- cs[cs[,1] == names(ns)[x],2]
    possbins <- seq(0, max, 40000)
    rand.s <- sample(possbins, ns[x])
    rand.e <- sample(possbins, ns[x])
    rdf <- rbind(rdf, cbind(V1=names(ns)[x],
                            V2=rand.s,
                            V3=rand.e))
  }
  rdf[,2] <- as.numeric(as.character(rdf[,2]))
  rdf[,3] <- as.numeric(as.character(rdf[,3]))
  return(rdf)
}

rand.comps <- function(intb){
  #null comp bound dist
  cs <- csizes
  ns <- table(intb[,1]) 
  rdf <- data.frame(V1=character(),
                    V2=numeric(),
                    V3=numeric())
  for ( x in 1:length(ns) ){
    max <- cs[cs[,1] == names(ns)[x],2]
    possbins <- seq(0, max, 100000)
    # get same number of start points
    bound <- sample(possbins, ns[x])
    rand.s <- bound - 50000
    # or just random as with starts
    rand.e <- bound + 50000
    rdf <- rbind(rdf, cbind(V1=names(ns)[x],
                            V2=rand.s,
                            V3=rand.e))
  }
  rdf[,2] <- as.numeric(as.character(rdf[,2]))
  rdf[,3] <- as.numeric(as.character(rdf[,3]))
  return(rdf)
}

tad.dist <- function(tb, k, g, r){
  ## calculate min distance to another tad boundary 
  # tb1 == tad boundaries you're comparing against
  # tb2, tb3 == the other two
  omat <- tb
  omat$kdiff.s <- NA
  omat$kdiff.e <- NA
  omat$gdiff.s <- NA
  omat$gdiff.e <- NA
  omat$rand.s <- NA
  omat$rand.e <- NA
  
  for( i in 1:nrow(tb)){
    chr <- tb[i,1]
    b1 <- tb[i,2]
    b2 <- tb[i,3]
    kcmp <- c(t(k[k[,1] == chr,2:3]))
    kdiff <- c(min(abs(b1-kcmp)), min(abs(b2-kcmp)))
    omat[i,c("kdiff.s", "kdiff.e")] <- kdiff
    gcmp <- c(t(g[g[,1] == chr,2:3]))
    gdiff <- c(min(abs(b1-gcmp)), min(abs(b2-gcmp)))
    omat[i,c("gdiff.s", "gdiff.e")] <- gdiff
    rcmp <- c(t(r[r[,1] == chr,2:3]))
    rdiff <- c(min(abs(b1-rcmp)), min(abs(b2-rcmp)))
    omat[i,c("rand.s", "rand.e")] <- rdiff
  }
  outdf <- data.frame(GM12878=c(omat$gdiff.s, omat$gdiff.e),
                      K562=c(omat$kdiff.s, omat$kdiff.e),
                      Null=c(omat$rand.s, omat$rand.e))
  return(outdf)
}

comp.dist <- function(tb, k, g, r){
  ## calculate min distance to another compartment boundary 
  # tb1 == tad boundaries you're comparing against
  # tb2, tb3 == the other two
  omat <- tb
  for( i in 1:nrow(tb) ){
    chr <- tb[i,1]
    b <- mean(as.numeric(tb[i,2:3]))
    # nb uses both start and end -- not useful for boundaries 
    # themselves
    kcmp <- apply(k[k[,1] == chr, 2:3], 1, mean)
    omat$K562[i] <- min(abs(b-kcmp))
    gcmp <- apply(g[g[,1] == chr, 2:3], 1, mean)
    omat$GM12878[i] <- min(abs(b-gcmp))
    # Do want to use fake comps, both sides
    rcmp <- apply(r[r[,1] == chr, 2:3], 1, mean)
    omat$Null[i] <- min(abs(b-rcmp))
  }
  return(as.data.frame(omat[,c("K562", "GM12878", "Null")]))
}

## 1) TADs:
## These boundaries are generated in 7_*.R
htb <- read.table("data/bedfiles/h1_tadbounds.bed")
ktb <- read.table("data/bedfiles/k5_tadbounds.bed")
gtb <- read.table("data/bedfiles/gm_tadbounds.bed")

rtb <- rand.tads(htb)
tdf <- tad.dist(htb, ktb, gtb, rtb)
taddf <- melt(tdf)
taddf <- cbind(taddf, bound="TADs")

# 2) Compartments:
## Original bounds (100 kb but according to 1 Mb bins)
hcb <- read.table("data/bedfiles/h1_compartmentbounds.bed")
kcb <- read.table("data/bedfiles/k5_compartmentbounds.bed")
gcb <- read.table("data/bedfiles/gm_compartmentbounds.bed")

rcb <- rand.comps(hcb)
cdf <- comp.dist(hcb, kcb, gcb, rcb)
compdf <- melt(cdf)
compdf <- cbind(compdf, bound="Compartments")
#saveRDS(compdf, "data/rds/compECDF_dat.rds")

# 3) Graph both:
bothdf <- rbind(taddf, compdf)
options(scipen=999)
cat("Drawing: figures/suppl/s2a_boundsEcdf.pdf\n")
pdf("figures/suppl/s2a_boundsEcdf.pdf", 6, 3.5)
ggplot(bothdf, aes(x=value+1, col=variable)) +
  stat_ecdf(geom="line", size=1.1) + scale_x_log10() +
  facet_grid(.~bound, scales="free_x") + theme_bw() +
  theme(legend.position=c(.15,.8), legend.background=element_blank()) +
  labs(col="", y="ECDF", x="Distance to nearest H1 boundary (bp)") +
  scale_color_brewer(type="qual", palette=3)
dev.off()

# 4) Misc., statistical tests etc.
median(tdf$K562)    # median: 300kb apart
median(tdf$GM12878) # median: 180 kb

# TADs: Prop of H1 bounds w/ ct bound <= 40 kb
tad.g <- 100 * nrow(tdf[tdf$GM12878 <= 40000,]) / nrow(tdf) # 33.38 %
tad.k <- 100 * nrow(tdf[tdf$K562 <= 40000,]) / nrow(tdf)    # 31.30 %
tad.null <- 100 * nrow(tdf[tdf$Null <= 40000,]) / nrow(tdf) # 17.84 %

# Compartments: Prop of H1 bounds w/ ct bound <= 100 kb 
comp.g <- 100 * nrow(cdf[cdf$GM12878 <= 100000,]) / nrow(cdf) # 36.90 %
comp.k <- 100 * nrow(cdf[cdf$K562 <= 100000,]) / nrow(cdf)    # 35.12 %
comp.null <- 100 * nrow(cdf[cdf$Null <= 100000,]) / nrow(cdf) # 7.49  %

ks.test(tdf$K562, tdf$Null)
ks.test(cdf$K562, cdf$Null)

## Fig 2b: corrgram of eigenvector correlations:
cat("Drawing: figures/suppl/s2b_compartmentCorrgram.pdf\n")
pdf("figures/suppl/s2b_compartmentCorrgram.pdf", 6, 6)
corrgram(data.frame(GM12878  = g.dat$eigen,
                   `H1 hESC` = h.dat$eigen,
                    K562     = k.dat$eigen),
         lower.panel = panel.pts2, pch = 20,
         upper.panel = panel.conf,
         diag.panel  = panel.density)
dev.off()
