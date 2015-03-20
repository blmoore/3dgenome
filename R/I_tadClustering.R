library("dplyr")
library("ggplot2")
library("gplots")
library("gridExtra")
library("RColorBrewer")
library("reshape2")
library("blmR")
set.seed(42)

g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

h1_tads <- add_bed_id("data/text/h1.tads")
write_bed(h1_tads, "data/bedfiles/h1_tads.bed")

gm_tads <- add_bed_id("data/text/gm.tads")
write_bed(gm_tads, "data/bedfiles/gm_tads.bed")

k5_tads <- add_bed_id("data/text/k5.tads")
write_bed(k5_tads, "data/bedfiles/k5_tads.bed")


buildTadDat <- function(ct=c("H1hesc", "Gm12878", "K562")){

  dir = "data/nonrepo/tadagg/"
  ct_short <- tolower(substr(ct, 1, 2))
  fl <- paste0(dir, list.files(dir, pattern=paste0("^", ct_short, ".*", ct,".*")))
  looper <- fl
  
  # initiate all.dat df
  a <- read.delim(fl[1], header=F)
  all.dat <- data.frame(row.names=a$V1)
  
  vars <- c("Atf3", "Cebp", "Chd1", "Chd2", "Cmyc",
            "Ctcf", "Dnase", "Egr1", "Ezh2", "Gabp",
            "Jund", "Max", "H2az", "H3k27ac", "H3k27me3",
            "H3k36me3", "H3k4me1", "H3k4me2", "H3k4me3", "H3k79me2", 
            "H3k9ac", "H3k9me3", "H4k20me1", "Mxi1", "Nrsf", 
            "P300", "Pol2", "Rad21", "Six5", "Sp1", 
            "Taf1", "Tbp", "Yy1", "Znf143")
  
  matches <- lapply(vars, function(x) grep(x, fl)[1])
  m <- fl[unlist(matches)]
  m <- na.omit(m)
  m <- unique(m)
  
  # Inspect matches:
  # cbind(gsub(".*\\/", "", m), vars)
  looper <- m
  
  # Here we know length of files, i.e. outfilename, check if already written
  # and don't remake if it is.
  outfilename <- paste0("data/rds/tadagg_", ct, "_", length(vars), "Vars.rds")
  
  
  cat("Reading in files... \n\n")
  for (i in looper){
    mark <- gsub(".*binnedBigWigs/(.*?)\\..*","\\1", i)
    print(mark)
    b <- read.delim(i, header=F)
    # NB: b$5 = mean0, missing vals = 0, b$6 = mean of nonzero values
    all.dat <- cbind(all.dat, new = b[,5])
    colnames(all.dat)[ncol(all.dat)] <- mark
    rm(b)
  }
  
  # Double-check column names 
  #cat("\n\nCheck these names match !\n")
  #print(cbind(colnames(all.dat), vars))
  colnames(all.dat) <- vars
  
  saveRDS(all.dat, outfilename)
  return(all.dat)
}

ktads <- buildTadDat("K562")
gtads <- buildTadDat("Gm12878")
htads <- buildTadDat("H1hesc")

## Parallel pv clustering, slow ##
# cli <- makeCluster(6, type="MPI")
# kpv <- parPvclust(cli, t(ktads), nboot=100, r=seq(.6, 1.2, by=.2))
# hpv <- parPvclust(cli, t(htads))
# gpv <- parPvclust(cli, t(gtads))
# stopCluster(cli)

get_comps <- function(dat, write=F){
  ct <- substr(deparse(substitute(dat)), 0, 1)
  
  comp <- data.frame(chr=as.character(gsub("-.*", "", rownames(dat))),
                     pos=as.numeric(gsub(".*-", "", rownames(dat))), 
                     comp=callStates(dat$eigen)$state)
  
  comp <- droplevels(subset(comp, chr != "chrX"))
  if(write == T)
    write_bed(comp, paste0("data/bedfiles/", ct, "_comps.bed"))
  comp
}

hcomp <- get_comps(h.dat)
gcomp <- get_comps(g.dat)
kcomp <- get_comps(k.dat)

call_comp <- function(j, tads, compartments){
  # 1) get compartments on chromosome
  rcb <- compartments[compartments$chr == tads[j,1],]
  
  to_cb <- tads[j,2]-rcb[,2]
  dist <- to_cb[which.min(to_cb[to_cb > 0])]
  if(length(dist) == 0) return(NA) # doesn't fit into any bin
  
  if(tads[j,3] - tads[j,2] < 1e6 & dist > 2e5)
    # size < 1 Mb fits into a single bin; 200 kb from boundary means > 80%
    if(!all(to_cb < 0))
      return(rcb[which.min(to_cb[to_cb > 0]), "comp"])
    else
      return(rcb[nrow(rcb), "comp"])
  else{
    # overlaps multiple bins, more tricky
    # get first bin, last bin calc % of states
    
    start <- rcb[which.min(to_cb[to_cb > 0]),]
    to_ce <- tads[j,3]-rcb[,2]
    end <- rcb[which.min(to_ce[to_ce > 0]),]
    comps <- compartments[rownames(compartments) %in% 
                            rownames(start):rownames(end),]$comp
    
    if(length(unique(comps)) == 1) return(unique(comps)) else {
      return(NA)
    }
  }
  # 2 == active; 1 == inactive; NA == mixed
  invisible()
}

comp_wrapper <- function(tads, comps)
  sapply(1:nrow(tads), call_comp, tads=tads, compartments=comps)

hc <- comp_wrapper(h1_tads, hcomp)
gc <- comp_wrapper(gm_tads, gcomp)
kc <- comp_wrapper(k5_tads, kcomp)

# heatmap with multiple row colourbars
devtools::source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

tad_heatmap2 <- function(tad_dat, compartments, k=3, km=F){
  set.seed(42)
  mat <- as.matrix(log2(tad_dat + 1))

  if(km == T){
    #clust <- kmeans(mat, centers=k)
    clust <- kmeans(mat[,c("Dnase", "H3k27me3", "H3k9me3")], centers=k, nstart=10)
    
    # non-consistent k-means colouring:
#     sidebars <- t(cbind(
#       RColorBrewer::brewer.pal(k+6, "Paired")[-c(1:6)][clust$cluster],
#       RColorBrewer::brewer.pal(k+2, "Paired")[-c(1:2)][factor(compartments)]))
#     
    # hand-hacked consistent colours (k=3):
    # red == null (bpal[1]); blue == active (bpal[2]); orange == PcG (bpal[3])
    bpal <- RColorBrewer::brewer.pal(7, "Set3")[4:7]
    ct <- substr(deparse(substitute(tad_dat)), 0, 1)
    
    bpal <- if(ct == "h") bpal[c(3,1,2)] else 
      if(ct == "k") bpal[c(2,1,3)] else bpal[c(1,3,2)]
    
    sidebars <- t(cbind(bpal[clust$cluster],
      RColorBrewer::brewer.pal(k+2, "Paired")[-c(1:2)][factor(compartments)]))
    
    rownames(sidebars) <- c("Clusters","Compartment")
    rsize <- 2
  } else {
    sidebars <- matrix(RColorBrewer::brewer.pal(k+2, "Paired")[-c(1:2)][factor(compartments)], 
                       nrow=1)
    rownames(sidebars) <- c("Compartment")
    rsize <- 1.5
  }

  sidebars[is.na(sidebars)] <- "#ffffff"
  hclust.fun <- function(...) hclust(..., method="complete")
  
  heatmap.3(mat, trace="none", labRow="",
            hclustfun=hclust.fun, key=F,
            #density.info="density",
            col=colorRampPalette(c("white", "steelblue4", "navy"))(256), 
            RowSideColors=sidebars, RowSideColorsSize=rsize,
            lhei=c(.1, .9), margins=c(6.5, 1),
            useRaster=F)
  invisible()
}


pdf("figures/inkscape/h1_tad_hm_3.pdf", 5, 5)
tad_heatmap2(htads, hc, k=2, km=T)
dev.off()

pdf("figures/inkscape/gm_tad_hm_3.pdf", 5, 5)
tad_heatmap2(gtads, gc, k=2, km=T)
dev.off()

pdf("figures/inkscape/k5_tad_hm_3.pdf", 5, 5)
tad_heatmap2(ktads, kc, k=2, km=T)
dev.off()

tad_tile <- function(tads, compartments, k=3){
  htp <- tads
  htp$block <- rownames(tads)
  htp$comp <- factor(compartments)
  htp$kmeans <- factor(kmeans(as.matrix(log2(tads[,c("Dnase", "H3k27me3", "H3k9me3")])+1)
                              , centers=k)$cluster)
  htp <- reshape2::melt(htp, id.vars=c("block", "comp", "kmeans"))
  
  hpp <- htp %>% group_by(comp, variable) %>% summarise(mean=mean(value))
  hp2 <- htp %>% group_by(kmeans, variable) %>% summarise(mean=mean(value))
  
  # reverse factor levels
  hpp$variable <- factor(hpp$variable, levels=rev(levels(hpp$variable)))
  hp2$variable <- factor(hp2$variable, levels=rev(levels(hp2$variable)))
  
  #pdf("figures/suppl/compartment_tile.pdf", 4.5, 7)
  grid.arrange(
    ggplot(subset(hpp, comp %in% c(1, 2)),
           aes(x=variable, y=comp, fill=log2(mean+1))) + 
      geom_tile() + coord_flip() + theme_minimal() +
      scale_fill_gradient(low="white", high=scales::muted("red"),
                          guide=guide_colorbar(title.position="top",
                                               title.hjust=.5)) +
      theme(legend.position="top") + 
      labs(x="", fill=expression(Mean~signal~(log[2])),
           y="Compartment"),
    ggplot(hp2, aes(x=variable, y=kmeans, fill=log2(mean+1))) + 
      geom_tile() + coord_flip() + theme_minimal() +
      scale_fill_gradient(low="white", high=scales::muted("blue"),
                          guide=guide_colorbar(title.position="top",
                                               title.hjust=.5)) +
      theme(legend.position="top") + 
      theme(axis.text.y=element_text(hjust=.5, vjust=.5))+
      labs(x="", fill=expression(Mean~signal~(log[2])),
           y="Clusters"),
    nrow=1
  )
  #dev.off()
  invisible()
}

svg("figures/inkscape//gm_tad_tiles.svg", 4.5, 7)
tad_tile(gtads, gc)
dev.off()

svg("figures/inkscape//h1_tad_tiles.svg", 4.5, 7)
tad_tile(htads, hc)
dev.off()

svg("figures/inkscape//k5_tad_tiles.svg", 4.5, 7)
tad_tile(ktads, kc)
dev.off()


## cluster on e.g. Pol2 (active) H3k27me3 (PcG) and H3k9me3 (null)?

## optimal number of clusters::
library("NbClust")
nb <- NbClust(ktads, min.nc=2, max.nc=10, method="kmeans")
# 8 proposed 2 (7 proposed 3)
nh <- NbClust(htads, min.nc=2, max.nc=10, method="kmeans")
# 10 proposed 2
ng <- NbClust(gtads, min.nc=2, max.nc=10, method="kmeans")
# 12 proposed 2
ncol(nb$Best.nc)


