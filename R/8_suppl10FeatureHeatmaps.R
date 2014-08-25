######## Cluster input fatures, draw corr. heatmap ##########
# Look at chromHMM / SegWay combined annotations in flipped #
# open, flipped closed compartments and test for enrichment #
# or depletion. For some states compare regions that are    #
# shared between cell types and those that are specific to  #
# one (cell type-specfic vs. shared; e.g. enhancers).       #
#############################################################
library("gplots")
library("pvclust")
library("RColorBrewer")
library("snow")
library("blmR")
set.seed(42)

g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

## Bootstrapped heirarchical clustering to estimate
## the significance (according to persistance across
## subsamples) of observed clusters.
parClust <- function(cluster, dat)
  parPvclust(cluster, cor(dat[, -1]), 
             nboot=5000)

caps <- function(c){
  c <- c("eigen", toupper(c[-1]))
  gsub("H2AZ", "H2A.Z",
       gsub("DNASE", "DNase", 
            gsub("AC", "ac", 
                 gsub("ME", "me", c))))
}

colnames(h.dat) <- colnames(g.dat) <- 
  colnames(k.dat) <- 
  caps(colnames(h.dat))


cli <- makeCluster(4, type="MPI")
kpv <- parClust(cli, k.dat)
hpv <- parClust(cli, h.dat)
gpv <- parClust(cli, g.dat)
stopCluster(cli)

clus.heat <- function(dat, clus, main){
  #dat = k.dat
  #clus = kpv
  vars <- colnames(dat)[-1][order.dendrogram(as.dendrogram(clus$hclust))]
  cs <- pvpick(clus)$clusters
  fcs <- sapply(1:length(cs), function(i) data.frame(name=c(cs[i]), 
                                 clus.n=rep(i, length(cs[i]))), simplify=F)
  fcs <- lapply(fcs, "names<-", c("var", "clus"))
  fcs <- do.call(rbind, fcs)
  all <- tryCatch({
    all <- data.frame(var=colnames(dat)[-1][!colnames(dat)[-1] %in% fcs[,1]], clus=0)
    merge(fcs, all, all=T)
    }, error=function(e) {
      # If all features fall into clusters, all will be empty
      return(fcs)
    })
  
  c.ind <- all[match(colnames(dat)[-1], all[,1]),]
  cols <- c("grey90", brewer.pal(6, "Pastel1"))
  
  heatmap.2(cor(dat[, -1]), 
            dendrogram="row",
            Rowv=as.dendrogram(clus$hclust),
            Colv=as.dendrogram(clus$hclust), 
            trace="none", main=main,
            ColSideColors=cols[c.ind[,2]+1],
            RowSideColors=cols[c.ind[,2]+1],
            col=colorRampPalette(brewer.pal(9, "RdBu"))(64),
            keysize=1, symkey=T,
            margins=c(7,7),
            cexRow = 1.3, cexCol = 1.3,
            breaks=seq(-1,1, length.out=65),
            zlim=c(-1,1))
}

pdf("figures/suppl/s4_gfeatmap.pdf", 7, 7)
clus.heat(g.dat, gpv, "GM12878")
dev.off()

pdf("figures/suppl/s4_hfeatmap.pdf", 7, 7)
clus.heat(h.dat, hpv, "H1 hESC")
dev.off()

pdf("figures/suppl/s4_kfeatmap.pdf", 7, 7)
clus.heat(k.dat, kpv, "K562")
dev.off()

