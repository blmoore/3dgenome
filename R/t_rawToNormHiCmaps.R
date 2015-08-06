## contact maps (split?) of raw counts -> normalised maps

library("ggplot2")
library("gridExtra")
library("readr")
library("reshape2")

melt_hic <- function(table, norm=F){
  # convert square interaction matrix to pairwise counts/IFs
  rownames(table) <- colnames(table) <- paste0("f", 1:nrow(table))
  table <- if(norm) 100 * table / max(table) else table
  table <- melt(as.matrix(table))
  table[,1:2] <- apply(table[,1:2], 2, function(i) as.integer(gsub("f", "", i)))
  table
}

plot_raw <- function(df){
  p <- ggplot(df, aes(x=Var1*.04, y=Var2*.04, fill=value+1)) +
      geom_raster() + 
      scale_fill_gradient2(name="IF", trans=scales::log10_trans(),
        low="white", high="steelblue4", na.value="grey70", 
        breaks=c(1, 10, 100, floor(max(df[,3]+1)))) +
      theme_minimal() +
      labs(x="Position along chromosome 18 (Mb)", y="", 
        fill="Contacts") +
    coord_fixed(xlim=c(35, 45), ylim=c(35, 45)) +
    scale_x_continuous(breaks=c(35, 40, 45)) +
    scale_y_continuous(breaks=c(35, 40, 45))
  return(p)
}

plot_hic <- function(df){
  p <- ggplot(df, aes(x=Var1*.04, y=Var2*.04, fill=value+1)) +
    geom_raster() + 
    scale_fill_gradient2(name="IF", trans=scales::log10_trans(),
      low="white", high="steelblue4", na.value="grey70", 
      breaks=c(1, 10, 100)) +
    theme_minimal() +
    labs(x="Position along chromosome 18 (Mb)", y="", fill="IF") +
    coord_fixed(xlim=c(35, 45), ylim=c(35, 45)) +
  scale_x_continuous(breaks=c(35, 40, 45)) +
  scale_y_continuous(breaks=c(35, 40, 45)) 
  return(p)
}

# raw counts
# graw <- melt_hic(read_tsv("/Volumes/BFD/hic_heatmaps/gm12878/all/40k/chr18_chr18.tsv",
#    col_names = F))
# gnorm <- melt_hic(read_csv("/Volumes/BFD/hic_heatmaps/icmats/hsa/gm12878/gm12878-all-40k-chr18.csv",
#   col_names = F), norm=T)

## compare: h1
hraw <- melt_hic(read_tsv("/Volumes/BFD/hic_heatmaps/h1/all/40k/chr18_chr18.tsv",
  col_names = F))
hnorm <- melt_hic(read_csv("/Volumes/BFD/hic_heatmaps/icmats/hsa/h1/h1-all-40k-chr18.csv",
  col_names=F), norm=T)

## compare: k562
kraw <- melt_hic(read_tsv("/Volumes/BFD/hic_heatmaps/k562/all/40k/chr18_chr18.tsv",
  col_names = F))
knorm <- melt_hic(read_csv("/Volumes/BFD/hic_heatmaps/icmats/hsa/k562/k562-all-40k-chr18.csv",
  col_names = F), norm=T)

g1 <- plot_raw(graw)
g2 <- plot_hic(gnorm)

h1 <- plot_raw(hraw)
h2 <- plot_hic(hnorm)

k1 <- plot_raw(kraw)
k2 <- plot_hic(knorm)

svg("~/hvl/thesis_plots/normalisedhic.svg", 7, 8)
grid.arrange(h1, h2, k1, k2, ncol=2)
dev.off()
