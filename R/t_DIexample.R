## t_TADschematic.R

library("R.matlab")
library("RColorBrewer")

hm <- readMat("~/hvl/hiclib/hm/gm/heatmap.mat")$heatmap
hm <- readMat("~/hvl/hiclib/hm/h140kb/(13, 13).mat")$`(13, 13)`

region <- 675:775
pal <- colorRampPalette(c("white", "white", brewer.pal(8, "Reds")))(256)
#image(log10(as.matrix(hm[region,region])+1), col=pal)

di <- read.table("~/hvl/tads/chr14.di")

# match co-ords between heatmap and directionality index
start <- seq(from=0, by=4e4, length.out = nrow(hm))
bed <- data.frame(c="chr14", s=start, e=start+4e5)

# region of interest is:
s <- min(bed[region,]$s)
# chr13: 269.6 - 310 Mb
e <- max(bed[region,]$e)

di_region <- c(which(di$V2 == s):which(di$V3 == e))
di <- di[di_region,]

par(xaxs="i")

svg("~/hvl/thesis_plots/di_example.svg", 5, 7)
layout(matrix(c(1,1,1,2), ncol=1))
image(log10(as.matrix(hm[region,region])+1), col=pal, 
  axes=F, useRaster=F, asp=1)
plot(di$V4, type="l", axes=F, xlab="", ylab="")

dev.off()
