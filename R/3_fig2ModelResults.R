############## Fig. 2 Model results #####################
# Plot a figure showing: (a) a scatterplot of predicted #
# vs. empirical values for Mb eignevectors along with   #
# various metrics to evaluate model performance; (b) a  #
# bar chart ranking variable importance (with estimate  #
# of variability estimated by K-fold cross validation.  #
#########################################################
library("randomForest")
library("blmR")

# We've already calculated + plotted OOB prediction results 
# in the script: 0_buildDatfiles.R, first look at OOB var imp
g.mod <- readRDS("data/rds/Gm12878_35Vars_RFmod.rds")
h.mod <- readRDS("data/rds/H1hesc_35Vars_RFmod.rds")
k.mod <- readRDS("data/rds/K562_35Vars_RFmod.rds")

## Fig 2a, model predictive accuracy per cell type. Note 
## these can't be combined on plot device with mfrow etc.
## due to the complicated layout() used to add side densities.
pdf("figures/f2_gmRes.pdf", 6, 6)
plotPredRes(x=g.mod$predicted, y=g.dat$eigen,
                col="blue", ct="GM12878", scale.factor=.8)
dev.off()

pdf("figures/f2_h1Res.pdf", 6, 6)
plotPredRes(x=h.mod$predicted, y=h.dat$eigen,
                col="orange", ct="H1 hESC", scale.factor=.7)
dev.off()

pdf("figures/f2_k5Res.pdf", 6, 6)
plotPredRes(x=k.mod$predicted, y=k.dat$eigen,
                col="red", ct="K562", scale.factor=.6)
dev.off()

impBars <- function(mod, ...){
  #imp <- mod$importance[,1]
  imp <- importance(mod, type=1)
  rownames(imp) <-  c("ATF3", "CEBP", "CHD1", "CHD2", 
                   "MYC", "CTCF", "DNase", "EGR1", 
                   "EZH2", "GABP", "JUND", "MAX",
                   "H2A.Z", "H3K27ac", "H3K27me3", 
                   "H3K36me3", "H3K4me1", "H3K4me2", 
                   "H3K4me3", "H3K79me2", "H3K9ac", 
                   "H3K9me3", "H4K20me1", "MXI1", 
                   "NRSF", "P300", "POL2", "RAD21", 
                   "SIX5", "SP1", "TAF1", "TBP", 
                   "YY1", "ZNF143", "GC")
  imp <- imp[order(imp[,1], decreasing=T),]
  barplot(rev(imp[1:10]), horiz=T, las=1, cex.names=1.5, ...)
}

pdf("figures/f2b_varImpPerModel_v2.pdf", 9, 3.5)
par(mfrow=c(1,3), mar=c(3,7,1.5,0.5), oma=c(2,0,0,0), mgp=c(0,.5,0))
impBars(g.mod, main="GM12878", col="#0000ff92", border=NA)
impBars(h.mod, main="H1 hESC", col="#FFA50092", border=NA)  
impBars(k.mod, main="K562", col="#ff000092", border=NA)
mtext(1, outer=T, text="Variable importance (% increase in MSE when permuted)")
dev.off()

## Combine A + B in inkscape to get figure 2