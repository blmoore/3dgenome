############# Fig. 3 Cross application ################
# Plot a figure showing: (a) the two best examples of #
# cross-applied models and (b) summary results of all #
# crosses as a bar chart. Optionally run via k-fold   #
# to get standard errors.                             #
#######################################################
library("calibrate")
library("ggplot2")
library("gridExtra")
library("infotheo")
library("randomForest")
library("RColorBrewer")
library("reshape2")
devtools::load_all("blmR")
source("~/blmR/R/misc_functions.R")

hm2 <- modelEigens.all(h.dat)
plotPredRes.ice(x=hm2$predicted, y=h.dat$eigen)

h2 <- h.mod
k2 <- k.mod
g2 <- g.mod

hh <- h.mod$predicted
hk <- predict(h.mod, k.dat[,-1])
hg <- predict(h.mod, g.dat[,-1])

kk <- k.mod$predicted
kh <- predict(k.mod, h.dat[,-1])
kg <- predict(k.mod, g.dat[,-1])

gg <- g.mod$predicted
gh <- predict(g.mod, h.dat[,-1])
gk <- predict(g.mod, k.dat[,-1])

## Set up cor mat, 1 = g, 2 = h, 3 = k. Rows == mod, cols == applied
cmat <- matrix(c(cor(gg, g.dat$eigen), 
                 cor(gh, h.dat$eigen),
                 cor(gk, k.dat$eigen), ## row 1 GM ->
                 
                 cor(hg, g.dat$eigen),
                 cor(hh, h.dat$eigen),
                 cor(hk, k.dat$eigen), ## row 2 H1 ->
                 
                 cor(kg, g.dat$eigen),
                 cor(kh, h.dat$eigen),
                 cor(kk, k.dat$eigen)), # row 3: k562 ->
               byrow=T, nrow=3)
rownames(cmat) <- colnames(cmat) <- c("GM12878", "H1 hESC", "K562")

plotPredRes.homer.custom <- function (modelOut=NA, x=NA, y=NA, col="blue", ct="H1"){
  if(!is.na(x) & !is.na(y)){
    d1 <- x
    d2 <- y
    modelOut$preds <- cbind(x, y)
  } else {
    d1 <- modelOut$preds[,1]
    d2 <- modelOut$preds[,2]
  }
  colour  <- "#0000ff42"
  if ( col != "blue" ) {
    ifelse(col == "red", colour <- "#ff000042", colour  <- "#FFA50062")
  }
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  xhist <- hist(d1, plot=FALSE)
  yhist <- hist(d2, breaks=100, plot=F)
  top <- max(c(xhist$counts, yhist$counts))
  nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(6,1), c(1,6), TRUE)
  par(mar=c(5,5,0,0), mgp=c(1.8,.4,0))
  max <- max(abs(c(d1, d2))) * 1.05
  plot(d1, d2, xlab="Predicted eig", ylab="Empirical eig", col=colour, 
       pch=16, type='n', xlim=c(-max*.9,max*.9), ylim=c(-max/2,max/2))
  abline(h=0,v=0, lty=2)
  points(d1, d2, col=colour, pch=16)
  textlabs <- c(-max/1.5, max*.9, # ct x, y --TOP
                -max/1.5, max*.8, # pcc x, y
                -max/1.5, max*.7, # RMSE
                max/1.5, -max*.9, # AUROC --bottom
                max/1.5, -max*.8) # Acc.
  text(textlabs[1], textlabs[2], ct, font=2, cex=1.2)
  text(textlabs[3], textlabs[4], paste("PCC = ", signif(cor(d1,d2), 3), sep=""), col="navy")
  text(textlabs[5], textlabs[6], paste("RMSE = ", signif(rmse(d1,d2), 3), sep = ""), col="navy")
  text(textlabs[9], textlabs[10], paste("Acc. = ",signif(100*(sum(apply(modelOut$preds,1,
                                                                        function(x) if(all(x > 0) || all(x<0)) 1 else 0))
                                                              /nrow(modelOut$preds)),4), sep=""), col="darkgreen")
  text(textlabs[7], textlabs[8], paste("AUROC = ", signif(getAUC.gen(modelOut),3), sep=""), col="darkgreen") 
 
  par(mar=c(0,5,1,1)) # RHS
  plot(density(d1, from=-max, to=max), 
       frame=F, axes=F, type="n", col=colour,
       xlab="", ylab="", main="", xlim=c(-max*.9,max*.9))
  polygon(density(d1,from=-max, to=max), 
          col=colour, border=colour, lwd=4) 
  par(mar=c(5,0,1,1))
  plot(density(d2,from=-max, to=max)$y, density(d2,from=-max, to=max)$x, 
       frame=F, axes=F, type="n", col=colour,
       xlab="", ylab="", main="", ylim=c(-max,max))
  polygon(density(d2, from=-max, to=max)$y, density(d2, from=-max/2, to=max/2)$x, 
          col=colour, border=colour, lwd=4)
  par(def.par)
}

cmat
100 * (cmat[1,1] - cmat[1,2]) / cmat[1,1] # 20.2% drop GM -> H1
100 * (cmat[1,1] - cmat[1,3]) / cmat[1,1] # 7.8% drop GM -> K5
100 * (cmat[2,2] - cmat[2,1]) / cmat[2,2] # 1.7% drop H1 -> GM
100 * (cmat[2,2] - cmat[2,3]) / cmat[2,2] # 21.0% drop H1 -> K5
100 * (cmat[3,3] - cmat[3,1]) / cmat[3,3] # 5.9% drop K5 -> GM
100 * (cmat[3,3] - cmat[3,2]) / cmat[3,3] # 20.4% drop K5 -> H1

xdat <- melt(t(cmat))
colnames(xdat) <- c("Applied to", "Trained on", "Corr")

xdat <- xdat[c(1:3, 5,4,6, 9,7,8),]

pdf(6, 6, file="~/hvl/ice/plots/f3b_crossApplyBars.pdf")
par(mgp=c(1.5,.2,0))
plot(1:13, seq(0,1,length=13), type="n", axes=F, xlab="", ylab="PCC", ylim=c(0,.85))
barpos <- c(2:4,6:8,10:12)
xl <- barpos-.4
xr <- barpos+.4

#         blue, gm     orange, h1    red, k562
pal <- c("#0000ff98", "#FFA50098", "#ff000098")
topcols <- c(pal[1], pal[2], pal[3],
             pal[2], pal[1], pal[3],
             pal[3], pal[1], pal[2])

for(x in 1:nrow(xdat)){
  rect(xl[x], 0, xr[x], xdat[x,3], col=topcols[x])
}

axis(2, cex.axis=.7, tck=.02)
axis(1, at=c(1,3,7,11,13), labels=c("","GM12878", "H1 hESC", "K562",""),
     tick=F, cex.axis=.8, font=2)
mtext("Trained in cell type",1, line=2)
legend("topright", fill=c(pal[1], pal[2], pal[3]), xpd=T, bty="n",
       legend=c("GM12878", "H1 hESC", "K562"), cex=.7, title="Applied to cell type")
dev.off()

svg(6, 6, file="~/hvl/ice/plots/inkscape/gmToK5_s.svg")
plotPredRes.ice(x=gk, y=k.dat$eigen, scale.factor=.6, ct="GM12878 → K562")
dev.off()
svg(6, 6, file="~/hvl/ice/plots/inkscape/k5ToGm_s.svg")
plotPredRes.ice(x=kg, y=g.dat$eigen, scale.factor=.7, ct="K562 → GM12878", col="red")
dev.off()

