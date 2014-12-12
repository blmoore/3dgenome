############# Fig. 3 Cross application ################
# Plot a figure showing: (a) the two best examples of #
# cross-applied models and (b) summary results of all #
# crosses as a bar chart. Optionally run via k-fold   #
# to get standard errors.                             #
#######################################################
library("randomForest")
library("reshape2")
library("blmR")

# load previously-built models:
g.mod <- readRDS("data/rds/Gm12878_35Vars_RFmod.rds")
h.mod <- readRDS("data/rds/H1hesc_35Vars_RFmod.rds")
k.mod <- readRDS("data/rds/K562_35Vars_RFmod.rds")
# and raw data:
g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

# Cross-apply models trained in once cell type to data
# from the other two
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

# Print cross-application results
cmat
# 100 * (cmat[1,1] - cmat[1,2]) / cmat[1,1] # 20.2% drop GM -> H1
# 100 * (cmat[1,1] - cmat[1,3]) / cmat[1,1] # 7.8% drop GM -> K5
# 100 * (cmat[2,2] - cmat[2,1]) / cmat[2,2] # 1.7% drop H1 -> GM
# 100 * (cmat[2,2] - cmat[2,3]) / cmat[2,2] # 21.0% drop H1 -> K5
# 100 * (cmat[3,3] - cmat[3,1]) / cmat[3,3] # 5.9% drop K5 -> GM
# 100 * (cmat[3,3] - cmat[3,2]) / cmat[3,3] # 20.4% drop K5 -> H1

xdat <- melt(t(cmat))
colnames(xdat) <- c("Applied to", "Trained on", "Corr")

xdat <- xdat[c(1:3, 5,4,6, 9,7,8),]
saveRDS(xdat, "data/rds/rf_xdat.rds")

pdf(6, 6, file="figures/f3b_crossApplyBars.pdf")
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

pdf(6, 6, file="figures/f3ai_gmXapplyK5.pdf")
plotPredRes(x=gk, y=k.dat$eigen, scale.factor=.6, ct="GM12878 to K562")
dev.off()
pdf(6, 6, file="figures/f3aii_k5XapplyGm.pdf")
plotPredRes(x=kg, y=g.dat$eigen, scale.factor=.7, ct="K562 to GM12878", col="red")
dev.off()

