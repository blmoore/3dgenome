library("randomForest")
library("RRF")

# load data, 35 + full datasets
hdat <- readRDS("data/rds/H1hesc_35Vars.rds")
hfull <- readRDS("data/rds/H1hesc_71Vars.rds")
gdat <- readRDS("data/rds/Gm12878_35Vars.rds")
gfull <- readRDS("data/rds/Gm12878_115Vars.rds")
kdat <- readRDS("data/rds/K562_35Vars.rds")
kfull <- readRDS("data/rds/K562_187Vars.rds")

h <- rrfcv(hdat[,-1], trainy=hdat[,"eigen"], cv.fold=5)
g <- rrfcv(gdat[,-1], trainy=gdat[,"eigen"], cv.fold = 5)
k <- rrfcv(kdat[,-1], trainy=kdat[,"eigen"], cv.fold=5)

hf <- rrfcv(hfull[,-1], trainy=hfull$eigen, cv.fold=5)
gf <- rrfcv(gfull[,-1], trainy=gfull$eigen, cv.fold = 5)
kf <- rrfcv(kfull[,-1], trainy=kfull$eigen, cv.fold = 5)

rrf_scree <- function(rrf, ...)
  with(rrf, plot(n.var, error.cv, log="x", type="o", lwd=2, ...))

par(mfrow=c(2,3))
rrf_scree(g, main="GM12878")
rrf_scree(h, main="H1 hESC")
rrf_scree(k, main="K562")

rrf_scree(hf, main="H1 hESC (71 vars)")
rrf_scree(gf, main="GM12878 (115 vars)")
rrf_scree(kf, main="K562 (187 vars)")



str(h)
dev.off()

rf <- RRF(hdat[,colnames(hdat) != "eigen"], hdat$eigen, flagReg = 0, importance=T)
rf$importance

imrf <- importance(rf, type=2)
imrf/max(imrf)


hrrf <- RRF(hdat[,-1], hdat[,1], flagReg=0, importance=T)
plot(hrrf)
hrrf$feaSet

class(hrrf) <- "randomForest"
top_ten(hrrf)

randomForest(hdat[,colnames(hdat) != "eigen"], hdat$eigen)


# "Guided" regularised random forest?
# ad-hoc implementation here: http://stats.stackexchange.com/questions/152162/how-to-sign-and-score-risk-factors-in-a-guided-regularized-random-forest
# see desktop PDF presentation for details

  
# original (old) hdat ::
oldh <- readRDS("~/Downloads/H1hesc_35Vars.rds")
orf <- randomForest(y=oldh$eigen, oldh[,-1], importance=T)
top_ten(orf)

all(oldh == hdat)

plot(oldh$eigen, hdat$eigen)

importance(hrrf)
