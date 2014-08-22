## build a general model

devtools::load_all("blmR")

all.dat <- rbind(transform(g.dat, ct="gm"), 
                 transform(h.dat, ct="h1"), 
                 transform(k.dat, ct="k5"))


all.mod <- modelEigens.all(all.dat)
saveRDS(all.mod, file="~/hvl/ice/rds/allmod.rds")

plotPredRes.ice(y=all.dat$eigen, x=all.mod$predicted, scale.factor=.6)

par(mar=c(5,8,3,3))
i <-  importance(all.mod, type=1)
i <- i[order(i[,1], decreasing=F),]
barplot(i, beside=T, horiz=T, las=1)
