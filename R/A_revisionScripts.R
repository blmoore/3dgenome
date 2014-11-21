library("DAAG")
library("plsdepot")
library("RColorBrewer")
library("reshape2")
library("blmR")

setwd("~/hvl/production/3dgenome/")
g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

### Linear regression::

dat <- list(g.dat, h.dat, k.dat)
labs <- c("GM12878", "H1 hESC", "K562")
res <- c()

for(i in 1:3){
  others <- setdiff(1:3, i)
  .env <- environment()
  form <- formula(eigen~., env=.env)
  clm <- cv.lm(df=dat[[i]], form.lm=form, m=10, seed=42, 
               plotit=F)
  self <- cor(clm$Predicted, dat[[i]]$eigen)
  
  model <- lm(form, data=dat[[i]])  
  x1 <- cor(dat[[others[1]]]$eigen, predict(model, dat[[others[1]]][,-1]))
  x2 <- cor(dat[[others[2]]]$eigen, predict(model, dat[[others[2]]][,-1]))
  
  out <- c(self, x1, x2)
  names(out) <- c(paste0(labs[i], ".", labs[i]), 
                  paste0(labs[i], ".", labs[others[1]]),
                  paste0(labs[i], ".", labs[others[2]]))
  res <- c(res, out)
}


## Some -ve correlations (?) e.g.::
# hlm <- lm(eigen~., data=h.dat)
# plotPredRes(x=k.dat$eigen, y=predict(hlm, k.dat[,-1]))

lmdf <- data.frame(do.call(rbind, strsplit(names(res), "\\.")), res, row.names=NULL)
colnames(lmdf)[1:2] <- c("trained", "applied")
lm.hm <- dcast(lmdf, trained~applied)

### PLS regression ::

# choose number of PCs via) cross-validation
g.pls <- plsreg1(g.dat[,-1], g.dat$eigen, crosval=T, comps=NULL)
h.pls <- plsreg1(h.dat[,-1], h.dat$eigen, crosval=T, comps=NULL)
k.pls <- plsreg1(k.dat[,-1], k.dat$eigen, crosval=T, comps=NULL)

plot(k.pls)
plot(h.pls)
plot(g.pls)

g.pls$R2
g.pls$Q2

plot(g.pls, what="observations", comps=c(1,3))

# labelled reponses
plot(g.dat$eigen, g.pls$y.pred, type = "n", xlab = "Original",
     ylab = "Predicted")
title("Comparison of responses", cex.main = 0.9)
abline(a = 0, b = 1, col = "gray85", lwd = 2)
text(g.dat$eigen, g.pls$y.pred, col = "#5592e3")

plotPredRes(y=g.dat$eigen, x=g.pls$y.pred)

predict(g.pls, h.dat[,-1])

## apply model to new data
pls.reg <- function(dat, mod)
  cor(dat[,1], as.matrix(dat[,-1]) %*% mod$reg.coefs[-1] + mod$reg.coefs[1])

gh.pls <- pls.reg(h.dat, g.pls)
gk.pls <- pls.reg(k.dat, g.pls)
gg.pls <- cor(g.dat$eigen, g.pls$y.pred)

hg.pls <- pls.reg(g.dat, h.pls)
hk.pls <- pls.reg(k.dat, h.pls)
hh.pls <- cor(h.dat$eigen, h.pls$y.pred)

kg.pls <- pls.reg(g.dat, k.pls)
kh.pls <- pls.reg(h.dat, k.pls)
kk.pls <- cor(k.dat$eigen, k.pls$y.pred)

# trained, applied, cor
pls.df <- data.frame(rbind(
                 c(t="GM12878", a="GM12878", cor=gg.pls),
                 c(t="GM12878", a="H1 hESC", cor=gh.pls),
                 c(t="GM12878", a="K562",    cor=gk.pls),
                 c(t="H1 hESC", a="H1 hESC", cor=hh.pls),
                 c(t="H1 hESC", a="GM12878", cor=hg.pls),
                 c(t="H1 hESC", a="K562",    cor=hk.pls),
                 c(t="K562",    a="K562",    cor=kk.pls),
                 c(t="K562",    a="H1 hESC", cor=kh.pls),
                 c(t="K562",    a="GM12878", cor=kg.pls)))

pls.df$cor <- as.numeric(as.character(pls.df$cor))
str(pls.df)

### Plotting ::

p <- colorRampPalette(brewer.pal(9, "RdBu"))(256)

pdf("figures/suppl/method_compare.pdf", 9, 3)
par(mar=c(4,4,2,2), mfrow=c(1,3), mgp=c(1,0,0))

image(x=1:3, y=1:3, z=as.matrix(lm.hm[,-1]), col=p, zlim=c(-1,1), 
      asp=1, frame=F, axes=F, xlab="Training cell type", ylab="Test cell type", 
      main="Linear regression")
axis(1, at=1:3, labels=colnames(lm.hm)[-1], line=NA, tick=F)
axis(2, at=1:3, labels=colnames(lm.hm)[-1], line=NA, tick=F)


## RandomForest (xdat from 4_fig3CrossApplication.R)
#xdat
rf.hm <- dcast( `Trained on` ~ `Applied to` , data=xdat)
image(x=1:3, y=1:3, z=as.matrix(rf.hm[,-1]), col=p, zlim=c(-1,1), 
      asp=1, frame=F, axes=F, xlab="Training cell type", ylab="Test cell type", 
      main="Random Forest")
axis(1, at=1:3, labels=colnames(lm.hm)[-1], line=NA, tick=F)
axis(2, at=1:3, labels=colnames(lm.hm)[-1], line=NA, tick=F)

## PLS
# pls.df
pls.hm <- dcast(t~a, data=pls.df)
image(x=1:3, y=1:3, z=as.matrix(pls.hm[,-1]), col=p, zlim=c(-1,1), 
      asp=1, frame=F, axes=F, xlab="Training cell type", ylab="Test cell type", 
      main="PLS regression")
axis(1, at=1:3, labels=colnames(lm.hm)[-1], line=NA, tick=F)
axis(2, at=1:3, labels=colnames(lm.hm)[-1], line=NA, tick=F)

dev.off()

## Cell type performance means:
mean.diag <- function(df)
  mean(diag(as.matrix(df[,-1])))

sapply(list(lm.hm, rf.hm, pls.hm), mean.diag)
# [1] 0.787 0.790 0.750

mean.nondiag <- function(df){
  mat <- as.matrix(df[,-1])
  diag(mat) <- NA
  mean(mat, na.rm=T)
}

sapply(list(lm.hm, rf.hm, pls.hm), mean.nondiag)
# [1] 0.139 0.689 0.641

pdf("figures/suppl/clbar.pdf", 4, 2)
image.scale(z=as.matrix(lm.hm[,-1]), zlim=c(-1,1), col=p)
dev.off()

image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  ## from :: http://menugget.blogspot.co.uk/2011/08/adding-scale-to-image-plot.html#more
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}



