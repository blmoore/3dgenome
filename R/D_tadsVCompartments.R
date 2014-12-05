

# how do the tads look in this region?
# first: drawTadss from 2_fig1*.R

drawTads <- function(di.file, tad.file, dat.file, 
                     state.file=NULL, ytop=200, chr="chr2", 
                     range=0:1000, ctcf=NA, ...){
  # build composite plot of DI, tads, A/B
  all.di <- read.table(di.file)
  #all.di <- all.di[order(all.di[,1]),]
  c1.di <- all.di[all.di[,1] == gsub("chr","",chr),]
  colnames(c1.di) <- c("chr", "start", "end", "di")
  tads <- read.table(tad.file)
  t1 <- tads[tads$V1 == chr,]
  xcords <- apply(c1.di[,2:3], 1, mean) / 1e6
  
  d <- c1.di[,4]
  cols <- ifelse(d < -max(d)/50, scales::muted("red"),
                 ifelse(d > max(d)/50, scales::muted("green"), "grey40"))
  
  plot(xcords[range], d[range], type="n", xlab=paste0("Position along ", chr, " (Mb)"),
       ylab="Directionality index", frame=F, ylim=c(-ytop, ytop), axes=F, ...)
  axis(2, at=c(-(ytop*.9), 0, ytop*.9), cex.axis=.8,
       labels=c("Upstream", "0", "Downstream"))
  axis(1, at=seq(round(min(xcords)), round(max(xcords)), by=5), 
       labels=seq(round(min(xcords)), round(max(xcords)), by=5))
  abline(h=0, col="darkgrey")
  gradient.rect(xleft=t1$V2/1e6, ybottom=-ytop*1.2, ytop=ytop*1.2, xright=t1$V3/1e6,
                col=smoothColors(muted("green"), 10,
                                 "white", 10,
                                 muted("red"), alpha=90), border=NA)
  s <- seq(length(xcords)-1)
  segments(x0=xcords[s], x1=xcords[s+1], y0=d[s], y1=d[s+1], col=cols, lwd=1)
  if(!is.null(state.file)){
    #debugging
    sf <- read.table(state.file)
    points(t1[,2]/1e6, rep(30, nrow(t1)), pch="x", col="red")
    points(xcords[range], (sf[range,4]-2)*20)
  }  
  # now get A/B compartments for overlay:
  heigs <- read.table(dat.file)
  heigs <- heigs[heigs$V2 == chr,]
  heigs$V6 <- scale(heigs$V6) * (ytop/4)
  points(((heigs$V3 + heigs$V4)/2)/1e6, smooth(heigs$V6), type="l", col="blue")
  
  if(!is.na(ctcf)){
    allc <- read.table(ctcf)
    allc <- allc[allc[,1] == chr,]
    ycord <- replicate(nrow(allc), -(ytop*.9)+rnorm(1, ytop/40, ytop/40))
    points(((allc[,2] + allc[,3])/2)/1e6,  ycord, 
           pch=15, col=rgb(0,0,0,.3), cex=.4)
    rm(allc)
  }
}

library("plotrix")
library("scales")

drawTads("data/text/h1.di",
         "data/text/h1.tads",
         "data/bedfiles/h100.bed",
         ytop=400, main="H1 hESC", chr=csome, range=region)

nrow(read.table("data/text/h1.di"))

## 1) GM open 
csome = "chr5"
region = 3850:4150
##
par(mfrow=c(3,1), mar=c(2,5,2,2), xaxs="i")
drawTads("data/text/h1.di",
         "data/text/h1.tads",
         "data/bedfiles/h100.bed",
         ytop=400, main="H1 hESC", chr=csome, range=region)
drawTads("data/text/gm.di",
         "data/text/gm.tads",
         "data/bedfiles/g100.bed",
         ytop=30, main="GM12878", chr=csome, range=region)
drawTads("data/text/k5.di",
         "data/text/k5.tads",
         "data/bedfiles/k100.bed",
         ytop=30, main="K5 (new)", chr=csome, range=region)

## 2) H1 open 
csome = "chr20"
region = 450:650
##
par(mfrow=c(3,1), mar=c(2,5,2,2), xaxs="i")
drawTads("data/text/h1.di",
         "data/text/h1.tads",
         "data/bedfiles/h100.bed",
         ytop=400, main="H1 hESC", chr=csome, range=region)
drawTads("data/text/gm.di",
         "data/text/gm.tads",
         "data/bedfiles/g100.bed",
         ytop=30, main="GM12878", chr=csome, range=region)
drawTads("data/text/k5.di",
         "data/text/k5.tads",
         "data/bedfiles/k100.bed",
         ytop=30, main="K5 (new)", chr=csome, range=region)

## 3) K5 open 
csome = "chr8"
region = 2550:2850
##
par(mfrow=c(3,1), mar=c(2,5,2,2), xaxs="i")
drawTads("data/text/h1.di",
         "data/text/h1.tads",
         "data/bedfiles/h100.bed",
         ytop=400, main="H1 hESC", chr=csome, range=region)
drawTads("data/text/gm.di",
         "data/text/gm.tads",
         "data/bedfiles/g100.bed",
         ytop=30, main="GM12878", chr=csome, range=region)
drawTads("data/text/k5.di",
         "data/text/k5.tads",
         "data/bedfiles/k100.bed",
         ytop=30, main="K5 (new)", chr=csome, range=region)


## example, gm12878
g.open <- read.table("data/bedfiles/g.open.bed")
g.tads <- read.table("data/bedfiles/gm_tadbounds.bed")

none <- read.table("data/bedfiles/none.bed")
g.closed <- read.table("data/bedfiles/g.closed.bed")

# g.open / closed may have separate 1Mb bins for e.g. 3 Mb region

mergeAdjacent <- function(bedlike){
  # empty dataframe
  df <- bedlike[0,]
  
  mergeDf <- function(d){
    # d must be bed -like (i.e. chr, start, end, id, ...)
    i <- 1
    while(i < nrow(d)){
      # if region > 1 Mb, 
      if(d[i+1, 2] == d[i,3]){
        row <- d[i,]
        row[,3] <- d[i+1, 3]
        i <- i+2
      } else {
        row <- d[i,]
        i <- i+1
      }
      df <- rbind(df, row)
    }
    # don't forget last row
    df <- rbind(df, d[nrow(d),])
    df
  }
  # recursively merge:
  v1 <- mergeDf(bedlike)
  while(nrow(v2) != nrow(v1)){
    v2 <- v1
    v1 <- mergeDf(v2) 
  } 
  v1
}

gom <- mergeAdjacent(g.open)
gcm <- mergeAdjacent(g.closed)
nm <- mergeAdjacent(none)

nearest <- function(n, tads)
  min(abs(c(tads[,2], tads[,3]) - n))

tadDistances <- function(bed, tads=g.tads){
  out <- vector(mode="numeric", length=0)
  for(c in paste0("chr", 1:22)){
    flip <- bed[bed[,1] == c,]
    tad <- tads[tads[,1] == c,]
    s <- sapply(flip[,2], nearest, tads=tads)
    e <- sapply(flip[,3], nearest, tads=tads)
    out <- c(out, s, e)
  }
  message(head(out))
  do.call(rbind, as.list(out))
}

open.tads <- tadDistances(gom)
closed.tads <- tadDistances(gcm)
none.tads <- tadDistances(nm)

open.tads

dev.off()
plot(ecdf(open.tads), col="darkgreen", log="x", xlim=c(1e4, 1e7), verticals=T)
lines(ecdf(none.tads))
lines(ecdf(closed.tads), col="darkred")

mean(go2)
mean(go.none)
mean(go.closed)
median(go2)
median(go.none)
