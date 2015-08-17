### Test repeat classes over boundaries
# RepeatMasked tracks from UCSC table browser 30 June
library("dplyr")
library("ggplot2")
library("readr")
library("reshape2")
library("blmR")

setwd("~/hvl/repeats/")

hrep <- read_tsv("~/hvl/repeats/hg19_ucscRepeatMasker.txt")
mrep <- read_tsv("~/hvl/repeats/mm10_ucscRepeatMasker.txt")

get_classes <- function(df, col){
  # return levels of repClass or repFamily
  out <- unique(df[[col]])
  out <- out[!grepl("\\?", out)]
  out
}
  
split_bed <- function(colvar = c("repFamily", "repClass")){
  colvar <- match.arg(colvar)
  
  per_var <- function(rep_class, df=df, genome=genome){
    #df <- subset(df, col_var == rep_class)
    df <- df[df[[colvar]] == rep_class,]
    # don't want unknowns with ?
    stopifnot(!grepl("\\?", df[[colvar]]))
    bed <- df[,c("genoName", "genoStart", "genoEnd", "repName")]
    outname <- paste0("~/hvl/repeats/", genome, "_", colvar, "_", rep_class, ".bed")
    message("Writing to ", outname)
    invisible(write_bed(bed, fname = outname))
  }
  
  classes <- get_classes(hrep, colvar)
    
  sapply(classes, per_var, df=hrep, genome="hg19")
  message("hg19 done")
  
  unique(mrep$repClass)[(which(!unique(mrep$repClass) %in% classes))]
  
  mclass <- get_classes(mrep, colvar)
    
  # repeat classes are the same in humand and mouse data (families not?)
  # stopifnot(all(classes %in% mclass) & all(mclass %in% classes))
  sapply(mclass, per_var, df=mrep, genome="mm10")
  
  message("mm10 done")
}


#split_bed("repFamily")
#split_bed("repClass")

## -------------------- sh intersect.sh ----------------------------- ##
options(warn=2)

buildBoundariesDat <- function(ct=c("h1", "k5", "gm"), type=c("c", "t"),
                              col=c("repFamily", "repClass"), summary=T, steps=nsteps){

  ct <- match.arg(ct)
  type <- match.arg(type)
  col <- match.arg(col)
  clab <- if(ct == "h1") "H1 hESC" else if(ct == "gm") "Gm12878" else "K562"
  
  if(type == "c"){
    # compartments
    pattern <- paste0("^", ct, "cb_.*?", col, ".*\\.bed\\.out")
    type <- "Compartments"
    #steps <<- 39 # 1 Mb, 40 steps
    #steps <- 19
  } else {
    # TADs
    pattern <- paste0("^", ct, "tb_.*?", col, ".*\\.bed\\.out")
    type <- "TADs"
    #steps <<- 39 # 1 Mb, 40 steps
    # steps <- 19 # 500 kb, 20 steps
  }
  
  fl <- list.files("~/hvl/repeats/", pattern=pattern, full.names = T)
  looper <- fl  
  print(col)
  vars <- get_classes(hrep, col)

  # get vars and filelist into same ordering
  vars <- sort(vars)
  m <- fl[order(gsub("\\..*", "", fl))]
  stopifnot(length(m) == length(vars))
  
  gdf <- data.frame(name=character(), 
    pos=numeric(),
    mean=numeric(), 
    mean.l=numeric(), 
    mean.u=numeric())
  bins <- matrix(nrow=0, ncol=steps)
  if(summary){
    for(f in seq_along(m)){
      message(m[f])
      feat <- read.table(m[f])  
      ac <- matrix(feat$V5, ncol=steps, byrow=T)
      ac[is.nan(ac)] <- 0
      
      means <- colMeans(ac)
      ## NB Confidence intervals ? *1.96, Std error? as-is
      se <- apply(ac, 2, plotrix::std.error)
      outr <- data.frame(name=vars[f], pos=1:steps, 
        mean=means, means.l=means-se, means.u=means+se)
      gdf <- rbind(gdf, outr)
    }
    gdf$feat <- rep(gsub("_", " ", vars), each=steps)
    gdf$ct <- clab
    gdf$type <- type
    return(gdf)
  } else {
    for(i in seq_along(m)){
      f <- m[i]
      message(f)
      var <- vars[i]
      feat <- read.table(f)  
      ac <- matrix(feat$V5, ncol=steps, byrow=T)
      ac[is.nan(ac)] <- 0
      outr <- data.frame(name=f, bound=1:nrow(ac), 
        feat=rep(gsub("_", " ", var), nrow(ac)), ac)
      gdf <- rbind(gdf, outr)
    }
    gdf$ct <- clab
    gdf$type <- type
    return(gdf)
  }
}

blob <- function(each){
  #nsteps == ncol(each) - x (three? four?)
  library("lazyeval")
  edges <- paste0("X", c(1:3, nsteps-2, nsteps-1, nsteps))
  mid <- paste0("X", round(nsteps/2))
  
  # use dplyr with standard evaluation via lazyeval
  return(
    each %>% group_by(ct, type, feat) %>%
    summarise_(
      bound.mean = interp(~mean(mid), mid=as.name(mid)),
      edge.mean  = interp(~mean(edges), edges=as.name(edges)),
      p          = interp(~wilcox.test(mid, edges)$p.value, mid=as.name(mid), edges=as.name(edges)))
  )
}
  
# nsteps: number of steps in binAroundBed.py -1 !
nsteps <- 41

fam_sum <- Map(buildBoundariesDat, 
  c("h1", "gm", "k5"), rep(c("c", "t"), 3), 
  col="repFamily", summary=T, steps=nsteps)
fam_sum <- do.call(rbind, fam_sum)
fam_each <- Map(buildBoundariesDat, 
  c("h1", "gm", "k5"), rep(c("c", "t"), 3), 
  col="repFamily", summary=F, steps=nsteps)
fam_each <- do.call(rbind, fam_each)

class_sum <- Map(buildBoundariesDat, 
  c("h1", "gm", "k5"), rep(c("c", "t"), 3), 
  col="repClass", summary=T, steps=nsteps)
class_sum <- do.call(rbind, class_sum)
class_each <- Map(buildBoundariesDat, 
  c("h1", "gm", "k5"), rep(c("c", "t"), 3), 
  col="repClass", summary=F, steps=nsteps)
class_each <- do.call(rbind, class_each)

saveRDS(fam_each, "~/hvl/production/3dgenome/data/rds/rpt_fam.rds")
saveRDS(class_each, "~/hvl/production/3dgenome/data/rds/rpt_class.rds")

head(class_each)
unique(fam_each$feat)

tad_sine <- subset(fam_each, feat == "Alu" & type == "TADs")
head(tad_sine)

# code matches F_predictTADbounds.R
ts <- melt(tad_sine, id.vars = c("name", "bound", "feat", "ct", "type"))
ts <- subset(ts, variable %in% c("X1", "X12"))
ts$variable <- ifelse(ts$variable == "X12", T, F)
ts$name <- NULL

rptdf <- dcast(ts, bound + type + ct + variable ~ feat, value.var="value")
ggplot(rptdf, aes(x=variable, y=Alu)) + geom_violin()
head(rptdf)

  
## 1) profiles:
pdf("plots/repFamily_1Mb.pdf", 6, 28)
ggplot(fam_sum, aes(x=pos, y=mean, col=type, fill=type)) + 
  facet_grid(feat~ct, scales="free") + 
  geom_vline(xintercept=ceiling(nsteps/2), lwd=4, col=I(rgb(.8,.8,.8, .5))) +
  geom_ribbon(aes(ymin=means.l, ymax=means.u), alpha=I(.3), col=I(NULL)) +
  geom_line() + theme_bw() +
  scale_color_brewer(type = "qual", palette = 2) + 
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_x_continuous(breaks=c(0, ceiling(nsteps/2), nsteps), 
    labels = c("-1 Mb", "Boundary", "+1 Mb")) +
  scale_y_continuous(breaks=scales::pretty_breaks(3)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=I(8)),
    panel.grid=element_blank(), strip.text.y=element_text(size=I(6))) +
  labs(x="", y="Mean number of annotated repeats per 50 kb", col="", fill="")
dev.off()

pdf("plots/repClass_1Mb.pdf", 6, 14)
ggplot(class_sum, aes(x=pos, y=mean, col=type, fill=type)) + 
  facet_grid(feat~ct, scales="free") + 
  geom_vline(xintercept=ceiling(nsteps/2), lwd=4, col=I(rgb(.8,.8,.8, .5))) +
  geom_ribbon(aes(ymin=means.l, ymax=means.u), alpha=I(.3), col=I(NULL)) +
  geom_line() + theme_bw() +
  scale_color_brewer(type = "qual", palette = 2) + 
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_x_continuous(breaks=c(0, ceiling(nsteps/2), nsteps), 
    labels = c("-1 Mb", "Boundary", "+1 Mb")) +
  scale_y_continuous(breaks=scales::pretty_breaks(3)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=I(8)),
    panel.grid=element_blank(), strip.text.y=element_text(size=I(7))) +
  labs(x="", y="Mean number of annotated repeats per 50 kb", col="", fill="")
dev.off()

options(warn=1)
## 2) bubble significance plots
pdf("plots/repClass_bubble_1Mb.pdf", 6, 5)
ggplot(blob(class_each), aes(x=feat, y=-log10(p), col=type,
  size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~.) + 
  geom_point(position=position_dodge(.15)) + theme_bw() +
  labs(list(size="Absolute difference\nat boundary",
    y=expression(-log[10](italic(p))),
    x="", col="Boundary type")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(range=c(2,7)) +
  guides(col=guide_legend(override.aes=list(size=4))) +
  geom_hline(yintercept=-log10(.05 / 264), linetype="dashed") +
  scale_colour_brewer(type="qual", palette=6) +
  ggtitle("Repeat classes over boundaries")
dev.off()

pdf("plots/repFamily_bubble_1Mb.pdf", 12, 5)
ggplot(blob(fam_each), aes(x=feat, y=-log10(p), col=type,
  size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~.) + 
  geom_point(position=position_dodge(.15)) + theme_bw() +
  labs(list(size="Absolute difference\nat boundary",
    y=expression(-log[10](italic(p))),
    x="", col="Boundary type")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(range=c(2,7)) +
  guides(col=guide_legend(override.aes=list(size=4))) +
  geom_hline(yintercept=-log10(.05 / 264), linetype="dashed") +
  scale_colour_brewer(type="qual", palette=6) +
  ggtitle("Repeat families over boundaries")
dev.off()


## ---- bp per bin ---- ##
  