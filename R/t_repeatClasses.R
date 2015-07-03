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
  
  classes <- get_classes(hrep[[colvar]])
    
  sapply(classes, per_var, df=hrep, genome="hg19")
  message("hg19 done")
  
  unique(mrep$repClass)[(which(!unique(mrep$repClass) %in% classes))]
  
  mclass <- get_classes(mrep[[colvar]])
    
  # repeat classes are the same in humand and mouse data (families not?)
  # stopifnot(all(classes %in% mclass) & all(mclass %in% classes))
  sapply(mclass, per_var, df=mrep, genome="mm10")
  
  message("mm10 done")
}


split_bed("repFamily")
split_bed("repClass")

## -------------------- sh intersect.sh ----------------------------- ##

buildBoundariesDat <- function(ct=c("h1", "k5", "gm"), type=c("c", "t"),
                              col=c("repFamily", "repClass"), summary=T){

  ct <- match.arg(ct)
  type <- match.arg(type)
  col <- match.arg(col)

  if(type == "c"){
    # compartments
    pattern <- paste0("^", ct, "cb_.*?", col, ".*\\.bed\\.out")
    #steps <- 30
    steps <- 19
  } else {
    # TADs
    pattern <- paste0("^", ct, "tb_.*?", col, ".*\\.bed\\.out")
    steps <- 19
  }
  
  fl <- list.files("~/hvl/repeats/", pattern=pattern, full.names = T)
  looper <- fl  
  print(col)
  vars <- get_classes(hrep, col)
  
  matches <- lapply(vars, function(x) grep(x, fl)[1])
  m <- fl[unlist(matches)]
  m <- na.omit(m)
  m <- unique(m)
  
  gdf <- data.frame(name=character(), 
    pos=numeric(),
    mean=numeric(), 
    mean.l=numeric(), 
    mean.u=numeric())
  bins <- matrix(nrow=0, ncol=steps)
  if(summary){
    for(f in seq_along(looper)){
      message(looper[f])
      feat <- read.table(looper[f])  
      ac <- matrix(feat$V5, ncol=steps, byrow=T)
      ac[is.nan(ac)] <- 0
      
      means <- colMeans(ac)
      ## NB Confidence intervals ? *1.96, Std error? as-is
      se <- apply(ac, 2, plotrix::std.error)
      outr <- data.frame(name=vars[f], pos=1:steps, 
        mean=means, means.l=means-se, means.u=means+se)
      gdf <- rbind(gdf, outr)
    }
    gdf$feat <-  rep(vars, each=steps)
    gdf$ct <- ct
    gdf$type <- type
    return(gdf)
  } else {
    for(i in seq_along(looper)){
      f <- looper[i]
      message(f)
      var <- vars[i]
      feat <- read.table(f)  
      ac <- matrix(feat$V5, ncol=steps, byrow=T)
      ac[is.nan(ac)] <- 0
      outr <- data.frame(name=f, bound=1:nrow(ac), 
        feat=rep(var, nrow(ac)), ac)
      gdf <- rbind(gdf, outr)
    }
    gdf$ct <- ct
    gdf$type <- type
    return(gdf)
  }
}

blob <- function(each){
  each %>% group_by(ct, type, feat) %>% 
    summarise(bound.mean = mean(X10),
      edge.mean  = mean(c(X1, X2, X3, X17, X18, X19)),
      p = wilcox.test(X10, c(X1, X2, X3,  X17, X18, X19))$p.value)
}

fam_sum <- Map(buildBoundariesDat, c("h1", "gm", "k5"), rep(c("c", "t"), 3), col="repFamily", summary=T)
fam_sum <- do.call(rbind, fam_sum)
fam_each <- Map(buildBoundariesDat, c("h1", "gm", "k5"), rep(c("c", "t"), 3), col="repFamily", summary=F)
fam_each <- do.call(rbind, fam_each)

class_sum <- Map(buildBoundariesDat, c("h1", "gm", "k5"), rep(c("c", "t"), 3), col="repClass", summary=T)
class_sum <- do.call(rbind, class_sum)
class_each <- Map(buildBoundariesDat, c("h1", "gm", "k5"), rep(c("c", "t"), 3), col="repClass", summary=F)
class_each <- do.call(rbind, class_each)

## 1) profiles:
pdf("plots/repFamily_v2.pdf", 6, 14)
ggplot(fam_sum, aes(x=pos, y=mean, col=type, fill=type)) + 
  facet_grid(feat~ct, scales="free") + 
  geom_vline(xintercept=10, lwd=4, col=I(rgb(.8,.8,.8, .5))) +
  geom_ribbon(aes(ymin=means.l, ymax=means.u), alpha=I(.3), col=I(NULL)) +
  geom_line() + theme_bw() +
  scale_color_brewer(type = "qual", palette = 2) + 
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_x_continuous(breaks=c(0, 10, 19), labels = c("-450 kb", "boundary", "+450 kb")) +
  labs(x="", y="Mean annotated repeats per bin", col="", fill="")
dev.off()

pdf("plots/repClass_v2.pdf", 6, 14)
ggplot(class_sum, aes(x=pos, y=mean, col=type, fill=type)) + 
  facet_grid(feat~ct, scales="free") + 
  geom_vline(xintercept=10, lwd=4, col=I(rgb(.8,.8,.8, .5))) +
  geom_ribbon(aes(ymin=means.l, ymax=means.u), alpha=I(.3), col=I(NULL)) +
  geom_line() + theme_bw() +
  scale_color_brewer(type = "qual", palette = 2) + 
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_x_continuous(breaks=c(0, 10, 19), labels = c("-450 kb", "boundary", "+450 kb")) +
  labs(x="", y="Mean annotated repeats per bin", col="", fill="")
dev.off()


## 2) bubble significance plots
pdf("plots/repClass_bubble_v2.pdf", 9, 5)
ggplot(blob(class_each), aes(x=feat, y=-log10(p), col=type,
  size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~., scales="free_y") + 
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

pdf("plots/repFamily_bubble_v2.pdf", 9, 5)
ggplot(blob(fam_each), aes(x=feat, y=-log10(p), col=type,
  size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~., scales="free_y") + 
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
files <- list.files("~/hvl/repeats/", pattern = ".*?\\.reptable", full.names = T)

bins <- read_tsv("~/hvl/repeats/h1_tb.bed", col_names = F)
head(bins)
nrow(bins)
# source of bins: 
s <- read_tsv("~/hvl/production/3dgenome/data/bedfiles/h1_tadbounds.bed")
nrow(s)



for( f in files ){
  # debugging
  f = files[[104]]
  
  f <- read_tsv(f, col_names = c("bin", "rep", "bp"))
  
  head(f)
  library("reshape2")
  rep <- dcast(bin ~ rep, value.var = "bp", data=f, fill=0)
  
  
  
  
  
  