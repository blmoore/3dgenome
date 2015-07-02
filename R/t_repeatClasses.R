### Test repeat classes over boundaries
# RepeatMasked tracks from UCSC table browser 30 June
library("blmR")
library("readr")

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
    steps <- 25
  } else {
    # TADs
    pattern <- paste0("^", ct, "tb_.*?", col, ".*\\.bed\\.out")
    steps <- 25
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

hc <- buildBoundariesDat("h1", type = "c", col="repFamily")
gc <- buildBoundariesDat("gm", type = "c", col="repFamily")
kc <- buildBoundariesDat("k5", type = "c", col="repFamily")

ht <- buildBoundariesDat("h1", type = "t", col="repFamily")
kt <- buildBoundariesDat("k5", type = "t", col="repFamily")
gt <- buildBoundariesDat("gm", type = "t", col="repFamily")


pdf <- rbind(hc, gc, kc, ht, kt, gt)

library("ggplot2")
pdf("~/hvl/repeats/repFamily.pdf", 6, 14)
ggplot(pdf, aes(x=pos, y=mean, col=type, fill=type)) + 
  facet_grid(feat~ct, scales="free") + 
  geom_vline(xintercept=12.5, lwd=4, col=I(rgb(.8,.8,.8, .5))) +
  geom_ribbon(aes(ymin=means.l, ymax=means.u), alpha=I(.3), col=I(NULL)) +
  geom_line() + theme_bw() +
  scale_color_brewer(type="qual", palette = 2) + 
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_x_continuous(breaks=c(0, 12.5, 25), labels = c("-450 kb", "boundary", "+450 kb")) +
  labs(x="", y="Mean annotated repeats per bin", col="", fill="")
dev.off()


hcf <- buildBoundariesDat("h1", type = "c", col="repClass", F)
gcf <- buildBoundariesDat("gm", type = "c", col="repClass", F)
kcf <- buildBoundariesDat("k5", type = "c", col="repClass", F)

htf <- buildBoundariesDat("h1", type = "t", col="repClass", F)
ktf <- buildBoundariesDat("k5", type = "t", col="repClass", F)
gtf <- buildBoundariesDat("gm", type = "t", col="repClass", F)

pdff <- rbind(hcf, gcf, kcf, htf, ktf, gtf)

library("dplyr")
blobs <- group_by(pdff, ct, type, feat) %>% 
  summarise(bound.mean = mean(X12),
    edge.mean  = mean(c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21)),
    p = wilcox.test(X12, c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21))$p.value)
head(blobs)

pdf("~/hvl/repeats/rptclass_bubble.pdf", 9, 5)
ggplot(blobs, aes(x=feat, y=-log10(p), col=type,
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
