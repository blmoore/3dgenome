### Test repeat classes over boundaries
# RepeatMasked tracks from UCSC table browser 30 June
library("blmR")
library("readr")

hrep <- read_tsv("~/hvl/repeats/hg19_ucscRepeatMasker.txt")
mrep <- read_tsv("~/hvl/repeats/mm10_ucscRepeatMasker.txt")

split_bed <- function(rep_class, df, genome = "hg19"){
    df <- subset(df, repClass == rep_class)
    # don't want unknowns with ?
    stopifnot(!grepl("\\?", df$repClass))
    bed <- df[,c("genoName", "genoStart", "genoEnd", "repName")]
    outname <- paste0("~/hvl/repeats/", genome, "_", rep_class, ".bed")
    message("Writing to ", outname)
    invisible(write_bed(bed, fname = outname))
}

classes <- unique(hrep$repClass)
classes <- classes[!grepl("\\?", classes)]

sapply(classes, split_bed, df=hrep)

unique(mrep$repClass)[(which(!unique(mrep$repClass) %in% classes))]

mclass <- unique(mrep$repClass)
mclass <- mclass[!grepl("\\?", mclass)]

# repeat classes are the same in humand and mouse data
stopifnot(all(classes %in% mclass) & all(mclass %in% classes))
sapply(mclass, split_bed, df=mrep, genome="mm10")

## ------------------------------------------------------------------------- ##

# read results on intersect.sh (bedtools intersect boundary bins vs. repeats)
buildBoundariesDat <- function(ct=c("h1", "k5", "gm"), type=c("c", "t")){

  ct <- match.arg(ct)
  type <- match.arg(type)
#   
#   ct="h1"
#   type="c"
#   
  if(type == "c"){
    # compartments
    pattern <- paste0("^", ct, "cb_.*\\.bed\\.out")
    #steps <- 30
    steps <- 25
  } else {
    # TADs
    pattern <- paste0("^", ct, "tb_.*\\.bed\\.out")
    steps <- 25
  }
  
  fl <- list.files("~/hvl/repeats/", pattern=pattern, full.names = T)
  looper <- fl  
  vars <- classes # classes, from above
  
  matches <- lapply(vars, function(x) grep(x, fl)[1])
  m <- fl[unlist(matches)]
  m <- na.omit(m)
  m <- unique(m)
  
  # Inspect matches:
  # cbind(gsub(".*\\/", "", m), vars)
  
  gdf <- data.frame(name=character(), 
    pos=numeric(),
    mean=numeric(), 
    mean.l=numeric(), 
    mean.u=numeric())
  bins <- matrix(nrow=0, ncol=steps)
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
}

hc <- buildBoundariesDat("h1", type = "c")
gc <- buildBoundariesDat("gm", type = "c")
kc <- buildBoundariesDat("k5", type = "c")

ht <- buildBoundariesDat("h1", type = "t")
kt <- buildBoundariesDat("k5", type = "t")
gt <- buildBoundariesDat("gm", type = "t")


pdf <- rbind(hc, gc, kc, ht, kt, gt)

library("ggplot2")
pdf("~/hvl/repeats/ctVrna.pdf", 6, 14)
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


