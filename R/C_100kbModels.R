############### build .dat files #######################
# Adaptation of 0_buildDatFiles.R for 100kb resolution #
########################################################
library("blmR")
set.seed(42)

read.pc <- function(fn)
  as.numeric(t(read.delim(fn, sep=" ", header=F))[,1])

# 100 kb eigs:
gpc <- read.pc("~/hvl/hiclib/hm/gm_corrected_hm_100kb.hdf5_PCA.out")
hpc <- read.pc("~/hvl/hiclib/hm/h1_corrected_hm_100kb.hdf5_PCA.out")
kpc <- read.pc("~/hvl/hiclib/hm/k5_corrected_hm_100kb.hdf5_PCA.out")

calcGC <- function(bedFile){
  # Get GC content for given bed file
  library("BSgenome.Hsapiens.UCSC.hg19")
  bed <- read.delim(bedFile, header=F)
  colnames(bed) <- c("chromo", "start", "end", "id")
  seqs <- getSeq(Hsapiens, bed$chromo, start=bed$start+1, end=bed$end, width=NA, 
                 as.character=F) 
  nucFreqs <- alphabetFrequency(seqs, as.prob=T, collapse=F, baseOnly=T)
  gc <- (nucFreqs[,2] + nucFreqs[,3]) 
  gc
}

buildAllDat <- function(ct=c("H1hesc", "Gm12878", "K562"), all=F, force=F){
  # Builds the feature set data frames per cell type. Run with 
  # all == TRUE, uses all available variables, when FALSE it 
  # only uses those features available for all cell types (35).
  # Force will rebuild df and rf model, overwriting existing.
  fl <- paste0("~/hvl/100/res/", list.files("~/hvl/100/res/", 
                                                 pattern=paste0(".*", ct,".*")))
  looper <- fl
  
  # initiate all.dat df
  a <- read.delim(fl[1], header=F)
  all.dat <- data.frame(row.names=a$V1)
  
  cat("\nUse all variables is set to: ", all, "\n\n")
  if(all == FALSE){
    vars <- c("Atf3", "Cebp", "Chd1", "Chd2", "Cmyc",
              "Ctcf", "Dnase", "Egr1", "Ezh2", "Gabp",
              "Jund", "Max", "H2az", "H3k27ac", "H3k27me3",
              "H3k36me3", "H3k4me1", "H3k4me2", "H3k4me3", "H3k79me2", 
              "H3k9ac", "H3k9me3", "H4k20me1", "Mxi1", "Nrsf", 
              "P300", "Pol2", "Rad21", "Six5", "Sp1", 
              "Taf1", "Tbp", "Yy1", "Znf143") # Mafk - not in gm12878
    
    matches <- lapply(vars, function(x) grep(x, fl)[1])
    m <- fl[unlist(matches)]
    m <- na.omit(m)
    m <- unique(m)
    
    # Inspect matches:
    cbind(gsub(".*\\/", "", m), vars)
    looper <- m
  } else {
    vars <- gsub(".*binnedBW/(.*?)\\..*","\\1", fl)
  }
  
  # Here we know length of files, i.e. outfilename, check if already written
  # and don't remake if it is. (+1 is for GC, not added yet.)
  outfilename <- paste0("~/hvl/ice/rds/100kb_", ct, "_", length(vars)+1, "Vars.rds")
  
  if(!file.exists(outfilename) | force == TRUE){
    cat("Reading in files... \n\n")
    for (i in looper){
      mark <- gsub(".*res/(.*?)\\..*","\\1", i)
      print(mark)
      b <- read.delim(i, header=F)
      # NB: b$5 = mean0, missing vals = 0, b$6 = mean of nonzero values
      all.dat <- cbind(all.dat, new = b[,5])
      colnames(all.dat)[ncol(all.dat)] <- mark
      rm(b)
    }
    
    if(!file.exists("~/hvl/100/gc.vec")){
      ## write "trimmed" bed file with actual bin ends
      library("BSgenome")
      library("BSgenome.Hsapiens.UCSC.hg19")
      library("dplyr")
      options(scipen=99)
      df <- data.frame(chr=gsub("-.*", "", rownames(all.dat)),
                       start=as.numeric(gsub(".*-", "", rownames(all.dat))))
      max <- group_by(df, chr) %.% summarise(start=max(start))
      max$end <- max$start + 1e6-1
      binrefs <- group_by(df, chr) %.% summarise(binRefs=which.max(start))
      
      seqs <- as.data.frame(Hsapiens@seqinfo)[1:23,]
      seqs <- seqs[match(max$chr, rownames(seqs)),]
      max$cend <- seqs$seqlengths
      max$binRefs <- binrefs$binRefs
      # which is smaller: Mb bin or chromosome end?
      max$binEnds <- apply(max[,3:4], 1, min)
      max$id <- paste0(max$chr, "-", max$start)
      
      bins <- read.table("~/hvl/100/bins.bed")
      bins$new <- bins$V3
      
      max <- max[order(match(as.character(max$chr), unique(as.character(bins$V1)))),]
      
      bins[bins$V4 %in% max$id,]$new <- max$binEnds
      bins$V3 <- bins$new
      bins$new <- NULL
      write.table(bins, "~/hvl/100/bins_trimmed.bed", sep="\t", quote=F,
                  row.names=F, col.names=F)
      
      # max length of hg19 chromosome length(Hsapines$chrN)
      
      gc <- calcGC("~/hvl/100/bins_trimmed.bed")
      write.table(gc, "~/hvl/100/gc.vec", quote=F, row.names=F,
                  col.names=F)
    }
    
    # calc GC and append (i.e. bedtools nuc -fi hg19.fa -bed 1mb_2811.bed > 1mb_gc.out)
    gc <- read.table("~/hvl/100/gc.vec")
    all.dat$GC <- gc[,1]
    
    # Double-check column names 
    cat("\n\nEnsure these names match !\n")
    eig <- if(ct == "H1hesc") hpc else 
      if(ct == "K562") kpc else gpc
    all.dat <- data.frame(eigen=eig, all.dat)
    print(cbind(colnames(all.dat), c("eigen", vars, "GC")))
    colnames(all.dat)[-1] <- c(vars, "GC")
    
    cat("\nBuilding Random Forest...\n\n")
    print(colnames(all.dat))
    saveRDS(all.dat, outfilename)
  } else {
    cat("Dataframe already exists! Reading...\n\n")
    all.dat <- readRDS(outfilename)
  }
  
  rf.file <- paste0("~/hvl/ice/rds/100kb_", ct, "_", ncol(all.dat[,-1]), "Vars_RFmod.rds")
  if(!file.exists(rf.file) | force == TRUE){
    res <- modelEigens.all(all.dat, n=200)
    saveRDS(res, rf.file)
  } else {
    res  <- readRDS(rf.file)
  }
  
  # Plot model results
  cat("Plotting ...\n\n")
  matched.col <- if(ct == "H1hesc") "orange" else 
    if(ct == "K562") "red" else "blue"
  pdf(paste0("~/hvl/ice/plots/100kb_", ct, "_", ncol(all.dat[,-1]), "Vars_RFres.pdf"),
      6, 6)
  plotPredRes.100kb(x=res$predicted, y=all.dat$eigen, 
                  ct=ct, col=matched.col)
  dev.off()
}

# Build 35 feature data frames:
buildAllDat("H1hesc", F, force=F)
buildAllDat("K562", F, force=F)
buildAllDat("Gm12878", F, force=F)

g.kbdat <- readRDS("~/hvl/ice/rds/100kb_Gm12878_35Vars.rds")
h.kbdat <- readRDS("~/hvl/ice/rds/100kb_H1hesc_35Vars.rds")
k.kbdat <- readRDS("~/hvl/ice/rds/100kb_K562_35Vars.rds")

## Custom plotPredRes for these smaller eigenvector values::
plotPredRes.100kb <- function (modelOut = NA, x = NA, y = NA, col = "blue", ct = "H1", 
          scale.factor = 0.7) 
{
  require("blmR")
  require("calibrate")
  if (!is.na(x) & !is.na(y)) {
    d1 <- x
    d2 <- y
  }
  else {
    stop("Must supply x and y.")
  }
  xy <- cbind(d1, d2)
  colour <- "#0000ff42"
  if (col != "blue") {
    ifelse(col == "red", colour <- "#ff000042", colour <- "#FFA50062")
  }
  def.par <- par(no.readonly = TRUE)
  xhist <- hist(d1, plot = FALSE)
  yhist <- hist(d2, breaks = 100, plot = F)
  top <- max(c(xhist$counts, yhist$counts))
  nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(6, 
                                                            1), c(1, 6), TRUE)
  par(mar = c(5, 5, 0, 0), mgp = c(1.8, 0.4, 0))
  #max <- max(abs(c(d1, d2))) * 1.02
  max=.18
  xmin <- -max *scale.factor
  xmax <- max *scale.factor
  plot(d1, d2, xlab = "Predicted eig", ylab = "Empirical eig", 
       col = colour, pch = 16, type = "n", xlim = c(xmin, xmax), 
       ylim = c(-max, max))
  abline(h = 0, v = 0, lty = 2)
  set.seed(42)
  # sample to reduce overplotting
  subsamp <- sample(1:length(d1), 5e3)
  points(d1[subsamp], d2[subsamp], col = colour, pch = 16)
  options(max.contour.segments = 50000)
  contour(kde2d(d1, d2, h = 0.1), nlevels = 8, drawlabels = F, 
          add = T)
  pos <- c(-max/2, max * 0.9, -max/2, max * 0.8, -max/2, max * 
             0.7, max/2, -max * 0.9, max/2, -max * 0.8)
  text(pos[1], pos[2], ct, font = 2, cex = 1.2)
  text(pos[3], pos[4], paste("PCC = ", signif(cor(d1, d2), 
                                              3), sep = ""), col = "navy")
  text(pos[5], pos[6], paste("RMSE = ", signif(blmR:::rmse(d1, d2), 
                                               3), sep = ""), col = "navy")
  text(pos[9], pos[10], paste("Acc. = ", 
                              signif(100 * (sum(apply(xy, 
               1, function(x) if (all(x > 0) || all(x < 0)) 1 else 0))/nrow(xy)), 
                                                4), sep = ""), col = "darkgreen")
  text(pos[7], pos[8], paste("AUROC = ", signif(blmR:::getAUC.gen(xy), 
                                                3), sep = ""), col = "darkgreen")
  par(mar = c(0, 5, 1, 1))
  plot(density(d1, from = xmin, to = xmax), frame = F, axes = F, 
       type = "n", col = colour, xlab = "", ylab = "", main = "", 
       xlim = c(xmin, xmax))
  polygon(density(d1, from = -max, to = max), col = colour, 
          border = colour, lwd = 4)
  par(mar = c(5, 0, 1, 1))
  plot(density(d2, from = -max, to = max)$y, density(d2, from = -max, 
                                                     to = max)$x, frame = F, axes = F, type = "n", col = colour, 
       xlab = "", ylab = "", main = "", ylim = c(-max, max))
  polygon(density(d2, from = -max, to = max)$y, density(d2, 
                                                        from = -max, to = max)$x, col = colour, border = colour, 
          lwd = 4)
  par(def.par)
}

g.kbmod <- readRDS("~/hvl/ice/rds/100kb_Gm12878_35Vars_RFmod.rds")
h.kbmod <- readRDS("~/hvl/ice/rds/100kb_H1hesc_35Vars_RFmod.rds")
k.kbmod <- readRDS("~/hvl/ice/rds/100kb_K562_35Vars_RFmod.rds")

impBars <- function(mod, ...){
  #imp <- mod$importance[,1]
  imp <- importance(mod, type=1)
  rownames(imp) <-  c("ATF3", "CEBP", "CHD1", "CHD2", 
                      "MYC", "CTCF", "DNase", "EGR1", 
                      "EZH2", "GABP", "JUND", "MAX",
                      "H2A.Z", "H3K27ac", "H3K27me3", 
                      "H3K36me3", "H3K4me1", "H3K4me2", 
                      "H3K4me3", "H3K79me2", "H3K9ac", 
                      "H3K9me3", "H4K20me1", "MXI1", 
                      "NRSF", "P300", "POL2", "RAD21", 
                      "SIX5", "SP1", "TAF1", "TBP", 
                      "YY1", "ZNF143", "GC")
  imp <- imp[order(imp[,1], decreasing=T),]
  barplot(rev(imp[1:10]), horiz=T, las=1, cex.names=1.5, ...)
}

pdf("~/hvl/ice/plots/f2b_100kb_varImpPerModel.pdf", 9, 3.5)
par(mfrow=c(1,3), mar=c(3,7,1.5,0.5), oma=c(2,0,0,0), mgp=c(0,.5,0))
impBars(g.kbmod, main="GM12878", col="#0000ff92", border=NA)
impBars(h.kbmod, main="H1 hESC", col="#FFA50092", border=NA)  
impBars(k.kbmod, main="K562", col="#ff000092", border=NA)
mtext(1, outer=T, text="Variable importance (% increase in MSE when permuted)")
dev.off()


## what about applying old models at new higher res??

g.mod <- readRDS("~/hvl/ice/rds/Gm12878_35Vars_RFmod.rds")
gkp <- predict(g.mod, g.kbdat)
gkkp <- predict(g.mod, k.kbdat)

plotPredRes.100kb(x=gkp, y=g.kbdat$eigen, scale.factor=1.2, 
                  ct="H1 hESC", col="orange")
plotPredRes.100kb(x=gkkp, y=k.kbdat$eigen, ct="GM12878 -> K562",
                  col="red", scale.factor=1.2)

dev.off()
