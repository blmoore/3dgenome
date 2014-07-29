############### build .dat files #####################
# Build key dataframes containing a) all features    #
# available for each cell type b) features available # 
# for all (~34). Also add in Dnase, GC to the MACS2  #
# ENCODE datasets.                                   #
######################################################
source("~/hvl/R/hvl.R")
set.seed(42)

read.pc <- function(fn)
  as.numeric(t(read.delim(fn, sep=" ", header=F))[,1])

# 1 MB eigs:
gpc <- read.pc("~/hvl/hiclib/hm/gm_corrected_hm_1Mb.hdf5_PCA.out")
hpc <- read.pc("~/hvl/hiclib/hm/h1_corrected_hm_1Mb.hdf5_PCA.out")
kpc <- read.pc("~/hvl/hiclib/hm/k5_corrected_hm_1Mb.hdf5_PCA.out")

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
  fl <- paste0("~/hvl/ice/binnedBW/", list.files("~/hvl/ice/binnedBW/", 
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
  outfilename <- paste0("~/hvl/ice/rds/", ct, "_", length(vars)+1, "Vars.rds")
  
  if(!file.exists(outfilename) | force == TRUE){
    cat("Reading in files... \n\n")
    for (i in looper){
      mark <- gsub(".*binnedBW/(.*?)\\..*","\\1", i)
      print(mark)
      b <- read.delim(i, header=F)
      # NB: b$5 = mean0, missing vals = 0, b$6 = mean of nonzero values
      all.dat <- cbind(all.dat, new = b[,5])
      colnames(all.dat)[ncol(all.dat)] <- mark
      rm(b)
    }
    
    if(!file.exists("~/hvl/ice//gc.vec")){
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
      
      bins <- read.table("~/hvl/ice/bins.bed")
      bins$new <- bins$V3
      
      max <- max[order(match(as.character(max$chr), unique(as.character(bins$V1)))),]
      
      bins[bins$V4 %in% max$id,]$new <- max$binEnds
      bins$V3 <- bins$new
      bins$new <- NULL
      write.table(bins, "~/hvl/ice/bins_trimmed.bed", sep="\t", quote=F,
                  row.names=F, col.names=F)
      
      # max length of hg19 chromosome length(Hsapines$chrN)
      
      gc <- calcGC("~/hvl/ice//bins_trimmed.bed")
      write.table(gc, "~/hvl/ice/gc.vec", quote=F, row.names=F,
                  col.names=F)
    }
    
    # calc GC and append (i.e. bedtools nuc -fi hg19.fa -bed 1mb_2811.bed > 1mb_gc.out)
    gc <- read.table("~/hvl/ice//gc.vec")
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
  
  rf.file <- paste0("~/hvl/ice/rds/", ct, "_", ncol(all.dat[,-1]), "Vars_RFmod.rds")
  if(!file.exists(rf.file) | force == TRUE){
    res <- modelEigens.all(all.dat, n=500)
    saveRDS(res, rf.file)
  } else {
    res  <- readRDS(rf.file)
  }
  
  # Plot model results
  cat("Plotting ...\n\n")
  matched.col <- if(ct == "H1hesc") "orange" else 
    if(ct == "K562") "red" else "blue"
  pdf(paste0("~/hvl/ice/plots/", ct, "_", ncol(all.dat[,-1]), "Vars_RFres.pdf"),
      6, 6)
  plotPredRes.ice(x=res$predicted, y=all.dat$eigen, 
                  ct=ct, col=matched.col)
  dev.off()
}

# Build 35 feature data frames:
buildAllDat("H1hesc", F, force=T)
buildAllDat("K562", F, force=T)
buildAllDat("Gm12878", F, force=T)

# Build all (~70-150?):
buildAllDat("H1hesc", T)
buildAllDat("K562", T)
buildAllDat("Gm12878", T, force=T)

## Write bed files of the new eigenvector data:
writeEigBed <- function(dat, fn){
  options(scipen=99)
  df <- data.frame(id     = rownames(dat),
                   chr    = gsub("-.*", "", rownames(dat)),
                   start  = as.numeric(gsub(".*-", "", rownames(dat))),
                   end    = as.numeric(gsub(".*-", "", rownames(dat))) + 1e6,
                   strand = "+",
                   eig    = dat$eigen)
  write.table(df, file=fn, quote=F, row.names=F, col.names=F, sep="\t")
}

writeEigBed(h.dat, "~/hvl/ice/h1_eigs.bed")
writeEigBed(g.dat, "~/hvl/ice/gm_eigs.bed")
writeEigBed(k.dat, "~/hvl/ice/k5_eigs.bed")
