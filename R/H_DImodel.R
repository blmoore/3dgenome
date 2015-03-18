## predict abs(Directionality index)

di_bedfile <- function(in_di, out_bed){
  di <- read.table(in_di)
  
  # proper chromosome formatting
  di[,1] <- gsub("chr23", "chrX", paste0("chr", di[,1]))
  
  # id column
  di[,4] <- gsub(" ", "", apply(di[,1:2], 1, paste0, collapse="-"))
  
  blmR::write_bed(di, out_bed)
  invisible()
}

di_bedfile("data/text/h1.di", "data/bedfiles/h1_40kb.bins")
di_bedfile("data/text/gm.di", "data/bedfiles/gm_40kb.bins")
di_bedfile("data/text/k5.di", "data/bedfiles/k5_40kb.bins")


h1_di <- read.table("data/text/h1.di")
h1_tads <- read.table("data/text/h1.tads")

c12 <- h1_di[h1_di[,1] == 12,]

par(mfrow=c(4,1))
plot(c12[100:1000,2], log2(abs(c12[100:1000,4])+1), type="l")
abline(v=h1_tads[h1_tads[,1] == "chr12",2], col="blue")

plot(c12[100:1000,2], loess.smooth(log2(abs(c12[100:1000,4])+1), c12[100:1000,2])$y)

loess.smooth(log2(abs(c12[,4]))+1)

head(h1_tads)
## Analyse

g_ctcf <- read.table("data/nonrepo/all40kb/gm_40kb.bins_BroadHistoneGm12878CtcfStdAlnRep0.bam_VS_BroadHistoneGm12878ControlStdAlnRep0.bam.fc.signal.bw")
g_di <- read.table("data/text/gm.di")


cor(g_ctcf[,5], abs(g_di[,4]), method="spearman")

## build df of all features w/ TAD calls (adapted from 0_buildDatFiles.R)::
g_di <- read.table("data/text/gm.di")
h_di <- read.table("data/text/h1.di")
k_di <- read.table("data/text/k5.di")

buildAllDat <- function(ct=c("H1hesc", "Gm12878", "K562"), all=F, force=F){
  #ct = "H1hesc"
  dir = "data/nonrepo/all40kb/"
  fl <- paste0(dir, list.files(dir, pattern=paste0(".*", ct,".*")))
  looper <- fl
  
  # initiate all.dat df
  a <- read.delim(fl[1], header=F)
  all.dat <- data.frame(row.names=a$V1)
  
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
  # cbind(gsub(".*\\/", "", m), vars)
  looper <- m
  
  # Here we know length of files, i.e. outfilename, check if already written
  # and don't remake if it is. (+1 is for GC, not added yet.)
  outfilename <- paste0("data/rds/all40kb_", ct, "_", length(vars)+1, "Vars.rds")
  
  if(!file.exists(outfilename) | force == TRUE){
    cat("Reading in files... \n\n")
    for (i in looper){
      mark <- gsub(".*binnedBigWigs/(.*?)\\..*","\\1", i)
      print(mark)
      b <- read.delim(i, header=F)
      # NB: b$5 = mean0, missing vals = 0, b$6 = mean of nonzero values
      all.dat <- cbind(all.dat, new = b[,5])
      colnames(all.dat)[ncol(all.dat)] <- mark
      rm(b)
    }
    
    
    eig <- if(ct == "H1hesc") h_di[,4] else 
      if(ct == "K562") k_di[,4] else g_di[,4]
    all.dat <- data.frame(eigen=eig, all.dat)
    # Double-check column names 
    cat("\n\nCheck these names match !\n")
    print(cbind(colnames(all.dat), c("eigen", vars)))
    colnames(all.dat)[-1] <- vars
    
    cat("\nBuilding Random Forest...\n\n")
    saveRDS(all.dat, outfilename)
  } else {
    cat("Dataframe already exists! Reading...\n\n")
    all.dat <- readRDS(outfilename)
  }
  
  rf.file <- paste0("data/rds/di_", ct, "_", ncol(all.dat[,-1]), "Vars_RFmod.rds")
  if(!file.exists(rf.file) | force == TRUE){
    res <- blmR::modelEigens.all(all.dat, n=200)
    saveRDS(res, rf.file)
  } else {
    res  <- readRDS(rf.file)
  }
}

buildAllDat("H1hesc")

h_40 <- readRDS("data/rds/all40kb_H1hesc_35Vars.rds")
str(h_40)

rfmod <- randomForest::randomForest(eigen~., data=h_40[sample(1:nrow(h_40), 5e3),], ntrees=200)
cor(predict(rfmod, h_40[1:100,]), h_40[1:100,"eigen"])

# try supervised classification, not regression
# is boundary: T/F
h_bounds <- read.table("data/bedfiles/h1_tadbounds.bed")
h_40$is_bound <- factor(ifelse(rownames(h_40) %in% h_bounds[,4], T, F))

## 1133 boundaries out of 72598 bins (1.56%)
t_samp <- sample(which(h_40$is_bound == T), 200)
f_samp <- sample(which(h_40$is_bound == F), 1000)

cmod <- randomForest::randomForest(is_bound~., data=h_40[c(t_samp, f_samp),-1], ntrees=500, nodesize=20)
cmod

# try SVM, supposedly less sensitive to class imbalance
library("e1071")
smod <- svm(is_bound~., data=h_40[c(t_samp, f_samp),-1])
plot(smod, h_40[c(t_samp, f_samp),-1], Ctcf~H3k36me3)

test <- (1:nrow(h_40))[!1:nrow(h_40) %in% c(t_samp, f_samp)]
t <- h_40[test,][sample(1:length(test), 500),]
class <- predict(cmod, t)

100 * length(which(t$is_bound == class)) / nrow(t)
length(which(t$is_bound == F & class == F))
length(which(t$is_bound == T & class == T))

pdf("~/Desktop/h1_tad_bounds_hm.pdf", 5, 6)
par(mar=c(3,3,3,2))
gplots::heatmap.2(log2(as.matrix(h_40[c(t_samp, f_samp),!colnames(h_40) %in% c("is_bound", "hmmclass", "eigen")])
                       +1), trace="none", RowSideColors=c("white", "black")[h_40[c(t_samp, f_samp),]$is_bound],
                  labRow="", keysize=1, key.title=NA, dendrogram="column", Rowv=F)
dev.off()

cmod
par(mar=c(4,8,3,3))
randomForest::varImpPlot(cmod)

## what about Bing Ren's HMM calls, can we predict those?
hmmo <- read.table("~/hvl/ice/tads/hmmout/h1_hmm.out", skip=2)
nrow(hmmo)
nrow(h_40)


h_40$hmmclass <- factor(hmmo[,1])

hmod <- randomForest::randomForest(hmmclass~., data=h_40[,!colnames(h_40) %in% c("eigen", "is_bound")], ntrees=500)
hmod


## average epigenetic signals over TADs (e.g. BWAOB tad, epi files)
## cluster