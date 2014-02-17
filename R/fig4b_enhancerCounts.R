############### homerFlipped.R ######################
# Consider the regions called by homer as "flipped" #
# and test for chromatin state enrichments between  #
# cell types.                                       #
#####################################################

f4 <- function(){
  library("RHmm")
  library("ggplot2")
  library("gridExtra")
  source("~/hvl/R/hvl.R")
  
  callStates <- function(indat){
    set.seed(42)
    hmm <- HMMFit(indat, dis="NORMAL", nStates=2, control=list(init="RANDOM", iter=5000))
    calls <- viterbi(hmm, indat)
    out <- data.frame(eigen=indat, state=calls$states)
    return(out)
  }
  
  # *.dat opened via hvl.R
  h.calls <- callStates(h.dat$eigen)
  k.calls <- callStates(k.dat$eigen)
  g.calls <- callStates(g.dat$eigen)
  
  all.calls <- data.frame(h1=h.calls[,2], k5=k.calls[,2], gm=g.calls[,2])
  rownames(all.calls) <- rownames(h.dat)
  
  ## These RDS are built using homerSigInt_BG.R with the output
  ## from countStates.py -- tables of C/S states per MB
  gcounts <- readRDS("~/hvl/chromhmm/gm_1mb_ecounts.rds")
  hcounts <- readRDS("~/hvl/chromhmm/h1_1mb_ecounts.rds")
  kcounts <- readRDS("~/hvl/chromhmm/k5_1mb_ecounts.rds")
  
  ensEnrich <- function(state=c(1,2,3), counts){
    # state 2 == open/A, 1 == closed/B
    others <- setdiff(1:3, state)
    ovec <- apply(all.calls[,others], 1, sum)
    closed <- all.calls[all.calls[,state] == 1 & ovec == 4,]
    open <- all.calls[all.calls[,state] == 2 & ovec == 2,]
    
    # signif diff dist?
    open <- counts[counts$bin %in% rownames(open),]
    close <- counts[counts$bin %in% rownames(closed),]
    #print(wilcox.test(counts$E.TS, open$E.TS))
    
    out <- rbind(
      cbind("shared", "open", open$E.CN),
      cbind("shared", "all", counts$E.CN),
      cbind("shared", "closed", close$E.CN), 
      cbind("tissue-specific", "open", open$E.TS),
      cbind("tissue-specific", "all", counts$E.TS),
      cbind("tissue-specific", "closed", close$E.TS))
    return(out)
  }
  
  hsum <- ensEnrich(1, hcounts)
  ksum <- ensEnrich(2, kcounts)
  gsum <- ensEnrich(3, gcounts) 
  asum <- rbind(cbind(hsum, "H1"),
                cbind(ksum, "K562"),
                cbind(gsum, "GM12878"))
  
  asum <- as.data.frame(asum)
  colnames(asum) <- c("etype", "compartment", "num", "ct")
  asum$num <- as.numeric(as.character(asum$num))
  asum$compartment <- factor(asum$compartment, levels=c("closed", "all", "open"))
  
  shared <- subset(asum, asum$etype == "shared")
  tissue <- subset(asum, asum$etype == "tissue-specific")
  grid.arrange(
    ggplot(shared, aes(x=compartment, y=num, fill=compartment)) + 
      coord_cartesian(ylim = quantile(shared$num, c(0, 0.95))) +
      facet_grid(etype~ct, scale="free_y", space="free_y") +
      geom_boxplot(notch=T, outlier.size=0) + theme_bw() + 
      scale_fill_brewer(palette="Blues") + 
      labs(y="Number of enhancers per Mb", x="", fill="Relative\ncompartment\nstate"),
    ggplot(tissue, aes(x=compartment, y=num, fill=compartment)) + 
      coord_cartesian(ylim = quantile(tissue$num, c(0, 0.95))) +
      facet_grid(etype~ct, scale="free_y", space="free_y") +
      #geom_violin() + theme_bw() +
      geom_boxplot(notch=T, outlier.size=0) + theme_bw() + 
      scale_fill_brewer(palette="Blues") + 
      labs(y="Number of enhancers per Mb", x="", fill=""), ncol=1
  )
}
  
tsignif <- function(comp1 = c("open", "all", "closed"), 
                    comp2 = c("open", "all", "closed"),
                    ct = c("K562", "GM12878", "H1")){
  print(wilcox.test(shared[shared$compartment == comp1 & shared$ct == ct, "num"], 
                    shared[shared$compartment == comp2 & shared$ct == ct, "num"],
                    alternative="two.sided"))
  print(wilcox.test(tissue[tissue$compartment == comp1 & tissue$ct == ct, "num"], 
                    tissue[tissue$compartment == comp2 & tissue$ct == ct, "num"],
                    alternative="two.sided"))
}

f4.p <- function(){
  tsignif("open", "closed", "K562")    # shared: .02, tissue: .58
  tsignif("open", "closed", "GM12878") # shared: .76, tissue: 2e-9
  tsignif("open", "closed", "H1")      # shared: .04, tissue: .83
  
  tsignif("open", "all", "K562")    # shared: .72, tissue: .005
  tsignif("open", "all", "GM12878") # shared: .5, tissue: 2e-10
  tsignif("open", "all", "H1")      # shared: 2e-5, tissue: .29
  
  tsignif("closed", "all", "K562")    # shared: 4e-6, tissue: 1e-4
  tsignif("closed", "all", "GM12878") # shared: .66, tissue: .22
  tsignif("closed", "all", "H1")      # shared: .93, tissue: .65
}