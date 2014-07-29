########### Additional figures collection ############
# A series of plots used in supplementary material,  #
# currently includes:                                #
#                                                    #
# * S3: Variable importance and mutual information   # 
# * S4: Boxplot of all Open vs. Closed features      #
######################################################
library("ggplot2")
library("gridExtra")
library("infotheo")
library("randomForest")
library("RColorBrewer")
library("reshape2")
devtools::load_all("blmR")
source("~/blmR/R/misc_functions.R")

############################################### S3 ###

## 1) Variable importance
tt <- data.frame(ct    = c(rep("GM12878", 35),
                           rep("H1 hESC", 35),
                           rep("K562", 35)),
                 feat  = rep(rownames(importance(g.mod)), 3),
                 value = c(importance(g.mod, type=1), 
                           importance(h.mod, type=1),
                           importance(k.mod, type=1)),
                 measure = as.character("Variable importance"))


## 2) Mutual information
calcsLong <- function(dat, cellt){
  g.disc <- discretize(dat,n=10)
  g.mi <- mutinformation(g.disc)[1,][-1]
  g <- data.frame(ct   = rep(cellt, length(g.mi)), 
                  feat  = names(g.mi), 
                  value = g.mi, 
                  measure = as.character("Mutual information"))
  g
}

gb <- calcsLong(g.dat, "GM12878")
hb <- calcsLong(h.dat, "H1 hESC")
kb <- calcsLong(k.dat, "K562")
ab <- rbind(gb, hb, kb)

bothdf <- rbind(tt, ab)

## Sort by mean Variable importance (it)
bothdf$feat <- with(bothdf, reorder(feat, value))
bothdf$ct <- factor(iraw$ct, levels=c("GM12878", "H1", "K562"))

# Proper nomenclature, caps:
caps <- toupper(levels(bothdf$feat))
caps <- gsub("H2AZ", "H2A.Z",
             gsub("DNASE", "DNase", 
                  gsub("AC", "ac", 
                       gsub("ME", "me", caps))))
levels(bothdf$feat) <- caps

pdf("~/hvl/ice/plots/supplementary/s3_MIandVarImp_v2.pdf", 9, 10)
grid.arrange(
  ggplot(subset(bothdf, measure!="Variable importance"), aes(x=feat, y=value, col=ct)) + 
    geom_point(size=3) + coord_flip() + facet_wrap(~measure) +
    theme_bw() + 
    ylab("Mutual information") + xlab("") + #ylim(0,.4) +
    scale_colour_manual(values=c("#0000ff98", "#ff000098", "#FFA50098"), guide=F),
  ggplot(subset(bothdf, measure=="Variable importance"), 
         aes(x=feat, y=value, col=NULL, fill=ct)) + 
    geom_bar(position=position_dodge(.9), stat="identity") + coord_flip() +
    facet_wrap(~measure) + theme_bw() + 
    labs(list(fill="", x="", y="Increase in RMSE when permuted (% RMSE)")) +
    theme(legend.position = c(0.75, 0.15), legend.background=element_blank()) +
    scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) 
  , ncol=2 
)
dev.off()


############################################### S4 ###

# Chrom features boxplots
all.dat <- rbind(
  transform(melt(g.dat[,-1]), ct = "GM12878", state=callStates(g.dat$eigen)$state),
  transform(melt(h.dat[,-1]), ct = "H1 hESC", state=callStates(h.dat$eigen)$state),
  transform(melt(k.dat[,-1]), ct = "K562", state=callStates(k.dat$eigen)$state))

all.dat$classify <- ifelse(grepl("^H", all.dat$variable), "Histone related",
                           "DNA binding proteins")
all.dat$state <- factor(ifelse(all.dat$state == 1, "Closed", "Open"),
                        levels=c("Open", "Closed"))

# Proper nomenclature, caps:
caps <- toupper(levels(all.dat$variable))
caps <- gsub("H2AZ", "H2A.Z",
             gsub("DNASE", "DNase", 
                  gsub("AC", "ac", 
                       gsub("ME", "me", caps))))
levels(all.dat$variable) <- caps

pdf("~/hvl/ice/plots/supplementary/s4_featureBoxplots_v2.pdf", 11, 9)
ggplot(subset(all.dat, variable != "GC"), 
       aes(x=variable, y=value, fill=state)) + 
  facet_grid(ct~classify, scales="free", space="free") + 
  xlab("") + theme_bw() + labs(fill="") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("grey90","grey50")) +
  ylab("Fold-change signal relative to control") + 
  geom_boxplot(notch=T, outlier.size=.7) +
  coord_cartesian(ylim=c(0,2.6))
dev.off()
