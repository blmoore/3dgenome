########### Additional figures collection ############
# A series of plots used in supplementary material,  #
# currently includes:                                #
#                                                    #
# * S3: Variable importance and mutual information   # 
# * S5: Boxplot of all Open vs. Closed features      #
######################################################
library("dplyr")
library("ggplot2")
library("gridExtra")
library("infotheo")
library("randomForest")
library("reshape2")
library("blmR")

# load previously-built models:
g.mod <- readRDS("data/rds/Gm12878_35Vars_RFmod.rds")
h.mod <- readRDS("data/rds/H1hesc_35Vars_RFmod.rds")
k.mod <- readRDS("data/rds/K562_35Vars_RFmod.rds")
# and raw data:
g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

############################################### S4 ###

## 1) Variable importance
tt <- data.frame(ct    = c(rep("GM12878", 35),
                           rep("H1 hESC", 35),
                           rep("K562", 35)),
                 feat  = rep(colnames(g.dat)[-1], 3),
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

#bothdf <- rbind(tt, ab)
## no longer using M.I.:
bothdf <- tt

ordering <- group_by(bothdf, measure, feat) %>% 
    dplyr::summarise(tot=sum(value))
ordering <- ordering[order(ordering$tot, decreasing=T),]

bothdf$feat <- factor(bothdf$feat, 
                      rev(levels(bothdf$feat)[match(as.character(ordering$feat), levels(bothdf$feat))]))

# Proper nomenclature, caps:
caps <- toupper(levels(bothdf$feat))
caps <- gsub("H2AZ", "H2A.Z",
             gsub("DNASE", "DNase", 
                  gsub("AC", "ac", 
                       gsub("ME", "me", caps))))
levels(bothdf$feat) <- caps

v.i <- subset(bothdf, measure=="Variable importance")
v.i <- dcast(v.i, feat~ct)
v.i$gh <- with(v.i, abs(GM12878 - `H1 hESC`))
v.i$hk <- with(v.i, abs(`H1 hESC` - K562))
v.i$kg <- with(v.i, abs(K562 - GM12878))
v.i <- v.i[,c(1,5:7)]

v.i <- melt(v.i)
v.i$feat <- factor(v.i$feat, levels=rev(unique(as.character(v.i[order(v.i$value, decreasing=T),]$feat))))

v.i$variable <- plyr::revalue(v.i$variable, c("gh" = "GM12878 vs. H1",
                                        "hk" = "H1 vs. K562",
                                        "kg" = "K562 vs. GM12878"))

## Combine above two:
pdf("figures/suppl/s3_VarImpDifferences.pdf", 9, 9)
grid.arrange(  
  ggplot(subset(bothdf, measure=="Variable importance"), 
         aes(x=feat, y=value, col=NULL, fill=ct)) + 
    geom_bar(position=position_dodge(.9), stat="identity") + coord_flip() +
    facet_wrap(~measure) + theme_bw() + 
    labs(list(fill="Cell type", x="", y="Increase in MSE when permuted (% MSE)")) +
    # formerly .08
    theme(legend.position = c(0.82, 0.08), legend.background=element_blank()) +
    scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")),
  ggplot(cbind(v.i, title="Difference between cell types"), 
         aes(x=feat, y=value, fill=variable, col=variable)) + 
    geom_point(size=I(3.4), shape=21) + coord_flip() + theme_bw() +
    facet_wrap(~title) +
    ## Inside of points (ct #1)
    scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) +
    labs(y="Absolute difference in variable importance",
         x="", col="Comparison", fill="Comparison") +
    theme(legend.position = c(0.7, 0.08), legend.background=element_blank()) +
    ## Outline of points (ct #2)
    scale_color_manual(values=c( "#FFA50098", "#ff000098", "#0000ff98")),
  ncol=2
)
dev.off()

############################################### S5 ###

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

pdf("figures/suppl/s5_featureBoxplots.pdf", 11, 9)
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

