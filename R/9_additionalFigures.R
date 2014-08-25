########### Additional figures collection ############
# A series of plots used in supplementary material,  #
# currently includes:                                #
#                                                    #
# * S3: Variable importance and mutual information   # 
# * S4: Boxplot of all Open vs. Closed features      #
######################################################
library("dplyr")
library("ggplot2")
library("gridExtra")
library("infotheo")
#library("plyr")
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

## Sort by mean Variable importance (it)
#bothdf$feat <- with(bothdf, reorder(feat, ordering$tot))
#bothdf$feat <- with(bothdf, reorder(feat, value))
#bothdf$ct <- factor(bothdf$ct, levels=c("GM12878", "H1", "K562"))
# 
# # Proper nomenclature, caps:
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
    scale_colour_manual(values=c("#0000ff98", "#FFA50098", "#ff000098"), guide=F),
  ggplot(subset(bothdf, measure=="Variable importance"), 
         aes(x=feat, y=value, col=NULL, fill=ct)) + 
    geom_bar(position=position_dodge(.9), stat="identity") + coord_flip() +
    facet_wrap(~measure) + theme_bw() + 
    labs(list(fill="", x="", y="Increase in MSE when permuted (% MSE)")) +
    theme(legend.position = c(0.75, 0.15), legend.background=element_blank()) +
    scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) 
  , ncol=2 
)
dev.off()

pdf("~/hvl/ice/plots/supplementary/s3_MIandVarImp_v4.pdf", 4.5, 8.5)
ggplot(subset(bothdf, measure=="Variable importance"), 
       aes(x=feat, y=value, col=NULL, fill=ct)) + 
  geom_bar(position=position_dodge(.9), stat="identity") + coord_flip() +
  facet_wrap(~measure) + theme_bw() + 
  labs(list(fill="", x="", y="Increase in MSE when permuted (% MSE)")) +
  # formerly .08
  theme(legend.position = c(0.82, 0.18), legend.background=element_blank()) +
  scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) 
dev.off()

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

ggplot(v.i, aes(x=feat, y=value, fill=variable, col=variable)) + 
  geom_point(size=I(3.4), shape=21) + coord_flip() + theme_bw() +
  #facet_grid(.~variable) +
  ## Inside of points (ct #1)
  scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) +
  labs(y="Absolute difference in variable importance (Percentage points)",
       x="", col="Comparison", fill="Comparison") +
  theme(legend.position = c(0.82, 0.18), legend.background=element_blank()) +
  ## Outline of points (ct #2)
  scale_color_manual(values=c( "#FFA50098", "#ff000098", "#0000ff98")) 

## Combine above two:
pdf("~/hvl/ice/plots/supplementary/s3_MIandVarImp_v5.pdf", 9, 9)
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

mvi <- group_by(bothdf, feat) %>%
  summarise(mean.varimp=mean(value))
mvv <- group_by(v.i, feat) %>%
  summarise(mean.variability=mean(value))

mv <- merge(mvi, mvv, by="feat")
library("MASS")
ggplot(mv, aes(x=mean.varimp, y=mean.variability)) +
  geom_smooth(method="rlm") +
  geom_text(aes(label=feat)) 

## Issues: 
#   not percentages of the same abslolutes
#   interpretability

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

##### Full models variable importance
nameTidy <- function(string, ct=c("h", "k", "g")){
  string  <- rownames(gimp)
  ## 1) dump tail-end
  string <- gsub(pattern="Aln.*|Std.*|Ucd.*|V04.*|Control.*|Igg.*|sc\\d+.*|Ifn.*|Pcr.*", "", x=string)
  ## 2) bin front end
  if(ct == "k"){
    string <- gsub("^.*?K562", "", string)
    string[178] <- "DNase"
  } else {
    if(ct == "h"){
      string <- gsub("^.*H1hesc", "", string)
      string[70] <- "DNase"
      } else {
        string <- gsub("^.*Gm12878", "", string)
        string[111] <- "DNase"
      }
  }
  string
}

kvars <- readRDS("~/hvl/ice/rds/K562_187Vars.rds")
kfull <- randomForest(eigen~., data=kvars, importance=T)
# saveRDS(kfull, "~/hvl/ice/rds/K562_187Vars_RFmod_imp.rds")
kimp <- importance(kfull, type=1)

hvars <- readRDS("~/hvl/ice/rds/H1hesc_71Vars.rds")
hfull <- randomForest(eigen~., data=hvars, importance=T)
#saveRDS(hfull, "~/hvl/ice/rds/H1hesc_35Vars_RFmod_imp.rds")
himp <- importance(hfull, type=1)

gvars <- readRDS("~/hvl/ice/rds/Gm12878_115Vars.rds")
gfull <- randomForest(eigen~., data=gvars, importance=T)
#saveRDS(gfull, "~/hvl/ice/rds/Gm12878_115Vars_RFmod_imp.rds")
gimp <- importance(gfull, type=1)

head(kimp)
rownames(gimp) <- nameTidy(rownames(gimp), "g")
rownames(himp) <- nameTidy(rownames(himp), "h")
rownames(kimp) <- nameTidy(rownames(kimp), "k")


par(mar=c(4,6,3,3), mfrow=c(1,3))
barplot(gimp[order(gimp[,1], decreasing=T),][20:1], horiz=T, las=1)
barplot(himp[order(himp[,1], decreasing=T),][20:1], horiz=T, las=1)
barplot(kimp[order(kimp[,1], decreasing=T),][20:1], horiz=T, las=1)


cor(hfull$y, h.dat$eigen)
