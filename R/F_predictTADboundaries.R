library("ggplot2")
library("reshape2")

# from 7_fig5boundaryEnrichments.R
all.t <- readRDS("data/rds/tad_boundary_features.rds")

# ctcf <- subset(all.t, feat == "Ctcf")
# plot(ctcf$X12, pch=20, cex=.1, col="darkred")
# points(ctcf$X1, pch=20, cex=.1, col="darkblue")
# 
# mctcf <- melt(ctcf, id.vars=c("name", "bound", "feat", "ct", "type"))
# 
# pdf("~/hvl/ice/plots/ctcf_over_tad_bounds.pdf", 8, 3)
# ggplot(#subset(mctcf, variable %in% c("X1", "X12")),
#   mctcf, aes(x=variable, y=value, group=bound)) + #geom_violin() +
#   geom_line(alpha=.025) + xlab("") + ylab("ctcf") +
#   facet_wrap(~ct) + theme_classic() +
#   theme(axis.ticks.x=element_blank(),
#         axis.text.x=element_blank())
# dev.off()
# 
# hist(ctcf$X12 - ctcf$X1)
# 
# ggplot(mctcf[mctcf$bound %in% sample(1:2809, 20),], 
#   aes(x=variable, y=value, group=bound)) + #geom_violin() +
#   geom_line(alpha=.5) + xlab("") + ylab("ctcf") +
#   facet_wrap(~ct) + theme_classic() +
#   theme(axis.ticks.x=element_blank(),
#     axis.text.x=element_blank()) +
#   geom_vline(xintercept=15)


feats <- melt(all.t, id.vars=c("name", "bound", "feat", "ct", "type"))
teg <- subset(feats, variable %in% c("X1", "X12"))
teg$variable <- ifelse(teg$variable == "X12", T, F)
teg$name <- NULL
#table(teg$variable)

taddf <- dcast(teg, bound + type + ct + variable ~ feat, value.var="value")
head(taddf)


taddf$variable <- as.factor(ifelse(taddf$variable == T, 1, 0))

library("AUCRF")
library("ROCR")
set.seed(42)

ctrf_auroc <- function(ct = c("Gm12878", "H1hesc", "K562")){
  ct <- match.arg(ct)
  h1df <- subset(taddf, ct == ct)[,-c(1:3)]
  train <- sample(1:nrow(h1df), nrow(h1df)/2)
  message("training : ", length(train), " out of ", nrow(h1df), " total rows")
  
  auc_class <- AUCRF(variable ~ ., data=h1df[train,], ntrees=200, pdel=0.1)
  opt_rf <- randomForest(as.formula(paste0("variable~", paste(as.vector(OptimalSet(auc_class)$Name), collapse="+"))),
    data=h1df[train,])

  preds <- predict(opt_rf, data=h1df[-train,], type="prob")
  message("predictions : ", length(preds))
  
  #print(head(preds))
  
  #perf <- performance(  prediction(opt_rf$votes[,2], h1df$variable), "tpr", "fpr")
  perf <- performance(prediction(preds[,2], h1df[-train,]$variable), "tpr", "fpr")
  df <- data.frame(fp=unlist(perf@x.values), tp=unlist(perf@y.values), ct=ct)
  obj <- list()
  
  obj$data <- h1df
  obj$train_ind <- train
  obj$predictions <- preds
  obj$auroc <- df
  obj$aurcrf <- auc_class
  obj$rf <- opt_rf
  return(obj)
}

h1auc <- ctrf_auroc("H1hesc")
h1auc$auroc$ct <- "H1 hESC"
gmauc <- ctrf_auroc("Gm12878")
gmauc$auroc$ct <- "GM12878"
k5auc <- ctrf_auroc("K562")

# saveRDS(h1auc, "data/rds/h1_tadpred.rds")
# saveRDS(gmauc, "data/rds/gm_tadpred.rds")
# saveRDS(k5auc, "data/rds/k5_tadpred.rds")
gmauc <- readRDS("data/rds/gm_tadpred.rds")
h1auc <- readRDS("data/rds/h1_tadpred.rds")
k5auc <- readRDS("data/rds/k5_tadpred.rds")

auc <- rbind(h1auc$auroc, gmauc$auroc, k5auc$auroc)

pdf("~/hvl/thesis_plots/tad_bounds_auroc_v2.pdf", 4.4, 4.4)
ggplot(auc, aes(x=fp, y=tp, col=ct)) + geom_line() +
  annotate("segment", x=0, xend=1, y=0, yend=1, linetype="dashed", col=I("grey20")) +
  scale_x_continuous(expand=c(0,0)) + coord_fixed() +
  scale_y_continuous(expand=c(0,0)) +
  theme_minimal() + theme(legend.position=c(.85,.15)) +
  labs(col="Cell type", x="False positive rate", y="True positive rate")
dev.off()

  
get_imp <- function(rf, top=5){
  i <- importance(rf)
  i <- data.frame(imp=i, feat=rownames(i))
  i <- i[order(i[,1], decreasing=T),][1:top,]
  rownames(i) <- NULL
  i
}

df <- get_imp(h1auc$rf)

impdf <- rbind(data.frame(get_imp(gmauc$rf), ct="GM12878"),
  data.frame(get_imp(h1auc$rf), ct="H1 hESC"),
  data.frame(get_imp(k5auc$rf), ct="K562"))

ggplot(impdf,aes(x=feat, y=MeanDecreaseGini, fill=ct)) + 
  geom_bar(stat="identity") + #facet_wrap(~ct) + #, scales="free") +
  coord_flip()

impdf

