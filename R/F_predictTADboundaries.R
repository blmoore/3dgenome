library("ggplot2")
library("reshape2")

rm(list=ls())
# from 7_fig5boundaryEnrichments.R
all.t <- readRDS("data/rds/tad_boundary_features.rds")

# ctcf <- subset(all.t, feat == "Ctcf")
# plot(ctcf$X12, pch=20, cex=.1, col="darkred")
# points(ctcf$X1, pch=20, cex=.1, col="darkblue")
# 
# mctcf <- melt(ctcf, id.vars=c("name", "bound", "feat", "ct", "type"))
# 
# pdf("~/hvl/ice/plots/ctcf_eover_tad_bounds.pdf", 8, 3)
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

taddf$variable <- as.factor(ifelse(taddf$variable == T, 1, 0))

# Import ALU counts from t_repeatClasses.R script
taddf$Alu <- rptdf$Alu

#saveRDS(taddf, "data/rds/tadpred.rds")


library("AUCRF")
library("ROCR")
set.seed(42)

ctrf_auroc <- function(cell, data=taddf){
  #cell <- match.arg(cell)
  h1df <- subset(taddf, ct == cell)[,-c(1:3)]
  train <- sample(1:nrow(h1df), (nrow(h1df)/10) * 8)
  message("training : ", length(train), " out of ", nrow(h1df), " total rows")
  
  auc_class <- AUCRF(variable ~ ., data=h1df[train,], ntrees=200, pdel=0.1)
  opt_rf <- randomForest(as.formula(paste0("variable~", paste(as.vector(OptimalSet(auc_class)$Name), collapse="+"))),
    data=h1df[train,])

  preds <- predict(opt_rf, newdata=h1df[-train,], type="prob")
  message("predictions : ", length(preds[,2]))
  
  print(head(preds))
  
  #perf <- performance(  prediction(opt_rf$votes[,2], h1df$variable), "tpr", "fpr")
  perf <- performance(prediction(preds[,2], h1df[-train,]$variable), "tpr", "fpr")
  df <- data.frame(fp=unlist(perf@x.values), tp=unlist(perf@y.values), ct=cell)
  obj <- list()
  
  obj$data <- h1df
  obj$train_ind <- train
  obj$predictions <- preds
  obj$auroc <- df
  obj$aurcrf <- auc_class
  obj$rf <- opt_rf
  return(obj)
}


h1auc <- ctrf_auroc(cell="H1hesc")
h1auc$auroc$ct <- "H1 hESC"
gmauc <- ctrf_auroc(cell="Gm12878")
gmauc$auroc$ct <- "GM12878"
k5auc <- ctrf_auroc(cell="K562")

saveRDS(h1auc, "data/rds/h1_tadpred.rds")
saveRDS(gmauc, "data/rds/gm_tadpred.rds")
saveRDS(k5auc, "data/rds/k5_tadpred.rds")
gmauc <- readRDS("data/rds/gm_tadpred.rds")
h1auc <- readRDS("data/rds/h1_tadpred.rds")
k5auc <- readRDS("data/rds/k5_tadpred.rds")

auc <- rbind(h1auc$auroc, gmauc$auroc, k5auc$auroc)

# actual AUROC values
get_auroc <- function(aucobj)
  signif(unlist(performance(prediction(aucobj$predictions[,2], aucobj$data[-aucobj$train_ind,]$variable), "auc")@y.values),3)

pdf("~/hvl/thesis_plots/tad_bounds_auroc_v3.pdf", 4.4, 4.4)
ggplot(auc, aes(x=fp, y=tp, col=ct)) + 
  annotate("segment", x=0, xend=1, y=0, yend=1, linetype="dashed", col=I("grey30")) +
  geom_line(size=I(1.1)) +
  scale_x_continuous(expand=c(0,0)) + coord_fixed() +
  scale_y_continuous(expand=c(0,0)) +
  theme_minimal() + theme(legend.position=c(.15,.86)) +
  labs(col="Cell type", x="False positive rate", y="True positive rate") +
  scale_color_manual(values=c("#0000ff98", "#FFA50098", "#ff000098")) +
  # add AUC values
  annotate("text", y = .24, x=.88, label="AUC", fontface=2, size=4.5) +
  annotate("text", y = .06, x=.88, label=get_auroc(gmauc), colour="#0000ff", size=4.5) +
  annotate("text", y = .18, x=.88, label=get_auroc(h1auc), colour="#FFA500", size=4.5) +
  annotate("text", y = .12, x=.88, label=get_auroc(k5auc), colour="#ff0000", size=4.5)
dev.off()

get_auroc(gmauc)
get_auroc(h1auc)
get_auroc(k5auc)
  
# Variable importance rankings
get_imp <- function(rf, top=10){
  i <- importance(rf)
  i <- data.frame(imp=i, feat=rownames(i))
  i <- i[order(i[,1], decreasing=T),][1:top,]
  rownames(i) <- NULL
  i
}

gp <- get_imp(gmauc$rf)
hp <- get_imp(h1auc$rf)
kp <- get_imp(k5auc$rf)

plot_imp <- function(imp, header){
  plot(y=1:length(imp$feat), x=rev(imp$MeanDecreaseGini), axes=F,
    xlab="", ylab="", type="n", main=header)
  axis(1)
  axis(2, labels = rev(imp$feat), at=1:length(imp$feat), las=1, tick=F)
  abline(h=1:length(imp$feat), lty=2, col="grey70")
  points(y=1:length(imp$feat), x=rev(imp$MeanDecreaseGini), pch=20, cex=1.2)
}

dev.off()
pdf("~/hvl/thesis_plots/tadpred_varimp.pdf", 7, 3)
#plot.new()
par(mar=c(3,6,2,2), mfrow=c(1,3), oma=c(1,0,0,0), mgp=c(.3,.4,0))
plot_imp(gp, header="GM12878")
plot_imp(hp, header="H1 hESC")
plot_imp(kp, header = "K562")
mtext("Variable importance (mean decrease in Gini coefficient)", side=1, outer=T)
dev.off()

