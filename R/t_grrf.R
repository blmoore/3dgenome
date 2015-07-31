library("DAAG")
library("ggplot2")
library("randomForest")
library("RRF")

set.seed(42)
if(getwd() != "/Users/benmoore/hvl/production/3dgenome/") setwd("/Users/benmoore/hvl/production/3dgenome/")

# load data, 35 + full datasets
hdat <- readRDS("data/rds/H1hesc_35Vars.rds")
hfull <- readRDS("data/rds/H1hesc_71Vars.rds")
gdat <- readRDS("data/rds/Gm12878_35Vars.rds")
gfull <- readRDS("data/rds/Gm12878_115Vars.rds")
kdat <- readRDS("data/rds/K562_35Vars.rds")
kfull <- readRDS("data/rds/K562_187Vars.rds")

### simple BIC-penalised lm

# rows for picking subset
train <- sample(1:nrow(gfull), floor(nrow(gfull)/2))
test <- (1:nrow(gfull))[! 1:nrow(gfull) %in% train]
stopifnot(all(!train %in% test))

# Full linear models
gflm <- lm(eigen~., data = gfull[train,])
gcv_full <- cv.lm(gfull, gflm, m=10)
hflm <- lm(eigen~., data = hfull[train,])
hcv_full <- cv.lm(hfull, hflm, m=10)
kflm <- lm(eigen~., data = kfull[train,])
kcv_full <- cv.lm(kfull, kflm, m=10)

# 35 linear models
gcv_35 <- cv.lm(gdat, lm(eigen~., data=gdat), m=10)
hcv_35 <- cv.lm(hdat, lm(eigen~., data=hdat), m=10)
kcv_35 <- cv.lm(kdat, lm(eigen~., data=kdat), m=10)

# Regularised linear models
greg <- stepAIC(gflm, k=log(nrow(gfull)), direction="both")
gcv <- cv.lm(gfull[test,], greg, m=10)
hreg <- stepAIC(hflm, k=log(nrow(hfull)), direction="both")
hcv <- cv.lm(hfull[test,], hreg, m=10)
kreg <- stepAIC(kflm, k=log(nrow(kfull)), direction="both")
kcv <- cv.lm(kfull[test,], kreg, m=10)

message("GM :: from ", ncol(gfull)-1, " to ", length(attr(greg$terms, "term.labels")))
message("H1 :: from ", ncol(hfull)-1, " to ", length(attr(hreg$terms, "term.labels")))
message("K5 :: from ", ncol(kfull)-1, " to ", length(attr(kreg$terms, "term.labels")))

saveRDS(greg, "data/rds/Gm12878_reg.rds")
saveRDS(hreg, "data/rds/H1hesc_reg.rds")
saveRDS(kreg, "data/rds/K562_reg.rds")

greg <- readRDS("data/rds/Gm12878_reg.rds")
hreg <- readRDS("data/rds/H1hesc_reg.rds")
kreg <- readRDS("data/rds/K562_reg.rds")


## full models 
rgf <- readRDS("data/rds/Gm12878_115Vars_RFmod.rds")
rhf <- readRDS("data/rds/H1hesc_71Vars_RFmod.rds")
rkf <- readRDS("data/rds/K562_187Vars_RFmod.rds")

## 35 models 
gmod <- readRDS("data/rds/Gm12878_35Vars_RFmod.rds")
hmod <- readRDS("data/rds/H1hesc_35Vars_RFmod.rds")
kmod <- readRDS("data/rds/K562_35Vars_RFmod.rds")

reg_rf <- function(fulldf, reglm){
  terms <- attr(reglm$terms, "term.labels")
  rf <- randomForest(as.formula(paste0("eigen~", 
          paste(terms, collapse = "+"))), data=fulldf)
  stopifnot(length(terms) == length(attr(rf$terms, "term.labels")))
  cor(rf$predicted, fulldf$eigen)
}


gg <- reg_rf(gfull, greg)
hh <- reg_rf(hfull, hreg)
kk <- reg_rf(kfull, kreg)


## ct (plot)
#
# model    |
#  perf    |    
#          |
#          |
#          |
#          ___________________
#                 nvars
# ( col :: ct; point :: rf|lm )

km <- data.frame(ct=rep("K562", 6),
  model=c(rep("rf", 3), rep("lm", 3)),
  nvars=c(
    c(ncol(kfull)-1, 35, length(attr(kreg$terms, "term.labels"))), # rf
    c(ncol(kfull)-1, 35, length(attr(kreg$terms, "term.labels")))), # lm
  pcc=c(cor(rkf$predicted, kfull$eigen),# RF full
    cor(kmod$predicted, kdat$eigen),    # normal RF
    kk,                                 # regularised RF
    cor(kcv_full$cvpred, kfull$eigen),  # lm full
    cor(kcv_35$Predicted, kdat$eigen),  # lm normal (35)
    cor(kcv$cvpred, kfull$eigen)        # lm reg
  )
)

hm <- data.frame(ct=rep("H1 hESC", 6),
  model=c(rep("rf", 3), rep("lm", 3)),
  nvars=c(
    c(ncol(hfull)-1, 35, length(attr(hreg$terms, "term.labels"))), # rf
    c(ncol(hfull)-1, 35, length(attr(hreg$terms, "term.labels")))), # lm
  pcc=c(cor(rhf$predicted, hfull$eigen),# RF full
    cor(hmod$predicted, hdat$eigen),    # normal RF
    hh,                                 # regularised RF
    cor(hcv_full$cvpred, hfull$eigen),  # lm full
    cor(hcv_35$Predicted, hdat$eigen),  # lm normal (35)
    cor(hcv$cvpred, hfull$eigen)        # lm reg
  )
)

gm <- data.frame(ct=rep("GM12878", 6),
  model=c(rep("rf", 3), rep("lm", 3)),
  nvars=c(
    c(ncol(gfull)-1, 35, length(attr(greg$terms, "term.labels"))), # rf
    c(ncol(gfull)-1, 35, length(attr(greg$terms, "term.labels")))), # lm
  pcc=c(cor(rgf$predicted, gfull$eigen), # RF full
    cor(gmod$predicted, gdat$eigen),    # normal RF
    gg,                                 # regularised RF
    cor(gcv_full$cvpred, gfull$eigen),  # lm full
    cor(gcv_35$Predicted, gdat$eigen),  # lm normal (35)
    cor(gcv$cvpred, gfull$eigen)        # lm reg
  )
)

df <- rbind(km, hm, gm)

df


# RF uses lm BIC vars

ggplot(df, aes(x=nvars, y=pcc, col=ct)) +
  geom_point(aes(shape=model)) + geom_line(aes(group=interaction(ct, model))) +
  scale_y_continuous(limits=c(.7, .9))

ggplot(df, aes(x=as.factor(nvars), y=pcc, fill=model, group=ct)) +
  geom_bar(stat="identity", position = position_dodge(.8))

df


# LASSO regression
library("glmnet")
# reset seed
set.seed(42)

lasso <- function(full){
  netcv <- cv.glmnet(as.matrix(full[train,-1]), 
                     y=full$eigen[train], alpha=1)
  model <- glmnet(as.matrix(full[train, -1]), y=full[train,]$eigen, 
    alpha=1, lambda=netcv$lambda.1se)
  pred <- predict(model, newx=as.matrix(full[test, -1]),
                  s=netcv$lambda.1se, type="response")
  message("nvars : ", length(which(coef(model)[,1] != 0))-1)
  message("lasso : ", cor(pred, full$eigen[test]))
  # use terms in RF too
  terms <- names(coef(model)[which(coef(model) != 0),])[-1]
  rf <- randomForest(as.formula(paste0("eigen~", paste(terms, collapse="+"))),
    data=full)
  message("   rf : ", cor(rf$predicted, full$eigen))
}

# record 
lasso(gfull)
lasso(hfull)
lasso(kfull)

ggcv <- cv.glmnet(as.matrix(hfull[train,-1]), y=hfull[train,]$eigen, alpha=1)
par(mfrow=c(2,1))
plot(ggg, xvar="lambda", label=T)
plot(ggcv)
ggcv$lambda.1se

gg2 <- glmnet(as.matrix(hfull[train,-1]), y=hfull[train,]$eigen, alpha=1, 
  family="gaussian", lambda=ggcv$lambda.1se)

p <- predict(gg2, newx=as.matrix(hfull[test,-1]), s=ggcv$lambda.1se,
  type="response")
cor(p, hfull$eigen[test])

cs <- as.data.frame(t(as.matrix(ggg$beta)))
cs$lambda <- ggg$lambda
head(cs)
cs <- melt(cs, id.vars =  "lambda")

#kept / unkept var
chosen <- cs$lambda[which.min(abs(unique(cs$lambda) - ggcv$lambda.1se))]
kept <- as.character(unique(cs$variable)[which(subset(cs, lambda == chosen)$value != 0)])
cs$kept <- ifelse(cs$variable %in% kept, T, F)

labdf <- subset(cs, lambda == min(cs$lambda) & kept == T)
labdf$labs <- gsub("V04\\d*|Std0|Iggrab0|Ucd0", "",
  gsub("HaibTfbs|SydhTfbs|OpenChromChip", "" , 
     gsub("AlnRep?", "", gsub("BroadHistone", "", 
          gsub("H1hesc|Gm12878|K562", "", labdf$variable)))))

pdf("~/hvl/thesis_plots/gfull_coefTrace.pdf", 4, 4)
ggplot(cs, aes(x=-log10(lambda), y=value, col=kept, group=variable)) +
  scale_x_reverse(expand=c(0,0)) +
  scale_color_manual(values=c("grey70", "black")) +
  theme_minimal() + 
  geom_vline(xintercept=-log10(ggcv$lambda.1se),
    col=I("grey50"), linetype="dashed") +
  geom_line(aes(alpha=kept)) + 
  theme(legend.position="none") +
  scale_alpha_manual(values=c(.5, 1)) +
  labs(y="Standardised coefficient", 
    x=expression(-log[10](lambda))) +
  scale_y_continuous(breaks=c(.5, 0, -.5),
    labels=c(.5, "", -.5)) +
  geom_text(data=labdf, inherit.aes=F, hjust=0,
    aes(label = labs, x=5.5, y=value), size=I(2))
dev.off()

library("rCharts")

cs$log <- -log10(cs$lambda)

p <- nPlot(value~log, data=cs,
  type="lineChart", group="variable")

p
