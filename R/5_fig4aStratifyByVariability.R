################## stratVar.R #######################
# Stratify blocks by amount of variation in eigs    #
# and model the subsets, compare results between    #
# cell-types.                                       #
#####################################################
library("caret")
library("randomForest")
library("plotrix")
library("blmR")
options(scipen=999)

# load previously-built models:
g.mod <- readRDS("data/rds/Gm12878_35Vars_RFmod.rds")
h.mod <- readRDS("data/rds/H1hesc_35Vars_RFmod.rds")
k.mod <- readRDS("data/rds/K562_35Vars_RFmod.rds")
# load featureset data.frames
g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

justEigs <- cbind(h.dat$eigen, g.dat$eigen, k.dat$eigen)
justEigs <- cbind(justEigs, apply(justEigs,1,sd))
rownames(justEigs) <- rownames(h.dat)
justEigs[,4] <- apply(justEigs[,1:3], 1, function(x) mad(x))
orderedEigs <- justEigs[order(justEigs[,4], decreasing=T),]

upper <- orderedEigs[1:round(quantile(1:nrow(justEigs), .33)),]
mid <- orderedEigs[round(quantile(1:nrow(justEigs), .34)):round(quantile(1:nrow(justEigs), .66)),]
low <- orderedEigs[round(quantile(1:nrow(justEigs), .67)):nrow(justEigs),]

# Build random forest models for each subset of the genome.
# This can take a while to run, could be parallelised (either
# with sprint's prandomforest / per model OR embarassingly
# parallel out each blmR::modelEigens call.
# Alternatively, skip these lines and call:
#   sm <- readRDS("data/rds/stratvar_df.rds")
# before plotting.

k <- k.dat
k5.all <- readRDS("data/rds/K562_35Vars_RFmod.rds")
k5.upper <- modelEigens(k[rownames(k) %in% rownames(upper),])
k5.mid <- modelEigens(k[rownames(k) %in% rownames(mid),])
k5.low <- modelEigens(k[rownames(k) %in% rownames(low),])

g <- g.dat
gm.all <- readRDS("data/rds/Gm12878_35Vars_RFmod.rds")
gm.upper <- modelEigens(g[rownames(g) %in% rownames(upper),])
gm.mid <- modelEigens(g[rownames(g) %in% rownames(mid),])
gm.low <- modelEigens(g[rownames(g) %in% rownames(low),])

h <- h.dat
h1.all <- readRDS("data/rds/H1hesc_35Vars_RFmod.rds")
h1.upper <- modelEigens(h[rownames(h) %in% rownames(upper),])
h1.mid <- modelEigens(h[rownames(h) %in% rownames(mid),])
h1.low <- modelEigens(h[rownames(h) %in% rownames(low),])

df <- rbind(cbind("H1", "low", h1.low$f.cors),
            cbind("H1", "mid", h1.mid$f.cors),
            cbind("H1", "high", h1.upper$f.cors),
            cbind("GM12878", "low", gm.low$f.cors),
            cbind("GM12878", "mid", gm.mid$f.cors),
            cbind("GM12878", "high", gm.upper$f.cors),
            cbind("K562", "low", k5.low$f.cors),
            cbind("K562", "mid", k5.mid$f.cors),
            cbind("K562", "high", k5.upper$f.cors))

colnames(df) <- c("ct", "variability", "PCC")
df <- as.data.frame(df)
df$ct <- as.factor(df$ct)
df$variability <- factor(df$var, levels=c("high", "mid", "low"))
df$PCC <- as.numeric(as.character(df$PCC))
sm <- summarySE(df, measure="PCC", groupvars=c("ct", "variability"))

saveRDS(sm, "data/rds/stratvar_df.rds")
sm <- readRDS("data/rds/stratvar_df.rds")

meandf <- data.frame(ct  = unique(sm$ct),
                     PCC = c(cor(g.mod$predicted, g.dat$eigen),
                             cor(h.mod$predicted, h.dat$eigen),
                             cor(k.mod$predicted, k.dat$eigen)))

pdf("figures/f4a_stratByVar.pdf", 6, 4.5)
ggplot(sm, aes(x=ct, y=PCC, group=variability)) + 
  geom_bar(stat="identity", position="dodge", aes(fill=variability)) +
  scale_fill_brewer(palette="Oranges") + theme_bw() + 
  geom_errorbar(aes(ymin=PCC-ci, ymax=PCC+ci),
                width=.2, position=position_dodge(.9)) +
  scale_y_continuous(breaks=seq(0,1, by=.1), limits=c(0,.925)) +
  labs(list(y="Model accuracy (PCC)", x="", fill="Variability")) +
  # geom_segment for each overall model acc
  geom_segment(aes(x=.5, xend=1.5, y=cor(g.mod$predicted, g.dat$eigen), 
                   yend=cor(g.mod$predicted, g.dat$eigen)), linetype=3) +
  geom_segment(aes(x=1.5, xend=2.5, y=cor(h.mod$predicted, h.dat$eigen), 
                   yend=cor(h.mod$predicted, h.dat$eigen)), linetype=3) +
  geom_segment(aes(x=2.5, xend=3.5, y=cor(k.mod$predicted, k.dat$eigen), 
                   yend=cor(k.mod$predicted, k.dat$eigen)), linetype=3)
dev.off()

sciNotation <- function(x, digits = 1) { 
  # R mailing list: http://r.789695.n4.nabble.com/Scientific-notation-in-plots-td791499.html
  if (length(x) > 1) { 
    return(append(sciNotation(x[1]), sciNotation(x[-1]))) 
  } 
  if (!x) return(0) 
  exponent <- floor(log10(x)) 
  base <- round(x / 10^exponent, digits) 
  as.expression(substitute(base %*% 10^exponent, 
                           list(base = base, exponent = exponent))) 
} 

sciNotation(t.test(gm.upper$f.cors, gm.low$f.cors)$p.value)
sciNotation(t.test(h1.upper$f.cors, h1.low$f.cors)$p.value)
sciNotation(t.test(k5.upper$f.cors, k5.low$f.cors)$p.value)

### END