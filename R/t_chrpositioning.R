library("ggplot2")
library("blmR")
## nuclear positioning (Boyle et al., 2001)
## Figure 3

edge <- c(7, 3, 13, 2, 8, 18, 11, 4)
near_edge <- c(12, 9, "X", "Y", 20)
no_bias <- c(5, 21, 15, 10, 6, 14)
central <- c(1, 17, 16, 22, 19)

pos <- data.frame(loc=c(
  rep("Edge", length(edge)),
  rep("near_edge", length(near_edge)),
  rep("no_bias", length(no_bias)),
  rep("Central", length(central))),
  chr=c(edge, near_edge, no_bias, central))
  
hdat <- readRDS("data/rds/H1hesc_35Vars.rds")
kdat <- readRDS("data/rds/K562_35Vars.rds")
gdat <- readRDS("data/rds/Gm12878_35Vars.rds")

chrs <- gsub("chr(.*?)-.*", "\\1", rownames(hdat))
locales <- data.frame(pos=rep(pos[match(chrs, pos$chr),1], 3),
  ct=c(
      rep("GM12878", length(chrs)), 
      rep("H1 hESC", length(chrs)), 
      rep("K562", length(chrs))
    ),
  eig=c(gdat$eigen, hdat$eigen, kdat$eigen))
  

locales$pos <- factor(locales$pos, levels=c("Edge", "no_bias", "Central"))

pdf("~/hvl/thesis_plots/nucpos.pdf", 6, 2)
ggplot(subset(locales, pos %in% c("Edge", "Central")),
  aes(x=eig, fill=pos)) +
  facet_wrap(~ct) + geom_density(alpha=I(.4), position = position_dodge()) +
  scale_fill_brewer(type="qual") +
  scale_y_continuous(expand=c(0,0)) +
  theme_minimal() +
  theme(axis.ticks.y=element_blank(), 
    axis.text.y=element_blank(),
    panel.grid=element_blank(),
    strip.text.x=element_text(face=2),
    legend.key=element_rect(colour="black")) +
  labs(y=NULL, fill="Nuclear\nposition", 
    x="Compartment eigenvector") +
  guides(fill = guide_legend(override.aes = list(colour = NULL)))
dev.off()

locales$comp <- as.factor(c(
  callStates(gdat$eigen)$state,
  callStates(hdat$eigen)$state,
  callStates(kdat$eigen)$state))

ggplot(subset(locales, pos %in% c("edge", "central")), aes(x=comp, fill=pos)) + 
  geom_histogram(position = position_fill()) + facet_wrap(~ct)

# h1
ks.test(
  subset(locales, ct == "H1 hESC" & pos == "Edge")$eig,
  subset(locales, ct == "H1 hESC" & pos == "Central")$eig)

# gm12878  
ks.test(
  subset(locales, ct == "GM12878" & pos == "Edge")$eig,
  subset(locales, ct == "GM12878" & pos == "Central")$eig)

# k562
ks.test(
  subset(locales, ct == "K562" & pos == "Edge")$eig,
  subset(locales, ct == "K562" & pos == "Central")$eig)
