library("ggplot2")
library("reshape2")

hd <- readRDS("data/rds/H1hesc_35Vars.rds")
#hdat[,-1] <- scale(hdat[,-1])
hmelt <- melt(hd, id.var="eigen")
head(hmelt)

cors <- apply(hdat[,-1], 2, function(k) cor(k, hdat$eigen))
cors <- ifelse(cors > .1, "+ve", ifelse(cors < -.1, "-ve", "Weak"))

hmelt$cor <- cors[match(hmelt$variable, names(cors))]

pdf("~/hvl/thesis_plots/h1_feats.pdf", 8, 8)
ggplot(hmelt, aes(x=eigen, y=value)) +
  facet_wrap(~variable, scales="free_y") + 
  theme_minimal() +
  stat_density2d(geom="polygon", aes(fill = ..level..)) +
  scale_fill_gradient(low=rgb(254,240,217, max=255), 
    high=rgb(153,0,0, max=255)) +
  theme(strip.text.x=element_text(face=2),
    legend.position="none", panel.grid.minor=element_blank()) + 
  scale_x_continuous(limits=c(-.25, .25),
    breaks=c(-.15, 0, .15),
    labels=c("B", "0", "A"), expand=c(0,0)) +
  labs(x="Compartment eigenvector", y="Signal")
dev.off()

pdf("~/hvl/thesis_plots/h1_feats.pdf", 7, 7)
ggplot(hmelt, aes(x=eigen, y=value, fill=cor)) +
  facet_wrap(~variable, scales="free_y") + 
  theme_minimal() +
  stat_density2d(geom="polygon", aes(alpha = ..level..)) +
  scale_alpha(guide="none") +
  theme(strip.text.x=element_text(face=2),
    legend.title=element_text(size=9, face=1),
    legend.position=c(.92, .075), panel.grid=element_blank()) + 
  scale_x_continuous(limits=c(-.25, .25),
    breaks=c(-.15, 0, .15),
    labels=c("B", "0", "A"), expand=c(0,0)) +
  labs(x="Compartment eigenvector", y="Signal",
    fill="Correlation \nwith eigenvector")
dev.off()

## base version, colour +vely correlated on colour, negatve another



