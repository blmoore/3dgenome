library("ggplot2")
library("reshape2")

# from 7_fig5boundaryEnrichments.R
all.t <- readRDS("data/rds/tad_boundary_features.rds")

ctcf <- subset(all.t, feat == "Ctcf")
plot(ctcf$X12, pch=20, cex=.1, col="darkred")
points(ctcf$X1, pch=20, cex=.1, col="darkblue")

mctcf <- melt(ctcf, id.vars=c("name", "bound", "feat", "ct", "type"))

pdf("~/hvl/ice/plots/ctcf_over_tad_bounds.pdf", 8, 3)
ggplot(#subset(mctcf, variable %in% c("X1", "X12")),
  mctcf, aes(x=variable, y=value, group=bound)) + #geom_violin() +
  geom_line(alpha=.025) + xlab("") + ylab("ctcf") +
  facet_wrap(~ct) + theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())
dev.off()

hist(ctcf$X12 - ctcf$X1)
