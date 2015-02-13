## split ctcf and yy1, i.e. is ctcf * yy1 more
## significantly enriched at boundaries than
## ctcf or yy1

all.t <- readRDS("data/rds/tad_boundary_features.rds")
# avoid melt
mall <- melt(all.t, id.vars=c("name", "bound", "feat", "ct", "type"))

mall[mall$feat == "Atf3",]$value <- subset(mall, feat == "Ctcf")$value * subset(mall, feat == "Yy1")$value

mall

all.t[all.t$feat=="Atf3", c(4:28)] <- all.t[all.t$feat=="Ctcf", c(4:28)] * all.t[all.t$feat=="Yy1", c(4:28)]

comp.t <- group_by(all.t, ct, type, feat) %>% 
  summarise(bound.mean = mean(X12),
            edge.mean  = mean(c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21)),
            p = wilcox.test(X12, c(X1, X2, X3, X4, X5, X25, X24, X23, X22, X21))$p.value)


ggplot(comp.t, aes(x=feat, y=-log10(p), col=type,
                 size=abs(bound.mean - edge.mean))) +
  facet_grid(ct~., scales="free_y") + 
  geom_point(position=position_dodge(.15)) + theme_bw() +
  labs(list(size="Absolute difference\nat boundary",
            y=expression(-log[10](italic(p))),
            x="", col="Boundary type")) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(range=c(2,7)) +
  guides(col=guide_legend(override.aes=list(size=4))) +
  geom_hline(yintercept=-log10(.01 / 297), linetype="dashed") +
  scale_colour_brewer(type="qual", palette=6) 