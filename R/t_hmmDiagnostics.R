library("blmR")

h.dat <- readRDS("data/rds/H1hesc_35Vars.rds")
g.dat <- readRDS("data/rds/Gm12878_35Vars.rds")
k.dat <- readRDS("data/rds/K562_35Vars.rds")

df <- data.frame(chr=factor(gsub("-.*", "", rownames(h.dat))), 
  pos=as.numeric(gsub(".*-", "", rownames(h.dat))),
  eig=h.dat$eigen, 
  states=factor(blmR::callStates(h.dat$eigen)$state))

library("ggplot2")

pdf("~/hvl/thesis_plots/hmm_calls.pdf", 9, 3)
ggplot(subset(df, chr == "chr8"),
  aes(x=pos/1e6, y=eig, group=NA)) +
  #facet_wrap(~chr, ncol=4, scales="free") + 
  geom_hline(yintercept=0, col=I("grey80")) + theme_minimal() +
  geom_line(col=I("grey50")) + geom_point(size=I(2.5), 
    aes(col=states)) + 
  labs(x="Genomic position on chromosome 8 (Mb)", 
    y="Compartment eigenvector",
    col="HMM state") +
  scale_color_brewer(type="qual", palette=4) +
  theme(legend.position="top")
dev.off()
