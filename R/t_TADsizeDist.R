library("ggplot2")
library("dplyr")

gt <- read.table("data/bedfiles/gm_tads.bed")
ht <- read.table("data/bedfiles/h1_tads.bed")
kt <- read.table("data/bedfiles/k5_tads.bed")

cts <- c("Gm12878", "H1 hESC", "K562")

lengths <- sapply(list(gt, ht, kt), function(b) with(b, V3 - V2))
names(lengths) <- cts

lengths <- sapply(cts, function(c) 
  t(rbind(`[[`(lengths, c), c))
)

df <- as.data.frame(do.call(rbind, lengths))
colnames(df) <- c("length", "ct")
df$length <- as.numeric(as.character(df$length))

mids <- df %>% group_by(ct) %>% 
  summarise(mid=median(length)/1e6)

pdf("figures/t_tadsizes.pdf", 5, 3)
ggplot(df, aes(x=length/1e6, fill=ct)) +
  geom_density(col=NA) +
  scale_x_log10() +
  theme_minimal() +
  geom_vline(data=mids, aes(xintercept=mid, col=ct)) +
  scale_fill_manual(values=c("#0000ff88", "#FFA50088", "#ff000088")) +
  scale_color_manual(values=c("#0000ff", "#FFA500", "#ff0000")) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position=c(.83,.8)) +
  labs(x="TAD size (Mb)", y="Frequency density",
    fill="Cell type")
dev.off()
