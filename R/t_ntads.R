library("dplyr")
library("ggplot2")

ht <- read.table("data/bedfiles/h1_tads.bed")
gt <- read.table("data/bedfiles/gm_tads.bed")
kt <- read.table("data/bedfiles/k5_tads.bed")

ht$ct <- "H1 hESC"
gt$ct <- "GM12878"
kt$ct <- "K562"

df <- rbind(ht, gt, kt)
df <- group_by(df, ct, V1) %>% tally

colnames(df) <- c("ct", "chr", "ntads")
df$chr <- factor(gsub("chr", "", as.character(df$chr)), levels=1:22)

pdf("~/hvl/thesis_plots/ntads.pdf", 8, 4)
ggplot(df, aes(x=chr, y=ntads, fill=ct)) +
  geom_bar(stat="identity", position = position_dodge()) +
  scale_fill_manual(values=c("#0000ff98", "#FFA50098", "#ff000098"))  + 
  theme_minimal() +
  scale_y_continuous(expand=c(0,0)) +
  labs(fill="Cell type", y="Number of TADs", x="Chromosome") +
  theme(legend.position=c(.93, .85))
dev.off()
  
df <- rbind(ht, gt, kt)
group_by(df, ct) %>% tally
