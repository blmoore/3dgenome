library("ggplot2")
library("readr")
library("reshape2")

# 40kb boundaries written by score.sh

bounds <- read.table("/Volumes/BFD/hic_heatmaps/insulation/h1/all/chr_18.is520001.ids320001.insulation.boundaries",
  header=T)
mat <- read_tsv("/Volumes/BFD/hic_heatmaps/h1/all/40k/chr18_chr18.ic.tsv", col_names = F)

mat <- as.matrix(mat)

rownames(mat) <- colnames(mat)
contacts <- melt(mat)
head(contacts)

# bins as factors, now converting to int drops missing levels implicitly
contacts$n1 <- as.integer(contacts$Var1)
contacts$n2 <- as.integer(contacts$Var2)

cubroot_trans <- function() 
  scales::trans_new('cubroot', transform = function(x) x^(1/3), inverse = function(x) x^3 )

bounds$midpos <- rowMeans(bounds[,c("start", "end")])

pdf("figures/t_insScore.pdf", 8, 8)
p <- ggplot(contacts, aes(x=n1*.04, y=n2*.04, fill=value+1)) +
  geom_raster() + 
  scale_fill_gradient2(name="IF", trans=scales::log10_trans(),
    low="white", high="steelblue4", na.value="grey70", 
    breaks=c(1, 10, 100, floor(max(contacts$value+1)))) +
  theme_minimal() +
  labs(x="Position along chromosome 18 (Mb)", y="", fill="") +
  coord_fixed(xlim=c(20,55), ylim=c(20,55))

p + geom_point(data=bounds, inherit.aes=F,
  aes(x=midpos/1e6, y=midpos/1e6, size=insulationScore),
  col=I(rgb(.2,.2,.2,.5))) +
  scale_size_continuous(range=c(1, 5))
dev.off()

hist(bounds$insulationScore)
