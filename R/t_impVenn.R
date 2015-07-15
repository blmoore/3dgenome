## venn of varimp + probabilities
library("randomForest")
library("VennDiagram")

gmod <- readRDS("~/hvl/production/3dgenome/data/rds/Gm12878_35Vars_RFmod.rds")
hmod <- readRDS("~/hvl/production/3dgenome/data/rds/H1hesc_35Vars_RFmod.rds")
kmod <- readRDS("~/hvl/production/3dgenome/data/rds/K562_35Vars_RFmod.rds")

top_ten <- function(mod){
  imp <- importance(mod, type=1)
  names(imp[order(imp, decreasing=T),])[1:10]
}

g10 <- top_ten(gmod)
h10 <- top_ten(hmod)
k10 <- top_ten(kmod)

rownames(importance(hmod, type=1)) == rownames(importance(gmod, type=1)) == rownames(importance(kmod, type=1))

v <- venn.diagram(x=list(K562=k10, GM12878=g10, H1=h10), 
  filename=NULL, type = "ChowRuskey", doWeights=F,
  alpha=0.4, col=c("red","blue","orange"),
  fill=c("red","blue","orange"))

#svg("~/hvl/thesis_plots/top10_venn.svg", 4, 4)
dev.new()
grid::grid.draw(v)
dev.off()

# permutations, how many do we expect in each?

all <- rownames(hmod$importance)

## all 3 intersection
all3 <- replicate(10000, length(intersect(sample(all, 10),
  intersect(sample(all, 10), sample(all, 10)))))
counts <- table(all3)

1 - sum(counts[rownames(counts) == 0]) / sum(counts)

## any two (three ways of doing this)
any2 <- replicate(10000, {
          a <- sample(all, 10)
          b <- sample(all, 10)
          c <- sample(all, 10)
          aVbc <- length(intersect(a, union(b, c))) - length(intersect(a, intersect(b, c)))
          bVc <- length(intersect(b, c)) - length(intersect(a, intersect(b, c)))
          aVbc + bVc
})
   
counts2 <- table(res)
sum(counts2[rownames(counts2) >= 7]) / sum(counts2)
       

library("BiasedUrn")
vignette("UrnTheory")

