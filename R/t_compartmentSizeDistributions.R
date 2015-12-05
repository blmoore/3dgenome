############### blockCharacterise.R ##################
# Calculate and visualise things like average length #
# of +ve / -ve blocks, compare between cell types.   #
######################################################

# Here be dragons: old code ressurected from 2012
library(plotrix)
library(sfsmisc)
library(TeachingDemos)
library("blmR")
library("RColorBrewer")

hdat <- readRDS("data/rds/H1hesc_35Vars.rds")
kdat <- readRDS("data/rds/K562_35Vars.rds")
gdat <- readRDS("data/rds/Gm12878_35Vars.rds")

eigs <- data.frame(`H1 hESC`=hdat$eigen,
                  GM12878=gdat$eigen,
                  K562=kdat$eigen)

states <- apply(eigs, 2, function(c) callStates(c)$state)
# 2 == A, 1 == B

cols <- brewer.pal(n=10, "Paired")
h.rle <- rle(states[,"H1.hESC"])
k.rle <- rle(states[,"K562"])
g.rle <- rle(states[,"GM12878"])

plotSizeDists <- function(){
  par(lwd=2, mfrow=c(1,3), mar=c(4,3,2,1), mgp=c(2,.6,0))

  ## GM12878
  plot(density(g.rle$lengths[g.rle$values==1], from=0), type="l", col=cols[1],
    xlab="",  main="GM12878", lwd=2, ylim=c(0, .2), xlim=c(0,50))
  points(density(g.rle$lengths[g.rle$values==2], from=0), type="l", col=cols[2])#, lty=2)
  text(27,.174, paste("median size = ",signif(median(g.rle$lengths),3), 
    " Mb\n mean = ",signif(mean(g.rle$lengths),3), 
    " ± ",signif(1.96*std.error(g.rle$lengths),2), " Mb", sep=""), col=cols[2])
  legend("top", horiz=T, col=c(cols[1], cols[2]), bty="n", box.lwd=.5, cex=.7, lwd=2,
    legend=c(paste("Open (n = ",length(g.rle$lengths[g.rle$values==1]), ")", sep=""),
      paste("Closed (n = ",length(g.rle$lengths[g.rle$values==2]), ")", sep="")))
    
  ## H1
  plot(density(h.rle$lengths[h.rle$values==1], from=0), type="l",
       xlab="Compartment size (Mb)", main="H1 hESC", col=cols[7], lwd=2, 
    ylim=c(0, .115), xlim=c(0, 50))
  points(density(h.rle$lengths[h.rle$values==2], from=0), 
         type="l", col=cols[8])
  text(27,.10, paste("median size = ",signif(median(h.rle$lengths),3), 
                     " Mb\n mean = ",signif(mean(h.rle$lengths),3), 
                     " ± ",signif(1.96*std.error(h.rle$lengths),2), " Mb", sep=""), col=cols[8])
  legend("top", horiz=T, col=c(cols[7], cols[8]), bty="n", box.lwd=.5, cex=.7, lwd=2,
         legend=c(paste("Open (n = ",length(h.rle$lengths[h.rle$values==1]), ")", sep=""),
                  paste("Closed (n = ",length(h.rle$lengths[h.rle$values==2]), ")", 
                        sep="")))
    
  ## K562
  plot(density(k.rle$lengths[k.rle$values==1], from=0), type="l", 
       col=cols[5], xlab="",  main="K562", lwd=2, ylim=c(0, .15),
    xlim=c(0,50))
  points(density(k.rle$lengths[k.rle$values==2], from=0), type="l", col=cols[6])
  text(27,.13, paste("median size = ",signif(median(k.rle$lengths),3), 
                     " Mb\n mean = ",signif(mean(k.rle$lengths),3), 
                     " ± ",signif(1.96*std.error(k.rle$lengths),2), " Mb", sep=""), col=cols[6])
  legend("top", horiz=T, col=c(cols[5], cols[6]), bty="n", box.lwd=.5, cex=.7, lwd=2,
         legend=c(paste("Open (n = ",length(k.rle$lengths[k.rle$values==1]), ")", sep=""),
                  paste("Closed (n = ",length(k.rle$lengths[k.rle$values==2]), ")", 
                        sep="")))

}

pdf("~/hvl/thesis_plots/comp_sizedists.pdf", 8, 3)
  plotSizeDists()
dev.off()

## TAD size distributions
read_tad <- function(file, celltype){
  tads <- read.table(file)
  tads$size <- tads$V3 - tads$V2
  tads$celltype <- celltype
  tads
}


gtads <- read_tad("data/bedfiles/gm_tads.bed", "GM12878")
htads <- read_tad("data/bedfiles/h1_tads.bed", "H1 hESC")
ktads <- read_tad("data/bedfiles/k5_tads.bed", "K562")

tads <- rbind(gtads, htads, ktads)

pdf("~/hvl/thesis_plots/tad_size_dist.pdf", 4.5, 3)
ggplot(tads, aes(x=size/1e6, fill=celltype)) + 
  geom_density(col=NA) +
  coord_cartesian(xlim=c(0, 5)) +
  theme_minimal() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c("#0000ff88", "#FFA50088", "#ff000088")) +
  labs(x="TAD size (Mb)", y="Density", fill="Cell type") +
  theme(legend.position=c(.85,.8))
dev.off()
  
