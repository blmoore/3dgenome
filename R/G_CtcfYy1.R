## split ctcf and yy1, i.e. is ctcf * yy1 more
## significantly enriched at boundaries than
## ctcf or yy1

all.t <- readRDS("data/rds/tad_boundary_features.rds")
# avoid melt
mall <- melt(all.t, id.vars=c("name", "bound", "feat", "ct", "type"))

mall[mall$feat == "Atf3",]$value <- subset(mall, feat == "Ctcf")$value * subset(mall, feat == "Yy1")$value

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

library("GenomicRanges")
library("rtracklayer")

read_peaks <- function(fname){
 # fname = "data/peaks//wgEncodeBroadHistoneGm12878CtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bam.regionPeak"
  tmp <- read.table(fname)
  with(tmp, GRanges(seqnames=V1, ranges=IRanges(start=V2, end=V3)))
}
  
## read peak calls, intersect, then boundary profile each
ctcf_g <- read_peaks("data/peaks/wgEncodeBroadHistoneGm12878CtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bam.regionPeak")
yy1_g <- read_peaks("data/peaks/wgEncodeHaibTfbsGm12891Yy1sc281V0416101AlnRep0.bam_VS_wgEncodeHaibTfbsGm12891RxlchPcr1xAlnRep0.bam.regionPeak")

both_g <- subsetByOverlaps(ctcf_g, yy1_g)

# now intersect with TAD boundaries:
tad_g <- read_peaks("data/bedfiles/gm_tadbounds.bed")
marked_tads_g <- subsetByOverlaps(tad_g, both_g)
# ~12% of TAD boundaries contain both::
# 100*length(marked_tads_g) / length(tad_g)

ctcf_tads_g <- subsetByOverlaps(tad_g, ctcf_g)
# ~61% of TAD boundaries (40kb) have signif CTCF peak::
# 100 * length(ctcf_tads_g) / length(tad_g)


ctcf_h <- read.table("data/peaks/wgEncodeBroadHistoneH1hescCtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneH1hescControlStdAlnRep0.bam.regionPeak")
ctcf_k <- read.table("data/peaks/wgEncodeBroadHistoneK562CtcfStdAlnRep0.bam_VS_wgEncodeBroadHistoneK562ControlStdAlnRep1.bam.regionPeak")
