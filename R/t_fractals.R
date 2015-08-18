# Compare hamiltonian path vs. hilbert curve 
# (i.e. illustration for equilibrium vs. fractalglobule in 2D)
library("TSP")
library("plotrix")
library("RColorBrewer")
set.seed(42)

## 1) Equilibrium globule
matsize <- 50 # i.e. 10 x 10

coords <- data.frame(x=rep(1:matsize, each=matsize), y=rep(1:matsize, matsize))
d <- dist(coords)
atsp <- as.ATSP(d)
t <- solve_TSP(atsp, method="nn")

t2 <- solve_TSP(atsp, method="two_opt", control = list(tour=t))
#t3 <- solve_TSP(atsp, method="two_opt", control = list(tour=t2))


## 2) Hilbert curve
library("HilbertVis")
h <- hilbertCurve(6)

pdf("~/hvl/thesis_plots/fractals.pdf", 4, 2)

par(mar=rep(0, 4), mfrow=c(1,2), oma=rep(.5, 4), col="white")
plot(coords, type="n", axes=F, frame=F, main="", asp=1)
color.scale.lines(coords[t2,],  col=colorRampPalette(brewer.pal(11, "Spectral"))(length(t)), lwd=5)

plot(h, type="n", axes=F, frame=F, main="", asp=1)
color.scale.lines(h,  col=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(h)), lwd=4.5)

dev.off()
