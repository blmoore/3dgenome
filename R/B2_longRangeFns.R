library("reshape2")
library("magrittr")
set.seed(42)

# symmetric matrix
a <- matrix(rpois(1e4, 2), 100)
a[upper.tri(a)] <- t(a)[upper.tri(a)]

# draw
par(mar=rep(1,4))
image(x=1:100, y=1:100, a, asp=1, frame=F, axes=F)

breaks <- c(1, 12, 14, 25, 60, 71, 89, 100)

abline(h=breaks-.5, lwd=2, col="white")
abline(v=breaks-.5, lwd=2, col="white")

b <- melt(a)

# generate groups from breaks
bb <- with(b, tapply(value, list(
  y=cut(Var1, breaks=breaks, include.lowest=T),
  x=cut(Var2, breaks=breaks, include.lowest=T)),
  sum)
)

bsq <- melt(bb)

getNum <- . %>%
  # rm brackets
  gsub("\\[|\\(|\\]|\\)", "", .) %>%
  # split digits and convert
  strsplit(",") %>%
  #xleft xright (?)
  unlist %>% as.numeric


y <- t(sapply(bsq[,1], getNum))
x <- t(sapply(bsq[,2], getNum))

bsq$size <- (y[,2] - y[,1]) * (x[,2] - x[,1])
bsq$norm <- bsq$value / bsq$size

plot(1:100, 1:100, type="n", frame=F, axes=F)
rect(ybottom=y[,1], ytop=y[,2],
     xleft=x[,1], xright=x[,2], 
     col=rgb(colorRamp(c("white", "steelblue4"))(bsq$norm / max(bsq$norm)), 
             alpha=255*(bsq$norm / max(bsq$norm)), max=255),
     border="white")

