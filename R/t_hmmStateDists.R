### Draw HMM state distributions, aking to HmmStatesGraphicDiagnostics.pdf
###   as produced by hvl/R/hmmModel.R
library("ggplot2")  
library("RColorBrewer")

# blmR::callStates edit
state_hmm <- function(indat){
  require("RHmm")
  set.seed(42)
  hmm <- HMMFit(indat, dis = "NORMAL", nStates = 2, control = list(init = "RANDOM", 
    iter = 5000))
  hmm
}

gdat <- readRDS("data/rds/Gm12878_35Vars.rds")
hdat <- readRDS("data/rds/H1hesc_35Vars.rds")  
kdat <- readRDS("data/rds/K562_35Vars.rds")

hhmm <- state_hmm(hdat$eigen)
ghmm <- state_hmm(gdat$eigen)
khmm <- state_hmm(kdat$eigen)

sim_state <- function(state, hmm, n=5e5){
  s <- list(mean=hmm$HMM$distribution$mean[[state]],
    var=hmm$HMM$distribution$var[[state]])
  rnorm(n, s$mean, sqrt(s$var))
}

theoretical <- function(hmm, ct){
  s1 <- sim_state(1, hmm)
  s2 <- sim_state(2, hmm)
  rbind(data.frame(value=s1, type="theo", state="1", cell=ct),
    data.frame(value=s2, type="theo", state="2", cell=ct))
}

g <- theoretical(ghmm, "GM12878")
h <- theoretical(hhmm, "H1 hESC")
k <- theoretical(khmm, "K562")  

alldf <- rbind(g, h, k)  
  
pal <- brewer.pal(9, "Paired")

pdf("~/hvl/thesis_plots/hmmDists.pdf", 3, 3.5)
ggplot(alldf, aes(x=value, fill=interaction(state, cell), group=state)) + 
  geom_density(alpha=I(.6)) +
  facet_grid(cell~.) + #, scales="free_x") +
  scale_fill_manual(values=pal[c(1,2,7,8,5,6)]) +
  theme_minimal() +
  theme(legend.position="none",
    axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
  labs(x="Compartment eigenvector",
    y="Density", fill="HMM state") #+
  #scale_y_continuous(expand=c(0,0))
dev.off()


