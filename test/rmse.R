rm(list = ls())

library(rjags)
library(ggplot2)
library(LaplacesDemon)
library(gridExtra)
library(tidyverse)

source("../src/functions.R")
source("../src/sim.111.mcyc.R")

# e_priorseq = c(0.45, 0.4, 0.3, 0.25, 0.2)
# nrep = 25
e_priorseq = c(0.45, 0.4, 0.3, 0.25, 0.2)
nrep = 50
dev_eL <- vector(mode = "list", length = length(e_priorseq))
dev_eLF <- vector(mode = "list", length = length(e_priorseq))
dev_Q.prime <- vector(mode = "list", length = length(e_priorseq))

for(i in 1:length(e_priorseq)){
  
  ##### Prior parameters #####
  e_prior = e_priorseq[i]
  
  band = 5                                # tuning width of uncertainty band
  G.prior = list.prior.data$G.prior       # uniform prior
  Q.prior = list.prior.data$Q.prior       # uniform prior
  m1.prior = c(0, 1E-2)                   # normal prior (mean, precision)
  tau1.prior = c(2, var(list.dat$y))      # IG prior for precision of log(e1)
  tau2.prior = c(band*2, band*var(list.dat$y))  # IG prior for precision of e2
  QL.prior = c(2, 10)
  QR.prior = c(2, 10)
  eL.prior = c(max(0, 0.5 - e_prior), min(1, 0.5 + e_prior))
  eLF.prior = c(max(0, 0.5 - e_prior), min(1, 0.5 + e_prior))
  eRF.prior = c(0.6, 1)
  
  for(j in 1:nrep){
    
    ########## jags ###########
    
    source("post_learn_111.R")
    
    ###########################
    
    dev_eL[[i]][j] <- abs(true.params$e.L - median(post.samp[, "eL"]))
    dev_eLF[[i]][j] <- abs(true.params$e.LF - median(post.samp[, "eLF"]))
    dev_Q.prime[[i]][j] <- abs(true.params$Q.prime - median(post.samp[, "Q.prime"]))
  }
  
}

dev_eL2 <- do.call(cbind, dev_eL)
colnames(dev_eL2) <- paste("p", 1:length(e_priorseq), sep = "")
dev_eL_long <- reshape2::melt(dev_eL2)
dev_eL_long <- dev_eL_long[, -1]
names(dev_eL_long) <- c("prior", "value")
dev_eL_long$par <- "eL"

dev_eLF2 <- do.call(cbind, dev_eLF)
colnames(dev_eLF2) <- paste("p", 1:length(e_priorseq), sep = "")
dev_eLF_long <- reshape2::melt(dev_eLF2)
dev_eLF_long <- dev_eLF_long[, -1]
names(dev_eLF_long) <- c("prior", "value")
dev_eLF_long$par <- "eLF"

dev_Q.prime2 <- do.call(cbind, dev_Q.prime)
colnames(dev_Q.prime2) <- paste("p", 1:length(e_priorseq), sep = "")
dev_Q.prime_long <- reshape2::melt(dev_Q.prime2)
dev_Q.prime_long <- dev_Q.prime_long[, -1]
names(dev_Q.prime_long) <- c("prior", "value")
dev_Q.prime_long$par <- "Q.prime"

dev_long <- rbind(dev_Q.prime_long, dev_eLF_long)

par_names <- list(
  # 'eL' = latex2exp::TeX("$\\epsilon_{L}$"),
  'Q.prime' = latex2exp::TeX("$Q'$"),
  'eLF' = latex2exp::TeX("$\\epsilon_{LF}$")
)

par_labeller <- function(variable,value){
  return(par_names[value])
}

plot1 <- ggplot(dev_long, aes(x = prior, y = value)) +
  geom_boxplot(color = "lightskyblue4", fill = "#B0E2FF",
               outlier.shape = NA) +
  facet_wrap(~ par,  ncol = 1, strip.position = "right",
             labeller = par_labeller, scales = "free") +
  ylab("Deviation from true value") +
  xlab("Choice of prior") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size = 12),
        strip.text.y.right = element_text(angle = 0))

source("post.learn.101.mcyc.R")

ggarrange(histplot, plot1, nrow = 1,
          labels = c("(a)", "(b)"), font.label = list(face = "plain"),
          vjust = 1.0, hjust = 0.0,
          widths = c(1, 1.5))
ggsave("../fig/fig3.jpg", units = "in",
       width = 6.7, height = 3.35, dpi=600)
