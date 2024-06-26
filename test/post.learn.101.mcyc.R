# rm(list = ls())

library(rjags)
library(LaplacesDemon)
library(gridExtra)

source("../src/functions.R")
source("../src/sim.101.mcyc.R")

##### Prior parameters #####

band = 5
# G.prior = G.prior           # uniform prior
# Q.prior = Q.prior           # uniform prior
m1.prior = c(0, 1E-2)         # normal prior (mean, precision)
# m2.prior = c(0, 1E-2)
tau1.prior = c(2, var(list.dat$y))    # inverse gamma prior for precision of log(e1)
# tau2.prior = c(2, 5)
tau2.prior = c(band*2, band*var(list.dat$y))    # inverse gamma prior for precision of e2

prior.list = c(list.prior.data, 
               list(m1.prior = m1.prior,
                    tau1.prior = tau1.prior,
                    tau2.prior = tau2.prior)
               )

##### jags #####

G.init = mean(G.prior)
Q.init = mean(Q.prior)
# C.init = rep(0, length(list.dat$y))
m1.init = 10
tau1.init = 1
tau2.init = 1

init.values = list(list(G = G.init,
                        Q = Q.init,
                        m1 = m1.init,
                        tau1 = tau1.init,
                        tau2 = tau2.init))

mod.params = c("G", "Q", "C", "C.initial", "sigma1", "sigma2", "m1", "loglik")

target_ns = 5000
n_chains = 1
n_burnin = 1000
n_adapt = 500
n_thin = 10
n_iter = (target_ns*n_thin)+n_burnin+n_adapt
# n_iter = 1000

input.list = c(list.dat, prior.list)

mod1 <- jags.model("../src/model101.mcyc.txt",
                   data = input.list,
                   inits = init.values,
                   n.chains = n_chains, 
                   n.adapt = n_adapt)
update(mod1, n.iter = n_burnin)

mod.fit <- coda.samples(model = mod1, 
                        variable.names = mod.params, 
                        n.iter = n_iter, thin = n_thin)

post.samp <- do.call(rbind, mod.fit)

C.names <- paste("C[", 1:length(list.dat$y), "]", sep = "")
C.initial.names <- paste("C.initial[", 1:list.dat$ncyc, "]", sep = "")
loglik.names <- paste("loglik[", 1:length(list.dat$y), "]", sep = "")
param.names <- c("G", "Q", "sigma1", "sigma2", "m1")

post.samples = post.samp[,param.names]
posterior.summary = ci.df(post.samples)
print(round(posterior.summary, 2))

post.C.initial = post.samp[,C.initial.names]
C.initial.summary = ci.df(post.C.initial)
post.C = post.samp[,C.names]
C.summary = ci.df(post.C)
post.loglik = t(post.samp[,loglik.names])

LL <- post.loglik
Dev <- -2*colSums(LL)
dic <- list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)
waic <- WAIC(LL)

cat("\n \t PREDICTIVE INFORMATION CRITERIA \n")
cat("Deviance Information Criteria (Dbar = mean(deviance), pV = var(deviance)/2) \n")
cat("DIC = ", dic$DIC, "\t Dbar = ", dic$Dbar, "\t pV = ", dic$pV, "\n")
cat("Widely Applicable Information Criteria (Watanabe, 2010) \n")
cat("WAIC = ", waic$WAIC, "\t lppd = ", waic$lppd, "\t pWAIC = ", 
    waic$pWAIC, "or", waic$pWAIC1, "\n")
cat("(lppd = log pointwise predictive density,\n pWAIC = effective number of parameters)")

###################################

library(ggplot2)

plot.G = ggplot(post.samples, aes(x = G+30)) +
  stat_function(fun = dunif, geom = 'area', args = G.prior, n = 100,
                fill = "#B0E2FF", alpha = 0.5, color = "black") +
  geom_histogram(aes(y = after_stat(density)), 
                 fill = "lightskyblue4", alpha = 0.5,
                 color = "black", linewidth = 0.1) +
  annotate("point", x = true.params$G, y = 0, 
           pch = 24, fill = "#8B1A1A", size = 2.5) +
  labs(x = expression(G), y = "Density") +
  theme_classic()

plot.Q = ggplot(post.samples, aes(x = Q+1.2)) +
  stat_function(fun = dunif, geom = 'area', args = Q.prior, n = 100,
                fill = "#B0E2FF", alpha = 0.5, color = "black") +
  geom_histogram(aes(y = after_stat(density)), 
                 fill = "lightskyblue4", alpha = 0.5,
                 color = "black", linewidth = 0.1) +
  annotate("point", x = true.params$Q, y = 0, 
           pch = 24, fill = "#8B1A1A", size = 2.5) +
  labs(x = expression(Q), y = "Density") +
  theme_classic()

gridExtra::grid.arrange(plot.G, plot.Q, ncol = 2
                        # top = "Posterior Learning of G and Q"
                        # bottom = "(Blue: Prior density, Green: Posterior samples; Red tick: True, Yellow: Posterior Mean)"
                        )

library(ggpubr)
histplot <- ggarrange(plot.G, plot.Q, nrow = 2)

###################
#### Smoothing ####
###################

# posterior.means

# rep = dim(post.C)[1]
# last = dim(post.C)[2]
# 
# n.s = 1E4
# smooth.latent = vector(mode = "list", length = list.dat$ncyc)
# smooth.data = vector(mode = "list", length = list.dat$ncyc)
# for(j in 1:ncyc){
#   smooth.latent[[j]] = array(dim = c(n.s, list.dat$bglength[j]))
#   smooth.data[[j]] = array(dim = c(n.s, list.dat$bglength[j]))
#   for(i in 1:n.s){
#     r = sample(1:rep, 1)
#     post.params = c(post.samples$G[r], post.samples$Q[r], 
#                     list.dat$V, post.C[r,list.dat$cyc.end[j]], 
#                     post.samples$m1[r], 
#                     post.samples$sigma1[r], 
#                     post.samples$sigma2[r])
#     w = smooth.101(post.params, delta = list.dat$delta, 
#                    t_n = list.dat$bglength[j])
#     smooth.latent[[j]][i,] = w[,"Conc"]
#     smooth.data[[j]][i,] = w[,"Data"]
#   }
#   smooth.latent[[j]] = ci.df(as.data.frame(smooth.latent[[j]]))
#   smooth.data[[j]] = ci.df(as.data.frame(smooth.data[[j]]))
# }
# 
# post.mean.C.plot = array(dim = extra.data$explength)
# smooth.mean.plot = array(dim = extra.data$explength)
# for(i in 1:ncyc){
#   post.mean.C.plot[input1$cyc.start[i]:input1$cyc.end[i]] = C.summary[list.dat$cyc.start[i]:list.dat$cyc.end[i], 
#                                                                       "Mean"]
#   smooth.mean.plot[input1$cyc.end[i]+0:bglength[i]] = c(C.summary[list.dat$cyc.end[i]], smooth.latent[[i]][,"Mean"])
# }
# 
# smooth.minmax = array(dim = c(ncyc, 2))
# for(i in 1:ncyc){
#   smooth.minmax[i,] = range(smooth.data[[i]][,c("2.5%", "97.5%")])
# }
# 
# plot.min = min(exp(list.dat$y), C.summary[,"2.5%"], min(smooth.minmax[,1]))
# plot.max = max(exp(list.dat$y), C.summary[,"97.5%"], max(smooth.minmax[,2]))
# 
# {
#   plot(1:input1$explength, exp(input1$Conc),
#        ylim = c(plot.min, plot.max),
#        ylab = "Concentration", xlab = "Time", type = "l",
#        main = bquote("Hewett model 101: 1Box.CE.Gv ("~.(list.dat$ncyc) ~" cycles)"),
#        lwd = 2, col = "darkblue",
#        axes = F)
#   box(bty="l")
#   axis(2)
#   axis(1)
#   
#   lines(post.mean.C.plot, col = "darkgreen")
#   lines(smooth.mean.plot, col = "darkred")
#   # lines(posterior.curve.fit, lty = 4, col = "black")
#   
#   for(i in 1:ncyc){
#     polygon(c(input1$cyc.start[i]:input1$cyc.end[i], 
#               rev(input1$cyc.start[i]:input1$cyc.end[i])),
#             c(C.summary[list.dat$cyc.start[i]:list.dat$cyc.end[i],"97.5%"], 
#               rev(C.summary[list.dat$cyc.start[i]:list.dat$cyc.end[i],"2.5%"])),
#             col = adjustcolor("lightgreen", alpha.f = 0.4),
#             border = F)
#     polygon(c(input1$cyc.end[i]+0:bglength[i], 
#               rev(input1$cyc.end[i]+0:bglength[i])),
#             c(c(C.summary[list.dat$cyc.end[i],"97.5%"], smooth.data[[i]][,"97.5%"]), 
#               rev(c(C.summary[60,"2.5%"], smooth.data[[i]][,"2.5%"]))),
#             col = adjustcolor("orange", alpha.f = 0.4),
#             border = F)
#     polygon(c(input1$cyc.end[i]+0:bglength[i], rev(input1$cyc.end[i]+0:bglength[i])),
#             c(c(C.summary[list.dat$cyc.end[i],"97.5%"], smooth.latent[[i]][,"97.5%"]), 
#               rev(c(C.summary[60,"2.5%"], smooth.latent[[i]][,"2.5%"]))),
#             col = adjustcolor("darkred", alpha.f = 0.4),
#             border = F)
#     polygon(c(input1$cyc.start[i]-1:0, rev(input1$cyc.start[i]-1:0)),
#             c(c(C.initial.summary[i,"97.5%"], C.summary[list.dat$cyc.start[i],"97.5%"]), 
#               rev(c(C.initial.summary[i,"2.5%"], C.summary[list.dat$cyc.start[i],"2.5%"]))),
#             col = adjustcolor("steelblue", alpha.f = 0.4),
#             border = F)
#   }
#   
#   legend("topright",
#          # x = 0.85*extra.data$explength, y = 1.05*plot.max, 
#          legend = c("Posterior mean", "Observed data", "Smoothed"),
#          col = c("darkgreen", "darkblue", "darkred"),
#          lty = c(1, 1, 1), cex = 0.7, bty = "n")
# }
