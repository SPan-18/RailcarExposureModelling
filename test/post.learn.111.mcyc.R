rm(list = ls())

library(rjags)
library(ggplot2)
library(LaplacesDemon)
library(gridExtra)

source("../src/functions.R")
source("../src/sim.111.mcyc.R")

##### Prior parameters #####

band = 5                                # tuning width of uncertainty band
G.prior = list.prior.data$G.prior       # uniform prior
Q.prior = list.prior.data$Q.prior       # uniform prior
m1.prior = c(0, 1E-2)                   # normal prior (mean, precision)
tau1.prior = c(2, var(list.dat$y))      # IG prior for precision of log(e1)
tau2.prior = c(band*2, band*var(list.dat$y))  # IG prior for precision of e2
QL.prior = c(2,10)
QR.prior = c(2,10)
eL.prior = c(0.3, 0.7)
eLF.prior = c(0.3, 0.7)
eRF.prior = c(0.6, 1)

prior.list = c(list.prior.data, 
               list(m1.prior = m1.prior,
                    tau1.prior = tau1.prior,
                    tau2.prior = tau2.prior,
                    QL.prior = QL.prior,
                    QR.prior = QR.prior,
                    eL.prior = eL.prior,
                    eLF.prior = eLF.prior,
                    eRF.prior = eRF.prior))

##### jags #####

G.init = mean(G.prior)
Q.init = mean(Q.prior)
QL.init = mean(QL.prior)
QR.init = mean(QL.prior)
eL.init = mean(eL.prior)
eLF.init = mean(eLF.prior)
eRF.init = mean(eRF.prior)
m1.init = 0
tau1.init = 1
tau2.init = 1

init.values = list(list(G = G.init,
                        Q = Q.init,
                        QL = QL.init,
                        QR = QR.init,
                        eL = eL.init,
                        eLF = eLF.init,
                        eRF = eRF.init,
                        m1 = m1.init,
                        tau1 = tau1.init,
                        tau2 = tau2.init))

mod.params = c("G", "Q", "QL", "QR", "eL", "eLF", "eRF", "G.prime", "Q.prime",
               "C", "C.initial", "sigma1", "sigma2", "m1", "loglik")

target_ns = 5000
n_chains = 1
n_burnin = 5000
n_adapt = 500
n_thin = 10
n_iter = (target_ns*n_thin)+n_burnin+n_adapt

input.list = c(list.dat, prior.list)

mod1 <- jags.model("../src/model111.mcyc.txt",
                   data = input.list,
                   inits = init.values,
                   n.chains = n_chains, 
                   n.adapt = n_adapt)

mod.fit <- coda.samples(model = mod1, 
                        variable.names = mod.params, 
                        n.iter = n_iter, thin = n_thin)

post.samp <- as.data.frame(window(mod.fit[[1]], start = n_burnin+1))

C.names <- paste("C[", 1:length(list.dat$y), "]", sep = "")
C.initial.names <- paste("C.initial[", 1:list.dat$ncyc, "]", sep = "")
loglik.names <- paste("loglik[", 1:length(list.dat$y), "]", sep = "")
param.names <- c("G", "Q", "QL", "QR", "eL", "eLF", "eRF", 
                 "sigma1", "sigma2", "m1", "G.prime", "Q.prime")

post.samples = post.samp[,param.names]
posterior.summary = ci.df(post.samples)
print(round(posterior.summary, 2))

post.C.initial = post.samp[,C.initial.names]
C.initial.summary = ci.df(post.C.initial)
post.C = post.samp[,C.names]
C.summary = ci.df(post.C)
post.loglik = t(post.samp[,loglik.names])

source("../src/printwaic.R")

######### Make Posterior Learning Plots ###############

source("../src/plot_postlearn_H111.R")

gridExtra::grid.arrange(plot.G, plot.Q, plot.QL, plot.QR,
                        plot.eL, plot.eLF, plot.eRF, plot.Q.prime,
                        ncol = 2)
# bottom = "(Blue: Prior density, Green: Posterior samples; Red tick: True, Yellow: Posterior Mean)"

###################
#### Smoothing ####
###################

# posterior.means

rep = dim(post.C)[1]
last = dim(post.C)[2]

n.s = 1E4
smooth.latent = vector(mode = "list", length = list.dat$ncyc)
smooth.data = vector(mode = "list", length = list.dat$ncyc)
for(j in 1:ncyc){
  smooth.latent[[j]] = array(dim = c(n.s, list.dat$bglength[j]))
  smooth.data[[j]] = array(dim = c(n.s, list.dat$bglength[j]))
  for(i in 1:n.s){
    r = sample(1:rep, 1)
    post.params = c(as.numeric(post.samples[r,1:10]),
                    list.dat$V, post.C[r,list.dat$cyc.end[j]])
    w = smooth.111(post.params, delta = list.dat$delta, 
                   t_n = list.dat$bglength[j])
    smooth.latent[[j]][i,] = w[,"Conc"]
    smooth.data[[j]][i,] = w[,"Data"]
  }
  smooth.latent[[j]] = ci.df(as.data.frame(smooth.latent[[j]]))
  smooth.data[[j]] = ci.df(as.data.frame(smooth.data[[j]]))
}

post.mean.C.plot = array(dim = extra.data$explength)
smooth.mean.plot = array(dim = extra.data$explength)
for(i in 1:ncyc){
  post.mean.C.plot[input1$cyc.start[i]:input1$cyc.end[i]] = C.summary[list.dat$cyc.start[i]:list.dat$cyc.end[i], 
                                                                      "Mean"]
  smooth.mean.plot[input1$cyc.end[i]+0:bglength[i]] = c(C.summary[list.dat$cyc.end[i]], smooth.latent[[i]][,"Mean"])
}

smooth.minmax = array(dim = c(ncyc, 2))
for(i in 1:ncyc){
  smooth.minmax[i,] = range(smooth.data[[i]][,c("2.5%", "97.5%")])
}

plot.min = min(exp(list.dat$y), C.summary[,"2.5%"], min(smooth.minmax[,1]))
plot.max = max(exp(list.dat$y), C.summary[,"97.5%"], max(smooth.minmax[,2]))

# Make plot
{
  plot(1:input1$explength, exp(input1$Conc),
       ylim = c(plot.min, plot.max),
       ylab = "Concentration", xlab = "Time", type = "l",
       # main = bquote("Hewett model 111: 1Box.CE.LevR.GvR ("~.(list.dat$ncyc) ~" cycles)"),
       col = "darkblue", lwd = 2,
       axes = F)
  box(bty="l")
  axis(2)
  axis(1)
  
  lines(post.mean.C.plot, col = "darkgreen", lty = 1, lwd = 1)
  lines(smooth.mean.plot, col = "darkred", lty = 1, lwd = 1)
  
  for(i in 1:ncyc){
    polygon(c(input1$cyc.start[i]:input1$cyc.end[i], 
              rev(input1$cyc.start[i]:input1$cyc.end[i])),
            c(C.summary[list.dat$cyc.start[i]:list.dat$cyc.end[i],"97.5%"], 
              rev(C.summary[list.dat$cyc.start[i]:list.dat$cyc.end[i],"2.5%"])),
            col = adjustcolor("lightgreen", alpha.f = 0.4),
            border = F)
    polygon(c(input1$cyc.end[i]+0:bglength[i], 
              rev(input1$cyc.end[i]+0:bglength[i])),
            c(c(C.summary[list.dat$cyc.end[i],"97.5%"], smooth.data[[i]][,"97.5%"]), 
              rev(c(C.summary[60,"2.5%"], smooth.data[[i]][,"2.5%"]))),
            col = adjustcolor("orange", alpha.f = 0.4),
            border = F)
    polygon(c(input1$cyc.end[i]+0:bglength[i], rev(input1$cyc.end[i]+0:bglength[i])),
            c(c(C.summary[list.dat$cyc.end[i],"97.5%"], smooth.latent[[i]][,"97.5%"]), 
              rev(c(C.summary[60,"2.5%"], smooth.latent[[i]][,"2.5%"]))),
            col = adjustcolor("darkred", alpha.f = 0.4),
            border = F)
    polygon(c(input1$cyc.start[i]-1:0, rev(input1$cyc.start[i]-1:0)),
            c(c(C.initial.summary[i,"97.5%"], C.summary[list.dat$cyc.start[i],"97.5%"]), 
              rev(c(C.initial.summary[i,"2.5%"], C.summary[list.dat$cyc.start[i],"2.5%"]))),
            col = adjustcolor("steelblue", alpha.f = 0.4),
            border = F)
  }
  
  legend("topright",
         # x = 0.85*extra.data$explength, y = 1.05*plot.max, 
         legend = c("Observed data", "Posterior mean", "Smoothing"),
         col = c("darkblue", "darkgreen", "darkred"), 
         lty = c(1, 1, 1), cex = 0.8, bty = "n", lwd = c(2, 2, 2),
         y.intersp = 0.5)
}
