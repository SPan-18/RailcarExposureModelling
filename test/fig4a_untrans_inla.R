rm(list = ls())

library(INLA)
library(splines)
library(R2jags)

source("../src/functions.r")
source("../src/exactsolsim.101.R")

##### Prior parameters #####

band = 5
G.prior = G.prior           # uniform prior
Q.prior = Q.prior           # uniform prior
m1.prior = c(0, 1E-2)       # normal prior (mean, precision)
tau1.prior = c(2, var(list.dat$y))    # inverse gamma prior for precision of log(e1)
# tau2.prior = c(2, 5)
tau2.prior = c(band*2, band*var(list.dat$y))    # inverse gamma prior for precision of e2

prior.list = list(G.prior = G.prior, 
                  Q.prior = Q.prior,
                  m1.prior = m1.prior,
                  tau1.prior = tau1.prior,
                  tau2.prior = tau2.prior
)

##### jags #####

n_iter = 10000
n_chains = 3
burnin = 5000

mod.params = c("G", "Q", "C.0", "C", "sigma1", "sigma2", "m1", "loglik")

input.list = c(list.dat, prior.list)

mod.fit = jags(data = input.list,
               parameters.to.save = mod.params, 
               n.chains = n_chains, n.iter = n_iter, n.burnin = burnin, 
               model.file = "../src/model101.txt")

# print(mod.fit)     # posterior means and credible intervals

post.G = mod.fit$BUGSoutput$sims.list$G
post.Q = mod.fit$BUGSoutput$sims.list$Q
post.C0 = mod.fit$BUGSoutput$sims.list$C.0
post.sigma1 = mod.fit$BUGSoutput$sims.list$sigma1
post.sigma2 = mod.fit$BUGSoutput$sims.list$sigma2
post.m1 = mod.fit$BUGSoutput$sims.list$m1

post.samples = data.frame(G = post.G,
                          Q = post.Q,
                          C.0 = post.C0,
                          sigma1 = post.sigma1,
                          sigma2 = post.sigma2,
                          m1 = post.m1)

posterior.summary = ci.df(post.samples)
print(round(posterior.summary, 2))

post.mean.C = mod.fit$BUGSoutput$mean$C
post.C = as.data.frame(mod.fit$BUGSoutput$sims.list$C)
C.summary = ci.df(post.C)
post.loglik = t(mod.fit$BUGSoutput$sims.list$loglik)

source("../src/printwaic.R")

###################
#### Smoothing ####
###################

rep = dim(post.C)[1]
last = dim(post.C)[2]

n.s = 1E4
smooth1 = array(dim = c(n.s, T.end - list.dat$T.obs))
smooth.dat.1 = array(dim = c(n.s, T.end - list.dat$T.obs))

for(i in 1:n.s){
  r = sample(1:rep, 1)
  post.params1 = c(post.G[r], post.Q[r], list.dat$V, post.C[r,list.dat$T.obs],
                   post.m1[r], post.sigma1[r], post.sigma2[r])
  w = smooth.101(post.params1, delta = list.dat$delta, 
                 t_n = T.end - list.dat$T.obs)
  smooth1[i,] = w[,"Conc"]
  smooth.dat.1[i,] = w[,"Data"]
}

smooth1.df = ci.df(as.data.frame(smooth1))
smooth.dat.1.df = ci.df(as.data.frame(smooth.dat.1))


##### INLA Smoothing #####

x.time <- 1:T.end
# C.dat <- c(list.dat$y, rep(NA, T.end - T.obs))
C.dat <- c(exp(list.dat$y), rep(NA, T.end - T.obs))
new.data <- data.frame(time = x.time, C = C.dat)
m.crw2 <- inla(C ~ -1 + f(time, model = "crw2", constr = FALSE),
               data = new.data, 
               control.predictor = list(compute = TRUE),
               control.compute = list(dic=TRUE, waic = TRUE))

plot.min = -50
plot.max = max(exp(list.dat$y), C.summary[,"97.5%"], 
               smooth1.df[,"97.5%"], smooth.dat.1.df[,"97.5%"])

jpeg("../fig/fig4a.jpg", width = 6.7, height = 4.47, units = "in", res = 300)
{
par(mar = c(4, 4, 0.1, 0.1))
plot(1:40, c(exp(list.dat$y), rep(NA, T.end-list.dat$T.obs)), 
     ylim = c(plot.min, plot.max+10),   # -10
     ylab = "Concentration", xlab = "Time",
     # main = "Simulated Hewett model 101: 1Box.CE.Gv",
     pch = 16, col = "darkblue", axes = F)
abline(h = 0, lty = 2, col = "darkgrey")
box(bty="l")
axis(2)
axis(1)
lines(m.crw2$summary.fitted.values$'mean', col = "mediumpurple4")
polygon(c(1:T.end, rev(1:T.end)),
        c(m.crw2$summary.fitted.values$'0.975quant', 
          rev(m.crw2$summary.fitted.values$'0.025quant')),
        col = adjustcolor("mediumpurple1", alpha.f = 0.3),
        border = F)


polygon(c(0:1, rev(0:1)),
        c(c(quantile(post.C0, probs = 0.975), C.summary[1,"97.5%"]),
          rev(c(quantile(post.C0, probs = 0.025), C.summary[1,"2.5%"]))),
        col = adjustcolor("steelblue", alpha.f = 0.4),
        border = F)

polygon(c(1:list.dat$T.obs, rev(1:list.dat$T.obs)),
        c(C.summary[,"97.5%"], rev(C.summary[,"2.5%"])),
        col = adjustcolor("lightgreen", alpha.f = 0.4),
        border = F)
lines(c(post.mean.C, rep(NA, T.end - list.dat$T.obs)), col = "darkgreen")
lines(c(rep(NA, list.dat$T.obs-1), post.mean.C[list.dat$T.obs], 
        smooth1.df[,"Mean"]), col = "darkred")
# lines(posterior.curve.fit, lty = 4)

polygon(c(list.dat$T.obs:T.end, rev(list.dat$T.obs:T.end)),
        c(c(C.summary[list.dat$T.obs,"97.5%"], smooth.dat.1.df[,"97.5%"]), 
          rev(c(C.summary[list.dat$T.obs,"2.5%"], smooth.dat.1.df[,"2.5%"]))),
        col = adjustcolor("orange", alpha.f = 0.4),
        border = F)

polygon(c(list.dat$T.obs:T.end, rev(list.dat$T.obs:T.end)),
        c(c(C.summary[list.dat$T.obs,"97.5%"], smooth1.df[,"97.5%"]), 
          rev(c(C.summary[list.dat$T.obs,"2.5%"], smooth1.df[,"2.5%"]))),
        col = adjustcolor("darkred", alpha.f = 0.4),
        border = F)
points(1:40, c(rep(NA, T.obs), exp(obs1[(T.obs+1):T.end])), pch = 4, col = "yellow")
legend(27.5, 60, legend = c("Semi-parametric", "Proposed model"),
       col = c("purple", "darkred"),
       title = "Smoothing", title.font = 2,
       pch = c(15, 15), cex = 1, bty = "n",
       inset = c(0.06, -0.06), y.intersp = 0.75)
}
fig_label("(a)", font = 2)
dev.off()
