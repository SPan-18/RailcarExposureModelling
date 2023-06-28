rm(list = ls())

library(R2jags)
library(LaplacesDemon)

source("../src/functions.r")
source("../src/exactsolsim.111.R")

##### Prior parameters #####

band = 5
G.prior = G.prior                   # uniform prior
Q.prior = Q.prior                   # uniform prior
m1.prior = c(0, 1E-2)               # normal prior (mean, precision)
tau1.prior = c(2, var(list.dat$y))  # IG prior for precision of log(e1)
tau2.prior = c(band*2, band*var(list.dat$y))  # IG prior for precision of e2
QL.prior = c(2,20)
QR.prior = c(2,20)
eL.prior = c(0, 1)
eLF.prior = c(0, 1)
eRF.prior = c(0, 1)

prior.list = list(G.prior = G.prior, 
                  Q.prior = Q.prior,
                  m1.prior = m1.prior,
                  tau1.prior = tau1.prior,
                  tau2.prior = tau2.prior,
                  QL.prior = QL.prior,
                  QR.prior = QR.prior,
                  eL.prior = eL.prior,
                  eLF.prior = eLF.prior,
                  eRF.prior = eRF.prior
)

##### jags #####

n_iter = 10000
n_chains = 3
burnin = 1000

mod.params = c("G", "Q", "QL", "QR", "eL", "eLF", "eRF", 
               "C.0", "C", "sigma1", "sigma2", "m1", "loglik")

input.list = c(list.dat, prior.list)

mod.fit = jags(data = input.list,
               parameters.to.save = mod.params, 
               n.chains = n_chains, n.iter = n_iter, n.burnin = burnin, 
               model.file = "../src/model111.txt")

post.G = mod.fit$BUGSoutput$sims.list$G
post.Q = mod.fit$BUGSoutput$sims.list$Q
post.QL = mod.fit$BUGSoutput$sims.list$QL
post.QR = mod.fit$BUGSoutput$sims.list$QR
post.eL = mod.fit$BUGSoutput$sims.list$eL
post.eLF = mod.fit$BUGSoutput$sims.list$eLF
post.eRF = mod.fit$BUGSoutput$sims.list$eRF
post.C0 = mod.fit$BUGSoutput$sims.list$C.0
post.sigma1 = mod.fit$BUGSoutput$sims.list$sigma1
post.sigma2 = mod.fit$BUGSoutput$sims.list$sigma2
post.m1 = mod.fit$BUGSoutput$sims.list$m1

post.samples = data.frame(G = post.G,
                          Q = post.Q,
                          QL = post.QL,
                          QR = post.QR,
                          eL = post.eL,
                          eLF = post.eLF,
                          eRF = post.eRF,
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
