
########### list of priors ############
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

mod.params = c("G", "Q", "QL", "QR", "eL", "eLF", "eRF", "G.prime", "Q.prime")

# target_ns = 5000
n_chains = 1
n_burnin = 1000
n_adapt = 500
n_thin = 1
# n_iter = (target_ns*n_thin)+n_burnin+n_adapt
n_iter = 1000

input.list = c(list.dat, prior.list)

mod1 <- jags.model("../src/model111.mcyc.txt",
                   data = input.list,
                   inits = init.values,
                   n.chains = n_chains, 
                   n.adapt = n_adapt)
update(mod1, n.iter = n_burnin)

mod.fit <- coda.samples(model = mod1, 
                        variable.names = mod.params, 
                        n.iter = n_iter, thin = n_thin)

post.samp <- do.call(rbind, mod.fit)
