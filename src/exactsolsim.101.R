set.seed(1728)

# source("functions.r")

###### Simulate Data (One-Zone model: Two Rise-Decay cycles) ######

T0 = 15              # generator switch is toggled
T.obs = 20           # observed until t = T.obs
T.end = 40           # length of the experiment

###################

G = 1000             # generation rate (mg/min or count/min)
V = 100              # volume of the room (m^3)
Q = 20               # average ventilation rate (m^3/min)
# Q.L = 5            # local exhaust ventilation rate (m^3/min)
# Q.R = 5            # room recirculation system ventilation rate (m^3/min)
# e.L = 0.5          # fraction of the source emissions immediately 
#                    # captured by the local exhaust
# e.LF = 0.5         # local exhaust return filtration efficiency
# e.RF = 0.9         # general ventilation recirculation filtration efficiency
C0 = 10
delta = 1

sd.obs = 0.1

###################

toy.vals = c(G, Q, V, T0, C0)
out1 = sapply(1:T.end, function(x) model101(x, params = toy.vals))
obs1 = sapply(out1, function(x) log(x) + rnorm(1, 0, sd.obs))

# plot.min = min(out1[1:T.obs], exp(obs1[1:T.obs]))
# plot.max = max(out1[1:T.obs], exp(obs1[1:T.obs]))
# plot(out1[1:T.obs], type = "l", ylim = c(plot.min, plot.max))
# lines(exp(obs1[1:T.obs]), lty = 2, col = "red")

obsdf = data.frame("Conc" = obs1[1:T.obs],
                   "G.Ind" = rep(c(1,0), times = c(T0, T.obs - T0)),
                   "Time" = 1:T.obs)

# priorGQ = find.priors.101(df = obsdf, V = V)
# G.prior = priorGQ[1,]
# Q.prior = priorGQ[2,]
G.prior = c(100, 3000)
Q.prior = c(2, 40)

list.dat = list(y = obsdf$Conc,
                T.obs = T.obs, 
                V = V,
                C0 = C0,
                G.ind = obsdf$G.Ind, delta = 1)

true.params = list(G = G,
                   Q = Q,
                   C0 = C0,
                   sigma1 = sd.obs)
