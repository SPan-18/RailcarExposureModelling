rm(list = ls())
set.seed(1728)

source("../src/functions.r")

###### Simulate Data (One-Zone model: Two Rise-Decay cycles) ######

T0 = c(0, 15, 40, 55, 80, 95, 120)    # generator switch is toggled

# T.end = 35           # length of the experiment

###################

G = 1000           # generation rate (mg/min or count/min)
V = 100            # volume of the room (m^3)
Q = 20             # average ventilation rate (m^3/min)
Q.L = 5            # local exhaust ventilation rate (m^3/min)
Q.R = 5            # room recirculation system ventilation rate (m^3/min)
e.L = 0.6          # fraction of the source emissions immediately captured by the local exhaust
e.LF = 0.3         # local exhaust return filtration efficiency
e.RF = 0.9         # general ventilation recirculation filtration efficiency
C0 = 10
delta = 1

###################

par.vals = c(G, Q, V, T0[2], C0, Q.L, Q.R, e.L, e.LF, e.RF)
cyc1 = sapply(1:T0[3], function(x) model111(x, params = par.vals))
par.vals[4] = T0[4] - T0[3]
par.vals[5] = cyc1[T0[3]]
cyc2 = sapply(1:(T0[5] - T0[3]), function(x) model111(x, params = par.vals))
par.vals[4] = T0[6] - T0[5]
par.vals[5] = cyc2[T0[5] - T0[3]]
cyc3 = sapply(1:(T0[7] - T0[5]), function(x) model111(x, params = par.vals))
obs1 = sapply(cyc1, function(x) x + exp(rnorm(1, -2, 1.5)))
obs2 = sapply(cyc2, function(x) x + exp(rnorm(1, -2, 1.5)))
obs3 = sapply(cyc3, function(x) x + exp(rnorm(1, -2, 1.5)))
# out2 = transition.111(toy.vals, delta, T.end)
out = c(C0, cyc1, cyc2, cyc3)
obs = c(obs1, obs2, obs3)

plot.min = min(out)
plot.max = max(out, obs)+1

jpeg("../fig/fig2.jpg", width = 7, height = 4.5, units = "in", res = 300)
par(mar = c(4, 4, 0.1, 0.1))
plot(out, type = "l", col = "darkred", lwd = 2, 
     ylim = c(plot.min, plot.max),
     # main = "Concentration during exposure rise and decay",
     xlab = "Time",
     ylab = expression(paste(C(t))),
     lty = 1,
     panel.first = rect(c(1, 41, 81), -1e6, c(30, 70, 110), 1e6, 
                        col = adjustcolor("lightgrey", alpha.f = 0.3), 
                        border=NA))
# lines(obs, lty = 2, col = "blue")
# polygon(c(0, 0, 30, 30)+1, c(0, plot.max, plot.max, 0), 
#         col = adjustcolor("lightgrey", alpha.f = 0.3),
#         border = F)
# polygon(c(0, 0, 30, 30)+41, c(0, plot.max, plot.max, 0), 
#         col = adjustcolor("lightgrey", alpha.f = 0.3),
#         border = F)
# polygon(c(0, 0, 30, 30)+81, c(0, plot.max, plot.max, 0), 
#         col = adjustcolor("lightgrey", alpha.f = 0.3),
#         border = F)
dev.off()
