set.seed(1728)

# source("functions.R")

###### Simulate Data (One-Zone model: Two Rise-Decay cycles) ######

ncyc = 3                                          # number of cycles
cyclength = rep(30, ncyc)                         # length of each cycle
bglength = rep(10, ncyc)                          # background data after each cycle
cyc.tlength = cyclength + bglength                # total cycle length
cyc.start = c(0, cumsum(c(cyclength[-ncyc] 
                          + bglength[-ncyc])))+1  # indices of cycle start
cyc.end = cyc.start + cyclength-1                 # indices of cycle end
bg.end = cyc.end + bglength                       # indices of background end
explength = bg.end[ncyc]                          # length of the experiment
T.G = cyc.start + 15                              # indices at which G turned off
C.initial = array(dim = ncyc+1)                   # concentrations at cycle ends

###################

G = 1000             # generation rate (mg/min or count/min)
V = 100              # volume of the room (m^3)
Q = 20               # average ventilation rate (m^3/min)
Q.L = 5            # local exhaust ventilation rate (m^3/min)
Q.R = 5            # room recirculation system ventilation rate (m^3/min)
e.L = 0.5          # fraction of the source emissions immediately
                   # captured by the local exhaust
e.LF = 0.5         # local exhaust return filtration efficiency
e.RF = 0.9         # general ventilation recirculation filtration efficiency
C0 = 10
delta = 1

mu2 = 0              # observation lognormal mean
sd2 = 0.1              # observation 

###################

out = array(dim = explength)
G.Ind = array(dim = explength)
Cyc.Ind = array(dim = explength)
obs.vec = array(dim = explength)
C.initial[1] = C0
for(cyc in 1:ncyc){
  par.vals = c(G, Q, V, T.G[cyc] - cyc.start[cyc], C.initial[cyc],
               Q.L, Q.R, e.L, e.LF, e.RF)
  out[cyc.start[cyc]:bg.end[cyc]] = sapply(1:(bg.end[cyc] - cyc.start[cyc] + 1),
                                           function(x) model111(x, params = par.vals))
  C.initial[cyc+1] = out[bg.end[cyc]]
  
  G.Ind[cyc.start[cyc]:bg.end[cyc]] = rep(c(1,0), times = c(T.G[cyc] - cyc.start[cyc], 
                                                            cyc.tlength[cyc] - T.G[cyc] + cyc.start[cyc]))
  Cyc.Ind[cyc.start[cyc]:bg.end[cyc]] = cyc
  obs.vec[cyc.start[cyc]:bg.end[cyc]] = rep(c(1,0), times = c(cyclength[cyc], bglength[cyc]))
}

# log-scale obs
obs = sapply(out, function(x) log(x) + rnorm(1, mu2, sd2))  

# plot.min = min(out, exp(obs))
# plot.max = max(out, exp(obs))
# plot(out, type = "p", ylim = c(plot.min, plot.max), col = obs.vec)
# lines(exp(obs), lty = 2, col = "red")

obs.dat = array(dim = explength)
for(cyc in 1:ncyc){
  obs.dat[cyc.start[cyc]:cyc.end[cyc]] = obs[cyc.start[cyc]:cyc.end[cyc]]
}

# Processed Data Input

input1 = list("Conc" = obs.dat,
              "Time" = 1:explength,
              "G.Ind" = G.Ind,
              "Cyc.Ind" = Cyc.Ind,
              "explength" = explength,
              "ncyc" = ncyc,
              "V" = V,
              "C0" = C0,
              "delta" = delta,
              "cyc.start" = cyc.start,
              "cyc.end" = cyc.end,
              "T.G" = T.G,
              "bglength" = bglength)

# Post-Processing of input and prior formation

pp.cyclength = input1$cyc.end - input1$cyc.start + 1
pp.cyc.start = cumsum(c(0, pp.cyclength[-ncyc])) + input1$cyc.start[1]
pp.cyc.end = pp.cyc.start + pp.cyclength - 1
pp.bglength = input1$bglength
pp.T.G = pp.cyc.start + input1$T.G - input1$cyc.start - 1
pp.Conc = input1$Conc[!is.na(input1$Conc)]
pp.Cyc.Ind = input1$Cyc.Ind[!is.na(input1$Conc)]
pp.G.Ind = input1$G.Ind[!is.na(input1$Conc)]

# priors
G.prior = c(200, 1800)
Q.prior = c(3, 50)
# 
list.prior.data = list(G.prior = G.prior,
                       Q.prior = Q.prior)
# 
list.dat = list("y" = pp.Conc,
                "G.Ind" = pp.G.Ind,
                "ncyc" = input1$ncyc,
                "V" = input1$V,
                "C0" = input1$C0,
                "delta" = input1$delta,
                "cyc.start" = pp.cyc.start,
                "cyc.end" = pp.cyc.end,
                "bglength" = pp.bglength)

extra.data = list("T.G" = pp.T.G,
                  "explength" = input1$explength)
# 
true.params = list(G = G,
                   Q = Q,
                   Q.L = Q.L,
                   Q.R = Q.R,
                   e.L = e.L,
                   e.LF = e.LF,
                   e.RF = e.RF,
                   C0 = C0,
                   sigma1 = sd2)

bg.data = obs[is.na(obs.dat)]
bg.index = which(is.na(obs.dat))
