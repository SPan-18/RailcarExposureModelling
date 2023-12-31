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
C0 = 10
delta = 1

obs.sd = 0.1         # log observation error s.d. 

###################

out = array(dim = explength)
G.Ind = array(dim = explength)
Cyc.Ind = array(dim = explength)
obs.vec = array(dim = explength)
C.initial[1] = C0
for(cyc in 1:ncyc){
  par.vals = c(G, Q, V, T.G[cyc] - cyc.start[cyc], C.initial[cyc])
  out[cyc.start[cyc]:bg.end[cyc]] = sapply(1:(bg.end[cyc] - cyc.start[cyc] + 1),
                                           function(x) model101(x, params = par.vals))
  C.initial[cyc+1] = out[bg.end[cyc]]
  
  G.Ind[cyc.start[cyc]:bg.end[cyc]] = rep(c(1,0), times = c(T.G[cyc] - cyc.start[cyc], 
                                                            cyc.tlength[cyc] - T.G[cyc] + cyc.start[cyc]))
  Cyc.Ind[cyc.start[cyc]:bg.end[cyc]] = cyc
  obs.vec[cyc.start[cyc]:bg.end[cyc]] = rep(c(1,0), times = c(cyclength[cyc], bglength[cyc]))
}

obs = sapply(out, function(x) log(x) + rnorm(1, 0, obs.sd))

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

# obsdf.raw = obsdf
# obsdf = obsdf[c(1:30, ((1:30)+40), ((1:30)+80)),]
# 
# 
# # priorGQ = find.priors.101(df = obsdf, V = V)
# # G.prior = priorGQ[1,]
# # Q.prior = priorGQ[2,]
# 
G.prior = c(200, 1800)
Q.prior = c(3.5, 46)
# 
list.prior.data = list(G.prior = c(80, 1500),
                       Q.prior = c(3.5, 46))
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
                   C0 = C0,
                   sigma1 = obs.sd)

bg.data = obs[is.na(obs.dat)]
bg.index = which(is.na(obs.dat))
