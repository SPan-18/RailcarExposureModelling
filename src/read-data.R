conc_dat <- read.table("../data/conc.txt")
G_Ind_dat <- read.table("../data/G_Ind.txt")
Cyc_Ind_dat <- read.table("../data/Cyc_Ind.txt")
C_0_dat <- read.table("../data/C_0.txt")
Cyc_start_dat <- read.table("../data/cyc_start.txt")
Cyc_end_dat <- read.table("../data/cyc_end.txt")
T.G_dat <- read.table("../data/T.G.txt")
bglength_dat <- read.table("../data/bglength.txt")

conc = as.vector(conc_dat$V1)
T.end = length(conc)
G_Ind = as.vector(G_Ind_dat$V1)
Cyc_Ind = as.vector(Cyc_Ind_dat$V1)
C_0 = as.vector(C_0_dat$V1)
ncyc = 3
Cyc_start = as.vector(Cyc_start_dat$V1)
Cyc_end = as.vector(Cyc_end_dat$V1)
T.G = as.vector(T.G_dat$V1)
bglength = as.vector(bglength_dat$V1)

# Input

input1 = list("Conc" = log(conc),
              "Time" = 1:T.end,
              "G.Ind" = G_Ind,
              "Cyc.Ind" = Cyc_Ind,
              "explength" = T.end,
              "ncyc" = ncyc,
              "V" = 150.5,
              "C0" = C_0,
              "delta" = 1,
              "cyc.start" = Cyc_start,
              "cyc.end" = Cyc_end,
              "T.G" = T.G,
              "bglength" = bglength)

# Post-Processing of input and prior formation

pp.cyclength = input1$cyc.end - input1$cyc.start + 1
pp.cyc.start = cumsum(c(0, pp.cyclength[-input1$ncyc])) + input1$cyc.start[1]
pp.cyc.end = pp.cyc.start + pp.cyclength - 1
pp.bglength = input1$bglength
pp.T.G = pp.cyc.start + input1$T.G - input1$cyc.start - 1
pp.Conc = input1$Conc[!is.na(input1$Conc)]
pp.Cyc.Ind = input1$Cyc.Ind[!is.na(input1$Conc)]
pp.G.Ind = input1$G.Ind[!is.na(input1$Conc)]

G.prior = c(100, 10000)
Q.prior = c(10, 150)

list.prior.data = list(G.prior = G.prior,
                       Q.prior = Q.prior)

######## Prepare data and its details #######

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

# bg.data = obs.dat.1[is.na(obs.dat.1)]
# bg.index = which(is.na(obs.dat.1))
# bglength = dat.1.bglength
# ncyc = 3