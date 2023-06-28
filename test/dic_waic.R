rm(list = ls())

library(INLA)
library(splines)

source("../src/functions.R")

ic_table_c1 <- function(list){
  temp <- data.frame("k" = as.numeric(), "dic" = as.numeric(), 
                     "waic" = as.numeric(), 
                     "df" = as.numeric(), "dic" = as.numeric(), 
                     "waic" = as.numeric(),
                     "RW" = as.character(), "dic" = as.numeric(), 
                     "waic" = as.numeric())
  x.time <- 1:T.end
  C.dat <- c(list$y, rep(NA, T.end - list$T.obs))
  new.data <- data.frame(time = x.time, C = C.dat)
  
  kseq = c(8, 5, 2)
  dfseq = c(3, 10, 20)
  for(k in 1:3){
    knots <- seq(0, T.end, by = kseq[k])
    cat('\r', "B-spline: knots every", kseq[k])
    m.bs3 <- inla(C ~ 1 + bs(time, knots = knots), 
                  data = new.data, 
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic=TRUE, waic = TRUE))
    flush.console()
    cat('\r', "Cubic spline: df =", dfseq[k], "\t\t")
    m.ns3 <- inla(C ~ 1 + ns(time, df = dfseq[k]),
                  data = new.data, 
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic=TRUE, waic = TRUE))
    temp[k, 1:3] = c(T.end/kseq[k], m.bs3$dic$dic, m.bs3$waic$waic)
    temp[k, 4:6] = c(dfseq[k], m.ns3$dic$dic, m.ns3$waic$waic)
  }
  flush.console()
  cat('\r', "RW2: Random walk of order 2")
  m.rw2 <- inla(C ~ -1 + f(time, model = "rw2", constr = FALSE),
                data = new.data, 
                control.predictor = list(compute = TRUE),
                control.compute = list(dic=TRUE, waic = TRUE))
  flush.console()
  cat('\r', "CRW2: Continuous RW2")
  m.crw2 <- inla(C ~ -1 + f(time, model = "crw2", constr = FALSE),
                 data = new.data, 
                 control.predictor = list(compute = TRUE),
                 control.compute = list(dic=TRUE, waic = TRUE))
  temp[1, 7:9] = c("RW2", m.rw2$dic$dic, m.rw2$waic$waic)
  temp[2, 7:9] = c("CRW2", m.crw2$dic$dic, m.crw2$waic$waic)
  flush.console()
  message('\r', "Predictive Information Criteria: Finished")
  return(temp)
}

ic_table_c3 <- function(input1.conc){
  temp <- data.frame("k" = as.numeric(), "dic" = as.numeric(), 
                     "waic" = as.numeric(), 
                     "df" = as.numeric(), "dic" = as.numeric(), 
                     "waic" = as.numeric(),
                     "RW" = as.character(), "dic" = as.numeric(), 
                     "waic" = as.numeric())
  T.end = length(input1.conc)
  x.time <- 1:T.end
  C.dat <- input1.conc
  new.data <- data.frame(time = x.time, C = C.dat)
  
  kseq = c(15, 10, 5)
  dfseq = c(7, 15, 20)
  for(k in 1:3){
    knots <- seq(0, T.end, by = kseq[k])
    cat('\r', "B-spline: knots every", kseq[k])
    m.bs3 <- inla(C ~ 1 + bs(time, knots = knots), 
                  data = new.data, 
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic=TRUE, waic = TRUE))
    flush.console()
    cat('\r', "Cubic spline: df =", dfseq[k], "\t\t")
    m.ns3 <- inla(C ~ 1 + ns(time, df = dfseq[k]),
                  data = new.data, 
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic=TRUE, waic = TRUE))
    temp[k, 1:3] = c(T.end/kseq[k], m.bs3$dic$dic, m.bs3$waic$waic)
    temp[k, 4:6] = c(dfseq[k], m.ns3$dic$dic, m.ns3$waic$waic)
  }
  flush.console()
  cat('\r', "RW2: Random walk of order 2")
  m.rw2 <- inla(C ~ -1 + f(time, model = "rw2", constr = FALSE),
                data = new.data, 
                control.predictor = list(compute = TRUE),
                control.compute = list(dic=TRUE, waic = TRUE))
  flush.console()
  cat('\r', "CRW2: Continuous RW2")
  m.crw2 <- inla(C ~ -1 + f(time, model = "crw2", constr = FALSE),
                 data = new.data, 
                 control.predictor = list(compute = TRUE),
                 control.compute = list(dic=TRUE, waic = TRUE))
  temp[1, 7:9] = c("RW2", m.rw2$dic$dic, m.rw2$waic$waic)
  temp[2, 7:9] = c("CRW2", m.crw2$dic$dic, m.crw2$waic$waic)
  flush.console()
  message('\r', "Predictive Information Criteria: Finished")
  return(temp)
}

source("../src/exactsolsim.101.R")
H101_1 <- ic_table_c1(list.dat)

source("../src/exactsolsim.111.R")
H111_1 <- ic_table_c1(list.dat)

source("../src/sim.101.mcyc.R")
H101_3 <- ic_table_c3(input1$Conc)

source("../src/sim.111.mcyc.R")
H111_3 <- ic_table_c3(input1$Conc)

cat("\t DIC and WAIC: Hewett 101 (1 cycle)\n")
print(H101_1, na.print = "")
cat("\n\n")
cat("\t DIC and WAIC: Hewett 101 (3 cycles)\n")
print(H101_3, na.print = "")
cat("\n\n")
cat("\t DIC and WAIC: Hewett 111 (1 cycle)\n")
print(H111_1, na.print = "")
cat("\n\n")
cat("\t DIC and WAIC: Hewett 111 (3 cycles)\n")
print(H111_3, na.print = "")
