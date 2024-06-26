# useful functions

# Get mean and variance of log X where
# X is univariate lognormal with some desired mean and variance

get.logn.par = function(des.mean, des.sd){
  den = des.mean^2 + des.sd^2
  mu = log(des.mean^2/sqrt(den))
  sd = sqrt(log(den/des.mean^2))
  return(c(mu, sd))
}


# indicator for generator G

G.indicator = function(n, stops){
  t0 = c(0, stops, n)
  vec = array(dim = n)
  flag = TRUE
  for(j in 1:(length(t0)-1)){
    vec[(t0[j]+1):t0[j+1]] = flag
    flag = !flag
  }
  return(as.numeric(vec))
}

# G.indicator(10, c(3,6,8))

# Mean and CI od posterior samples (inputs as columns of a data frame)

ci.df = function(df){
  rnames = colnames(df)
  out = array(dim = c(length(rnames), 5))
  for(i in 1:length(rnames)){
    out[i,1] = mean(df[,i])
    out[i,2] = sd(df[,i])
    out[i,3] = quantile(df[,i], 0.025)
    out[i,4] = quantile(df[,i], 0.5)
    out[i,5] = quantile(df[,i], 0.975)
  }
  colnames(out) = c("Mean", "SD", "2.5%", "50%", "97.5%")
  rownames(out) = rnames
  # cat("Summary of Posterior Samples \n")
  # print(out)
  return(out)
}

##########################################

model101 = function(x, params){
  G = params[1]
  Q = params[2]
  V = params[3]
  T0 = params[4]
  C0 = params[5]
  temp1 = Q/V
  temp2 = G/Q
  C.max = C0*exp(-T0*temp1) + temp2*(1 - exp(-T0*temp1))
  if(x <= T0) return(C0*exp(-x*temp1) + temp2*(1 - exp(-x*temp1)))
  else return(C.max*exp(-(x-T0)*temp1))
}

##########################################

model111 = function(x, params){
  G = params[1]
  Q = params[2]
  V = params[3]
  T0 = params[4]
  C0 = params[5]
  Q.L = params[6]
  Q.R = params[7]
  e.L = params[8]
  e.LF = params[9]
  e.RF = params[10]
  temp1 = (Q + e.LF*Q.L + e.RF*Q.R)/V
  temp2 = (1 - e.L*e.LF)*G/(Q + e.LF*Q.L + e.RF*Q.R)
  C.max = C0*exp(-T0*temp1) + temp2*(1 - exp(-T0*temp1))
  if(x <= T0) return(C0*exp(-x*temp1) + temp2*(1 - exp(-x*temp1)))
  else return(C.max*exp(-(x-T0)*temp1))
}

##########################################

transition.101 = function(params, delta, t_n){
  G = params[1]
  Q = params[2]
  V = params[3]
  T0 = params[4]
  C0 = params[5]
  temp1 = 1 - Q*delta/V
  temp2 = G*delta/V
  G.ind = rep(c(1,0), times = c(T0, t_n - T0))
  Conc = array(dim = t_n)
  Conc[1] = temp1*C0 + temp2*G.ind[1]
  for(i in 1:(t_n-1)){
    Conc[i+1] = temp1*Conc[i] + temp2*G.ind[i+1]
  }
  return(Conc)
}

##########################################

transition.111 = function(params, delta, t_n){
  G = params[1]
  Q = params[2]
  V = params[3]
  T0 = params[4]
  C0 = params[5]
  Q.L = params[6]
  Q.R = params[7]
  e.L = params[8]
  e.LF = params[9]
  e.RF = params[10]
  temp1 = 1 - (Q + e.LF*Q.L + e.RF*Q.R)*delta/V
  temp2 = (1 - e.L*e.LF)*G*delta/V
  G.ind = rep(c(1,0), times = c(T0, t_n - T0))
  Conc = array(dim = t_n)
  Conc[1] = temp1*C0 + temp2*G.ind[1]
  for(i in 1:(t_n-1)){
    Conc[i+1] = temp1*Conc[i] + temp2*G.ind[i+1]
  }
  return(Conc)
}

##################################

find.priors.101 = function(df, V){
  decay = subset(df, G.Ind == 0, select = c("Conc", "Time"))
  rise = subset(df, G.Ind == 1, select = c("Conc", "Time"))
  fit.decay = lm(log(Conc) ~ Time, decay)
  Q.hat1 = as.numeric((-fit.decay$coefficients[2])*V)
  fit.rise = lm(log(Conc) ~ Time, rise)
  Q.hat2 = as.numeric(fit.rise$coefficients[2]*V)
  G.hat1 = as.numeric(exp(fit.rise$coefficients[1])*Q.hat1)
  G.hat2 = G.hat1*Q.hat2/Q.hat1
  G.prior = c(0.6*min(G.hat1, G.hat2), 3.5*max(G.hat1, G.hat2))
  Q.prior = c(0.6*min(Q.hat1, Q.hat2), 2.5*max(Q.hat1, Q.hat2))
  return(rbind(G.prior, Q.prior))
}

###################################

smooth.101 = function(params, delta, t_n){
  G = params[1]
  Q = params[2]
  V = params[3]
  C0 = params[4]
  m1 = params[5]
  s1 = params[6]
  s2 = params[7]
  temp = 1 - Q*delta/V
  out = array(dim = c(t_n, 2))
  colnames(out) = c("Conc", "Data")
  out[1,1] = temp*C0 + exp(rnorm(1, mean = m1, sd = s1))
  for(i in 1:(t_n-1)){
    out[i+1,1] = temp*out[i,1] + exp(rnorm(1, mean = m1, sd = s1))
  }
  out[,2] = exp(log(out[,1]) + rnorm(n = t_n, mean = 0, sd = s2))
  return(out)
}

###################################

smooth.111 = function(params, delta, t_n){
  G = params[1]
  Q = params[2]
  Q.L = params[3]
  Q.R = params[4]
  e.L = params[5]
  e.LF = params[6]
  e.RF = params[7]
  s1 = params[8]
  s2 = params[9]
  m1 = params[10]
  # alpha = params[11]
  # beta = params[12]
  V = params[11]
  C0 = params[12]
  Q.prime = Q + e.LF*Q.L + e.RF*Q.R
  G.prime = (1 - e.L*e.LF)*G
  temp = 1 - Q.prime*delta/V
  out = array(dim = c(t_n, 2))
  err = array(dim = t_n)
  colnames(out) = c("Conc", "Data")
  out[1,1] = temp*C0 + exp(rnorm(1, mean = m1, sd = s1))
  for(i in 1:(t_n-1)){
    out[i+1,1] = temp*out[i,1] + exp(rnorm(1, mean = m1, sd = s1))
  }
  out[,2] = exp(log(out[,1]) + rnorm(n = t_n, mean = 0, sd = s2))
  return(out)
}

smooth.111.dynvar = function(params, delta, t_n){
  G = params[1]
  Q = params[2]
  Q.L = params[3]
  Q.R = params[4]
  e.L = params[5]
  e.LF = params[6]
  e.RF = params[7]
  s1 = params[8]
  s2 = params[9]
  # m1 = min(params[10], 1.5)
  m1 = params[10]
  alpha = params[11]
  beta = params[12]
  V = params[13]
  C0 = params[14]
  Q.prime = Q + e.LF*Q.L + e.RF*Q.R
  G.prime = (1 - e.L*e.LF)*G
  temp = 1 - Q.prime*delta/V
  out = array(dim = c(t_n, 2))
  err = array(dim = t_n)
  colnames(out) = c("Conc", "Data")
  out[1,1] = temp*C0 + beta*exp(rnorm(1, mean = -m1, sd = s1))
  for(i in 1:(t_n-1)){
    out[i+1,1] = temp*out[i,1] + (beta^(i+1))*exp(rnorm(1, mean = -m1, sd = s1))
  }
  out[,2] = exp(log(out[,1]) + rnorm(n = t_n, mean = 0, sd = s2))
  return(out)
}

# Source: https://www.r-bloggers.com/2017/03/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}
