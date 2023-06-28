# Model: Hewett 111 (Univariate)
# Use: Make posterior learning plot for simulated data

library(ggplot2)

plot.G = ggplot(post.samples, aes(x = G)) +
  stat_function(fun = dunif, geom = 'area', args = G.prior, n = 100,
                fill = "#619CFF", alpha = 0.5, color = "black") +
  geom_density(aes(y = after_stat(density)), fill = "green", alpha = 0.5,
               color = "black") +
  geom_point(aes(x = true.params$G, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(G), y = 0),
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(G), y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

plot.Q = ggplot(post.samples, aes(x = Q)) +
  stat_function(fun = dunif, geom = 'area', args = Q.prior, n = 100,
                fill = "#619CFF", alpha = 0.5, color = "black") +
  geom_density(aes(y = after_stat(density)), fill = "green", alpha = 0.5,
               color = "black") +
  geom_point(aes(x = true.params$Q, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(Q), y = 0),
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(Q), y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

plot.QL = ggplot(post.samples, aes(x = QL)) +
  stat_function(fun = dunif, geom = 'area', args = QL.prior, n = 100, 
                fill = "#619CFF", alpha = 0.5, color = "black") +
  geom_density(aes(y = ..density..), fill = "green", alpha = 0.5, 
               color = "black") +
  geom_point(aes(x = true.params$Q.L, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(QL), y = 0), 
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(Q[L]), y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

plot.QR = ggplot(post.samples, aes(x = QR)) +
  stat_function(fun = dunif, geom = 'area', args = QR.prior, n = 100, 
                fill = "#619CFF", alpha = 0.5, color = "black") +
  geom_density(aes(y = ..density..), fill = "green", alpha = 0.5, 
               color = "black") +
  geom_point(aes(x = true.params$Q.R, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(QR), y = 0), 
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(Q[R]), y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())


plot.eL = ggplot(post.samples, aes(x = eL)) +
  stat_function(fun = dunif, geom = 'area', args = eL.prior, n = 100, 
                fill = "#619CFF", alpha = 0.5, color = "black") +
  geom_density(aes(y = ..density..), fill = "green", alpha = 0.5, 
               color = "black") +
  geom_point(aes(x = true.params$e.L, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(eL), y = 0), 
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(e[L]), y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())


plot.eLF = ggplot(post.samples, aes(x = eLF)) +
  stat_function(fun = dunif, geom = 'area', args = eLF.prior, n = 100, 
                fill = "#619CFF", alpha = 0.5, color = "black") +
  geom_density(aes(y = ..density..), fill = "green", alpha = 0.5, 
               color = "black") +
  geom_point(aes(x = true.params$e.LF, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(eLF), y = 0), 
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(e[LF]), y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

plot.eRF = ggplot(post.samples, aes(x = eRF)) +
  stat_function(fun = dunif, geom = 'area', args = eRF.prior, n = 100, 
                fill = "#619CFF", alpha = 0.5, color = "black") +
  geom_density(aes(y = ..density..), fill = "green", alpha = 0.5, 
               color = "black") +
  geom_point(aes(x = true.params$e.RF, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(eRF), y = 0), 
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(e[RF]), y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

post.size = length(post.samples$Q)
Q.prime.prior.samples <- array(dim = post.size)
count = 0
count.in = 0
while(count.in < post.size){
  Q.prior.samples <- runif(n = 1, min = Q.prior[1], max = Q.prior[2])
  QL.prior.samples <- runif(n = 1, min = QL.prior[1], max = QL.prior[2])
  QR.prior.samples <- runif(n = 1, min = QR.prior[1], max = QR.prior[2])
  eLF.prior.samples <- runif(n = 1, min = eLF.prior[1], max = eLF.prior[2])
  eRF.prior.samples <- runif(n = 1, min = eRF.prior[1], max = eRF.prior[2])
  Q.prime.prior <- Q.prior.samples + 
    eLF.prior.samples*QL.prior.samples + 
    eRF.prior.samples*QR.prior.samples
  if(Q.prime.prior < (max(post.samples$Q.prime)+5) &&
     Q.prime.prior > min(post.samples$Q.prime)){
    Q.prime.prior.samples[count.in+1] = Q.prime.prior
    count.in = count.in + 1
    count = count + 1
  }else count = count + 1
}
prop.in = count.in/count
Q.prime.samples = data.frame(prior = Q.prime.prior.samples,
                             posterior = post.samples$Q.prime)
true.Qprime <- true.params$Q + true.params$e.LF*true.params$Q.L + 
  true.params$e.RF*true.params$Q.R 

plot.Q.prime <- ggplot(Q.prime.samples) +
  geom_density(aes(x = posterior, y = ..density..), 
               fill = "green", alpha = 0.5, 
               color = "black") +
  geom_density(aes(x = prior, y = ..density..*prop.in), 
               fill = "#619CFF", alpha = 0.5, 
               color = "black") +
  geom_point(aes(x = true.Qprime, y = 0), pch = 24, fill = "red", size = 5,
             inherit.aes = F) +
  geom_point(aes(x = mean(posterior), y = 0),
             pch = 24, fill = "yellow", size = 5,
             inherit.aes = F) +
  labs(x = expression(paste("Q +", "e"[RF], "Q"[R], "+", "e"[LF], "Q"[L])), 
       y = "Density") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())
