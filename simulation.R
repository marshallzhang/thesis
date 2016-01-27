source("/Users/marshall/Documents/senior/thesis/functions.R")

#
# Mimicking unconstrained.
#
N = 10000

# OU process.
theta = c(0, 1, 1)

o.sims = array(0, dim = c(N, 100))
o.sims.start = rnorm(N, mean = 0, sd = sqrt(1/2))
for (i in 1:N) {
  print(i)
  o.sims[i, ] = ou(o.sims.start[i], 100, theta)
}

cov(o.sims[,c(1,100)])

nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, exp(-1)/2, exp(-1)/2, 1/2), nrow = 2, ncol = 2))} 

nu.o.sims = gen.double.nu.bridge(ou, theta, 0, 0, 1000)
nu.o.sims[2, ] = rev(nu.o.sims[2,])
plot(nu.o.sims[1,], type = "l")
lines(nu.o.sims[2,])
sum(nu.o.sims[1,] - nu.o.sims[2,] > 0 | nu.o.sims[1,] - nu.o.sims[2,] < 0)

nu.o.sims = gen.nu.bridge(nu.o, ou, theta, 10000, 100)

# CIR process.
theta = c(1, 1, 1)
c.sims = array(0, dim = c(N, 100))
c.sims.start = rgamma(N, shape = 2, rate = 2)
for (i in 1:N) {
  print(i)
  c.sims[i, ] = cir(c.sims.start[i], 100, theta)
}
cov(c.sims[,c(1,100)])

nu.c = function(n) { mcmc(c(1, 1), 20 * N + 2000, 2000, everyother = 20, trace = 1000, d.posterior = joint, r.proposal = function(n, current) {rgamma(2, rate = 2, shape = 2)}) }

nu.c.sims = gen.nu.bridge(nu.c, cir, theta, 10000, 100)

# Plot.
th.q = qsOU(ppoints(o.sims[,50]), c(0, 1, 1))
data.q = sort(o.sims[,50])
my.data.q = sort(nu.o.sims[,50])
data = data.frame(dat = data.q, th = th.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
o.plot = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, color = variable)) + 
              xlab("Theoretical Quantiles") + 
              ylab("Sample Quantiles") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = -3:3) + 
              scale_x_continuous(breaks = -3:3) + 
              scale_color_manual(name = "", labels = c("Unconstrained diffusion       ", bquote(paste(nu, "-bridge"))), values = c(dat = "black", my.dat = "gray")) 

th.q = qsCIR(ppoints(c.sims[,50]), c(1,1,1))
data.q = sort(c.sims[,50])
my.data.q = sort(nu.c.sims[,50])
data = data.frame(dat = data.q, th = th.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
c.plot = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, color = variable)) + 
              xlab("Theoretical Quantiles") + 
              ylab("Sample Quantiles") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = 0:6) + 
              scale_x_continuous(breaks = 0:6) + 
              scale_shape_discrete(name = "Simulations", labels = c("CIR", bquote(paste(nu^C, "-bridge of diffusion")))) +
              scale_color_manual(name = "", values = c(dat = "black", my.dat = "gray")) 

grid_arrange_shared_legend(o.plot, c.plot)

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/invariant_density.pdf", width= 10, height = 5, #' see how it looks at this size
    useDingbats=F)
grid_arrange_shared_legend(o.plot, c.plot)
dev.off()

#
# Varying correlation.
#

N = 10000

# OU process.
o.theta = c(0, 1, 1)

nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, exp(-1)/2, exp(-1)/2, 1/2), nrow = 2, ncol = 2))} 
nu.o.sims.orig = gen.nu.bridge(nu.o, ou, o.theta, 10000, 100)
nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, 0, 0, 1/2), nrow = 2, ncol = 2))} 
nu.o.sims.zero= gen.nu.bridge(nu.o, ou, o.theta, 10000, 100)
nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, -0.4, -0.4, 1/2), nrow = 2, ncol = 2))} 
nu.o.sims.neg= gen.nu.bridge(nu.o, ou, o.theta, 10000, 100)
nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, 0.49, 0.49, 1/2), nrow = 2, ncol = 2))} 
nu.o.sims.pos= gen.nu.bridge(nu.o, ou, o.theta, 10000, 100)

# CIR process.
c.theta = c(1, 1, 1)

nu.c = function(n) { mcmc(c(1, 1), 20 * N + 2000, 2000, everyother = 20, trace = 1000, d.posterior = joint, r.proposal = function(n, current) {rgamma(2, rate = 2, shape = 2)}) }
nu.c.sims.orig = gen.nu.bridge(nu.c, cir, c.theta, 10000, 100)
nu.c = function(n) { gaussian.cop(10000, function(x) {qgamma(x, shape = 2, rate = 2)}, function(x) {qgamma(x, shape = 2, rate = 2)},  matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2))} 
nu.c.sims.zero = gen.nu.bridge(nu.c, cir, c.theta, 10000, 100)
nu.c = function(n) { gaussian.cop(10000, function(x) {qgamma(x, shape = 2, rate = 2)}, function(x) {qgamma(x, shape = 2, rate = 2)},  matrix(c(1, -0.8, -0.8, 1), nrow = 2, ncol = 2))} 
nu.c.sims.neg = gen.nu.bridge(nu.c, cir, c.theta, 10000, 100)
nu.c = function(n) { gaussian.cop(10000, function(x) {qgamma(x, shape = 2, rate = 2)}, function(x) {qgamma(x, shape = 2, rate = 2)},  matrix(c(1, .8, .8, 1), nrow = 2, ncol = 2))} 
nu.c.sims.pos = gen.nu.bridge(nu.c, cir, c.theta, 10000, 100)

# Plot.
i = 50
th.q = qsOU(ppoints(nu.o.sims.orig[,i]), c(0, 1, 1))
orig.q = sort(nu.o.sims.orig[,i])
zero.q = sort(nu.o.sims.zero[,i])
neg.q = sort(nu.o.sims.neg[,i])
pos.q = sort(nu.o.sims.pos[,i])
data = data.frame(orig = orig.q,  pos = pos.q, zero = zero.q, neg = neg.q)
data.melt = melt(data, id = "orig")
o.plot = pretty.plot(data.melt) + 
  geom_point(aes(x = orig, y = value, color = variable)) + 
  ylab("Sample Quantiles (varying covariance)") + 
  xlab("Sample Quantiles (approx. unconstrained diffusion)") + 
  geom_abline(slope = 1) +
  scale_y_continuous(breaks = -3:3) + 
  scale_x_continuous(breaks = -3:3) +
  scale_color_manual(name = "", 
                     labels = c(bquote(paste(rho, "=0.8   ")), 
                                bquote(paste(rho, "=0   ")),
                                bquote(paste(rho, "=-0.8   "))
                                ), 
                     values = c(pos = "grey66", 
                                zero = "grey33",
                                neg = "black")
                     ) 
th.q = qsCIR(ppoints(nu.c.sims.orig[,50]), c(1, 1, 1))
orig.q = sort(nu.c.sims.orig[,50])
zero.q = sort(nu.c.sims.zero[,50])
neg.q = sort(nu.c.sims.neg[,50])
pos.q = sort(nu.c.sims.pos[,50])
data = data.frame(orig = orig.q, zero = zero.q, neg = neg.q, pos = pos.q)
data.melt = melt(data, id = "orig")
c.plot = pretty.plot(data.melt) + 
  geom_point(aes(x = orig, y = value, color = variable)) + 
  ylab("Sample Quantiles (varying correlation)") + 
  xlab("Sample Quantiles (approx. unconstrained diffusion)") + 
  geom_abline(slope = 1) +
  scale_y_continuous(breaks = 0:6) + 
  scale_x_continuous(breaks = 0:6) +
  scale_color_manual(name = "",
                     values = c(pos = "grey66", 
                                zero = "grey33",
                                neg = "black")
                     ) 

grid_arrange_shared_legend(o.plot, c.plot)


pdf(file = "/Users/marshall/Documents/senior/thesis/figures/varying_correlation.pdf", width= 10, height = 5, #' see how it looks at this size
    useDingbats=F)
grid_arrange_shared_legend(o.plot, c.plot)
dev.off()

# Dirac delta intuition.
nu.o = function(n) { joint.dirac(n, -1, 1, 0.5, 0.5) }
nu.o.sims.even = gen.nu.bridge(nu.o, ou, o.theta, 10000, 100)
nu.o = function(n) { joint.dirac(n, -1, 1, 0.5, 1) }
nu.o.sims.same= gen.nu.bridge(nu.o, ou, o.theta, 10000, 100)
nu.o = function(n) { joint.dirac(n, -1, 1, 0.5, 0) }
nu.o.sims.diff = gen.nu.bridge(nu.o, ou, o.theta, 10000, 100)

nu.c = function(n) { joint.dirac(n, 1, 3, 0.5, 0.5) }
nu.c.sims.even = gen.nu.bridge(nu.c, cir, c.theta, 10000, 100)
nu.c = function(n) { joint.dirac(n, 1, 3, 0.5, 1) }
nu.c.sims.same= gen.nu.bridge(nu.c, cir, c.theta, 10000, 100)
nu.c = function(n) { joint.dirac(n, 1, 3, 0.5, 0) }
nu.c.sims.diff = gen.nu.bridge(nu.c, cir, c.theta, 10000, 100)

a = 0.0025

data = data.frame(t(nu.o.sims.even))
data$time<- 1:100 / 100
data.melt = melt(data, id = "time")
o.even.plot = pretty.plot(data.melt) + 
  geom_line(aes(x = time, y = value, group = variable), alpha = a) +
  ylab("Value") + 
  xlab("Time") +
  ylim(-3, 3)

data = data.frame(t(nu.o.sims.same))
data$time<- 1:100 / 100
data.melt = melt(data, id = "time")
o.same.plot = pretty.plot(data.melt) + 
  geom_line(aes(x = time, y = value, group = variable), alpha = a) +
  ylab("Value") + 
  xlab("Time") +
  ylim(-3, 3)

data = data.frame(t(nu.o.sims.diff))
data$time<- 1:100 / 100
data.melt = melt(data, id = "time")
o.diff.plot = pretty.plot(data.melt) + 
  geom_line(aes(x = time, y = value, group = variable), alpha = a) +
  ylab("Value") + 
  xlab("Time") +
  ylim(-3,3)

data = data.frame(t(nu.c.sims.even))
data$time<- 1:100 / 100
data.melt = melt(data, id = "time")
c.even.plot = pretty.plot(data.melt) + 
  geom_line(aes(x = time, y = value, group = variable), alpha = a) + 
  ylab("Value") + 
  xlab("Time") +
  ylim(0,5)

data = data.frame(t(nu.c.sims.same))
data$time<- 1:100 / 100
data.melt = melt(data, id = "time")
c.same.plot = pretty.plot(data.melt) + 
  geom_line(aes(x = time, y = value, group = variable), alpha = a) +
  ylab("Value") + 
  xlab("Time") +
  ylim(0,5)

data = data.frame(t(nu.c.sims.diff))
data$time<- 1:100 / 100
data.melt = melt(data, id = "time")
c.diff.plot = pretty.plot(data.melt) + 
  geom_line(aes(x = time, y = value, group = variable), alpha = a) +
  ylab("Value") + 
  xlab("Time") +
  ylim(0,5)

png(file = "/Users/marshall/Documents/senior/thesis/figures/dirac_delta.png", width= 10, height = 15, units = "in", res = 300)
multiplot(o.same.plot, o.even.plot, o.diff.plot, c.same.plot, c.even.plot, c.diff.plot, cols = 2)
dev.off()
