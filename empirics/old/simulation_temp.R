source("/Users/marshall/Documents/senior/thesis/empirics/simulation_functions.R")

#
# 
# STATIONARY SOLUTIONS.
#
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

nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, exp(-1)/2, exp(-1)/2, 1/2), nrow = 2, ncol = 2))} 

o.stationary.approx = timeit(nu.o.sims <- gen.nu.bridge(nu.o, ou, theta, N, 100))

o.stationary.exact = timeit(nu.o.sims.exact <- exact.bridges(1, nu.o, 10000, ou, theta, 100))

# CIR process.
theta = c(1, 1, 1)
c.sims = array(0, dim = c(N, 100))
c.sims.start = rgamma(N, shape = 2, rate = 2)
basic = timeit(for (i in 1:N) {
#   print(i)
  c.sims[i, ] = cir(c.sims.start[i], 100, theta)
})

nu.c = function(n) { mcmc(c(1, 1), 20 * N + 2000, 2000, everyother = 20, trace = 1000, d.posterior = cir.joint, r.proposal = function(n, current) {rgamma(2, rate = 2, shape = 2)}, special = "gamma") }
sims = timeit(c.sims2 <- nu.c(202000))
c.stationary.approx = timeit(nu.c.sims <- gen.nu.bridge(nu.c, cir, theta, 10000, 100))

c.start = nu.c(202000)
c.stationary.exact = timeit(nu.c.sims.exact <- exact.bridges(1, function(n) t(as.matrix(c.start[sample(1:10000, n, replace = T),])), 10000, cir, theta, 100))

# Plot.
th.q = qsOU(ppoints(1:10000), c(0, 1, 1))
data.q = sort(o.sims[1:10000,50])
my.data.q = sort(nu.o.sims[1:10000,50])
exact.data.q = sort(nu.o.sims.exact.cut[,50])
data = data.frame(th = data.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
o.plot.approx = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Approximate bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = -4:4) + 
              scale_x_continuous(breaks = -4:4) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-4,4), xlim=c(-4,4))+
              theme(legend.position = "none")

data = data.frame(th = data.q, my.dat = exact.data.q)
data.melt = melt(data, id = "th")
o.plot.exact = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Exact bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = -4:4) + 
              scale_x_continuous(breaks = -4:4) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-4,4), xlim=c(-4,4))+
              theme(legend.position = "none")

data.q = sort(c.sims[,50])
my.data.q = sort(nu.c.sims[,50])
exact.data.q = sort(nu.c.sims.exact[,50])
data = data.frame(th = data.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
c.plot.approx = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Approximate bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = 0:8) + 
              scale_x_continuous(breaks = 0:8) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(0,8), xlim=c(0,8))+
              theme(legend.position = "none")

data = data.frame(th = data.q, my.dat = exact.data.q)
data.melt = melt(data, id = "th")
c.plot.exact = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Exact bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = 0:8) + 
              scale_x_continuous(breaks = 0:8) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(0,8), xlim=c(0,8))+
              theme(legend.position = "none")

multiplot(o.plot.approx, o.plot.exact, c.plot.approx, c.plot.exact, cols = 2)

figure1 = list(real.ou = o.sims, 
               approx.ou = nu.o.sims,
               exact.ou = nu.o.sims.exact,
               real.cir = c.sims,
               approx.cir = nu.c.sims,
               exact.cir = nu.c.sims.exact)

save(figure1, file = "figure1data.Rdata")

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/invariant_densities.pdf", width= 10, height = 10, #' see how it looks at this size
    useDingbats=F)
multiplot(o.plot.approx, o.plot.exact, c.plot.approx, c.plot.exact, cols = 2)
dev.off()


#
#
# MORE EXTREME STARTING POINTS
#
#

# Mean 1 Normal distribution.
theta = c(0, 1, 1)

o.sims.m1 = array(0, dim = c(N, 100))
o.sims.start = rnorm(N, mean = 1, sd = sqrt(1/2))
for (i in 1:N) {
  print(i)
  o.sims.m1[i, ] = ou(o.sims.start[i], 100, theta)
}

o.sims2 = array(0, dim = c(N, 100))
o.sims.start2 = rnorm(N, mean = 1, sd = sqrt(1/2))
for (i in 1:N) {
#   print(i)
  o.sims2[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)])}
sims.m1.time = timeit(nu.o.sims.m1 <- gen.nu.bridge(nu.o, ou, theta, N, 100))

nu.o = function(n) {t(as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)]))}
sims.m1.exact.time = timeit(nu.o.sims.exact.m1 <- exact.bridges(5, nu.o, 10000, ou, theta, 100))

# Bernouilli starting distribution.
o.sims.binom = array(0, dim = c(N, 100))
o.sims.start = (rbinom(N, 1, 0.5) * 2) - 1
for (i in 1:N) {
  print(i)
  o.sims.binom[i, ] = ou(o.sims.start[i], 100, theta)
}

o.sims2 = array(0, dim = c(N, 100))
o.sims.start = (rbinom(N, 1, 0.5) * 2) - 1
for (i in 1:N) {
#   print(i)
  o.sims2[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)])}
sims.binom.time = timeit(nu.o.sims.binom <- gen.nu.bridge(nu.o, ou, theta, N, 100))

nu.o = function(n) {t(as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)]))}
sims.exact.binom.time = timeit(nu.o.sims.exact.rbinom <- exact.bridges(1, nu.o, 10000, ou, theta, 100))

# Expo starting distribution.
o.sims.exp = array(0, dim = c(N, 100))
o.sims.start = rexp(N, 2)
exp.time = timeit(
for (i in 1:N) {
#   print(i)
  o.sims.exp[i, ] = ou(o.sims.start[i], 100, theta)
})

o.sims2 = array(0, dim = c(N, 100))
o.sims.start = rexp(N, 2)
for (i in 1:N) {
#   print(i)
  o.sims2[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)])}
sims.exp.time = timeit(nu.o.sims.exp <- gen.nu.bridge(nu.o, ou, theta, N, 100))

nu.o = function(n) {t(as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)]))}
sims.exact.exp.time = timeit(nu.o.sims.exact.exp<- exact.bridges(1, nu.o, 10000, ou, theta, 100))

table1data = list(binom.real = o.sims.binom,
                  binom.approx = nu.o.sims.binom,
                  binom.exact = nu.o.sims.exact.rbinom,
                  m1.real = o.sims.m1,
                  m1.approx = nu.o.sims.m1,
                  m1.exact = nu.o.sims.exact.m1,
                  exp.real = o.sims.exp,
                  exp.approx = nu.o.sims.exp,
                  exp.exact = nu.o.sims.exact.exp)

save(table1data, file = "/Users/marshall/Documents/senior/thesis/figures/table1data.Rdata")

#
#
# EM WITH UNIT DIFFUSION
#
#

# Stationary.
theta = c(0, 1, 1)

o.sims = array(0, dim = c(N, 100))
o.sims.start = rnorm(N, mean = 0, sd = sqrt(1/2))
for (i in 1:N) {
  print(i)
  o.sims[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, exp(-1)/2, exp(-1)/2, 1/2), nrow = 2, ncol = 2))} 

DEoptim(f = function(theta) q.unit(theta, ou, c(2,2,1), nu.o, 10, 100, exact = F),
      lower = c(-1, 0),
      upper = c(5, 4),
      control = DEoptim.control(trace = 1))

oldpar = c(2,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = T),
        method = "CG",
        control = list(trace = 6))$par
  print(oldpar)
}

nu.c = function(n) { mcmc(c(1, 1), 20 * n + 2000, 2000, everyother = 20, trace = 1000, d.posterior = cir.joint, r.proposal = function(n, current) {rgamma(2, rate = 2, shape = 2)}, special = "gamma") }
c.start = nu.c(10000)
real.nu.c = function(n) t(as.matrix(c.start[sample(1:10000, n, replace = T),]))
oldpar = c(3,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, cir, c(oldpar,1), real.nu.c, round(i^1.5)*100, 100, exact = T),
        method = "CG",
        control = list(trace = 6))$par
  print(oldpar)
}

q.unit(c(1,1), ou, c(2,2,1), nu.o, 5000, 100, exact = F)
