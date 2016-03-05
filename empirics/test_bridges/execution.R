source("/Users/marshall/Documents/senior/thesis/empirics/general_functions.R")
source("/Users/marshall/Documents/senior/thesis/empirics/test_bridges/functions.R")
N = 10000
steps = 100

# QQ plot.

# OU process.
theta = c(0, 1, 1)

o.sims = array(0, dim = c(N, steps))
for (i in 1:N) {
  o.sims[i, ] = ou.fast(rnorm(1, theta[1] / theta[2], sqrt(theta[3]^2 / (2 * theta[2]))), steps, theta)
}

nu.o = make.joint(50000, theta, ou.fast, function(n) rnorm(n, theta[1] / theta[2], sqrt(theta[3]^2 / (2 * theta[2]))))

o.stationary.approx = timeit(nu.o.sims <- gen.nu.bridge(nu.o, ou.fast, theta, N, 100))

o.stationary.exact = timeit(nu.o.sims.exact <- exact.bridges(10, nu.o, N, ou.fast, theta, 100))

ks.test(nu.o.sims[,50], o.sims[,50])
ks.test(nu.o.sims.exact[,50], o.sims[,50])

# Plot.
edge = 3
skip = 1
data.q = sort(o.sims[1:N,50])
my.data.q = sort(nu.o.sims[1:N,50])
exact.data.q = sort(nu.o.sims.exact[,50])
data = data.frame(th = data.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
o.plot.approx = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Approximate bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge,skip)) + 
              scale_x_continuous(breaks = seq(-edge,edge,skip)) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge), xlim=c(-edge,edge))+
              theme(legend.position = "none")

data = data.frame(th = data.q, my.dat = exact.data.q)
data.melt = melt(data, id = "th")
o.plot.exact = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Exact bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge,skip)) + 
              scale_x_continuous(breaks = seq(-edge,edge,skip)) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge), xlim=c(-edge,edge))+
              theme(legend.position = "none")

multiplot(o.plot.approx, o.plot.exact, cols = 2)

save(figure1, file = "figure1data.Rdata")

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/invariant_densities.pdf", width= 10, height = 5, #' see how it looks at this size
    useDingbats=F)
multiplot(o.plot.approx, o.plot.exact, cols = 2)
dev.off()

# Alternative starting distributions.

# Expo
theta = c(0, 1, 1)

init = function(n) rexp(n)
o.sims = array(0, dim = c(N, steps))
for (i in 1:N) {
  o.sims[i, ] = ou.fast(init(1), steps, theta)
}

nu.o = make.joint(50000, theta, ou.fast, init)

sims.m1.time = timeit(nu.o.sims <- gen.nu.bridge(nu.o, ou, theta, N, 100))

sims.m1.exact.time = timeit(nu.o.sims.exact <- exact.bridges(10, nu.o, N, ou, theta, 100))

ks.test(nu.o.sims[,50], o.sims[,50])
ks.test(nu.o.sims.exact[,50], o.sims[,50])

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
sims.exact.binom.time = timeit(nu.o.sims.exact.rbinom <- exact.bridges(1, nu.o, N, ou, theta, 100))

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
sims.exact.exp.time = timeit(nu.o.sims.exact.exp<- exact.bridges(1, nu.o, N, ou, theta, 100))

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