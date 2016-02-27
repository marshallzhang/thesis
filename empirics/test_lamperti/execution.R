source("/Users/marshall/Documents/senior/thesis/empirics/general_functions.R")
source("/Users/marshall/Documents/senior/thesis/empirics/test_lamperti/functions.R")
N = 10000

# Beskos, known diffusion, approx.

c.sims2 = array(0, dim = c(N, 100))
c.sims.start = rgamma(N, shape = 1/2, rate = 1/2)
for (i in 1:N) {
  c.sims2[i, ] = cir.fast(c.sims.start[i], 100, c(3,1,2))
}

nu.c = function(n) {as.matrix(c.sims2[sample(1:N, replace = T),c(1,100)])}

oldpar = c(4,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.beskos.known(theta, cir.fast, c(oldpar,2), nu.c, round(i^1.5)*100, 100, exact = F),
        method = "CG", control = list(trace = 6, maxit = 30))$par
  print(oldpar)
}

# Beskos, known diffusion, exact.
oldpar = c(4,2)
c.start = nu.c(N)
nu.c.exact = function(n) t(as.matrix(c.start[sample(1:N, n, replace = T),]))
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.beskos.known(theta, cir.fast, c(oldpar,1), nu.c.exact, round(i^1.5)*100, 100, exact = T),
        method = "CG", control = list(trace = 6, maxit = 30))$par
  print(oldpar)
}


# Beskos, unknown diffusion approx.
oldpar = c(4,2,2)

c.sims2 = array(0, dim = c(N, 100))
c.sims.start = rgamma(N, shape = 1/2, rate = 1/2)
for (i in 1:N) {
  c.sims2[i, ] = cir.fast(c.sims.start[i], 100, c(3,1,2))
}

nu.c = function(n) {as.matrix(c.sims2[sample(1:100, n, replace = T),c(1,100)])}
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.beskos.unknown(theta, cir.fast, oldpar, nu.c, round(i^1.5)*100, 100, exact = F),
        method = "CG", control = list(trace = 6))$par
  print(oldpar)
}

# Beskos, unknown diffusion, exact.
oldpar = c(4,2,2)
c.start = nu.c(N)
nu.c.exact = function(n) t(as.matrix(c.start[sample(1:N, n, replace = T),]))
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.beskos.unknown(theta, cir.fast, oldpar, nu.c.exact, round(i^1.5)*100, 100, exact = T),
        method = "CG", control = list(trace = 6))$par
  print(oldpar)
}

c.start = nu.c(N)
nu.c.exact = function(n) t(as.matrix(c.start[sample(1:N, n, replace = T),]))

oldpar = c(1,1,1)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.bladt.unknown(theta, cir.fast, oldpar, nu.c.exact, round(i^1.5)*100, 100, exact = T),
        method = "CG", control = list(trace = 6))$par
  print(oldpar)
}