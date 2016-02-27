source("/Users/marshall/Documents/senior/thesis/empirics/general_functions.R")
source("/Users/marshall/Documents/senior/thesis/empirics/test_lamperti/functions.R")
N = 10000

# Beskos, unknown diffusion approx.
theta = c(0,1,2)
sims2 = array(0, dim = c(N, 100))
sims.start = rnorm(N, theta[1]/theta[2], sqrt((theta[3]^2) / (2 * theta[2])))
for (i in 1:N) {
  sims2[i, ] = ou.fast(sims.start[i], 100, theta)
}

nu = function(n) {(as.matrix(sims2[sample(1:N, n, replace = T),c(1,100)]))}

oldpar = c(2,2,2)
for (i in 1:10) {
  set.seed(2)
#   oldpar = optim(par = oldpar,
#         f = function(theta) q.beskos.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
#         method = "CG", control = list(trace = 6, maxit = 30))$par
  oldpar = DEoptim(fn = function(theta) q.beskos.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
                   lower = c(-2.5, 0, 0),
                   upper = c(2.5, 3.5, 4.5),
                   control = DEoptim.control(trace = TRUE, itermax = i*15, strategy= 6))[["optim"]]$bestmem
  print(oldpar)
}

# Bladt, unknown diffusion, exact.
theta = c(0,1,2)
sims2 = array(0, dim = c(N, 100))
sims.start = rnorm(N, theta[1]/theta[2], sqrt((theta[3]^2) / (2 * theta[2])))
for (i in 1:N) {
  sims2[i, ] = ou.fast(sims.start[i], 100, theta)
}

nu = function(n) {(as.matrix(sims2[sample(1:N, n, replace = T),c(1,100)]))}

oldpar = c(2,2,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.bladt.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
        method = "CG", control = list(trace = 6, maxit = 30))$par
  print(oldpar)
}

