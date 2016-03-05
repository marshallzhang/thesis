source("/Users/marshall/Documents/senior/thesis/empirics/general_functions.R")
source("/Users/marshall/Documents/senior/thesis/empirics/test_lamperti_ou/functions.R")
N = 10000

# Stationary.
theta = c(0,1,1)

init = function(n) rnorm(n, theta[1]/theta[2], sqrt((theta[3]^2) / (2 * theta[2])))
nu = make.joint(50000, theta, ou.fast, init, steps = 100)
  
oldpar = c(4,2,2)
set.seed(1)
for (i in 1:10) {
#   oldpar = optim(par = oldpar,
#                  f = function(theta) q.beskos.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
#                  method = "CG",
#                  control = list(trace = 6))$par
  oldpar = DEoptim(fn = function(theta) q.beskos.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
                   lower = c(-3.5, 0, 0),
                   upper = c(3.5, 3.5, 4.5),
                   control = DEoptim.control(trace = TRUE, itermax = 30 + i*2, strategy = 6))[["optim"]]$bestmem
  print(oldpar)
}

# 2Bern - 1.
theta = c(0,1,2)
sims2 = array(0, dim = c(N, 100))
sims.start = 2*rbinom(N, 1, 0.5) - 1
for (i in 1:N) {
  sims2[i, ] = ou.fast(sims.start[i], 100, theta)
}

nu = function(n) {(as.matrix(sims2[sample(1:N, n, replace = T),c(1,100)]))}

oldpar = c(2,2,2)
for (i in 1:10) {
  set.seed(2)
  oldpar = DEoptim(fn = function(theta) q.beskos.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
                   lower = c(-2.5, 0, 0),
                   upper = c(2.5, 3.5, 4.5),
                   control = DEoptim.control(trace = TRUE, itermax = i*15, strategy= 6))[["optim"]]$bestmem
  print(oldpar)
}

# Unif(0,1)
theta = c(0,1,2)
sims2 = array(0, dim = c(N, 100))
sims.start = runif(N)
for (i in 1:N) {
  sims2[i, ] = ou.fast(sims.start[i], 100, theta)
}

nu = function(n) {(as.matrix(sims2[sample(1:N, n, replace = T),c(1,100)]))}

oldpar = c(2,2,2)
for (i in 1:10) {
  set.seed(2)
  oldpar = DEoptim(fn = function(theta) q.beskos.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
                   lower = c(-2.5, 0, 0),
                   upper = c(2.5, 3.5, 4.5),
                   control = DEoptim.control(trace = TRUE, itermax = i*15, strategy= 6))[["optim"]]$bestmem
  print(oldpar)
}

# t
theta = c(0,1,2)
sims2 = array(0, dim = c(N, 100))
sims.start = theta[3] * rt(N, 5) / (sqrt(5/3) * sqrt(2 * theta[2]))
for (i in 1:N) {
  sims2[i, ] = ou.t(sims.start[i], 100, theta)
}

nu = function(n) {(as.matrix(sims2[sample(1:N, n, replace = T),c(1,100)]))}

oldpar = c(2,2,2)
for (i in 1:10) {
  set.seed(2)
  oldpar = DEoptim(fn = function(theta) q.beskos.unknown(theta, trans.ou, oldpar, nu, round(i^1.2)*100, 100, exact = F),
                   lower = c(-2.5, 0, 0),
                   upper = c(2.5, 3.5, 8),
                   control = DEoptim.control(trace = TRUE, itermax = i*15, strategy= 6))[["optim"]]$bestmem
  print(oldpar)
}