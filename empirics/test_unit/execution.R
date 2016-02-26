source("/Users/marshall/Documents/senior/thesis/empirics/general_functions.R")
source("/Users/marshall/Documents/senior/thesis/empirics/test_unit/functions.R")
N = 10000

# Stationary.
theta = c(0, 1, 1)

nu.o = function(n) {rmvnorm(n, c(0,0), matrix(c(1/2, exp(-1)/2, exp(-1)/2, 1/2), nrow = 2, ncol = 2))} 

oldpar = c(2,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = F),
        method = "CG")$par
  print(oldpar)
}

# N(1,1/2)
theta = c(0, 1, 1)

o.sims2 = array(0, dim = c(N, 100))
o.sims.start = rnorm(N, mean = 1, sd = sqrt(1/2))
for (i in 1:N) {
  o.sims2[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)])}

print("N(1,1/2) Approximate")
oldpar = c(2,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = F),
        method = "CG")$par
  print(oldpar)
}

print("N(1,1/2) Exact")
nu.o = function(n) {t(as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)]))}
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = T),
        method = "CG")$par
  print(oldpar)
}

# 2Bern - 1
o.sims2 = array(0, dim = c(N, 100))
o.sims.start = (rbinom(N, 1, 0.5) * 2) - 1
for (i in 1:N) {
#   print(i)
  o.sims2[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)])}

print("2Bern - 1 Approximate")
oldpar = c(2,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = F),
        method = "CG")$par
  print(oldpar)
}

nu.o = function(n) {t(as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)]))}
print("2Bern - 1 Exact")
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = T),
        method = "CG")$par
  print(oldpar)
}

# Expo(2)
o.sims2 = array(0, dim = c(N, 100))
o.sims.start = rexp(N, 2)
for (i in 1:N) {
#   print(i)
  o.sims2[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)])}

print("Expo(2) Approximate")
oldpar = c(2,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = F),
        method = "CG")$par
  print(oldpar)
}

print("Expo(2) Exact")
nu.o = function(n) {t(as.matrix(o.sims2[sample(1:N, n, replace = T),c(1,100)]))}
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = T),
        method = "CG")$par
  print(oldpar)
}

o.sims2 = array(0, dim = c(100, 100))
o.sims.start = (rnorm(N, 0, 0.5))
for (i in 1:100) {
#   print(i)
  o.sims2[i, ] = ou(o.sims.start[i], 100, theta)
}

nu.o = function(n) {as.matrix(o.sims2[sample(1:100, n, replace = T),c(1,100)])}

print("Small Sample Approximate")
oldpar = c(2,2)
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = F),
        method = "CG")$par
  print(oldpar)
}

nu.o = function(n) {t(as.matrix(o.sims2[sample(1:100, n, replace = T),c(1,100)]))}
print("Small Sample Exact")
for (i in 1:10) {
  oldpar = optim(par = oldpar,
        f = function(theta) q.unit(theta, ou, c(oldpar,1), nu.o, round(i^1.5)*100, 100, exact = T),
        method = "CG")$par
  print(oldpar)
}