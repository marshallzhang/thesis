# GENERAL
trans.ou = function(start, steps, theta, time = 1) {
  x = sde.sim(T = time, X0 = start, N = steps - 1,
              drift = eval(substitute(expression(
                (- lambda * ((x * sigma) - mu))/sigma
                ), list(lambda = theta[2], mu = theta[1]/theta[2], sigma = theta[3]))),
              drift.x = eval(substitute(expression(
                - lambda
                ), list(lambda = theta[2], mu = theta[1]/theta[2], sigma = theta[3]))),
              drift.xx = eval(substitute(expression(
                0
                ), list(lambda = theta[2], mu =theta[1]/theta[2], sigma = theta[3]))),
              drift.t = expression(0),
              sigma = expression(1),
              sigma.x = expression(0),
              method = "shoji")
  as.numeric(x)
}

gen.nu.bridge.special = function(r.nu, r.diffusion, theta.tilde, theta, M, steps, type) {
  if (type == "exact") {
    raw.bridges.0 = exact.bridges(1, function(n) eta(r.nu(n), theta.tilde), M, r.diffusion, theta.tilde, steps)
  } else if (type == "approx") {
    raw.bridges.0 = gen.nu.bridge(function(n) eta(r.nu(n), theta.tilde), r.diffusion, theta.tilde, M, steps)
  }
  X.dot = raw.bridges.0 - vec.seq(raw.bridges.0[,1], raw.bridges.0[,steps], length.out = steps)
  new.ends = eta(r.nu(M), theta)
  X.dot + vec.seq(new.ends[,1], new.ends[,2], length.out = steps) 
}

gen.nu.bridge.beskos = function(r.nu, r.diffusion, theta.tilde, theta, M, steps, type) {
  if (type == "exact") {
    raw.bridges.0 = exact.bridges(1, function(n) matrix(0, nrow = n, ncol = 2) ,M, r.diffusion, theta.tilde, steps)
  } else if (type == "approx") {
    raw.bridges.0 = gen.nu.bridge(function(n) matrix(0, nrow = n, ncol = 2), r.diffusion, theta.tilde, M, steps)
  }
  new.ends = eta((r.nu(M)), theta)
  difference = new.ends
  raw.bridges.0 + vec.seq(difference[,1], difference[,2], length.out = steps)
}

q.beskos.known = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(2)
  theta1 = theta[1]
  theta2 = theta[2]
  if ( theta2<0) {
    return(NA)
  }
  
  mu = theta1 / theta2
  lambda = theta2
    
  expectation = 0
  
  bridges = if (exact) {
    eta(exact.bridges(1, r.nu, M, r.diffusion, theta.tilde, steps), c(theta,2))
  } else {
    eta(gen.nu.bridge(r.nu, r.diffusion, theta.tilde, M, steps), c(theta, 2))
  }
  
  us = sample(1:steps, M, replace = T)
  
#   expectation2 = sum(rowSums(rowDiffs(bridges) * (alpha(bridges[,1:(steps-1)], c(theta,1))))) - sum((1/2) * (alpha(bridges[cbind(1:M, us)], c(theta,1)))^2)
  expectation = expectation + mean(A(bridges[,100], c(theta,1)) - A(bridges[,1], c(theta, 2)))
  expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], c(theta, 2)))
  
  if (is.nan(expectation)) {
    return(NA)
  }
  -expectation
}

q.beskos.unknown = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  if (theta[2]<0 |theta[3]<0) {
    return(NA)
  }
  
  bridges = if (exact) {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "exact")
  } else {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "approx")
  }
  
  us = sample(1:steps, M, replace = T)
  
  expectation = 0
  new.ends = eta(r.nu(M), theta)
  expectation = expectation + mean(A(new.ends[,2], theta) - A(new.ends[,1], theta))
  new.ends = eta(r.nu(M), theta)
  expectation = expectation + mean(dnorm(new.ends[,2] - new.ends[,1], 0, 1, log = T))
  new.ends = eta(r.nu(M), theta)
  expectation = expectation + mean(log(abs(eta.p(eta.inv(new.ends[,100], theta), theta))))
  expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
  
  if (is.nan(expectation)) {
    return(NA)
  }
  -expectation
}

q.bladt.unknown = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(3)
  if (theta[2]<0 |theta[3]<0) {
    return(NA)
  }
  
  bridges = if (exact) {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "exact")
  } else {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "approx")
  }
  
  us = sample(1:steps, M, replace = T)
  
  expectation = 0
  expectation = expectation + mean(A(bridges[,100], theta) - A(bridges[,1], theta))
  expectation = expectation + mean(dnorm(bridges[,100] - bridges[,1], 0, 1, log = T))
  expectation = expectation + mean(log(abs(eta.p(eta.inv(bridges[,100], theta), theta))))
  expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
  
#   expectation = 0
#   expectation = expectation + mean(g(eta.inv(bridges[,100], theta), theta) - g(eta.inv(bridges[,1], theta), theta))
#   expectation = expectation - (1/2) * mean((bridges[,100] - bridges[,1])^2)
#   expectation = expectation - log(theta[3])
#   expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
  
  if (is.nan(expectation)) {
    return(NA)
  }
  -expectation
}


eta.inv = function(x, theta) {
  x * theta[3]
}

eta = function(x, theta) {
  x / theta[3]
}

eta.p = function(x, theta) {
  1 / theta[3]
}

mu.ou = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  - lambda * (x - mu)
}

alpha = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  
  (- lambda * ((x * sigma) - mu))/sigma
}

a.2.a.prime = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  ((- lambda * ((x * sigma) - mu))/sigma)^2 - lambda
}

# BESKOS
A = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  
  (lambda / sigma) * (x * mu - ((sigma * x^2)/2))
}

# BLADT
s = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  
  (lambda / sigma^2) * (x * mu - ((x^2) /2))
}

g = function(x, theta) {
  sigma = theta[3]
  s(x,theta) * (1/2) * log(sigma)
}