# GENERAL
trans.ou.fast = function(start, steps, theta, time = 1) {
  Y <- numeric(steps)
  Y[1] = start
  Z = rnorm(steps - 1)
  Dt <- time / (steps - 1)
  for (i in 1:(steps - 1)) {
    Y[i+1] = Y[i] + ((theta[1] / theta[3]) - theta[2] * Y[i]) * Dt + sqrt(Dt) * Z[i]
  }
  as.numeric(Y)
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

q.beskos.unknown.single = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  if (theta[2]<0 |theta[3]<0) {
    return(NA)
  }
  
  expectation = 0
  new.ends = eta(r.nu(M), theta)
  expectation = expectation + mean(A(new.ends[,2], theta) - A(new.ends[,1], theta))
  new.ends = eta(r.nu(M), theta)
  expectation = expectation + mean(dnorm(new.ends[,2] - new.ends[,1], 0, 1, log = T))
  new.ends = eta(r.nu(M), theta)
  expectation = expectation - mean(log(theta[3]))
  
  us = sample(1:steps, M, replace = T)
  bridges = if (exact) {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "exact")
  } else {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "approx")
  }
  
  expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
  
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