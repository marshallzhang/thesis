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

gen.nu.bridge.special = function(r.nu, r.diffusion, theta.tilde, theta, M, steps, type, time = 1) {
  if (type == "exact") {
    raw.bridges.0 = exact.bridges(1, function(n) eta(r.nu(n), theta.tilde), M, r.diffusion, theta.tilde, steps, set.start = c(-Inf, Inf), time)
  } else if (type == "approx") {
    raw.bridges.0 = gen.nu.bridge(function(n) eta(r.nu(n), theta.tilde), r.diffusion, theta.tilde, M, steps, set.start = c(-Inf, Inf), time)
  }
  X.dot = raw.bridges.0 - vec.seq(raw.bridges.0[,1], raw.bridges.0[,steps], length.out = steps)
  new.ends = eta(r.nu(M), theta)
  X.dot + vec.seq(new.ends[,1], new.ends[,2], length.out = steps) 
}

q.beskos.unknown = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, horizons, tot.time, exact = F) {
  if (theta[2]<0 |theta[3]<0) {
    return(NA)
  }
  
  n = length(r.nu(1))

  expectation = 0
  new.ends = eta(r.nu(M), theta)
  expectation = expectation + mean(A(new.ends[,n], theta) - A(new.ends[,1], theta))
  
  new.ends = eta(r.nu(M), theta)
  mini.sum = 0
  
  for (i in 1:M) {
    this.new.end = new.ends[i, ]
    mini.sum = mini.sum + sum(dnorm(diff(this.new.end), 0, sqrt(diff(horizons) * tot.time), log = T))
  }
  
  expectation = expectation + (mini.sum / M)
  
  new.ends = eta(r.nu(M), theta)
  expectation = expectation - (n - 1) * mean(log(theta[3]))
  
  us = sample(1:steps, M, replace = T)
  samples = vector("list", length(horizons) - 1)
  for (i in 1:(length(horizons) - 1)) {
    samples[[i]] = array(NA, M)
  }
  counters = array(1, length(horizons) - 1)
  
  for (u in us) {
    nu.i = findInterval((u / steps), horizons)
    if (nu.i >= n) {
      samples[[nu.i-1]][counters[nu.i-1]] = a.2.a.prime(r.nu(1)[n], theta)
      counters[nu.i-1] = counters[nu.i-1] + 1
    } else {
      time.diff = horizons[nu.i + 1] - horizons[nu.i]
      bridge.steps = ceiling(time.diff * steps)
      bridge = if (exact) {
        gen.nu.bridge.special(function(n) t(r.nu(n)[nu.i:(nu.i+1)]), r.diffusion, theta.tilde, theta, 1, bridge.steps, "exact", time = tot.time * time.diff)
      } else {
        gen.nu.bridge.special(function(n) t(r.nu(n)[nu.i:(nu.i+1)]), r.diffusion, theta.tilde, theta, 1, bridge.steps, "approx", time = tot.time * time.diff)
      }
      ind = u - horizons[nu.i] * steps
      if (ind == 0) {
        ind = 1
      } else {
        ind = ind
      }
      samples[[nu.i]][counters[nu.i]] = a.2.a.prime(bridge[ind], theta)
      counters[nu.i] = counters[nu.i] + 1
    }
  }
  
  expectation = expectation - (1/2) * sum(sapply(samples, mean, na.rm = T) * diff(horizons) * tot.time)
  
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