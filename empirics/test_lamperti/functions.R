# GENERAL
trans.cir = function(start, steps, theta, time = 1) {
  x = sde.sim(T = time, X0 = start, N = steps - 1,
              drift = eval(substitute(expression((-(lambda * x) / 2) - 
                                   (1 / (2*x)) *(1- (4*lambda*mu)/(sigma^2))), list(lambda = theta[2], mu = theta[1]/theta[2], sigma = theta[3]))),
              drift.x = eval(substitute(expression(- 
                                                     (lambda / 2) +
                                                     (1/(2*x^2)) * (1 - (4*lambda*mu)/(sigma^2))), list(lambda = theta[2], mu = theta[1]/theta[2], sigma = theta[3]))),
              drift.xx = eval(substitute(expression((-1/(x^3)) * (1 - (4*lambda*mu)/(sigma^2))), list(lambda = theta[2], mu =theta[1]/theta[2], sigma = theta[3]))),
              drift.t = expression(0),
              sigma = expression(1),
              sigma.x = expression(0),
              method = "shoji")
  as.numeric(x)
}

gen.nu.bridge.special = function(r.nu, r.diffusion, theta.tilde, theta, M, steps, type) {
  if (type == "exact") {
    raw.bridges.0 = exact.bridges(1, r.nu ,M, r.diffusion, theta.tilde, steps)
  } else if (type == "approx") {
    raw.bridges.0 = gen.nu.bridge(r.nu, r.diffusion, theta.tilde, M, steps)
  }
  Z = eta(raw.bridges.0, theta.tilde)
  ends = Z[,c(1,steps)]
  new.ends = eta(raw.bridges.0[,c(1,steps)], theta)
  difference = new.ends - ends
  Z + vec.seq(difference[,1], difference[,2], length.out = steps) 
#   Z + cbind(vec.seq(difference[,1], rep(0, M), length.out = steps / 2),
#           vec.seq(rep(0, M), difference[,2], length.out = steps/ 2))
}

gen.nu.bridge.beskos = function(r.nu, r.diffusion, theta.tilde, theta, M, steps, type) {
  if (type == "exact") {
    raw.bridges.0 = exact.bridges(1, function(n) matrix(0.3, nrow = n, ncol = 2) ,M, r.diffusion, theta.tilde, steps)
  } else if (type == "approx") {
    raw.bridges.0 = gen.nu.bridge(function(n) matrix(0.3, nrow = n, ncol = 2), r.diffusion, theta.tilde, M, steps)
  }
  new.ends = eta(t(r.nu(M)), theta)
  difference = new.ends - 0.3
  raw.bridges.0 + vec.seq(difference[,1], difference[,2], length.out = steps)
#   Z + cbind(vec.seq(difference[,1], rep(0, M), length.out = steps / 2),
#           vec.seq(rep(0, M), difference[,2], length.out = steps/ 2))
}

q.beskos.known = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(2)
  theta1 = theta[1]
  theta2 = theta[2]
  if (2 * theta1 <= 1 | theta1<0 | theta2<0) {
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
  set.seed(2)
  theta1 = theta[1]
  theta2 = theta[2]
  theta3 = 2
  theta = c(theta1, theta2, 2)
  if (2 * theta1 <= theta3^2 | theta1<0 | theta2<0 |theta3<0) {
    return(NA)
  }
  
  mu = theta1 / theta2
  lambda = theta2
  sigma = theta3
    
  expectation = 0
  
  bridges = if (exact) {
    eta(exact.bridges(1, r.nu, M, r.diffusion, theta.tilde, steps), theta)
#     gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "exact")
  } else {
    eta(gen.nu.bridge(r.nu, r.diffusion, theta.tilde, M, steps), theta)
#     gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "approx")
  }
  
  us = sample(1:steps, M, replace = T)
  
  if (any(bridges[cbind(1:M, us)] < 0.05)) {
#     browser()
    print("Super small!")
  }
  expectation = expectation + mean(A(bridges[,100], theta) - A(bridges[,1], theta))
#   expectation = expectation + mean(dnorm(bridges[,100] - bridges[,1], 0, 1, log = T))
#   expectation = expectation + mean(log(abs(eta.p(eta.inv(bridges[,100], theta), theta))))
#   expectation = expectation + (sum(rowSums(rowDiffs(bridges) * (alpha(bridges[,1:(steps-1)], theta)))) - sum((1/2) * (alpha(bridges[cbind(1:M, us)], theta))^2)) / M
  expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
  
  if (is.nan(expectation)) {
    return(NA)
  } else if (expectation > 5) {
    browser()
  }
  -expectation
}


q.bladt.unknown = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(1)
  theta1 = theta[1]
  theta2 = theta[2]
  theta3 = theta[3]
  if (2 * theta1 <= theta3^2 | theta1<0 | theta2<0 |theta3<0) {
    return(NA)
  }
  
  mu = theta1 / theta2
  lambda = theta2
  sigma = theta3
    
  expectation = 0
  
  bridges = if (exact) {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "exact")
  } else {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "approx")
  }
  
  us = sample(1:steps, M, replace = T)
  
  expectation = expectation + mean(g(eta.inv(bridges[,100], theta), theta) - g(eta.inv(bridges[,1], theta), theta))
  expectation = expectation - (1/2) * mean((bridges[,100] - bridges[,1])^2)
  expectation = expectation - mean(log(sigma * sqrt(eta.inv(bridges[,100], theta))))
  expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
  
  if (is.nan(expectation)) {
    return(NA)
  } else if (expectation > 5) {
    browser()
  }
  -expectation
}

eta.inv = function(x, theta) {
  (theta[3]^2 / 4) * x^2
}

eta = function(x, theta) {
  2 * sqrt(x) / theta[3]
}

eta.p = function(x, theta) {
  ((1 / 2) * (2 * x^(-1/2))) / theta[3]
}

alpha = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  (-(lambda * x) / 2) - (1 / (2*x)) *(1- (4*lambda*mu)/(sigma^2))
}

a.2.a.prime = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  (
    (-(lambda * x) / 2) - ((1/(2*x)) * (1 - (4*lambda*mu)/(sigma^2)))
  )^2 - 
  (lambda / 2) +
  (1/(2*x^2)) * (1 - (4*lambda*mu)/(sigma^2))
}

# BESKOS
A = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  
  (1/2) *
    ( ((-(x^2) * lambda)/2) -
      log(x) +
      ((4 * lambda * mu * log(x)) / (sigma^2))
  )
}

# BLADT
s = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  
  (-lambda * (x - mu * log(x)))/(sigma^2)
}

g = function(x, theta) {
  sigma = theta[3]
  s(x,theta) * (1/2) * log(sigma * sqrt(x))
}