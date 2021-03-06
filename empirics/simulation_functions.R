

#
# GENERATE NU BRIDGE
#

gen.double.nu.bridge = function(r.diffusion, theta, start, end, steps) {
    y.1 = r.diffusion(start, steps, theta)
    y.2 = r.diffusion(end, steps, theta)
    
    y.1.1 = y.1
    y.2.1 = rev(y.2)
    difference.1 = y.1.1 - y.2.1
    
    y.1 = r.diffusion(end, steps, theta)
    y.2 = r.diffusion(start, steps, theta)
    
    y.1.2 = y.1
    y.2.2 = rev(y.2)
    difference.2 = y.1.2 - y.2.2
    
    counter = 0
    
    # If they don't intersect, re-generate.
    while (no.cross(y.1.1,y.2.1) | no.cross(y.1.2,y.2.2)){
      counter = counter + 1
      if (counter > 1000) {
        stop("Too many tries.")
      }

      y.1 = r.diffusion(start, steps, theta)
      y.2 = r.diffusion(end, steps, theta)
      
      y.1.1 = y.1
      y.2.1 = rev(y.2)
      difference.1 = y.1.1 - y.2.1
      
      y.1 = r.diffusion(end, steps, theta)
      y.2 = r.diffusion(start, steps, theta)
      
      y.1.2 = y.1
      y.2.2 = rev(y.2)
      difference.2 = y.1.2 - y.2.2
      
    }
      
    # Otherwise, find where they intersect and create the bridge.  
    if (difference.1[1] > 0) {
      iota = which(difference.1 < 0)[1] 
      b1 = c(y.1.1[1:(iota - 1)], y.2.1[iota:length(y.2.1)])
    } else {
      iota = which(difference.1 > 0)[1]
      b1 = c(y.1.1[1:(iota - 1)], y.2.1[iota:length(y.2.1)])
      }
    
    if (difference.2[1] > 0) {
      iota = which(difference.2 < 0)[1] 
      b2 = rev(c(y.1.2[1:(iota - 1)], y.2.2[iota:length(y.2.2)]))
    } else {
      iota = which(difference.2 > 0)[1]
      b2 = rev(c(y.1.2[1:(iota - 1)], y.2.2[iota:length(y.2.2)]))
      }
    
    rbind(b1, b2)
}

# Check crossing
no.cross = function(x, y) {
  if (length(x) == length(y)) {
    diff = x - y
    return(all(diff< 0) | all(diff> 0))
  } else {
    stop("Different lengths of paths.")
  }
}

# Helper to generate a bridge.
gen.nu.helper = function(r.diffusion, theta, start, end, steps, flip = F) {
    y.1 = r.diffusion(start, steps, theta)
    y.2 = rev(r.diffusion(end, steps, theta))

    counter = 0
    # If they don't intersect, re-generate.
    while (no.cross(y.1, y.2)) {
      counter = counter + 1
      if (counter > 5000) {
        stop("Too many tries.")
      }
      y.1 = r.diffusion(start, steps, theta)
      y.2 = rev(r.diffusion(end, steps, theta))
    }
      
    # Otherwise, find where they intersect and create the bridge.  
    difference = y.1 - y.2
    if ((difference[1] > 0)) {
      iota = which(difference < 0)[1] 
      c(y.1[1:(iota - 1)], y.2[iota:length(y.2)])
    } else {
      iota = which(difference > 0)[1]
      c(y.1[1:(iota - 1)], y.2[iota:length(y.2)])
      }
}

# Generating a nu bridge.
gen.nu.bridge = function(r.nu, r.diffusion, theta, N, steps, set.start = c(-Inf,-Inf)) {
  
  bridges = array(0, dim = c(N, steps))
  xy = r.nu(N)
  if (N == 1 & all(is.finite(set.start))) {
    xy = t(matrix(set.start))
  }
  
  for (i in 1:N){
  # Generate two independent diffusions, one starting at start and the other at end.
   tryCatch({
      bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy[i, 1], xy[i, 2], steps)
    }, error = function(e) {
      tryCatch({
#         browser()
        xy1 = r.nu(200)[200,]
        bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps)
      }, error = function(e) {
        tryCatch({
          xy1 = r.nu(200)[200,]
          bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps)
        }, error = function(e) {
          tryCatch({
            xy1 = r.nu(200)[200,]
            bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps)
          }, error = function(e) {
            tryCatch({
              xy1 = r.nu(200)[200,]
              bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps)
            }, error = function(e) {
               browser()
               print(c(i, xy[i,1], xy[i,2]))
            })
          })
        })
      })
    })
#     print(i)
  }
  bridges
}

#
#
#
# OU
#
#
#

# Generate OU process.
ou = function(start, steps, theta, time = 1) {
  x = sde.sim(T= time, X0 = start, N = steps - 1,
          model = "OU", theta = theta, method = "euler")
  x
}

ou.fast = function(start, steps, theta, time = 1) {
  x = sde.sim(T = time, X0 = start, N = steps - 1,
              drift = eval(substitute(expression(theta1 -  (-theta2) * x), list(theta1 = theta[1], theta2 = theta[2]))),
              drift.x = eval(substitute(expression(-theta2), list(theta2 = theta[2]))),
              sigma = eval(substitute(expression(theta3), list(theta3 = theta[3]))),
              sigma.x = expression(0),
              method = "milstein")
  as.numeric(x)
}

set.seed(123)
d <- expression(-5 * x)
s <- expression(3.5) 
sde.sim(X0=10,drift=d, sigma=s) -> X

#
#
#
# CIR
#
#
#


cir = function(start, steps, theta, time = 1) {
  x = sde.sim(T = time, X0 = start, N = steps - 1,
              model = "CIR", theta = theta )
  as.numeric(x)
}

cir.joint = function(current) {
  x = current[1]
  y = current[2]
  c = 2 / (1-exp(-1))
  q = 1
  u = c * x *exp(-1)
  v = c * y
  trans = c * exp(-u-v) * (v/u)^(q/2) * besselI(2 * sqrt(u * v), q)
  log(dgamma(x, shape = 2, rate = 2) * trans)
}


#
#
#
# EXACT SIMULATION
#
#
#

# hat rho
d.rho = function(x, M, r.nu, r.diffusion, theta, N, steps) {
  Ts = array(1, dim = M)
  for (i in 1:M) {
    init.1 = r.diffusion(x[1,steps], steps, theta)[steps]
    init.2 = r.diffusion(x[2,steps], steps, theta)[steps]
    hitting1 = as.vector(r.diffusion(init.1, steps, theta))
    hitting2 = as.vector(r.diffusion(init.2, steps, theta))
    while (no.cross(x[1,], hitting1) | no.cross(x[2,], hitting2)) {
      Ts[i] = Ts[i] + 1
      init.1 = r.diffusion(x[1,steps], steps, theta)[steps]
      init.2 = r.diffusion(x[2,steps], steps, theta)[steps]
      hitting1 = as.vector(r.diffusion(init.1, steps, theta))
      hitting2 = as.vector(r.diffusion(init.2, steps, theta))
    }
  }
  c(mean(Ts), !no.cross(x[1,], hitting1), !no.cross(x[2,], hitting2))
}

rho = function(x, M, r.diffusion, theta, steps) {
  Ts = array(1, dim = M)
  for (i in 1:M) {
    init.1 = r.diffusion(x[steps], steps, theta)[steps]
    hitting = as.vector(r.diffusion(init.1, steps, theta))
    while (no.cross(x, hitting)) {
      Ts[i] = Ts[i] + 1
      init.1 = r.diffusion(x[steps], steps, theta)[steps]
      hitting = as.vector(r.diffusion(init.1, steps, theta))
    }
  }
  mean(Ts)
}

exact.bridges = function(M, r.nu, samples, r.diffusion, theta, steps) {
  bridges = matrix(0, nrow = samples, ncol = 100)
  bridges[1, ] = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps)
  for (i in 2:samples) {
    tryCatch({
      start.end = r.nu(1)
      proposal.1 = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, set.start = start.end)
      proposal.2 = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, set.start = start.end)
    }, error = function(e) {
      bridges[i, ] = bridges[i-1,]
    })
    r = rho(proposal.2, M, r.diffusion, theta, steps) / rho(proposal.1, M, r.diffusion, theta, steps)
    alpha = min(1, r)
#     print(paste(i,alpha))
    if (runif(1) < alpha) {
      bridges[i,] = proposal.2
    } else {
      bridges[i,] = proposal.1
    }
  }
  bridges
}

exact.nu.double.bridge = function(M, r.nu, samples, r.diffusion, theta, steps) {
  bridges = vector("list", samples)
  xy = r.nu(1)
  bridges[[1]] = gen.double.nu.bridge(r.diffusion, theta, xy[1], xy[2], steps)
  for (i in 2:samples) {
    xy = r.nu(1)
    proposal = tryCatch({
      gen.double.nu.bridge(r.diffusion, theta, xy[1], xy[2], steps)
    }, error = function(e) {
     bridges[[i - 1]]})
    rho.im1 = d.rho(bridges[[i-1]], M, r.nu, r.diffusion, theta, N, steps)
    rho.i = d.rho(proposal, M, r.nu, r.diffusion, theta, N, steps)
    s.cov = matrix(c(1/2, exp(-1)/2, exp(-1)/2, 1/2), nrow = 2, ncol = 2)
    r = (rho.i[1] / rho.im1[1]) * (dmvnorm(proposal[,1], mean = c(0,0), sigma = s.cov) / dmvnorm(bridges[[i-1]][,1], mean = c(0,0), sigma = s.cov))
    alpha = min(1, r)
    print(paste(i,alpha))
    if (runif(1) < alpha) {
      bridges[[i]] = proposal
    } else {
#       xy = r.nu(1)
#       second.proposal = gen.double.nu.bridge(r.diffusion, theta, xy[1], xy[2], steps)
#       rho.ip1 = d.rho(second.proposal, M, r.nu, r.diffusion, theta, N, steps)
#       alpha.2 = (rho.ip1[1] / rho.im1[1]) * ((1 - min(1,(rho.i[1] / rho.ip1[1]))) / (1 - r))
#       print(paste(i, round(r, digits = 2), round(min(1,alpha.2), digits = 2)))
#       if (runif(1) < min(1,alpha.2)) {
#         bridges[[i]] = second.proposal
#       } else {
        bridges[[i]] = bridges[[i-1]]
#       }
    }
  }
  ret = matrix(0, nrow = samples, ncol = 100)
  for (i in 1:samples) {
    ret[i, ] = bridges[[i]][1, ]
  }
  ret
}

exact.nu.bridge = function(M, r.hitting.init, samples, ...) {
  bridges = array(0, dim = c(samples, 100))
  bridges[1,] = gen.nu.bridge(...)
  for (i in 2:samples) {
    proposal = gen.nu.bridge(...)
    rho.im1 = rho(bridges[i-1,], M, r.hitting.init, ...)
    rho.i = rho(proposal, M, r.hitting.init, ...)
    r = rho.i / rho.im1
    alpha = min(1, r)
    if (runif(1) < alpha) {
      bridges[i,] = proposal
    } else {
      second.proposal = gen.nu.bridge(...)
      rho.ip1 = rho(second.proposal, M, r.hitting.init, ...)
      alpha.2 = (rho.ip1 / rho.im1) * ((1 - min(1,(rho.i / rho.ip1))) / (1 - r))
      print(paste(i, round(r, digits = 2), round(min(1,alpha.2), digits = 2)))
      if (runif(1) < min(1,alpha.2)) {
        bridges[i, ] = second.proposal
      } else {
          bridges[i,] = bridges[i-1,]
      }
    }
  }
  bridges
}

#
#
#
# GENERAL FUNCTIONS LIKE PLOTTING AND MH
#
#
#

# Vectorized seq.

composite = function(f,g) function(...) f(g(...))
vec.seq = composite(t, Vectorize(seq.default, vectorize.args = c("from", "to")))

# Metropolis-Hastings with normal proposals.
mcmc.mh = function(current, d.posterior, r.proposal, special = "") {
  
    # Propose and find density according to random walk.
    proposal = r.proposal(length(current), current)
    if (special == "gamma") {
      alpha = min(1, exp(d.posterior(proposal) - d.posterior(current) + sum(dgamma(current, shape = 2, rate = 2, log =T )) - sum(dgamma(proposal, shape = 2, rate = 2, log = T))))
    } else if (special == "") {
      alpha = min(1, exp(d.posterior(proposal) - d.posterior(current)))
    }
    
    # Accept or reject.
    if (runif(1) < alpha) {
      proposal
    } else {
      current
    }
  
}

# MCMC wrapper function.
mcmc = function(start, iterations, burn, everyother, trace, ...) {
  
  accept = 0
  # Initialize storage.
  draws = array(dim = c(iterations + 1, length(start)))
  draws[1, ] = start
  
  # Iterate.
  for (i in 1:iterations) {
#     if (trace != 0) if (i %% trace == 0) print(paste("Iteration", i, "of", iterations))
    draws[i + 1, ] = mcmc.mh(draws[i, ], ...)
    if (all(draws[i+1, ] == draws[i, ])) {
      accept = accept + 1
    }
  }
  
#   print(accept / iterations)
  last = tail(draws, -burn)
  last[1:(dim(last)[1] / everyother) * everyother - 1, ]
  
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom", legend.key = element_blank(), legend.text.align = 0, legend.title.align = 0.5,legend.text=element_text(size=12)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    arrangeGrob(grobs = lapply(plots, function(x)
      x + theme(legend.position="none")), layout_matrix = t(matrix(c(1,2)))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight * 0.5, 0.5 * lheight))
}

gaussian.cop = function(n, m1, m2, cov) {
  S <- cov
  AB <- rmvnorm(mean=c(0,0),sig=S,n=n) #Our gaussian variables
  U <- pnorm(AB) #Now U is uniform - check using hist(U[,1]) or hist(U[,2])
  x <- m1(U[,1]) #x is gamma distributed
  y <- m2(U[,2]) #y is beta distributed
  cbind(x,y)
}

plot.many = function(data, N, ...){
  plot(data[1, ], type = "l", col = rgb(0,0,0,1/(N/50)), ...)
  for (i in 1:N) {
    lines(data[i, ], col = rgb(0,0,0,1/(N/75)))
  }
}

# Pretty plot.
pretty.plot = function(melted.data) {
  ggplot(data = data.melt) + theme_bw(base_size = 12, base_family = "Helvetica")
}

#
#
# MAXIMUM LAMBDA-NU ESTIMATION
#
#

ou.expected.ll = function(theta, x) {
  delta = 1/length(x)
  ll = 0
  for (i in 2:length(x[1,])) {
    ll = ll + log(dnorm(x[,i], x[,i-1] * exp(-delta), rep(1-exp(-2*delta) / 2, length(x[,1]))))
  }
  ll/length(x[,1])
}

q.unit = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(1)
  theta1 = theta[1]
  theta2 = theta[2]
  expectation = 0
  bridges = if (exact) {
    exact.bridges(1, r.nu, M, r.diffusion, theta.tilde, steps)
  } else {
    gen.nu.bridge(r.nu, r.diffusion, theta.tilde, M, steps)
  }
  us = sample(1:steps, M, replace = T)
  expectation = sum(rowSums(rowDiffs(bridges) * (theta1 - theta2 * bridges[,1:(steps-1)]))) - sum((1/2) * (theta1 - theta2 * bridges[cbind(1:M, us)])^2)
  -expectation / M
}

A = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  
  (-(lambda * x^2) / 4) - 
    (
      (1/2) - ((2*lambda*mu)/sigma^2)
    ) * log(x)
}

drift.cir = expression(theta1 - theta2 * x)

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

alpha = function(x, theta) {
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  (-(lambda * x) / 2) - (1 / (2*x)) *(1- (4*lambda*mu)/(sigma^2))
}


eta.inv = function(x, theta) {
  (theta[3]^2 / 4) * x^2
}

eta = function(x, theta) {
  2 * sqrt(x) / theta[3]
}

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

# N = 10000
# tran = matrix(0, nrow = N, ncol = 100)
# for (i in 1:10000) {
#   tran[i, ] = (trans.cir(sqrt(2),100,c(4,2,2)))
# }
# norm = matrix(0, nrow = N, ncol = 100)
# for (i in 1:10000) {
#   norm[i, ] = (eta(cir(2,100,c(4,2,2)), c(4,2,2)))
# }
# 
# plot(tran[1, ], ylim = c(0, 6), type = "l", col = rgb(0,0,0,0.05))
# for (i in 1:1000) {
#   lines(tran[i,], ylim = c(0, 6), type = "l", col = rgb(0,0,0,0.05))
# }
# 
# plot(norm[1, ], ylim = c(0, 6), type = "l", col = rgb(0,0,0,0.05))
# for (i in 1:1000) {
#   lines(norm[i,], ylim = c(0, 4), type = "l", col = rgb(0,0,0,0.05))
# }
# 
# ks.test(tran[,100], norm[,100])

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


# q = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
#   set.seed(1)
#   theta1 = theta[1]
#   theta2 = theta[2]
#   theta3 = theta[3]
#   if (2 * theta1 <= theta3^2 | theta1<0 | theta2<0 |theta3<0) {
#     return(NA)
#   }
#   
#   mu = theta1 / theta2
#   lambda = theta2
#   sigma = theta3
#     
#   expectation = 0
#   
#   r.nu.theta = function(n) {eta(r.nu(n), theta)}  
#   
#   bridges = if (exact) {
#     exact.bridges(1, r.nu.theta, M, r.diffusion, theta.tilde, steps)
#   } else {
#     gen.nu.bridge(r.nu.theta, r.diffusion, theta.tilde, M, steps, theta)
#   }
#   
#   us = sample(1:steps, M, replace = T)
#   
#   expectation = expectation + mean(A(bridges[,100], theta) - A(bridges[,1], theta))
#   
#   expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
#   
# #   print(c(theta, expectation))
#   
#   if (is.nan(expectation)) {
#     return(NA)
#   }
#   -expectation
# }
#  

gen.nu.bridge.special = function(r.nu, r.diffusion, theta.tilde, theta, M, steps, type) {
  if (type == "exact") {
    raw.bridges.0 = exact.bridges(1, function(n) {eta(r.nu(n), theta.tilde)},M, r.diffusion, theta.tilde, steps)
  } else if (type == "approx") {
    raw.bridges.0 = gen.nu.bridge(function(n) {eta(r.nu(n), theta.tilde)}, r.diffusion, theta.tilde, M, steps)
  }
#   if (type == "exact") {
#     raw.bridges.0 = exact.bridges(1, function(n) {matrix(0.5, nrow = n, ncol = 2)},M, r.diffusion, theta.tilde, steps)
#   } else if (type == "approx") {
#     raw.bridges.0 = gen.nu.bridge(function(n) {matrix(0.5, nrow = n, ncol = 2)}, r.diffusion, theta.tilde, M, steps)
#   }
  ends = raw.bridges.0[,c(1,steps)]
  new.ends = eta(eta.inv(ends, theta.tilde), theta)
#   new.ends = eta(matrix(1, nrow = M, ncol = 2), theta)
  difference = new.ends - ends
  raw.bridges.0 + vec.seq(difference[,1], difference[,2], length.out = steps)
}

q.1 = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(1)
  theta1 = theta[1]
  theta2 = theta[2]
  theta3 = theta[3]
  theta3 = 1
  if (2 * theta1 <= theta3^2 | theta1<0 | theta2<0 |theta3<0) {
    return(NA)
  }
  
  mu = theta1 / theta2
  lambda = theta2
  sigma = theta3
    
  expectation = 0
  
  r.nu.theta = function(n) {eta(r.nu(n), theta)}  
  
  bridges = if (exact) {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "exact")
  } else {
    gen.nu.bridge.special(r.nu, r.diffusion, theta.tilde, theta, M, steps, "approx")
  }
  
  us = sample(1:steps, M, replace = T)
  
  
  expectation = expectation + mean(g(eta.inv(bridges[,100], theta), theta) - g(eta.inv(bridges[,1], theta), theta))
  expectation = expectation - (1/2) * mean((bridges[,100] - bridges[,1])^2)
  expectation = expectation - mean(log(sigma * sqrt(eta.inv(bridges[,100], theta))))
#   expectation = sum(rowSums(rowDiffs(bridges) * (theta1 - theta2 * bridges[,1:(steps-1)]))) - sum((1/2) * (theta1 - theta2 * bridges[cbind(1:M, us)])^2)
#     expectation = expectation + mean(rowSums(rowDiffs(bridges) * alpha(bridges[,1:(steps-1)], theta)))
#   expectation = expectation - (1/2) * mean(alpha(bridges[cbind(1:M,us)], theta)^2)
  expectation = expectation - mean((1/2) * a.2.a.prime(bridges[cbind(1:M, us)], theta))
  
  if (is.nan(expectation)) {
    return(NA)
  } else if (expectation > 8){
    browser()
  }
  -expectation
}


q.2 = function(theta, r.diffusion, theta.tilde, r.nu, M, steps, exact = F) {
  set.seed(2)
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
  
  r.nu.theta = function(n) {eta(r.nu(n), theta)}  
  
  bridges = if (exact) {
    exact.bridges(1, r.nu.theta, M, r.diffusion, theta.tilde, steps)
  } else {
    gen.nu.bridge(r.nu.theta, r.diffusion, theta.tilde, M, steps, theta)
  }
  
  us = sample(1:steps, M, replace = T)
  
  expectation = sum(rowSums(rowDiffs(bridges) * (alpha(bridges[,1:(steps-1)], theta)))) - sum((1/2) * alpha(bridges[cbind(1:M, us)], theta)^2)
  
  if (is.nan(expectation)) {
    return(NA)
#   } else if (abs(expectation) > 100) {
#     browser()
  }
  -expectation / M
}






