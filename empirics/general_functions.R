# library(sde)
# library(lattice)
# library(mvtnorm)
# library(scales)
# library(grid)
# library(gridExtra)
# library(NCmisc)
# library(matrixStats)
# library(Ryacas)
# library(psych)
library(DEoptim)
# library(dfoptim)
# library(stabledist)
# library(VarianceGamma)
# library(ghyp)
# library(Matrix)
library(parallel)
# library(doParallel)
# library(somebm)

# APPROXIMATE BRIDGES

no.cross = function(x, y) {
  if (length(x) == length(y)) {
    diff = x - y
    return(all(diff< 0) | all(diff> 0))
  } else {
    stop("Different lengths of paths.")
  }
}

gen.nu.helper = function(r.diffusion, theta, start, end, steps, time = 1, flip = F) {
    y.1 = r.diffusion(start, steps, theta, time)
    y.2 = rev(r.diffusion(end, steps, theta, time))

    counter = 0
    # If they don't intersect, re-generate.
    while (no.cross(y.1, y.2)) {
      counter = counter + 1
      if (counter > 5000) {
        print("Too many tries.")
        stop("Too many tries.")
      }
      y.1 = r.diffusion(start, steps, theta, time)
      y.2 = rev(r.diffusion(end, steps, theta, time))
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

gen.nu.bridge = function(r.nu, r.diffusion, theta, N, steps, set.start = c(-Inf,-Inf), time = 1) {
  
  bridges = array(0, dim = c(N, steps))
  xy = r.nu(N)
  if (N == 1 & all(is.finite(set.start))) {
    xy = t((matrix(set.start)))
  }
#   stop(paste(as.character(xy), as.character(r.nu(1)), set.start))
  for (i in 1:N){
  # Generate two independent diffusions, one starting at start and the other at end.
   tryCatch({
      bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy[i, 1], xy[i, 2], steps, time)
    }, error = function(e) {
      tryCatch({
        xy1 = r.nu(1)
        bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps, time)
      }, error = function(e) {
        tryCatch({
          xy1 = r.nu(1)
          bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps, time)
        }, error = function(e) {
          tryCatch({
            xy1 = r.nu(1)
            bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps, time)
          }, error = function(e) {
            tryCatch({
              xy1 = r.nu(1)
              bridges[i, ] =gen.nu.helper(r.diffusion, theta, xy1[1], xy1[2], steps, time)
            }, error = function(e) {
               browser()
               print(c(i, xy[i,1], xy[i,2]))
            })
          })
        })
      })
    })
  }
  bridges
}

# EXACT BRIDGES

exact.bridges = function(M, r.nu, samples, r.diffusion, theta, steps, time = 1) {
  bridges = matrix(0, nrow = samples, ncol = 100)
  start.end = r.nu(1)
  bridges[1, ] = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, time, set.start = start.end)
  for (i in 2:samples) {
#       stop("IN THE FOR LOOP")
#     tryCatch({
      start.end = r.nu(1)
      proposal.1 = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, time, set.start = start.end)
      proposal.2 = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, time, set.start = start.end)
#     }, error = function(e) {
#       bridges[i, ] = bridges[i-1,]
#     })
    r = rho(proposal.2, M, r.diffusion, theta, steps, time) / rho(proposal.1, M, r.diffusion, theta, steps, time)
    alpha = min(1, r)
#     print(c(i, alpha))
#     print(paste(i,alpha))
    if (runif(1) < alpha) {
      bridges[i,] = proposal.2
    } else {
      bridges[i,] = proposal.1
    }
  }
  bridges
}

rho = function(x, M, r.diffusion, theta, steps, time = 1) {
  Ts = array(1, dim = M)
  for (i in 1:M) {
    init.1 = r.diffusion(x[steps], steps, theta, time)[steps]
    hitting = as.vector(r.diffusion(init.1, steps, theta, time))
    while (no.cross(x, hitting)) {
      if (Ts[i] > 5000) {
        Ts[i] = 5000
        break
      }
      Ts[i] = Ts[i] + 1
      init.1 = r.diffusion(x[steps], steps, theta, time)[steps]
      hitting = as.vector(r.diffusion(init.1, steps, theta, time))
    }
  }
  mean(Ts)
}

# DIFFUSIONS

sde.sim.t = function ( t0 =0 , T =1 , X0 =1 , N =100 , delta ,
                        drift , sigma , sigma.x, df){
  d1 <- function(t , x) eval( drift )
  s1 <- function(t , x) eval( sigma )
  sx <- function(t , x) eval( sigma.x)
  if ( t0 <0 | T <0)
    stop( " please use positive times !" )
  if ( missing( delta )){
    t <- seq(t0 ,T , length = N +1)
  } else {
    t <- c(t0 , t0 + cumsum ( rep ( delta , N )))
    T <- t[N +1]
    warning(" T set to = " ,T , "\ n" )
  }
  Z <- rt( N, df) / sqrt(df/(df-2))
  X <- numeric( N +1)
  Dt <- (T - t0 )/ N
  sDt <- sqrt( Dt )
  X [1] <- X0
  for (i in 2:( N +1)){
    X[ i] <- X[i -1] + d1( t[i -1] , X[i -1]) * Dt +
      s1(t[i -1] , X[i -1]) * sDt *Z[i -1] +
      0.5 * s1( t[i -1] , X [i -1]) * sx(t[i -1] , X[i -1]) *
      ( Dt *Z[i -1]^2 - Dt )
  }
  X <- ts(X , start = t0 , deltat = Dt )
  invisible(X)
}

ou = function(start, steps, theta, time = 1) {
  x = sde.sim(T= time, X0 = start, N = steps - 1,
              model = "OU", theta = theta)
  x
}

ou.fast = function(start, steps, theta, time = 1) {
  Y <- numeric(steps)
  Y[1] = start
  Z = rnorm(steps - 1)
  Dt <- time / (steps - 1)
  for (i in 1:(steps - 1)) {
    Y[i+1] = Y[i] + (theta[1] - theta[2] * Y[i]) * Dt + theta[3] * sqrt(Dt) * Z[i]
  }
  as.numeric(Y)
}

wiener = function(start, steps, theta, time = 1) {
  Y <- numeric(steps)
  Y[1] = start
  Z = rnorm(steps - 1)
  Dt <- time / (steps - 1)
  for (i in 1:(steps - 1)) {
    Y[i+1] = Y[i] + theta[3] * sqrt(Dt) * Z[i]
  }
  as.numeric(Y)
}

ou.stable = function(start, steps, theta, time = 1) {
  Y <- numeric(steps)
  Y[1] = start
  Dt <- time / (steps - 1)
  Z = rstable(steps - 1, alpha = 1.8, beta = 0, gamma = theta[3] * Dt^(1/1.8), delta = 0)
  for (i in 1:(steps - 1)) {
    Y[i+1] = Y[i] + (theta[1] - theta[2] * Y[i]) * Dt + Z[i]
  }
  as.numeric(Y)
}

ou.nig = function(start, steps, theta, time = 1) {
  Y <- numeric(steps)
  Y[1] = start
  Dt <- time / (steps - 1)
  Z = rghyp(steps - 1, NIG.ad(alpha = 1, delta = 1, beta = 0, mu = 0))
  for (i in 1:(steps - 1)) {
    Y[i+1] = Y[i] + (theta[1] - theta[2] * Y[i]) * Dt + sqrt(Dt) *  Z[i]
  }
  as.numeric(Y)
}

ou.vg = function(start, steps, theta, time = 1) {
  Y <- numeric(steps)
  Y[1] = start
  Dt <- time / (steps - 1)
  Z = rnorm(steps - 1)
  G = rgamma(steps - 1, shape = Dt / 0.1, scale = 1/0.1)
  for (i in 1:(steps - 1)) {
    Y[i+1] = Y[i] + (theta[1] - theta[2] * Y[i]) * Dt + sqrt(G[i]) *  Z[i]
  }
  as.numeric(Y)
}

hyp = function(start, steps, theta, time = 1) {
  Y <- numeric(steps)
  Y[1] = start
  Z = rnorm(steps - 1)
  Dt <- time / (steps - 1)
  for (i in 1:(steps - 1)) {
    Y[i+1] = Y[i] + 
      (- (theta[4]^2 / 2) * 
         (theta[2] * (Y[i] - theta[1])) / (sqrt(theta[3]^2 + (Y[i] - theta[1])^2))
       ) * Dt + 
      theta[4] * sqrt(Dt) *  Z[i]
  }
  as.numeric(Y)
}


make.joint = function(N, theta, r.diffusion, init, steps = 100, all = T) {
  sims2 = array(0, dim = c(N, steps))
  sims.start = init(N)
  for (i in 1:N) {
    sims2[i, ] = r.diffusion(sims.start[i], steps, theta)
  }
  if (all) {
    function(n) {(as.matrix(sims2[sample(1:N, n, replace = T),c(1,steps)]))}
  } else {
    function(n) {(as.matrix(sims2[sample(1:100, n, replace = T),c(1,steps)]))}
  }
}

make.joint.many = function(N, theta, r.diffusion, init, horizons, steps = 100, time = 1, all = T) {
  sims2 = array(0, dim = c(N, steps))
  sims.start = init(N)
  for (i in 1:N) {
    sims2[i, ] = r.diffusion(sims.start[i], steps, theta, time)
  }
  if (all) {
    function(n) {(as.matrix(sims2[sample(1:N, n, replace = T),c(1/steps,horizons,1) * steps]))}
  } else {
    function(n) {(as.matrix(sims2[sample(1:100, n, replace = T),c(1/steps,horizons,1) * steps]))}
  }
}

cir = function(start, steps, theta, time = 1) {
  x = sde.sim(T = time, X0 = start, N = steps - 1,
              model = "CIR", theta = theta )
  as.numeric(x)
}

cir.fast = function(start, steps, theta, time = 1) {
  x = sde.sim(T = time, X0 = start, N = steps - 1,
              drift = eval(substitute(expression(theta1 -  theta2 * x), list(theta1 = theta[1], theta2 = theta[2]))),
              drift.x = eval(substitute(expression(-theta2), list(theta2 = theta[2]))),
              sigma = eval(substitute(expression(theta3 * sqrt(x)), list(theta3 = theta[3]))),
              sigma.x = eval(substitute(expression(theta3 / (2 * sqrt(x))), list(theta3 = theta[3]))),
              method = "milstein")
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

# MCMC

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

# TOOLS

compare.diffusions = function(r.dif1, r.dif2, N, start1, start2, steps, theta1, theta2) {
  tran = matrix(0, nrow = N, ncol = steps)
  timeit(for (i in 1:N) {
    tran[i, ] = r.dif1(start1, steps, theta1)
  })[1,]
  norm = matrix(0, nrow = N, ncol = steps)
  timeit(
  for (i in 1:N) {
    norm[i, ] = r.dif2(start2, steps, theta2)
  })[1,]
  print(ks.test(tran[,steps], norm[,steps]))
}

composite = function(f,g) function(...) f(g(...))

vec.seq = composite(t, Vectorize(seq.default, vectorize.args = c("from", "to")))

trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

gen.many = function(M, r.diffusion, start, steps, theta, time = 1) {
  paths = matrix(0, nrow = M, ncol = steps)  
  for (i in 1:M) {
    paths[i, ] = r.diffusion(start, steps, theta, time)
  }
  paths
}