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

exact.bridges = function(M, r.nu, samples, r.diffusion, theta,steps, time = 1) {
  bridges = matrix(0, nrow = samples, ncol = 100)
  start.end = r.nu(1)
  bridges[1, ] = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, time, set.start = start.end)
  for (i in 2:samples) {
      start.end = r.nu(1)
      proposal.1 = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, time, set.start = start.end)
      proposal.2 = gen.nu.bridge(r.nu, r.diffusion, theta, 1, steps, time, set.start = start.end)
    r = rho(proposal.2, M, r.diffusion, theta, steps, time) / rho(proposal.1, M, r.diffusion, theta, steps, time)
    alpha = min(1, r)
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