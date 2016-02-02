library(Matrix)

# Simulate multivariate OU.
mvOU = function(mu, lambda, cov, start, steps) {
  n = length(mu)
  res = matrix(0, nrow = n, ncol = steps)
  res[,1] = start
  for (t in 2:steps) {
    el = exp(diag(lambda,3))
    el[el == 1] = 0
    mean = mu + el %*% (res[, t-1] - mu)
    res[,t] = rmvnorm(1, mean = mean, sigma = as.matrix(cov - el %*% cov %*% el))
  }
  res
}

# Helper to generate a bridge.
gen.mv.nu.helper = function(theta, start, end, steps) {
    co = gen.cov(theta[5:10])
    y.1 = mvOU(theta[1:3], theta[4], co, start, steps)
    y.2 = t(apply(mvOU(theta[1:3], theta[4], co, end, steps), 1, rev))
    difference = y.1 - y.2
    
    counter = 0
    # If they don't intersect, re-generate.
    while (any(apply(difference, 1, function(x) all(x < 0) | all(x > 0)))) {
      counter = counter + 1
      if (counter > 1500) {
        stop("Too many tries.")
      }
      y.1 = mvOU(theta[1:3], theta[4], co, start, steps)
      y.2 = t(apply(mvOU(theta[1:3], theta[4], co, end, steps), 1, rev))
      difference = y.1 - y.2
    }
      
    # Otherwise, find where they intersect and create the bridge.  
    bridge = matrix(0, nrow = length(start), ncol = steps)
    for (i in 1:length(start)) {
      if (difference[i,1] > 0) {
        iota = which(difference[i, ] < 0)[1] 
        bridge[i, ] = c(y.1[i, 1:(iota - 1)], y.2[i, iota:ncol(y.2)])
      } else {
        iota = which(difference[i, ] > 0)[1]
        bridge[i, ] = c(y.1[i, 1:(iota - 1)], y.2[i, iota:ncol(y.2)])
      }  
    }
  bridge  
}

# Generating an mv nu bridge.
gen.mv.nu.bridge = function(r.nu, theta, N, steps) {
  
  bridges = vector("list", N)
  xy = r.nu(N)
  
  for (i in 1:N){
  # Generate two independent diffusions, one starting at start and the other at end.
    bridges[[i]] = 
     tryCatch({
      gen.mv.nu.helper(theta, xy[i, 1:3], xy[i, 4:6], steps)
    }, error = function(e) {
      print("Too many tries.")
     bridges[[i - 1]]
    })
    print(i)
  }
  bridges
}

plot.many = function(series, alpha) {
  plot(series[1,], type = "l", col = rgb(0,0,0,alpha))
  for (i in 2:100)  {
    lines(infs[,i], col = rgb(0,0,0,0.1))
  }
}