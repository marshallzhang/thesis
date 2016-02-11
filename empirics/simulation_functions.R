library(sde)
library(lattice)
library(ggplot2)
library(reshape2)
library(mvtnorm)
library(scales)
library(grid)
library(gridExtra)
library(NCmisc)

#
# GENERATE NU BRIDGE
#

gen.double.nu.bridge = function(r.diffusion, theta, start, end, steps) {
    y.1 = r.diffusion(start, steps, theta)
    y.2 = r.diffusion(end, steps, theta)
  
    y.1.1 = y.1
    y.2.1 = rev(y.2)
    difference.1 = y.1.1 - y.2.1
    
    y.1.2 = rev(y.1)
    y.2.2 = y.2
    difference.2 = y.2.2 - y.1.2
    
    counter = 0
    
    # If they don't intersect, re-generate.
    while (all(difference.1 < 0) | all(difference.1 > 0) | all(difference.2 < 0) | all(difference.2 > 0)) {
      counter = counter + 1
      if (counter > 500) {
        stop("Too many tries.")
      }
      y.1 = r.diffusion(start, steps, theta)
      y.2 = r.diffusion(end, steps, theta)
    
      y.1.1 = y.1
      y.2.1 = rev(y.2)
      difference.1 = y.1.1 - y.2.1
      
      y.1.2 = rev(y.1)
      y.2.2 = y.2
      difference.2 = y.2.2 - y.1.2
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
      b2 = c(y.2.2[1:(iota - 1)], y.1.2[iota:length(y.1.2)])
    } else {
      iota = which(difference.2 > 0)[1]
      b2 = c(y.2.2[1:(iota - 1)], y.1.2[iota:length(y.1.2)])
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
      if (counter > 500) {
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
gen.nu.bridge = function(r.nu, r.diffusion, theta, N, steps) {
  
  bridges = array(0, dim = c(N, steps))
  xy = r.nu(N)
  
  for (i in 1:N){
  # Generate two independent diffusions, one starting at start and the other at end.
    bridges[i, ] = 
     tryCatch({
      gen.nu.helper(r.diffusion, theta, xy[i, 1], xy[i, 2], steps)
    }, error = function(e) {
     bridges[i - 1, ]
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
ou = function(start, steps, theta) {
  x = sde.sim(X0 = start, N = steps - 1,
          model = "OU", theta = theta, method = "euler")
  x
}


#
#
#
# CIR
#
#
#


cir = function(start, steps, theta) {
  x = sde.sim(X0 = start, N = steps - 1,
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
d.rho = function(x, M, r.hitting.init, r.nu, r.diffusion, theta, N, steps) {
  Ts = array(1, dim = M)
  for (i in 1:M) {
    hitting1 = as.vector(r.diffusion(rnorm(1, x[1,100]*exp(-1), sqrt((1-exp(-2))/2)), steps, theta))
    hitting2 = as.vector(r.diffusion(rnorm(1, x[2,100]*exp(-1), sqrt((1-exp(-2))/2)), steps, theta))
    while (no.cross(x[1,], hitting1) | no.cross(x[2,], hitting2)) {
      Ts[i] = Ts[i] + 1
      hitting1 = as.vector(r.diffusion(rnorm(1, x[1,100]*exp(-1), sqrt((1-exp(-2))/2)), steps, theta))
      hitting2 = as.vector(r.diffusion(rnorm(1, x[2,100]*exp(-1), sqrt((1-exp(-2))/2)), steps, theta))
    }
  }
  mean(Ts)
}

rho = function(x, M, r.hitting.init, r.nu, r.diffusion, theta, N, steps) {
  Ts = array(1, dim = M)
  for (i in 1:M) {
    hitting = as.vector(r.diffusion(rnorm(1, x[100]*exp(-1), sqrt((1-exp(-2))/2)), steps, theta))
    while (no.cross(x, hitting)) {
      Ts[i] = Ts[i] + 1
      hitting = as.vector(r.diffusion(rnorm(1, x[100]*exp(-1), sqrt((1-exp(-2))/2)), steps, theta))
    }
  }
  mean(Ts)
}

exact.nu.double.bridge = function(M, r.hitting.init, r.nu, samples, r.diffusion, theta, steps) {
  bridges = vector("list", samples)
  xy = r.nu(1)
  bridges[[1]] = gen.double.nu.bridge(r.diffusion, theta, xy[1], xy[2], steps)
  for (i in 2:samples) {
    xy = r.nu(1)
    proposal = gen.double.nu.bridge(r.diffusion, theta, xy[1], xy[2], steps)
    rho.im1 = d.rho(bridges[[i-1]], M, r.hitting.init, r.nu, r.diffusion, theta, N, steps)
    rho.i = d.rho(proposal, M, r.hitting.init, r.nu, r.diffusion, theta, N, steps)
    r = rho.i / rho.im1
    alpha = min(1, r)
    if (runif(1) < alpha) {
      bridges[[i]] = proposal
    } else {
      xy = r.nu(1)
      second.proposal = gen.double.nu.bridge(r.diffusion, theta, xy[1], xy[2], steps)
      rho.ip1 = d.rho(second.proposal, M, r.hitting.init, r.nu, r.diffusion, theta, N, steps)
      alpha.2 = (rho.ip1 / rho.im1) * ((1 - min(1,(rho.i / rho.ip1))) / (1 - r))
      print(paste(i, round(r, digits = 2), round(min(1,alpha.2), digits = 2)))
      if (runif(1) < min(1,alpha.2)) {
        bridges[[i]] = second.proposal
      } else {
        bridges[[i]] = bridges[[i-1]]
      }
    }
  }
  bridges
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
    if (trace != 0) if (i %% trace == 0) print(paste("Iteration", i, "of", iterations))
    draws[i + 1, ] = mcmc.mh(draws[i, ], ...)
    if (all(draws[i+1, ] == draws[i, ])) {
      accept = accept + 1
    }
  }
  
  print(accept / iterations)
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





















joint.dirac = function(N, a, b, lambda, kappa) {
  
  ret = array(0, dim = c(N, 2))
  
  for (i in 1:N) {
    if (runif(1) < lambda) {
      ret[i,1] = a
    } else {
      ret[i,1] = b
    }
    if (ret[i,1] == a) {
      if (runif(1) < kappa) {
        ret[i,2] = a
      } else {
        ret[i,2] = b
      }
    } else {
      if (runif(1) < kappa) {
        ret[i,2] = b
      } else {
        ret[i,2] = a
      }
    }
  }
  ret
}