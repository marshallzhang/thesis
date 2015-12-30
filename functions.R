library(sde)
library(lattice)
library(ggplot2)
library(reshape2)
library(mvtnorm)
library(scales)
library(grid)
library(gridExtra)

#
# GENERATE NU BRIDGE
#

# Helper to generate a bridge.
gen.nu.helper = function(r.diffusion, theta, start, end, steps) {
    y.1 = r.diffusion(start, steps, theta)
    y.2 = rev(r.diffusion(end, steps, theta))
    difference = y.1 - y.2
    
    counter = 0
    # If they don't intersect, re-generate.
    while (all(difference < 0) | all(difference > 0)) {
      counter = counter + 1
      if (counter > 500) {
        stop("Too many tries.")
      }
      y.1 = r.diffusion(start, steps, theta)
      y.2 = rev(r.diffusion(end, steps, theta))
      difference = y.1 - y.2
    }
      
    # Otherwise, find where they intersect and create the bridge.  
    if (difference[1] > 0) {
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
    print(i)
  }
  bridges
}

# MISC

# Generate OU process.
ou = function(start, steps, theta) {
  x = sde.sim(X0 = start, N = steps - 1,
          model = "OU", theta = theta )
  x
}

# Pretty plot.
pretty.plot = function(melted.data) {
  ggplot(data = data.melt) + theme_bw(base_size = 12, base_family = "Helvetica")
}

# CIR

cir = function(start, steps, theta) {
  x = sde.sim(X0 = start, N = steps - 1,
              model = "CIR", theta = theta )
  as.numeric(x)
}

joint = function(current) {
  x = current[1]
  y = current[2]
  c = 2 / (1-exp(-1))
  q = 1
  u = c * x *exp(-1)
  v = c * y
  trans = c * exp(-u-v) * (v/u)^(q/2) * besselI(2 * sqrt(u * v), q)
  log(dgamma(x, shape = 2, rate = 2) * trans)
}


# Metropolis-Hastings with normal proposals.
mcmc.mh = function(current, d.posterior, r.proposal) {
  
    # Propose and find density according to random walk.
    proposal = r.proposal(length(current), current)
    alpha = min(1, exp(d.posterior(proposal) - d.posterior(current) + sum(dgamma(current, shape = 2, rate = 2, log =T )) - sum(dgamma(proposal, shape = 2, rate = 2, log = T))))
    
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