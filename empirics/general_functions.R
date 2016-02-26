library(sde)
library(lattice)
library(ggplot2)
library(reshape2)
library(mvtnorm)
library(scales)
library(grid)
library(gridExtra)
library(NCmisc)
library(matrixStats)
library(Ryacas)
library(psych)

# APPROXIMATE BRIDGES

no.cross = function(x, y) {
  if (length(x) == length(y)) {
    diff = x - y
    return(all(diff< 0) | all(diff> 0))
  } else {
    stop("Different lengths of paths.")
  }
}

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
  }
  bridges
}

# EXACT BRIDGES

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

# DIFFUSIONS

ou = function(start, steps, theta, time = 1) {
  x = sde.sim(T= time, X0 = start, N = steps - 1,
          model = "OU", theta = theta)
  x
}

ou.fast = function(start, steps, theta, time = 1) {
  x = sde.sim(T = time, X0 = start, N = steps - 1,
              drift = eval(substitute(expression(theta1 -  theta2 * x), list(theta1 = theta[1], theta2 = theta[2]))),
              drift.x = eval(substitute(expression(-theta2), list(theta2 = theta[2]))),
              sigma = eval(substitute(expression(theta3), list(theta3 = theta[3]))),
              sigma.x = expression(0),
              method = "milstein")
  as.numeric(x)
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

multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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

grid_arrange_shared_legend = function(...) {
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

composite = function(f,g) function(...) f(g(...))

vec.seq = composite(t, Vectorize(seq.default, vectorize.args = c("from", "to")))

plot.many = function(data, N, ...){
  plot(data[1, ], type = "l", col = rgb(0,0,0,1/(N/50)), ...)
  for (i in 1:N) {
    lines(data[i, ], col = rgb(0,0,0,1/(N/75)))
  }
}

pretty.plot = function(melted.data) {
  ggplot(data = data.melt) + theme_bw(base_size = 12, base_family = "Helvetica")
}