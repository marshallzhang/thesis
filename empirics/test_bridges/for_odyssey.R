if (Sys.getenv("SLURM_JOB_ID") != "") {
  
  N = 10000
  steps = 100
  
  source("../general_functions.R")
  
  my_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  
  theta = c(0, 1, 1)

  inits = c(function(n) rnorm(n, 0, theta[3] / sqrt(2 * theta[2])),
            function(n) rexp(n),
            function(n) 2*rbinom(n, 1, 0.5))
  
  r.diffusions = c(ou.fast,
                  ou.stable,
                  wiener)
  
  init = inits[[ceiling(my_id / 3)]]
  r.diffusion = r.diffusions[[1 + ((my_id - 1) %% 3)]]
  
  print("Doing o.sims.")
  o.sims = array(0, dim = c(N, steps))
  timeit(for (i in 1:N) {
    o.sims[i, ] = r.diffusion(init(1), steps, theta)
  })
  
  nu.o = make.joint(50000, theta, r.diffusion, init)
  
  print("Doing nu.o.sims.")
  print(timeit(nu.o.sims <- gen.nu.bridge(nu.o, r.diffusion, theta, N, 100), suppressResult = T))
  
  print("Doing nu.o.sims.exact.")
  print(timeit(nu.o.sims.exact <- exact.bridges(1, nu.o, N, r.diffusion, theta, 100), suppressResult = T))
  
  print(ks.test(nu.o.sims[,50], o.sims[,50]))
  print(ks.test(nu.o.sims.exact[,50], o.sims[,50]))
  
  data = list(o.sims, nu.o.sims, nu.o.sims.exact)
  save(data, file = paste("dump/", my_id, ".Rdata", sep = ""))
  
}