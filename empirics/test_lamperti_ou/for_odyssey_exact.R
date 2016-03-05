if (Sys.getenv("SLURM_JOB_ID") != "") {
  
  steps = 100
  
  source("../general_functions.R")
  source("functions.R")
  
  my_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  
  theta = c(0, 1, 1)

  inits = c(function(n) rnorm(n, 0, theta[3] / sqrt(2 * theta[2])),
            function(n) rnorm(n, 0, theta[3] / sqrt(2 * theta[2])),
            function(n) 2*rbinom(n, 1, 0.5) - 1)
  
  r.diffusions = c(ou.fast,
                   ou.vg,
                   ou.nig
  )
  
  if (my_id > 9) {
    counter = 1 + (my_id - 1) %% 9
  } else {
    counter = my_id
  }
  
  init = inits[[ceiling(counter / 3)]]
  r.diffusion = r.diffusions[[1 + ((counter - 1) %% 3)]]
  
  print("Making nu.")
  if (ceiling(counter / 3) == 2) {
    nu = make.joint(1000000, theta, r.diffusion, init, steps = 100, all = F)
  } else {
    nu = make.joint(1000000, theta, r.diffusion, init, steps = 100)
  }
  
  print("Starting DEoptim.")
  oldpar = c(4,2,2)
  set.seed(my_id+2345234)
  for (i in 1:10) {
    oldpar = DEoptim(fn = function(theta) q.beskos.unknown(theta, trans.ou.fast, oldpar, nu, round(i^1.2)*100, steps, exact = T),
                     lower = c(-1.5, 0, 0),
                     upper = c(1.5, 2.5, 2.5),
                     control = DEoptim.control(trace = TRUE, strategy = 1 + (my_id - 1) %% 6))[["optim"]]$bestmem
    print(oldpar)
  }
}