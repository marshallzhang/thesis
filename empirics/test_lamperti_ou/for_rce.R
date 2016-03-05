source("../general_functions.R")
source("functions.R")

my_id = as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1

theta = c(0, 1, 1)

steps = 100 
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
  nu = make.joint(10000, theta, r.diffusion, init, steps = 100, all = F)
} else {
  nu = make.joint(10000, theta, r.diffusion, init, steps = 100)
}

print("Starting DEoptim.")
oldpar = c(4,2,2)
set.seed(3)
for (i in 1:10) {
  obj = function(theta) {
    source("to_source_parallel.R")
    q.beskos.unknown.single(theta, trans.ou.fast, oldpar, nu, round(i^1.2)*100, steps, exact = F)
  }
  oldpar = DEoptim(fn = obj,
                   lower = c(-1.5, 0, 0),
                   upper = c(1.5, 2.5, 2.5),
#                    control = DEoptim.control(trace = TRUE, strategy = 6))
                   control = DEoptim.control(trace = TRUE, strategy = 6, parallelType = 1,
                                             packages = c(),
                                             parVar = list("q.beskos.unknown.single", "eta", "nu", "i", "A",
                                                           "steps", "oldpar",
                                                           "vec.seq", "a.2.a.prime")
                                                        ))[["optim"]]$bestmem
  print(oldpar)
}