source("../../general_functions.R")
source("../mcem_functions.R")

process = as.integer(as.character(commandArgs(trailingOnly = TRUE)))

theta = c(0, 1, 1)

init = function(n) rnorm(n, 0, theta[3] / sqrt(2 * theta[2]))

r.diffusion = ou.fast

horizons = c(.5)
tot.time = 1
print("Making nu.")
nu = make.joint.many(10000, theta, r.diffusion, init, horizons = horizons, steps = 100, time = tot.time)

print("Starting DEoptim.")
oldpar = c(4,2,2)
steps = 100
# set.seed(process + 10)
for (i in 1:10) {
  obj = function(theta) {
    source("/Users/marshall/Documents/senior/thesis/empirics/ecb/for_rce/to_source_parallel.R")
    q.beskos.unknown(theta, trans.ou.fast, oldpar, nu, round(i^1.2)*100, steps, c(0, horizons, 1), tot.time, exact = F)
  }
  oldpar = DEoptim(fn = obj,
                   lower = c(-1.5, 0, 0),
                   upper = c(1.5, 2.5, 2.5),
                   control = DEoptim.control(trace = TRUE, strategy = 6, parallelType = 1,
                                             packages = c(),
                                             parVar = list("q.beskos.unknown", "eta", "nu", "i", "A",
                                                           "steps", "oldpar", "horizons", "tot.time",
                                                           "vec.seq", "a.2.a.prime")
                                                        ))[["optim"]]$bestmem
  print(oldpar)
}

horizons = c()
time = 5
nu = make.joint.many(10000, theta, ou, init, horizons, steps = 10, time = time)
hi = nu(10000)
colnames(hi) = c(0, time)

DEoptim(fn = function(theta) ou.ll(theta, fc(2011,3)),
                 lower = c(-1.5, 0, 0),
                 upper = c(1.5, 2, 2),
                 control = DEoptim.control(trace = TRUE, strategy = 6,itermax = 10000))$optim$bestmem

