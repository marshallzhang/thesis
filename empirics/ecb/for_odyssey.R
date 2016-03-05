if (Sys.getenv("SLURM_JOB_ID") != "") {
  
  source("../general_functions.R")
  source("functions.R")
  
  my_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  years = 1999:2015
  quarters = 1:4
  
  times = vector("list", length(years) * length(quarters))
  
  for (i in 1:length(years)) {
    for (j in 1:length(quarters)) {
      times[[(i - 1) * 4 + j]] = c(years[i], quarters[j])
    }
  }
  
  year = times[[my_id]][1]
  quarter = times[[my_id]][2]
  
  forecast = fc(year, quarter)
  
  print("Starting optim.")
  theta.star = DEoptim(fn = function(theta) ou.ll(theta,forecast),
                   lower = c(-9, 0, 0),
                   upper = c(9, 3, 3),
                   control = DEoptim.control(trace = TRUE, strategy = 6, itermax = 5000))[["optim"]]$bestmem
  print(paste(year, quarter, theta.star[1], theta.star[2], theta.star[3], sep = ","))
}