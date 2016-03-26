

fc = forecast.mats[["2015"]][["4"]]
      
N = 500

fixed = fc[sample(1:nrow(fc), N, replace = T), ]

fst = matrix(NA, nrow = N, ncol = 50)
for (i in 1:N) {
  fst[i, ] = gen.nu.bridge(function(n) rep(0,n), ou.fast, thetas[["2015"]][["4"]], 1, 50, set.start = fixed[i,1:2], time = 1)
}

snd = matrix(NA, nrow = N, ncol = 200)
for (i in 1:N) {
  snd[i, ] = gen.nu.bridge(function(n) rep(0,n), ou.fast, thetas[["2015"]][["4"]], 1, 200, set.start = fixed[i,2:3], time = 4)
}

third = matrix(NA, nrow = N, ncol = 200)
for (i in 1:N) {
  third[i, ] = gen.nu.bridge(function(n) rep(0,n), ou.fast, thetas[["2015"]][["4"]], 1, 200, set.start = fixed[i,3:4], time = 4)
}
fourth = matrix(NA, nrow = N, ncol = 600)
for (i in 1:N) {
  fourth[i, ] = gen.nu.bridge(function(n) rep(0,n), ou.fast, thetas[["2015"]][["4"]], 1, 600, set.start = fixed[i,4:5], time = 12)
}

bridges = cbind(fst, snd, third, fourth)
bridges = rbind(bridges, 21 * 1:1050 / 1050)
bridges = t(bridges)
bridges = as.data.frame(bridges)
melted = melt(bridges, id = "V501")

o = pretty.plot(melted) + geom_line(aes(x=V501, y = value, group = variable), alpha = 0.05) + 
              xlab("Quarters Elapsed") + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 12)) + 
              ylab("Inflation (%)") + 
              scale_y_continuous(breaks = seq(0,3,0.5)) + 
              theme(legend.position = "none")

fc = forecast.mats[["2015"]][["4"]]
      
N = 500

fixed = fc[sample(1:nrow(fc), N, replace = T), ]

fst = matrix(NA, nrow = N, ncol = 50)
for (i in 1:N) {
  fst[i, ] = gen.nu.bridge(function(n) rep(0,n), hyp.fast, c(2.00,7.14,0.03, 0.125), 1, 50, set.start = fixed[i,1:2], time = 1)
}

snd = matrix(NA, nrow = N, ncol = 200)
for (i in 1:N) {
  print(i)
  tryCatch({
    snd[i, ] = gen.nu.bridge(function(n) rep(0,n), hyp.fast, c(2.00,7.14,0.03, 0.125), 1, 200, set.start = fixed[i,2:3], time = 4)
  }, error = function(e) {
    snd[i, ] = snd[i-1,]
  }
  })
}

third = matrix(NA, nrow = N, ncol = 200)
for (i in 1:N) {
  third[i, ] = gen.nu.bridge(function(n) rep(0,n), hyp.fast, c(2.00,7.14,0.03, 0.125), 1, 200, set.start = fixed[i,3:4], time = 4)
}
fourth = matrix(NA, nrow = N, ncol = 600)
for (i in 1:N) {
  fourth[i, ] = gen.nu.bridge(function(n) rep(0,n), hyp.fast, c(2.00,7.14,0.03, 0.125), 1, 600, set.start = fixed[i,4:5], time = 12)
}

bridges = cbind(fst, snd, third, fourth)
bridges = rbind(bridges, 21 * 1:1050 / 1050)
bridges = t(bridges)
bridges = as.data.frame(bridges)
melted.hyp = melt(bridges, id = "V501")

h = pretty.plot(melted.hyp) + geom_line(aes(x=V501, y = value, group = variable), alpha = 0.05) + 
              xlab("Quarters Elapsed") + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 12)) + 
              ylab("Inflation (%)") + 
              scale_y_continuous(breaks = seq(0,3,0.5)) + 
              theme(legend.position = "none")

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/cont_path.pdf", width= 10, height = 7, #' see how it looks at this size
    useDingbats=F)
multiplot(o, h, col = 1)
dev.off()