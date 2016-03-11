

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

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/est_kl.pdf", width= 10, height = 5, #' see how it looks at this size
    useDingbats=F)
pretty.plot(melted) + geom_line(aes(x=V501, y = value, group = variable), alpha = 0.05) + 
              xlab("Quarters Elapsed") + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 12)) + 
              ylab("Inflation (%)") + 
              scale_y_continuous(breaks = seq(0,3,0.5)) + 
              theme(legend.position = "none") + 
            ggtitle("Imputed Path of Inflation Based on 2015Q4 SPF")

dev.off()