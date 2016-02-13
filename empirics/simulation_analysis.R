load("/Users/marshall/Documents/senior/thesis/figures/figure1data.Rdata")

approx.ou = figure1$approx.ou[,70]
exact.ou = figure1$exact.ou[,20]
real.ou = figure1$real.ou[,20]

ks.test(real.ou, approx.ou)
ks.test(real.ou, exact.ou)

approx.cir= figure1$approx.cir[,50]
exact.cir= figure1$exact.cir[,50]
real.cir= figure1$real.cir[,50]

ks.test(real.cir, approx.cir)
ks.test(real.cir, exact.cir)

load("/Users/marshall/Documents/senior/thesis/figures/table1data.Rdata")

i = 30
ks.test(table1data$binom.real[,i], table1data$binom.approx[,i])
ks.test(table1data$binom.real[,50], table1data$binom.exact[,50])

ks.test(table1data$m1.real[,i], table1data$m1.approx[,i])
ks.test(table1data$m1.real[,50], table1data$m1.exact[,50])

ks.test(table1data$exp.real[,i], table1data$exp.approx[,i])
ks.test(table1data$exp.real[,50], table1data$exp.exact[,50])
