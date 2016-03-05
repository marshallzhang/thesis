load("/Users/marshall/Documents/senior/thesis/empirics/test_bridges/1.Rdata")

# Plot.
edge = 3
skip = 1
data.q = sort(data[[1]][,50])
my.data.q = sort(data[[2]][,50])
exact.data.q = sort(data[[3]][,50])
data = data.frame(th = data.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
o.plot.approx = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Exact bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge,skip)) + 
              scale_x_continuous(breaks = seq(-edge,edge,skip)) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge), xlim=c(-edge,edge))+
              ggtitle("N(0,1/2)") + 
              theme(legend.position = "none",
                    text = element_text(size=16))

data = data.frame(th = data.q, my.dat = exact.data.q)
data.melt = melt(data, id = "th")
o.plot.exact = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Approximate bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge,skip)) + 
              scale_x_continuous(breaks = seq(-edge,edge,skip)) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge), xlim=c(-edge,edge))+
              ggtitle("N(0,1/2)") + 
              theme(legend.position = "none",
                    text = element_text(size=16))

load("/Users/marshall/Documents/senior/thesis/empirics/test_bridges/4.Rdata")

# Plot.
edge = 5
skip = 1
data.q = sort(data[[1]][,50])
my.data.q = sort(data[[2]][,50])
exact.data.q = sort(data[[3]][,50])
data = data.frame(th = data.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
exp.plot.approx = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Approximate bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge-1,skip)+2) + 
              scale_x_continuous(breaks = seq(-edge,edge-1,skip) +2) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge-1)+2, xlim=c(-edge,edge-1)+2)+
              ggtitle("Expo") + 
              theme(legend.position = "none",
                    text = element_text(size=16))

data = data.frame(th = data.q, my.dat = exact.data.q)
data.melt = melt(data, id = "th")
exp.plot.exact = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Exact bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge-1,skip)+2) + 
              scale_x_continuous(breaks = seq(-edge,edge-1,skip)+2) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge-1)+2, xlim=c(-edge,edge-1)+2)+
              ggtitle("Expo") + 
              theme(legend.position = "none",
                    text = element_text(size=16))

load("/Users/marshall/Documents/senior/thesis/empirics/test_bridges/7.Rdata")

# Plot.
edge = 3
skip = 1
data.q = sort(data[[1]][,50])
my.data.q = sort(data[[2]][,50])
exact.data.q = sort(data[[3]][,50])
data = data.frame(th = data.q, my.dat = my.data.q)
data.melt = melt(data, id = "th")
binom.plot.approx = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Approximate bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge,skip)) + 
              scale_x_continuous(breaks = seq(-edge,edge,skip)) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge), xlim=c(-edge,edge))+
              ggtitle("2Bern - 1") + 
              theme(legend.position = "none",
                    text = element_text(size=16))

data = data.frame(th = data.q, my.dat = exact.data.q)
data.melt = melt(data, id = "th")
binom.plot.exact = pretty.plot(data.melt) + 
              geom_point(aes(x = th, y = value, shape = variable), size = 3) + 
              xlab("Unconstrained diffusion") + 
              ylab("Exact bridge") + 
              geom_abline(slope = 1) +
              scale_y_continuous(breaks = seq(-edge,edge,skip)) + 
              scale_x_continuous(breaks = seq(-edge,edge,skip)) +
              scale_shape_manual(values = c(1)) +
              coord_cartesian(ylim=c(-edge,edge), xlim=c(-edge,edge))+
              theme(legend.position = "none",
                    text = element_text(size=16)) + 
              ggtitle("2Bern - 1")


multiplot(o.plot.approx, exp.plot.approx, binom.plot.approx,
          o.plot.exact, exp.plot.exact, binom.plot.exact,
          cols = 2)

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/approx-vs-exact.pdf", width= 15, height = 10, #' see how it looks at this size
    useDingbats=F)
multiplot(o.plot.exact, o.plot.approx,
          binom.plot.approx, binom.plot.exact, 
          exp.plot.approx, exp.plot.exact,
          cols = 3)
dev.off()