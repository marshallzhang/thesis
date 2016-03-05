library(ggplot2)
library(reshape2)

plot.many = function(data, N, ...){
  plot(data[1, ], type = "l", col = rgb(0,0,0,1/(N/50)), ...)
  for (i in 1:N) {
    lines(data[i, ], col = rgb(0,0,0,1/(N/75)))
  }
}

pretty.plot = function(melted.data) {
  ggplot(data = data.melt) + theme_bw(base_size = 12, base_family = "Helvetica")
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