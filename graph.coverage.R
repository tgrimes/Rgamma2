graph.overall <- function(coverage.list, tests = NULL, significance.level = 0.05) {
  description <- coverage.list[[1]]
  coverage.list[[1]] <- NULL
  
  n <- length(coverage.list)
  p <- ncol(coverage.list[[1]])
  plots <- vector("list", n)
  
  i <- 1
  for(data in coverage.list) {
    parameter.vals <- data[1, 4]
    if(p > 4) {
      for(j in 5:p) {
        parameter.vals <- c(parameter.vals, data[1, j])
      }
    }
    plots[[i]] <- get.coverage.plot(data, parameter.vals, tests, 
                                     significance.level)
    i <- i + 1
  }
  
  layout <- NULL
  for(divisor in c(3, 4, 5, 6, 7, 2)) {
    if((n %% divisor) == 0) {
      for(k in 1:(n/divisor)) {
        layout <- rbind(layout, ((divisor*k - (divisor - 1)):(divisor*k)))
      }
      layout <- rbind(layout, rep(n + 1, divisor))
      get.plot(plots, layout = layout)
      
      return(NULL)
    }
  }
}


get.coverage.plot <- function(coverage.data, parameter.vals, tests = NULL,
                               significance.level) {
  
  main <- paste("(", parameter.vals[1], sep = "")
  i <- 2
  while(i <= length(parameter.vals)) {
    main <- paste(main, ", ", parameter.vals[i])
    i <- i + 1
  }
  main <- paste(main, ")", sep = "")
  
  if(!is.null(tests)) {
    coverage.data <- subset(coverage.data, Test %in% tests)
  }
  
  plot <- ggplot(data = coverage.data, aes(x = Delta, y = Coverage, linetype = Test)) +
    geom_line() + 
    geom_point(size=2, shape=21, fill="white", alpha = 0.5) +
    geom_hline(yintercept = significance.level) +
    labs(title = main, y = "Power") +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size = 10),
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size = 10, face = "bold"))
}



get.plot <- function(plots, distribution, parameter.names, parameter.vals, 
                     layout = rbind(1:3, 4:6, 7:9, rep(10,3))) {
  
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  grid.arrange(grobs = c(gl, list(legend)), layout_matrix = layout,
               heights = grid::unit.c(rep((unit(1, "npc") - lheight) *
                                            (1/(nrow(layout)-1)),nrow(layout)-1), 
                                      lheight))
  
}

