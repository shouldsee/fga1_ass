###### Reference: ttps://stackoverflow.com/questions/35982640/
grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right"), top = NULL) {
  require(ggplot2)
  require(gridExtra)
  require(grid)
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  ltitle = 0.1
  usetitle = !is.null(top)
  ltitle = ltitle * (usetitle)
  # cat(ltitle)
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1.0 - ltitle, "npc") - lheight, lheight)
                                            ,top = top
                     ),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1.0, "npc") - lwidth, lwidth)
                                           ,top = top
                     )
                     
  )
  grid.newpage()
  grid.draw(combined)
}


