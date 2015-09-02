#' Display multiple plots created by ggplot2 in the same window.
#' @param ... One or more ggplot2 plots
#' @param plotlist One or more ggplot2 plots already grouped into ggplot2 plotlist structure
#' @param layout A matrix of integers that maps the location of each plot to the resulting plot grid
#' @param ncol The number of columns desired in the output grid. Only used if layout=NULL.
#' @author Winston Chang - Code was taken from "Cookbook for R"
#' @examples 
#' # This example uses the ChickWeight dataset, which comes with ggplot2
#' # First plot
#'  p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) +
#'    geom_line() +
#'    ggtitle("Growth curve for individual chicks")
#'
#' # Second plot
#' p2 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet)) +
#'   geom_point(alpha=.3) +
#'   geom_smooth(alpha=.2, size=1) +
#'   ggtitle("Fitted growth curve per diet")
#' 
#' # Third plot
#' p3 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, colour=Diet)) +
#'   geom_density() +
#'   ggtitle("Final weight, by diet")
#' 
#' # Fourth plot
#' p4 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, fill=Diet)) +
#'   geom_histogram(colour="black", binwidth=50) +
#'   facet_grid(Diet ~ .) +
#'   ggtitle("Final weight, by diet") +
#'   theme(legend.position="none")        # No legend (redundant in this graph) 
#'   
#' # A simple 2x2 grid of plots   
#' multiplot(p1, p2, p3, p4, ncol=2)
#' 
#' # Make plot1 take up 3 rows, and put the other three in column two
#' my.layout <- cbind( c(1,1,1), c(2,3,4))
#' my.layout
#' multiplot(p1, p2, p3, p4, layout=my.layout)
#' @export
multiplot <- function(..., plotlist=NULL, layout=NULL, ncol=1) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, ncol * ceiling(numPlots/ncol)),
                     ncol = ncol, nrow = ceiling(numPlots/ncol))
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
