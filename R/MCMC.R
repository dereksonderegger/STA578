#' Graph the first few steps of a MCMC chain
#' @param density A ggplot2 object that shows a density function
#' @param chain A vector of values of the chain
#' @param proposed A vector of proposed values for the chain 
#' @export
plot_chain <- function(density, chain, proposed){
  n <- length(chain)
  data <- data.frame(x=chain, x.star=proposed, y=-(1:n)/n*.2, group=1)
  data$group[n] <- 2
  data$group <- factor(data$group)
  out <- density + 
    geom_path(data=data) +
    geom_point(data=data, aes(x=x.star), color='blue') +
    geom_point(data=data, aes(color=group)) +
    scale_color_manual(values=c('black','red')) +
    theme(legend.position='none') +
    scale_y_continuous(
      breaks=seq(-.2,.3,by=.1),
      labels=c('','',0.0,0.1,0.2,0.3))
  return(out)
}