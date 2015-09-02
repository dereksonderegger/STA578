#' Graph the first few steps of a MCMC chain
#' @param density A ggplot2 object that shows a density function
#' @param chain A vector of values of the chain
#' @param proposed A vector of proposed values for the chain 
#' @examples 
#' # Target Distribution
#' df <- function(x){   
#'    return(.7*dnorm(x, mean=2, sd=1) + .3*dnorm(x, mean=5, sd=1)) 
#' }
#' x <- seq(-3,12, length=200) 
#' density.data <- data.frame(x=x, y=df(x)) 
#' density <- ggplot( density.data, aes(x=x, y=y)) +
#'    geom_line() + 
#'    labs(y='f(x)', x='x') 
#' density
#' 
#' # Define the proposal distribution
#' rproposal <- function(x.i){   
#'   out <- x.i + runif(1, -2, 2)  # x.i + 1 observation from Uniform(-2, 2)   
#'   return(out) 
#' } 
#' 
#' # First chain step
#' x <- 3;       # starting value for the chain 
#' x.star <- 3   # initialize the vector of proposal values
#' x.star[2] <- rproposal( x[1] ) 
#' x.star[2]
#' if( df(x.star[2]) / df(x[1])  >  runif(1) ){
#'   x[2] <- x.star[2]
#' }else{
#'   x[2] <- x[1] 
#' }
#'
#' # Show the first proposed value
#' plot_chain(density, x, x.star)
#' 
#' # Take a few more steps in the chain
#' for( i in 2:10 ){
#'   x.star[i+1] <- rproposal( x[i] )
#'   if( df(x.star[i+1]) / df(x[i]) > runif(1) ){
#'     x[i+1] <- x.star[i+1]
#'   }else{
#'     x[i+1] <- x[i]
#'   }
#' } 
#' 
#' # Show the first 10 steps of the chain
#' plot_chain(density, x, x.star)
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


#' Create a single Markov Chain Monte Carlo simulation
#' @param df A function that is the target distribution
#' @param start The starting point for the MCMC
#' @param rprop A function that takes the current chain step and proposes a new value
#' @param dprop The density function of the proposal distribution. If NULL, it is assumed to be symetric.
#' @param N How many chain steps to calculate.
#' @examples 
#' # Target Distribution
#' df <- function(x){   
#'    return(.7*dnorm(x, mean=2, sd=1) + .3*dnorm(x, mean=5, sd=1)) 
#' }
#' 
#' # Proposal distribution 
#' rprop <- function(x){ x + runif(1, -2, 2)}
#' dprop <- function(x.star, x){ dunif(x.star-x, -2, 2) }
#' 
#' # Create 1 chain
#' chain1 <- MCMC(df, start=4, rprop, dprop, N=1000)
#' trace_plot(chain1)
#' 
#' # Create several chains
#' chain2 <- MCMC(df, start= 8, rprop, dprop, N=1000)
#' chain3 <- MCMC(df, start=-3, rprop, dprop, N=1000)
#' chain4 <- MCMC(df, start=-6, rprop, dprop, N=1000)
#' chains <- cbind(chain1, chain2, chain3, chain4)
#' trace_plot(chains)
#' @export
MCMC <- function(df, start, rprop, dprop=NULL, N=1000){
  if(is.null(dprop)){
    dprop <- function(to,from){1}
  }
  chain <- rep(NA, N)
  chain[1] <- start
  for( i in 1:(N-1) ){
    x.star <- rprop( chain[i] )
    r1 <- df(x.star) / df(chain[i])
    r2 <- dprop(chain[i], x.star) / dprop(x.star, chain[i])
    if( r1*r2 > runif(1) ){
      chain[i+1] <- x.star
    }else{
      chain[i+1] <- chain[i]
    }
  }
  return(chain)
}

#' Create a traceplot from a MCMC chain or chains
#' @export
trace_plot <- function( chain ){
  if(is.null(ncol(chain))){
    out <- ggplot(data.frame(y=chain, x=seq(1:length(chain)))) +
      geom_line(aes(x=x, y=y)) +
      xlab("Iteration") + ylab("Value")
  }else{
    N <- nrow(chain)
    num.chains <- ncol(chain)
    data <- data.frame( Iteration = rep(1:N, times=num.chains),
                        chain     = factor(rep(1:num.chains, each=N)),
                        value     = c(chain))
    out <- ggplot(data, aes(x=Iteration, y=value, color=chain)) +
      geom_path() +
      xlab("Iteration") + ylab("Value")
  }  
  return(out)
}
