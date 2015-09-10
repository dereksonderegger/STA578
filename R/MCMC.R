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

#' Run multiple MCMC chains. 
#' 
#'  This function generalizes the MCMC function to allow for
#'  running multiple chains and alows for multi-diminsionality 
#'  in the distribution. 
#'  
#' @return An object of type \code{\link[coda]{mcmc.list}}. This type of object
#' can be used by functions in the \pkg{coda} package.
#' @param df A function that is the target distribution
#' @param start The starting point for the MCMC
#' @param rprop A function that takes the current chain step and proposes a new value
#' @param dprop The density function of the proposal distribution. If NULL, it is assumed to be symetric.
#' @param N How many chain steps to calculate.
#' @param num.chains How many chains to create.
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
#' # Create 4 chains
#' chains <- mMCMC(df, 4, rprop, dprop, N=1000, num.chains=4)
#' traceplot(chains)
#' 
#' ##############################
#' ##  A multi-variate example
#' ##############################
#' dtarget <- function(x){
#'   dmvnorm(x, mean=c(3,10), sigma=matrix(c(3,3,3,7), nrow=2))
#' }
#' 
#' x1 <- seq(-6,12,length=101)
#' x2 <- seq(-11,31, length=101)
#' contour.data <- expand.grid(x1=x1, x2=x2)
#' contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
#' 
#' target.map <- ggplot(contour.data, aes(x=x1, y=x2, z=Z)) +
#'   stat_contour()
#' target.map
#' 
#' rprop <- function(x){
#'   rmvnorm(1, mean=x, sigma=diag(c(1,1)))
#' }
#' start <- c(x1=0, x2=0)
#' 
#' chain <- mMCMC(df=dtarget, start, rprop=rprop, 
#'                N=20, num.chains=1)
#' target.map + 
#'   geom_path(data=data.frame(chain[[1]]), aes(z=0), 
#'             color='red', alpha=.6) +
#'   geom_point(data=data.frame(chain[[1]]), aes(z=0),
#'              color='red', alpha=.6, size=2)
#' 
#' chain <- mMCMC(df=dtarget, start, rprop=rprop, 
#'                N=100, num.chains=1)
#' target.map + 
#'   geom_path(data=data.frame(chain[[1]]), aes(z=0), 
#'             color='red', alpha=.6) +
#'   geom_point(data=data.frame(chain[[1]]), aes(z=0),
#'              color='red', alpha=.6, size=2)
#' 
#' chain <- mMCMC(df=dtarget, start, rprop=rprop, 
#'                N=1000, num.chains=1)
#' target.map + 
#'   geom_path(data=data.frame(chain[[1]]), aes(z=0), 
#'             color='red', alpha=.6) +
#'   geom_point(data=data.frame(chain[[1]]), aes(z=0),
#'              color='red', alpha=.6, size=2)
#' traceplot( chain )
#' plot(chain)
#' 
#' chains <- mMCMC(df=dtarget, start, rprop=rprop, 
#'                N=1000, num.chains=4)
#' plot(chains)
#' @export
mMCMC <- function(df, start, rprop, dprop=NULL, N=1000, num.chains=4){
  # If no proposal density, assume it is symetric
  if(is.null(dprop)){
    dprop <- function(to,from){1}
  }
  # If we have only one starting point, use the same for every chain
  if( !is.list(start) ){
    start.list <- list(start)
    for( i in 2:num.chains ){
      start.list[[i]] <- start
    }
    start <- start.list
  }
  # If there are no names associated with start, give them names
  for( m in 1:length(start)){
    if(is.null(names(start[[m]]))){
      names(start[[m]]) <- paste('X',1:length(start[[m]]),sep='')
    }
  }
  
  out <- list()
  for( j in 1:num.chains){
    chain <- matrix(NA, nrow=N, ncol=length(start[[j]]), 
                    dimnames=list(NULL, names(start[[j]])))
    chain[1,] <- start[[j]]
    for( i in 1:(N-1) ){
      x.star <- rprop( chain[i,] )
      r1 <- df(x.star) / df(chain[i,])
      r2 <- dprop(chain[i,], x.star) / dprop(x.star, chain[i,])
      if( r1*r2 > runif(1) ){
        chain[i+1,] <- x.star
      }else{
        chain[i+1,] <- chain[i,]
      }
    }
    out[[j]] <- as.mcmc(chain)
  }
  return(as.mcmc.list(out))
}



#' Hamiltonian Monte Carlo
#' 
#' @param dtarget The target density function
#' @param start The starting point. If a vector, all chains start at that point. If 
#' is a list, then each element should be a chain starting point.
#' @param Eps The stepsize epsilon
#' @param L The number of leapfrog steps to take for each MCMC step.
#' @param m The particle mass.
#' @param N The number of Chain steps to take.
#' @param num.chains The number of independent chains to create.
#' @examples 
#' dtarget <- function(x){
#'   dmvnorm(x, mean=c(3,10), sigma=matrix(c(1,0,0,1), nrow=2)) 
#' }
#' x1 <- seq(0,6,length=101) 
#' x2 <- seq(7,13, length=101) 
#' contour.data <- expand.grid(x1=x1, x2=x2) 
#' contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
#' target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
#'   stat_contour(aes(z=Z)) 
#' target.map 
#' 
#' # Poor behavior of the chain
#' chains <- HMC(dtarget, start=c(3,9), Eps=.5, L=20, N=1000, num.chains=4)
#' trace_plot(chains)
#' plot_2D_chains(target.map, chains)
#' # better
#' chains <- H.MCMC(dtarget, start=c(3,9), Eps=.5, L=10, N=1000, num.chains=4)
#' trace_plot(chains)
#' plot_2D_chains(target.map, chains)
#' # Poor again
#' chains <- HMC(dtarget, start=c(3,9), Eps=.5, L=7, N=1000, num.chains=4)
#' trace_plot(chains)
#' plot_2D_chains(target.map, chains)
#'           
#'           
#' # A slightly harder example
#' dtarget <- function(x){
#'   dmvnorm(x, mean=c(3,10), sigma=matrix(c(3,3,3,7), nrow=2)) 
#' }
#' x1 <- seq(-6,12,length=101) 
#' x2 <- seq(-11,31, length=101) 
#' contour.data <- expand.grid(x1=x1, x2=x2) 
#' contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
#' target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
#'   stat_contour(aes(z=Z)) 
#' target.map 
#' 
#' chains <- HMC(dtarget, start=c(2,10), Eps=.4, L=6, N=200, num.chains=4)
#' trace_plot(chains)
#' plot_2D_chains(target.map, chains)
#' 
#' chains <- HMC(dtarget, start=c(2,10), Eps=.4, L=10, N=200, num.chains=4)
#' trace_plot(chains)
#' plot_2D_chains(target.map, chains)
#'
#'
#' # Now a hard distribution!
#' dtarget <- function(x){
#'   B <- .05
#'   exp( -x[1]^2 / 200 - (1/2)*(x[2]+B*x[1]^2 -100*B)^2 )
#' }
#' x1 <- seq(-20,20, length=201)
#' x2 <- seq(-15,10, length=201)
#' contour.data <- expand.grid(x1=x1, x2=x2) 
#' contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
#' target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
#'   stat_contour(aes(z=Z)) 
#' target.map
#'  
#' chains <- HMC(dtarget, start=c(2,5), Eps=.5, L=10)
#' trace_plot(chains)
#' plot_2D_chains(target.map, chains)
#' 
#' chains <- HMC(dtarget, start=c(2,5), Eps=1, L=10)
#' trace_plot(chains)
#' plot_2D_chains(target.map, chains)
#' @export
HMC <- function(dtarget, start, Eps=.2, L=10, m=1, N=1000, num.chains=4){
  
  neg.log.dtarget <- function(x){
    return(-log(dtarget(x)))
  }
  chains <- list()
  
  # If we have only one starting point, use the same for every chain
  if( !is.list(start) ){
    start.list <- list(start)
    for( i in 2:num.chains ){
      start.list[[i]] <- start
    }
    start <- start.list
  }
  # If there are no names associated with start, give them names
  for( m in 1:length(start)){
    if(is.null(names(start[[m]]))){
      names(start[[m]]) <- paste('X',1:length(start[[m]]),sep='')
    }
  }
  
  for(m in 1:num.chains){
    chain <- matrix(NA, nrow=N, ncol=length(start[[m]]), dimnames = list(NULL, names(start[[m]])))
    chain[1, ] <- start[[m]]
    for( i in 2:N){
      chain[i,] <- HMC.one.step(neg.log.dtarget, chain[i-1,], Eps, L, m)$value
    }
    chains[[m]] <- as.mcmc( chain )
  }
  chains <- as.mcmc.list(chains)
  return(chains)
}

#' One step of Hamiltonian Monte Carlo
#' 
#' @param U A function taking a single argument, which is the position.
#' @param q The current position.
#' @param Eps The step size, epsilon.
#' @param L The number of leapfrog steps.
#' @param m The mass of the particle. 
#' @examples
#' dtarget <- function(x){
#'   dmvnorm(x, mean=c(3,10), sigma=matrix(c(1,0,0,1), nrow=2)) 
#' }
#' x1 <- seq(0,6,length=101) 
#' x2 <- seq(7,13, length=101) 
#' contour.data <- expand.grid(x1=x1, x2=x2) 
#' contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
#' target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
#'   stat_contour(aes(z=Z)) 
#' target.map 
#' 
#' U <- function(x){ return( -log(dtarget(x))) }
#' steps <- HMC.one.step(U, c(runif(1,1,5),runif(1,8,12)), Eps=.4, L=20, m=1)
#' plot_HMC_proposal_path(target.map, steps)
#'   
#'  
#' # A little harder...     
#' dtarget <- function(x){
#'   dmvnorm(x, mean=c(3,10), sigma=matrix(c(3,3,3,7), nrow=2)) 
#' }
#' x1 <- seq(-6,12,length=101) 
#' x2 <- seq(-11,31, length=101) 
#' contour.data <- expand.grid(x1=x1, x2=x2) 
#' contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
#' target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
#'   stat_contour(aes(z=Z)) 
#' target.map 
#' 
#' U <- function(x){ return( -log(dtarget(x))) }
#' steps <- HMC.one.step(U, c(runif(1,-2,8),runif(1,4,16)), Eps=.4, L=40, m=1)
#' plot_HMC_proposal_path(target.map, steps)
#' 
#' 
#' # Now a hard example!
#' dtarget <- function(x){
#'   B <- .05
#'   exp( -x[1]^2 / 200 - (1/2)*(x[2]+B*x[1]^2 -100*B)^2 )
#' }
#' x1 <- seq(-20,20, length=201)
#' x2 <- seq(-15,10, length=201)
#' contour.data <- expand.grid(x1=x1, x2=x2) 
#' contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
#' target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
#'   stat_contour(aes(z=Z)) 
#' target.map
#' 
#' U <- function(x){ return( -log(dtarget(x))) }
#' steps <- HMC.one.step(U, c(runif(1,-10,10),runif(1,0,10)), Eps=.15, L=1000, m=1)
#' plot_HMC_proposal_path(target.map, steps)
#' @export
HMC.one.step <- function(U, current_q, Eps, L, m=1){
  out <- list()
  p <- rmvnorm(1, rep(0, length(current_q)))
  current_p <- p
  q <- current_q
  
  out$q <- matrix(NA, nrow=L, ncol=length(current_q))
  out$p <- matrix(NA, nrow=L, ncol=length(current_q))
  out$q[1,] <- q
  out$p[1,] <- p
  
  p <- p - Eps*grad(U,q) / 2
  for(i in 2:L){
    q <- q + Eps*p/m
    out$q[i,] <- q
    if( i != L ){
      p <- p - Eps * grad(U,q)
      out$p[i,] <- p
    }
  }
  p <- p - Eps*grad(U,q) / 2
  out$p[L,] <- p
  p <- -p
  
  current_U <- U(current_q)
  current_K <- sum( current_p^2 )/2
  proposed_U <- U(q)
  proposed_K <- sum( p^2 )/2

  if( runif(1) < exp( current_U - proposed_U + current_K - proposed_K ) ){
    out$value <- q
  }else{
    out$value <- current_q
  }
  return(out)
}

#' Plot a Hamiltonian proposal path
#' @param target.map A density plot of a 2-D distribution
#' @param steps The output of HMC.one.step
#' @export
plot_HMC_proposal_path <- function(target.map, steps){
  path.data <- as.data.frame(steps$q)
  colnames(path.data) <- c('x','y')
  out <- target.map + 
    geom_path(data=path.data, aes(x=x, y=y), color='orange') + 
    geom_point(data=path.data, aes(x=x, y=y), color='orange') + 
    geom_point(data=path.data[1,], aes(x=x, y=y), color='red', size=4) +
    geom_point(data=path.data[nrow(path.data),], aes(x=x, y=y), color='blue', size=4)
  return(out)
}

#' Plot MCMC chains along a 2D density
#' @param target.map A density plot of a 2-D distribution
#' @param chains An mcmc.list
#' @export
plot_2D_chains <- function(target.map, chains){
  data <- ggs(chains) %>% 
      mutate( 
        Parameter = factor(Parameter, levels=unique(Parameter), labels=c('x','y')),
        Chain = factor(Chain)) %>%
      spread(Parameter, value)
  out <- target.map + geom_path(data=data, aes(x=x, y=y, color=Chain))
  return(out)
}


#' Create a traceplot from a MCMC chain or chains
#' @param chain Either a mcmc.list or a matrix or a vector
#' @export
trace_plot <- function( chain ){
  if( is.mcmc.list(chain) ){
    out <- ggs_traceplot( ggs(chain) )
  }else if(is.null(ncol(chain))){
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


#' Calculate the Gelman & Ruben Diagnostic Ratio
#' 
#' @param x Either a matrix or a mcmc.list
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
#' # Create chains
#' chain1 <- MCMC(df, start= 4, rprop, dprop, N=1000)
#' chain2 <- MCMC(df, start= 8, rprop, dprop, N=1000)
#' chain3 <- MCMC(df, start=-3, rprop, dprop, N=1000)
#' chain4 <- MCMC(df, start=-6, rprop, dprop, N=1000)
#' chains <- cbind(chain1, chain2, chain3, chain4)
#' Gelman(chains)
#' 
#' # Or using mMCMC
#' chains <- mMCMC(df, start=list(4,8,-3,-6), rprop, dprop, N=1000, num.chains=4)
#' traceplot(chains)
#' chains <- window(chains, 100, 1000)
#' Gelman(chains)
#' @export
Gelman <- function(x){
  if( is.mcmc.list(x) ){ 
    data <- ggs(x)
  }else{
    # covert matrix to data.frame with 4 columns - Chain, Iteration, Parameter, Value
    X <- as.data.frame(x)
    colnames(X) <- paste('Chain', 1:ncol(x), sep='')
    X$Iteration <- 1:nrow(X)
    X$Parameter <- factor('x1')
    data <- tidyr::gather_( X, key_col='Chain', value_col='value',
                            gather_cols=paste('Chain', 1:ncol(x), sep='') )
  }    
  # Calculate a few easy things.
  xbars <- data %>% group_by(Parameter) %>%
    summarise(xbar = mean(value))
  m <- length( unique(data$Chain) )

  # Use the package dplyr to summarise the data frame into W,B
  temp <- data %>% group_by(Parameter, Chain) %>%
      summarise( s2_i = var(value),
                 xbar_i = mean(value),
                 n = n()) %>%
      merge(xbars) %>%
      group_by(Parameter) %>%
      summarise( n = mean(n),
                 W = mean(s2_i),
                 B = n/(m-1) * sum( (xbar_i - xbar)^2 ) ) %>%
      mutate( sigma2.hat = (n-1)/n * W + (1/n)*B,
              R.hat      = sqrt(sigma2.hat/W),
              n.eff      = m*n * sigma2.hat / B ) %>%
      select( -n )
  return(temp)
}
