% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/MCMC.R
\name{HMC}
\alias{HMC}
\title{Hamiltonian Monte Carlo}
\usage{
HMC(dtarget, start, Eps = 0.2, L = 10, m = 1, N = 1000,
  num.chains = 4)
}
\arguments{
\item{dtarget}{The target density function}

\item{start}{The starting point. If a vector, all chains start at that point. If
is a list, then each element should be a chain starting point.}

\item{Eps}{The stepsize epsilon}

\item{L}{The number of leapfrog steps to take for each MCMC step.}

\item{m}{The particle mass.}

\item{N}{The number of Chain steps to take.}

\item{num.chains}{The number of independent chains to create.}
}
\description{
Hamiltonian Monte Carlo
}
\examples{
dtarget <- function(x){
  dmvnorm(x, mean=c(3,10), sigma=matrix(c(1,0,0,1), nrow=2))
}
x1 <- seq(0,6,length=101)
x2 <- seq(7,13, length=101)
contour.data <- expand.grid(x1=x1, x2=x2)
contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
  stat_contour(aes(z=Z))
target.map

# Poor behavior of the chain
chains <- HMC(dtarget, start=c(3,9), Eps=.5, L=20, N=1000, num.chains=4)
trace_plot(chains)
plot_2D_chains(target.map, chains)
# better
chains <- H.MCMC(dtarget, start=c(3,9), Eps=.5, L=10, N=1000, num.chains=4)
trace_plot(chains)
plot_2D_chains(target.map, chains)
# Poor again
chains <- HMC(dtarget, start=c(3,9), Eps=.5, L=7, N=1000, num.chains=4)
trace_plot(chains)
plot_2D_chains(target.map, chains)


# A slightly harder example
dtarget <- function(x){
  dmvnorm(x, mean=c(3,10), sigma=matrix(c(3,3,3,7), nrow=2))
}
x1 <- seq(-6,12,length=101)
x2 <- seq(-11,31, length=101)
contour.data <- expand.grid(x1=x1, x2=x2)
contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
  stat_contour(aes(z=Z))
target.map

chains <- HMC(dtarget, start=c(2,10), Eps=.4, L=6, N=200, num.chains=4)
trace_plot(chains)
plot_2D_chains(target.map, chains)

chains <- HMC(dtarget, start=c(2,10), Eps=.4, L=10, N=200, num.chains=4)
trace_plot(chains)
plot_2D_chains(target.map, chains)


# Now a hard distribution!
dtarget <- function(x){
  B <- .05
  exp( -x[1]^2 / 200 - (1/2)*(x[2]+B*x[1]^2 -100*B)^2 )
}
x1 <- seq(-20,20, length=201)
x2 <- seq(-15,10, length=201)
contour.data <- expand.grid(x1=x1, x2=x2)
contour.data$Z <- apply(contour.data, MARGIN=1, dtarget)
target.map <- ggplot(contour.data, aes(x=x1, y=x2)) +
  stat_contour(aes(z=Z))
target.map

chains <- HMC(dtarget, start=c(2,5), Eps=.5, L=10)
trace_plot(chains)
plot_2D_chains(target.map, chains)

chains <- HMC(dtarget, start=c(2,5), Eps=1, L=10)
trace_plot(chains)
plot_2D_chains(target.map, chains)
}

