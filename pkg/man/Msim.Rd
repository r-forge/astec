\name{Msim}
\alias{Msim}
\title{
  Metacommunity Simulator
}
\description{
  This function launches the metacommunity simulator.
}
\usage{
  Msim(parameters)
}
\arguments{
  \item{parameters}{
    parameters of the simulator.
  }
}
%\details{}
\value{
  The returned value is a list containing the following components: 
  \item{param}{
  	The model parameters used in the \code{model} simulations.}
  \item{stats}{
  	The summary statistics obtained at the end of the \code{model} simulations.}
  \item{weights}{
  	The weights of the different \code{model} simulations. In the standard rejection scheme, all \code{model} simulations have the same weights.}
  \item{stats_normalization}{
  	The standard deviation of the summary statistics across the \code{model} simulations.}
  \item{nsim}{
  	The number of \code{model} simulations performed.}
  \item{nrec}{
  	The number of retained simulations (if targeted summary statistics are provided).}
  \item{computime}{
  	The computing time to perform the simulations.}
}
\references{
  Pritchard, J.K., and M.T. Seielstad and A. Perez-Lezaun and
  M.W. Feldman (1999) Population growth of human Y chromosomes: a study
  of Y chromosome microsatellites. \emph{Molecular Biology and
  Evolution}, \bold{16}, 1791--1798.
}
\author{Franck Jabot}
\seealso{
}
\examples{
    ##### EXAMPLE 1 #####
    #####################
    set.seed(1)

    ## artificial example to show how to use the 'Msim' function.
    parameters=1:4
    stats_output=array(0,4)
    Msim(parameters,stats_output)$stats_output
}

