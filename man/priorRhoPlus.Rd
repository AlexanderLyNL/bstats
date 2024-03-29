% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pearson.R
\name{priorRhoPlus}
\alias{priorRhoPlus}
\title{Prior for Pearson's rho restricted to positive values}
\usage{
priorRhoPlus(rho, kappa = 1, betaA = NULL, betaB = NULL)
}
\arguments{
\item{rho}{numeric in (-1, 1) at which the prior needs to be evaluated}

\item{kappa}{numeric > 0 which provides the scale of}
}
\value{
numeric, the density value at rho
}
\description{
Prior for Pearson's rho restricted to positive values
}
\examples{
priorRhoPlus(0.4)
priorRhoPlus(-0.4)

rhoDomain <- seq(-0.99, 0.99, by=0.01)

# kappa <- 2, e.g., beta(1/2, 1/2), is needed for information consistency
yLine <- priorRhoPlus(rhoDomain, kappa=2)

graphics::plot(rhoDomain, yLine, type="l")
}
