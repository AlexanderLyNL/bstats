% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pearson.R
\name{bf10CorJeffreysIntegrate}
\alias{bf10CorJeffreysIntegrate}
\title{BF10 for Pearson's correlation based on Jeffreys's approximation to the reduced likelihood,
see Jeffreys (1961, pp. 289-292).}
\usage{
bf10CorJeffreysIntegrate(n, r, kappa = 1, h0 = 0)
}
\arguments{
\item{kappa}{numeric > 0 that specifies the}

\item{h0}{numeric between -1 and 1 that specifies the point null hypothesis}
}
\value{
Returns approximate Bayes factor BF10 in favour of the alternative over the null
}
\description{
BF10 for Pearson's correlation based on Jeffreys's approximation to the reduced likelihood,
see Jeffreys (1961, pp. 289-292).
}
