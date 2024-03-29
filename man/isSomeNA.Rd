% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0helpers.R
\name{isSomeNA}
\alias{isSomeNA}
\title{Checks whether some object is NA}
\usage{
isSomeNA(...)
}
\arguments{
\item{...}{objects that need testing}
}
\value{
Returns TRUE if there's some object that's an NA, FALSE when all objects are not NA.
}
\description{
Checks whether some object is NA
}
\examples{
x <- 1
y <- "a"
z <- NA

isSomeNA(x, y)
isSomeNA(x, y, z)
}
