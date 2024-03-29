% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pearson.R
\name{bf10CorExact}
\alias{bf10CorExact}
\title{Title}
\usage{
bf10CorExact(n, r, kappa = 1, h0 = 0, oneThreshold = 0.001)
}
\arguments{
\item{n}{integer representing sample size}

\item{r}{numeric in (-1, 1) sample Pearson correlation r}

\item{kappa}{numeric > 0 sample}

\item{h0}{numeric in (-1, 1) which represents the test point of H0}

\item{oneThreshold}{numeric > 0 such that if abs(1 - stat) < oneThreshold, then abs(stat) is viewed as one.}

\item{methodNumber}{integer either 1, or 2. Methodnumber 1 uses the exact results of Ly et.al (2018),
whereas Methodnumber 2 uses the exact integral based on Jeffreys's approximation to the reduced likelihood.}

\item{hyperGeoOverFlowThreshold}{integer set to 24. Some experiments show that when log(bf10) > 24, that the
one-sided Bayes factors become instable [hypergeo v. 1.2-13]. For instance,

myN <- 300
myR <- 0.42 # 0.415

(bf10 <- bf10CorExact(n=myN, stat=myR))
(bfPlus0 <- bf10CorExact(n=myN, stat=myR) + mPlusCorExact(n=myN, stat=myR))
(bfMin0 <- bf10CorExact(n=myN, stat=myR) + mPlusCorExact(n=myN, stat=-myR))

Changing myR <- 0.42. Shows that the result is okay}
}
\value{
Returns bf10, bfPlus0, bfMin0 and the stretched beta fit of the posterior.

Exact BF10 for Pearson's correlation based on the exact reduced likelihood, see Ly et al. (2018).

Returns Bayes factor BF10 in favour of the alternative over the null
}
\description{
#' Exact bf10, bfPlus0, bfMin0 based on the exact
}
\details{
for Pearson's correlation based on the exact reduced likelihood, see Ly et al. (2018).
}
