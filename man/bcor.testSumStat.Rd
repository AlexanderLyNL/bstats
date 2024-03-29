% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainWrappers.R
\name{bcor.testSumStat}
\alias{bcor.testSumStat}
\title{Summary stats version of "bcor.test()"}
\usage{
bcor.testSumStat(
  n,
  stat,
  alternative = c("two.sided", "less", "greater"),
  method = c("pearson", "kendall", "spearman"),
  ciValue = 0.95,
  h0 = 0,
  kappa = 1,
  hyperGeoOverFlowThreshold = 25,
  oneThreshold = 0.001,
  var = 1,
  k = 0
)
}
\arguments{
\item{n}{number of observations}

\item{stat}{sample correlation coefficient}

\item{alternative}{indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
"greater" corresponds to positive association, "less" to negative association.}

\item{method}{a character string indicating which correlation coefficient is to b used for the test. One of
"pearson", "kendall" or "spearman".}

\item{ciValue}{numeric in (0, 1) credible level for the returned credible interval.}

\item{h0}{numeric between -1 and 1 that specifies the point null hypothesis}

\item{kappa}{numeric > 0 that specifies the}

\item{hyperGeoOverFlowThreshold}{numeric > 0 such that if log(bf10) > hyperGeoOverFlowThreshold then the
Savage-Dickey adaptation is used for to compute the one-sided Bayes factors, instead of the analytical ones.}

\item{oneThreshold}{numeric > 0 such that if abs(1 - stat) < oneThreshold, then abs(stat) is viewed as one.}

\item{var}{numeric > 0 that specifies the asymptotic variance of the approximate likelihood for Kendall's tau}

\item{k}{number of controlling variables}
}
\value{
A list with class "btest" containing the following components:
}
\description{
Summary stats version of "bcor.test()"
}
\examples{
bcor.testSumStat(n=34, stat=0.4)
bcor.testSumStat(n=34, stat=0.4, method="kendall")
bcor.testSumStat(n=34, stat=0.45, k = 1)
}
