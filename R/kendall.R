# 1. Priors for Kendall's Tau -------------
#
stretchedBetaTau <- function(tauPop, alpha=1, beta=1) {
  logResult <- tryOrFailWithNA((-alpha-beta)*log(2) + (alpha-1)*log(1+sin(pi/2*tauPop))
                               + (beta-1)*log(1-sin(pi/2*tauPop)) + log(cos(pi/2*tauPop) - lbeta(alpha, beta))
  )

  if (is.na(logResult))
    result <- tryOrFailWithNA(
      pi * 2^(-alpha-beta)/beta(alpha, beta) * (1+sin(pi/2*tauPop))^(alpha-1) *
        (1-sin(pi/2*tauPop))^(beta-1) * cos(pi/2*tauPop)
    )
  else {
    result <- pi * exp(logResult)
  }
  return(result)
}

stretchedBetaTauSymmetric <- function(tauPop, alpha=1) {
  result <- ((pi*2^(-2*alpha))/beta(alpha, alpha))  * cos((pi*tauPop)/2)^(2*alpha-1)
  return(result)
}


priorTau <- function(tauPop, kappa=1, alternative="two.sided") {
  if (alternative == "two.sided") {
    priorLine <- stretchedBetaTauSymmetric(tauPop, alpha = 1/kappa)
  } else if (alternative == "greater") {
    priorLine <- priorTauPlus("tauPop"=tauPop, "kappa"=kappa)
  } else if (alternative == "less") {
    priorLine <- priorTauMin("tauPop"=tauPop, "kappa"=kappa)
  }
  return(priorLine)
}

priorTauPlus <- function(tauPop, kappa=1) {
  nonNegativeIndex <- tauPop >= 0
  lessThanOneIndex <- tauPop <= 1
  valueIndex <- as.logical(nonNegativeIndex*lessThanOneIndex)
  result <- tauPop*0
  result[valueIndex] <- 2*priorTau(tauPop[valueIndex], kappa)
  return(result)
}

priorTauMin <- function(tauPop, kappa=1) {
  negativeIndex <- tauPop <= 0
  greaterThanMinOneIndex <- tauPop >= -1
  valueIndex <- as.logical(negativeIndex*greaterThanMinOneIndex)
  result <- tauPop*0
  result[valueIndex] <- 2*priorTau(tauPop[valueIndex], kappa)
  return(result)
}

# 2. Bayes factors for Kendall's Tau -------------
# These are the functions used for to compute the Bayes factors
#
computeKendallBCor <- function(n, tauObs, h0=0, kappa=1, ciValue=0.95, var=1,
                                   methodNumber=1L, oneThreshold=1e-3) {
  #
  sidedResult <- list("n"=n, "stat"=tauObs, "bf"=NA, "tooPeaked"=NA,
                      "lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA)

  result <- list("two.sided"=sidedResult,
                 "less"=sidedResult,
                 "greater"=sidedResult,
                 "kappa"=kappa, "ciValue"=ciValue, "h0"=h0,
                 "methodNumber"=methodNumber, "var"=var, "call"=match.call()
  )

  failedSidedResult <- list("n"=n, "stat"=NaN, "bf"=NA, "tooPeaked"=TRUE,
                            "ciValue"=ciValue, "lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA)

  failedResult <- list("two.sided"=failedSidedResult,
                       "less"=failedSidedResult,
                       "greater"=failedSidedResult,
                       "kappa"=kappa, "ciValue"=ciValue, "acceptanceRate"=1,
                       "methodNumber"=6, "call"=match.call()
  )

  # When the prior is trivial (null is alternative) or when the data is predictively matched
  #
  predictiveMatchingList <- list("two.sided"=list("bf"=1, "tooPeaked"=FALSE),
                                 "greater"=list("bf"=1, "tooPeaked"=FALSE),
                                 "less"=list("bf"=1, "tooPeaked"=FALSE),
                                 "methodNumber"=0)

  # Information consistent result
  #
  plusSidedInfList <- list("two.sided"=list("bf"=Inf, "tooPeaked"=TRUE, "posteriorMedian"=tauObs),
                           "less"=list("bf"=0, "tooPeaked"=TRUE, "posteriorMedian"=tauObs),
                           "greater"=list("bf"=Inf, "tooPeaked"=TRUE)
  )

  minSidedInfList <- list("two.sided"=list("bf"=Inf, "tooPeaked"=TRUE, "posteriorMedian"=tauObs),
                          "less"=list("bf"=Inf, "tooPeaked"=TRUE),
                          "greater"=list("bf"=0, "tooPeaked"=TRUE, "posteriorMedian"=tauObs)
  )

  checkTypeInput <- !isEveryNumeric(tauObs, kappa, n)

  if (checkTypeInput) {
    result <- failedResult
    errorMessage <- "Input error: the sample size n, the summary statistic stat, or kappa are not numeric"
    result[["error"]] <- errorMessage
    return(result)
  }

  checkData <- failIfNot(abs(tauObs) <= 1, n >= 0, kappa > 0)

  if (!is.null(checkData)) {
    result <- failedResult
    result[["error"]] <- checkData
    return(result)
  }

  # Note: Data: OK
  # "No" prior, alternative model is the same as the null model
  # The bound kappa=0.002 is chosen arbitrarily I should choose this based on a trade off
  # between r and n, but it doesn't really matter.
  if (kappa <= 0.002 || n <= 2) {
    result <- modifyList(result, predictiveMatchingList)
    return(result)
  }

  checkTauObs <- (1 - abs(tauObs) < oneThreshold) # check whether 1 - |r| < oneThreshold

  # Information consistent result:
  if (kappa >= 1 && n > 2 && checkTauObs) {
    if (tauObs >= 0) {
      result <- modifyList(result, plusSidedInfList)
      return(result)
    } else if (tauObs < 0) {
      result <- modifyList(result, minSidedInfList)
      return(result)
    }
  }

  # TODO(Alexander): Check for information inconsistent kappas
  #
  if (n <= 2 || kappa==0) {
    result <- modifyList(result, predictiveMatchingList)
    return(result)
  } else if (kappa >= 1 && n > 2 && checkTauObs) {
    if (tauObs > 0) {
      result <- modifyList(result, plusSidedInfList)
      return(result)
    } else if (tauObs <= 0) {
      result <- modifyList(result, minSidedInfList)
      return(result)
    }
    result <- modifyList(result, infoConsistentList)
    return(result)
  }

  # Note(Alexander): Add stretched beta fitted posterior
  #
  for (alternative in c("two.sided", "greater", "less")) {
    result[[alternative]] <- bfKendallTauSavageDickey("n"=n, "tauObs"=tauObs, "kappa"=kappa, "var"=var,
                                                       "h0"=h0, "alternative"=alternative)
  }

  tempList <- computeKendallCredibleInterval("n"=n, "tauObs"=tauObs, "kappa"=kappa, "var"=var,
                                              "ciValue"=ciValue)
  result <- modifyList(result, tempList)

  return(result)
}

bfKendallTauSavageDickey <- function(n, tauObs, kappa=1, var=1, h0=0, alternative="two.sided") {
  result <- list("n"=n, "stat"=tauObs, "bf"=NA, "tooPeaked"=TRUE, "lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA)

  bf <- tryOrFailWithNA(
    priorTau("tauPop"=h0, "kappa"=kappa, "alternative"=alternative) /
      posteriorTau("n"=n, "tauObs"=tauObs, "tauPop"=h0, "kappa"=kappa, "var"=var, "alternative"=alternative)
  )

  if (is.finite(bf)) {
    result[["bf"]] <- bf
    result[["tooPeaked"]] <- FALSE
  }

  return(result)
}

# 3. Posteriors ---------
# These are the functions used for to the posteriors
#
posteriorTauU <- function(n, Tstar, tauPop, kappa=1, var=1, alternative="two.sided") {
  #
  result <- stats::dnorm("x"=Tstar, "mean"=(1.5*tauPop*sqrt(n)), "sd"=sqrt(var)) *
    priorTau("tauPop"=tauPop, "kappa"=kappa, "alternative"=alternative )
  return(result)
}

posteriorTau <- function(n, tauObs, tauPop, kappa=1, var=1, alternative="two.sided") {
  Tstar <- (tauObs * ((n*(n-1))/2))/sqrt(n*(n-1)*(2*n+5)/18)
  lims <- switch(alternative,
                 "two.sided"=c(-1, 1),
                 "greater"=c(0, 1),
                 "less"=c(-1, 0)
  )

  logicalCensor <- (tauPop >= lims[1] && tauPop <= lims[2])

  integrand <- function(x) {
    posteriorTauU("n"=n, "Tstar"=Tstar, "tauPop"=x, "kappa"=kappa, "var"=var, "alternative"=alternative)
  }
  normalisingConstant <- try(silent=TRUE, integrate(integrand, lims[1], lims[2])[["value"]])

  if (isTryError(normalisingConstant)) {
    result <- NA
  } else {
    # I could also do
    result <- integrand(tauPop)/normalisingConstant
  }
  return(result)
}


computeKendallCredibleInterval <- function(n, tauObs, kappa=1, var=1, ciValue=0.95, h0=0) {
  # Compute Kendall's correlation credible interval based on a sampling
  #
  check <- failIfNot(ciValue > 0, ciValue < 1, !isSomeNull(n, tauObs))

  sidedResult <- list("lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA, "ciValue"=NA)
  result <- list("two.sided"=sidedResult, "greater"=sidedResult, "less"=sidedResult, "ciValue"=ciValue)

  if (!is.null(check))
    return(result)

  for (alternative in c("two.sided", "greater", "less")) {
    result[[alternative]] <- credibleIntervalKendallTauSingle("n"=n, "tauObs"=tauObs, "kappa"=kappa, "var"=var,
                                                               "alternative"=alternative, "ciValue"=ciValue)
  }

  return(result)
}

# Compute credible intervals kendalls tau
credibleIntervalKendallTauSingle <- function(n, tauObs, kappa=1, var=1, alternative="two.sided", ciValue = 0.95, m=4000) {
  # TODO(Alexander): Interesting use-case: n=800, tauObs=-0.8, m=1000
  lowCI <- (1-ciValue)/2
  upCI <- (1+ciValue)/2
  tauPopDomain <- seq(-1, 1, length.out=(m-1))
  densVals <- posteriorTau("n"=n, "tauObs"=tauObs, "tauPop"=tauPopDomain, "kappa"=kappa, "var"=var,
                            "alternative"=alternative)
  cdfVals <- cumsum((densVals[1:(m-1)] + densVals[2:m]) * 0.5 * (tauPopDomain[2]-tauPopDomain[1]))

  lowerCi <- tauPopDomain[which(cdfVals>=lowCI)[1]]
  upperCi <- tauPopDomain[which(cdfVals>=upCI)[1]]
  posteriorMedian <- tauPopDomain[which(cdfVals>=0.5)[1]]

  result <- list("lowerCi"=lowerCi, "upperCi"=upperCi, "posteriorMedian"=posteriorMedian)

  if (purrr::some(result, is.na)) {
    result <- list("lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA, "tooPeaked"=TRUE,
                   "error"="Can't compute credible intervals")
    return(result)
  }

  if (abs(upperCi-lowerCi) <= .Machine$double.eps)
    result[["tooPeaked"]] <- TRUE

  return(result)
}
