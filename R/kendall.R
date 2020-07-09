# 1. Priors for Kendall's Tau -------------
#
stretchedBetaTau <- function(tauPop, betaA=1, betaB=1) {
  result <- stretchedBeta(sin(pi*tauPop/2), betaA=betaA, betaB=betaB)*cos(pi*tauPop/2)*pi/2
  return(result)
}

stretchedBetaTauSymmetric <- function(tauPop, betaA=1) {
  result <- ((pi*2^(-2*betaA))/beta(betaA, betaA)) * cos((pi*tauPop)/2)^(2*betaA-1)
  return(result)
}

priorTau <- function(tauPop, kappa=1, alternative="two.sided", betaA=NULL, betaB=NULL) {
  priorLine <- switch(alternative,
                      "two.sided"=priorTauTwoSided("tauPop"=tauPop, "kappa"=kappa, "betaA"=betaA, "betaB"=betaB),
                      "greater"=priorTauPlus("tauPop"=tauPop, "kappa"=kappa, "betaA"=betaA, "betaB"=betaB),
                      "less"=priorTauMin("tauPop"=tauPop, "kappa"=kappa, "betaA"=betaA, "betaB"=betaB))
  return(priorLine)
}

priorTauTwoSided <- function(tauPop, kappa=1, betaA=NULL, betaB=NULL) {
  if (is.null(betaA) || is.null(betaB))
    result <- stretchedBetaTauSymmetric("tauPop"=tauPop, "betaA" = 1/kappa)
  else
    result <- stretchedBetaTau("tauPop"=tauPop, "betaA"=betaA, "betaB"=betaB)

  return(result)
}

priorTauPlus <- function(tauPop, kappa=1, betaA=NULL, betaB=NULL) {
  nonNegativeIndex <- tauPop >= 0
  lessThanOneIndex <- tauPop <= 1
  valueIndex <- as.logical(nonNegativeIndex*lessThanOneIndex)
  result <- tauPop*0

  if (is.null(betaA) || is.null(betaB))
    normalisationConstant <- 1/2
  else
    normalisationConstant <- pbeta(0.5, "shape1"=betaA, "shape2"=betaB, "lower.tail"=FALSE)

  if (is.null(betaA) || is.null(betaB)) {
    result[valueIndex] <- stretchedBetaTauSymmetric("tauPop"=tauPop[valueIndex],
                                                    "betaA" = 1/kappa)/normalisationConstant
  } else {
    result[valueIndex] <- stretchedBetaTau("tauPop"=tauPop[valueIndex], "betaA"=betaA,
                                           "betaB"=betaB)/normalisationConstant
  }
  return(result)
}

priorTauMin <- function(tauPop, kappa=1, betaA=NULL, betaB=NULL) {
  negativeIndex <- tauPop <= 0
  greaterThanMinOneIndex <- tauPop >= -1
  valueIndex <- as.logical(negativeIndex*greaterThanMinOneIndex)
  result <- tauPop*0

  if (is.null(betaA) || is.null(betaB))
    normalisationConstant <- 1/2
  else
    normalisationConstant <- pbeta(0.5, "shape1"=betaA, "shape2"=betaB, "lower.tail"=TRUE)

  if (is.null(betaA) || is.null(betaB)) {
    result[valueIndex] <- stretchedBetaTauSymmetric("tauPop"=tauPop[valueIndex],
                                                    "betaA" = 1/kappa)/normalisationConstant
  } else {
    result[valueIndex] <- stretchedBetaTau("tauPop"=tauPop[valueIndex], "betaA"=betaA,
                                           "betaB"=betaB)/normalisationConstant
  }
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
# These are the functions used for the posteriors
#
posteriorTauU <- function(n, tStar, nu=NULL, kappa=1, betaA=NULL, betaB=NULL, var=1,
                          alternative="two.sided") {
  result <- function(tauPop) {
    tempResult <- 0

    for (i in seq_along(n))
      tempResult <- tempResult + stats::dnorm("x"=tStar[i], "mean"=(1.5*tauPop*sqrt(n[i])), "sd"=sqrt(var), log=TRUE)

    tempResult <- exp(tempResult)*priorTau("tauPop"=tauPop, "kappa"=kappa, "alternative"=alternative,
                                           "betaA"=betaA, "betaB"=betaB)
    return(tempResult)
  }
  return(result)
}

makePosteriorTauFunc <- function(n, tauObs, kappa=1, var=1, alternative="two.sided",
                             betaA=NULL, betaB=NULL) {

  tStar <- vector("numeric", length=length(n))

  for (i in seq_along(n))
    tStar[i] <- (tauObs[i] * ((n[i]*(n[i]-1))/2))/sqrt(n[i]*(n[i]-1)*(2*n[i]+5)/18)

  lims <- switch(alternative,
                 "two.sided"=c(-1, 1),
                 "greater"=c(0, 1),
                 "less"=c(-1, 0)
  )

  integrandF <- posteriorTauU("n"=n, "tStar"=tStar, "kappa"=kappa, "betaA"=betaA,
                              "betaB"=betaB, "var"=var, "alternative"=alternative)

  normalisingConstant <- try(silent=TRUE, integrate(integrandF, lims[1], lims[2])[["value"]])

  if (isTryError(normalisingConstant))
    result <- NA
  else
    result <- function(tauPop){integrandF(tauPop)/normalisingConstant}

  return(result)
}

posteriorTau <- function(n, tauObs, tauPop, kappa=1, var=1, alternative="two.sided",
                         betaA=NULL, betaB=NULL) {
  posteriorTauF <- makePosteriorTauFunc("n"=n, "tauObs"=tauObs, "kappa"=kappa, "var"=var,
                                        "alternative"=alternative, "betaA"=betaA,
                                        "betaB"=betaB)
  posteriorTauF(tauPop)
}

computeKendallCredibleInterval <- function(n, tauObs, kappa=1, var=1, ciValue=0.95, h0=0, betaA=NULL, betaB=NULL) {
  # Compute Kendall's correlation credible interval based on a sampling
  #
  check <- failIfNot(ciValue > 0, ciValue < 1, !isSomeNull(n, tauObs))

  sidedResult <- list("lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA, "ciValue"=NA)
  result <- list("two.sided"=sidedResult, "greater"=sidedResult, "less"=sidedResult, "ciValue"=ciValue)

  if (!is.null(check))
    return(result)

  for (alternative in c("two.sided", "greater", "less")) {
    result[[alternative]] <- credibleIntervalKendallTauSingle("n"=n, "tauObs"=tauObs, "kappa"=kappa, "var"=var,
                                                               "alternative"=alternative, "ciValue"=ciValue,
                                                              "betaA"=betaA, "betaB"=betaB)
  }

  return(result)
}

# SINGLE refers to the side
credibleIntervalKendallTauSingle <- function(n, tauObs, kappa=1, var=1, alternative="two.sided",
                                             ciValue = 0.95, m=4000, betaA=NULL, betaB=NULL) {
  # TODO(Alexander): Interesting use-case: n=800, tauObs=-0.8, m=1000
  posteriorTauFunc <- makePosteriorTauFunc("n"=n, "tauObs"=tauObs, "kappa"=kappa, "var"=var,
                                           "alternative"=alternative, "betaA"=betaA, "betaB"=betaB)

  lowCI <- (1-ciValue)/2
  upCI <- (1+ciValue)/2
  tauPopDomain <- seq(-1, 1, length.out=(m-1))

  densVals <- posteriorTauFunc(tauPopDomain)
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
