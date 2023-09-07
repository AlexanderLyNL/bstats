# 1. Main test function --------
#' Bayesian Test for Association/Correlation Between Paired Samples
#' This mimics cor.test
#'
#' @param x,y numeric vectors of data values. x and y must have the same length
#' @param z data frame of data values of the controlling variables.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' "greater" corresponds to positive association, "less" to negative association.
#' @param method a character string indicating which correlation coefficient is to b used for the test. One of
#' "pearson", "kendall" or "spearman".
#' @param ciValue numeric in (0, 1) credible level for the returned credible interval.
#' @param use
#' @param h0 numeric between -1 and 1 that specifies the point null hypothesis
#' @param kappa numeric > 0 that specifies the
#' @param oneThreshold numeric > 0 such that if abs(1 - stat) < oneThreshold, then abs(stat) is viewed as one.
#' @param var numeric > 0 that specifies the asymptotic variance of the approximate likelihood for Kendall's tau
#' @param hyperGeoOverFlowThreshold numeric > 0 such that if log(bf10) > hyperGeoOverFlowThreshold then the
#' Savage-Dickey adaptation is used for to compute the one-sided Bayes factors, instead of the analytical ones.
#'
#' @return A list with class "btest" containing the following components:
#'
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#'
#' bcor.test(x, y)
#' bcor.test(x, y, method="kendall")
#'
#' z <- data.frame(rnorm(100))
#' bcor.test(x, y, z)
bcor.test <- function(x, y, z = NULL, alternative=c("two.sided", "less", "greater"),
                      method=c("pearson", "kendall", "spearman"), ciValue=0.95,
                      use="pairwise.complete.obs", h0=0, kappa=1, hyperGeoOverFlowThreshold=25,
                      oneThreshold=0.001, var=1) {
  stopifnot(length(x) == length(y))
  stopifnot(is.null(z) || length(x) == nrow(z))

  method <- match.arg(method)
  if (is.null(method)) {
    result <- computePearsonBCor(NULL, NULL)
    result[["error"]] <- "No method selected"
    return(result)
  }

  if (is.null(z)) {
    k <- 0
    xyz <- cbind(x, y)
    stat <- tryOrFailWithNA(cor(x, y, use=use, method=method))
  } else {
    k <- ncol(z)
    xyz <- cbind(x, y, z)
    stat <- tryOrFailWithNA(sampleParCor(xyz, use=use, method=method))
  }

  n <- sum(complete.cases(xyz))

  result <- bcor.testSumStat(n=n, stat=stat, alternative=alternative,
                             method=method, ciValue=ciValue, h0=h0, kappa=kappa,
                             hyperGeoOverFlowThreshold=hyperGeoOverFlowThreshold,
                             oneThreshold=oneThreshold, var=var, k = k)
  result[["call"]] <- match.call()
  return(result)
}

#' Summary stats version of "bcor.test()"
#'
#' @param n number of observations
#' @param stat sample correlation coefficient
#' @param k number of controlling variables
#' @inherit bcor.test
#'
#' @export
#'
#' @examples
#' bcor.testSumStat(n=34, stat=0.4)
#' bcor.testSumStat(n=34, stat=0.4, method="kendall")
#' bcor.testSumStat(n=34, stat=0.45, k = 1)
bcor.testSumStat <- function(n, stat, alternative=c("two.sided", "less", "greater"),
                             method=c("pearson", "kendall", "spearman"), ciValue=0.95,
                             h0=0, kappa=1, hyperGeoOverFlowThreshold=25, oneThreshold=0.001, var=1, k = 0) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  if (method == "spearman")
    warning("`method`=='spearman' was called but Spearman method is not yet implemented!")

  if (is.na(stat) || is.na(n) || n <= k) {
    if (method == 'kendall') {
      result <- computeKendallBCor(NaN, NaN)
    } else {
      result <- computePearsonBCor(NaN, NaN)
    }
    result[["error"]] <- "Can't compute the correlation"
  } else {
    result <- switch(
      method,
      pearson=computePearsonBCor(n=n-k, r=stat, h0=h0, kappa=kappa, ciValue=ciValue,
                                 hyperGeoOverFlowThreshold=hyperGeoOverFlowThreshold,
                                 oneThreshold = oneThreshold),
      kendall=computeKendallBCor(n=n, tauObs=stat, h0=h0, kappa=kappa, ciValue=ciValue,
                                 oneThreshold=oneThreshold, var=var),
      spearman=computePearsonBCor(n=NA, r=NA, h0=h0, kappa=kappa, ciValue=ciValue,
                                  oneThreshold=oneThreshold)
      )
  }

  result[["alternative"]] <- alternative
  result[["method"]] <- method
  result[["k"]] <- k
  return(result)
}

#' Function to grab sided ("two.sided", "greater", "less") information of an bfObject with extra items
#'
#' @param bfObject
#' @param alternative
#' @param itemNames optional character strings, that specify additional items from bfObject
#'
#' @return a (sub)list of bfObject
#'
#' @examples
getSidedObject <- function(bfObject, alternative="two.sided", itemNames=NULL) {
  result <- modifyList(bfObject[[alternative]], bfObject[itemNames])
}


# 2 Compute posteriors --------
#' Title
#'
#' @param bfObject
#' @param method
#' @param alternative
#' @param minX
#' @param maxX
#'
#' @return
#' @export
#'
#' @examples
computeCorPosteriorLine <- function(bfObject, alternative="two.sided", minX=-0.985, maxX=0.985) {
  method <- bfObject[["method"]]

  xDomain <- seq(minX, maxX, length.out = 1001)

  if (alternative %in% c("two-sided", "two.sided")) {
    alternative <- "two.sided"
  } else if (alternative %in% c("greater", "right", "positive")) {
    alternative <- "greater"
  } else if (alternative %in% c("less", "left", "negative")) {
    alternative <- "less"
  }

  sidedObject <- getSidedObject("bfObject"=bfObject, "alternative"=alternative,
                                 itemNames=c("error", "h0", "betaA", "betaB", "ciValue"))

  # Note(Alexander): Don't compute if it's already computed
  #
  if (!is.null(sidedObject[["posteriorLine"]]))
    return(sidedObject)

  # Note(Alexander): Don't compute if there's an error
  #
  errorMessage <- sidedObject[["error"]]

  if (!is.null(errorMessage)) {
    sidedObject[["posteriorLine"]] <- errorMessage
    return(sidedObject)
  }

  # Note(Alexander): Don't compute if it's too peaked
  #
  if (sidedObject[["tooPeaked"]]) {
    sidedObject[["posteriorLine"]] <- "Posterior is too peaked"
    return(sidedObject)
  }

  if (method=="pearson") {
    if (isSomeNA(bfObject[["betaA"]], bfObject[["betaB"]])) {
      sidedObject[["posteriorLine"]] <- "Posterior is too peaked"
      return(sidedObject)
    }

    # TODO(Alexander): Derive this for shifted h0/Find master student
    #
    sidedObject[["priorAtH0"]] <- priorRho("rho"=sidedObject[["h0"]], "kappa"=bfObject[["kappa"]],
                                           "alternative"=alternative)
    sidedObject[["priorLine"]] <- priorRho("rho"=xDomain, "kappa"=bfObject[["kappa"]],
                                           "alternative"=alternative)

    subCounter <- 1

    while (subCounter <= 3) {
      if (subCounter==1) {
        posteriorAtH0 <- posteriorRho("bfObject"=bfObject, "rho"=sidedObject[["h0"]], "alternative"=alternative)
        posteriorLine <- posteriorRho("bfObject"=bfObject, "rho"=xDomain, "alternative"=alternative)
      } else if (subCounter==2) {
        sidedObject[["approximation"]] <- "beta"
        posteriorAtH0 <- posteriorRhoBetaApprox("rho"=sidedObject[["h0"]], "bfObject"=bfObject,
                                                 alternative=alternative)
        posteriorLine <- posteriorRhoBetaApprox("rho"=xDomain, "bfObject"=bfObject,
                                                 alternative=alternative)
      } else if (subCounter==3) {
        sidedObject[["approximation"]] <- "fisher"
        posteriorAtH0 <- posteriorRhoFisherApprox("bfObject"=bfObject, "rho"=sidedResult[["h0"]],
                                                   "alternative"=alternative)
        posteriorLine <- posteriorRhoFisherApprox("bfObject"=bfObject, "rho"=xDomain,
                                                   "alternative"=alternative)
      }

      if (isSomeNA(posteriorLine) || any(posteriorLine < 0) || isSomeInfinite(posteriorLine)) {
        subCounter <- subCounter +1
      } else {
        sidedObject[["posteriorAtH0"]] <- posteriorAtH0
        sidedObject[["posteriorLine"]] <- posteriorLine
        break()
      }
    }
  } else if (method=="kendall") {
    sidedObject[["priorAtH0"]] <- priorTau("tauPop"=sidedObject[["h0"]], "kappa"=bfObject[["kappa"]],
                                            "alternative"=alternative)
    sidedObject[["priorLine"]] <- priorTau("tauPop"=xDomain, "kappa"=bfObject[["kappa"]], "alternative"=alternative)

    subCounter <- 1

    while (subCounter <= 1) {
      if (subCounter==1) {
        posteriorAtH0 <- posteriorTau("tauPop"=sidedObject[["h0"]], "n"=sidedObject[["n"]],
                                       "tauObs"=sidedObject[["stat"]],"kappa"=bfObject[["kappa"]],
                                       "alternative"=alternative)
        posteriorLine <- posteriorTau("tauPop"=xDomain, "n"=sidedObject[["n"]], "tauObs"=sidedObject[["stat"]],
                                       "kappa"=bfObject[["kappa"]], "alternative"=alternative)
      } else if (subCounter==2) {
        # TODO(Alexander): Perhaps add a Fisher approximation to this
      }

      if (isSomeNA(posteriorLine) || any(posteriorLine < 0) || isSomeInfinite(posteriorLine)) {
        subCounter <- subCounter + 1
      } else {
        sidedObject[["posteriorAtH0"]] <- posteriorAtH0
        sidedObject[["posteriorLine"]] <- posteriorLine
        break()
      }
    }
  }

  if (isSomeNA(posteriorLine) || any(posteriorLine < 0) || isSomeInfinite(posteriorLine)){
    sidedObject[["posteriorLine"]] <- "Posterior is too peaked"
    return(sidedObject)
  }

  sidedObject[["yMax"]] <- max(sidedObject[["priorLine"]], sidedObject[["posteriorLine"]])
  sidedObject[["xDomain"]] <- xDomain
  return(sidedObject)
}

#' Title
#'
#' @param x
#' @param y
#' @param bfObject
#' @param method
#'
#' @return
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' bfObject <- bcor.test(x, y, method="kendall")
#'
#' computeCorSequentialLine <- computeCorSequentialLine(x, y, bfObject)
#' graphics::plot(computeCorSequentialLine$nDomain,
#'                log(computeCorSequentialLine$less$sequentialLine),
#'                type="l", xlab="n", ylab=expression(log("BF"[10])))
computeCorSequentialLine <- function(x, y, bfObject) {
  # sidedObject <- getSidedObject(bfObject, alternative=alternative)
  #
  method <- bfObject[["method"]]

  error <- bfObject[["error"]]
  bf10 <- bfObject[["two.sided"]][["bf"]]

  if (bfObject[["two.sided"]][["tooPeaked"]] || bfObject[["greater"]][["tooPeaked"]] |
      bfObject[["less"]][["tooPeaked"]]) {
    error <- "Posterior is too peaked"
  }

  if (!is.null(error) || is.na(bf10)) {
    if (is.null(error)) {
      # Note(Alexander): This means that there's no error message
      error <- "Could not compute"
    }

    sideError <- list("sequentialLine"=error)
    result <- list("two.sided"=sideError, "greater"=sideError, "less"=sideError)
    return(result)
  }

  n <- bfObject[[1]][["n"]]

  compN <- 3:n
  nDomain <- c(1:2, compN)

  calculateSequentialCor <- function(i, x, y, method) {
    return(tryOrFailWithNA(cor(x[1:i], y[1:i], use="pairwise.complete.obs", method=method)))
  }

  statSeq <- purrr::map_dbl(compN, calculateSequentialCor, "x"=x, "y"=y, "method"=method)

  if (sum(is.na(statSeq)) >= 1) {
    sideError <- list("sequentialLine"="Could not compute")
    result <- list("two.sided"=sideError, "greater"=sideError, "less"=sideError)
  }

  h0 <- bfObject[["h0"]]
  kappa <- bfObject[["kappa"]]

  methodNumber <- bfObject[["methodNumber"]]

  placeHolder <- vector("numeric", length=length(nDomain))
  placeHolder[1] <- placeHolder[2] <- 1
  sideResult <- list("sequentialLine"=placeHolder)

  result <- list("two.sided"=sideResult, "greater"=sideResult, "less"=sideResult)

  if (method=="pearson") {
    calculateSequentialBCorPearson <- function(n, r) {
      bfObject <- computePearsonBCor(n, r, "h0"=h0, "kappa"=kappa,
                                         methodNumber=methodNumber)
      list(bfObject[["two.sided"]][["bf"]],
           bfObject[["greater"]][["bf"]],
           bfObject[["less"]][["bf"]]
      )
    }

    allBfs <- purrr::map2(compN, statSeq, calculateSequentialBCorPearson)

    for (i in seq_along(allBfs)){
      result[[1]][[1]][i+2] <- allBfs[[i]][[1]]
      result[[2]][[1]][i+2] <- allBfs[[i]][[2]]
      result[[3]][[1]][i+2] <- allBfs[[i]][[3]]
    }

  } else if (method=="kendall") {
    calculateSequentialBCorTau <- function(n, tauObs, alternative) {
      bfObject <- unlist(bfKendallTauSavageDickey("n"=n, "tauObs"=tauObs, "kappa"=kappa,"h0"=h0, alternative=alternative)[["bf"]])
    }

    alternativeItems <- c("two.sided", "greater", "less")

    for (i in seq_along(alternativeItems)) {
      alternative <- alternativeItems[i]
      tempResult <- purrr::map2_dbl(compN, statSeq, calculateSequentialBCorTau, alternative=alternative)
      result[[alternative]][[1]] <- c(1, 1, tempResult)
    }
  } else if (method=="spearman") {
    # TODO(Johnny)
  }

  for (j in 1:3) {
    if (isSomeInfinite(result[[j]][[1]])) {
      result[[j]][[1]] <- "Bayes factor hits infinity"
    }

    if (sum(is.na(result[[j]][[1]])) >= 1) {
      result[[j]][[1]] <- "Some of the Bayes factors in the sequence could not be computed, likely caused by an extreme feature of the data, such as a large number of observations and observed correlation."
    }
  }
  result[["nDomain"]] <- nDomain
  return(result)
}

#' Title
#'
#' @param bfObject
#' @param method
#'
#' @return
#' @export
#'
#' @examples
#' bfObject <- bcor.testSumStat(n=34, stat=0.3)
#' result <- computeCorRobustnessLine(bfObject)
#'
#' xLine <- result$kappaDomain
#' yLine <- result$two.sided$robustnessLine
#' yMax <- result$two.sided$robustnessMaxBf
#' plot(xLine, yLine, ylim=c(0, yMax), type="l")
computeCorRobustnessLine <- function(bfObject) {
  method <- bfObject[["method"]]
  error <- bfObject[["error"]]
  bf10 <- bfObject[["two.sided"]][["bf"]]

  if (bfObject[["two.sided"]][["tooPeaked"]] || bfObject[["greater"]][["tooPeaked"]] |
      bfObject[["less"]][["tooPeaked"]]) {
    error <- "Posterior is too peaked"
  }

  if (!is.null(error) || is.na(bf10)) {
    if (is.null(error)) {
      # Note(Alexander): This means that there's no error message
      error <- "Could not compute"
    }

    sideError <- list("robustnessLine"=error)
    result <- list("two.sided"=sideError, "greater"=sideError, "less"=sideError)
    return(result)
  }

  n <- bfObject[["two.sided"]][["n"]]
  stat <- bfObject[["two.sided"]][["stat"]]
  h0 <- bfObject[["h0"]]
  methodNumber <- bfObject[["methodNumber"]]

  kappas <- makeKappas(50)
  compKappas <- kappas[3:50]

  kappaDomain <- c(kappas[1:2], compKappas)

  placeHolder <- c(1, 1, vector("numeric", length=48))
  sideResult <- list("robustnessLine"=placeHolder)
  result <- list("two.sided"=sideResult, "greater"=sideResult, "less"=sideResult)

  if (method=="pearson") {
    calculatePearsonRobustness <- function(kappa) {
      bfObject <- computePearsonBCor("n"=n, "r"=stat, "h0"=h0, "kappa"=kappa, "methodNumber"=methodNumber)
      list(bfObject[["two.sided"]][["bf"]],
           bfObject[["greater"]][["bf"]],
           bfObject[["less"]][["bf"]])
    }

    allBfs <- purrr::map(compKappas, .f=calculatePearsonRobustness)

    for (i in seq_len(48)){
      result[[1]][[1]][i+2] <- allBfs[[i]][[1]]
      result[[2]][[1]][i+2] <- allBfs[[i]][[2]]
      result[[3]][[1]][i+2] <- allBfs[[i]][[3]]
    }
  } else if (method=="kendall") {
    calculateKendallRobustness <- function(kappa, alternative) {
      bfKendallTauSavageDickey("n"=n, "tauObs"=stat, "kappa" = kappa, "h0"=0, "alternative"=alternative)[["bf"]]
    }

    alternativeItems <- c("two.sided", "greater", "less")

    for (i in seq_along(alternativeItems)) {
      alternative <- alternativeItems[i]
      tempResult <- purrr::map_dbl(compKappas, ".f"=calculateKendallRobustness, "alternative"=alternative)
      result[[alternative]][[1]] <- c(1, 1, tempResult)
    }
  }

  for (j in 1:3) {
    if (isSomeInfinite(result[[j]][[1]])) {
      result[[j]][[1]] <- "Bayes factor hits infinity"
    }

    if (sum(is.na(result[[j]][[1]])) >= 1) {
      result[[j]][[1]] <- "Could not compute"
    } else {
      robustnessMaxBf <- max(result[[j]][[1]])
      robustnessKappaOfMaxBf <- kappaDomain[which.max(result[[j]][[1]])]
      result[[j]] <- modifyList(result[[j]],
                                list("robustnessMaxBf"=robustnessMaxBf,
                                     "robustnessKappaOfMaxBf"=robustnessKappaOfMaxBf))
    }
  }
  result[["kappaDomain"]] <- kappaDomain

  return(result)
}
