
stretchedBeta <- function(rho, betaA, betaB) {
  result <- 1/2*dbeta((rho+1)/2, betaA, betaB)
  return(result)
}


#' Prior for Pearson's rho
#'
#' @param rho numeric in (-1, 1) at which the prior needs to be evaluated
#' @param kappa numeric > 0 which provides the scale of
#' @param alternative
#'
#' @return
#' @export
#'
#' @examples
priorRho <- function(rho, kappa=1, alternative="two.sided") {
  if (alternative == "two.sided") {
    priorLine <- stretchedBeta("rho"=rho, "betaA"=1/kappa, "betaB"=1/kappa)
  } else if (alternative == "greater") {
    priorLine <- priorRhoPlus("rho"=rho, "kappa"=kappa)
  } else if (alternative =="less" ) {
    priorLine <- priorRhoMin("rho"=rho, "kappa"=kappa)
  }
  return(priorLine)
}

priorRhoPlus <- function(rho, kappa=1) {
  nonNegativeIndex <- rho >=0
  lessThanOneIndex <- rho <=1
  valueIndex <- as.logical(nonNegativeIndex*lessThanOneIndex)
  result <- rho*0
  result[valueIndex] <- 2*priorRho(rho[valueIndex], kappa)
  return(result)
}

priorRhoMin <- function(rho, kappa=1) {
  negativeIndex <- rho <=0
  greaterThanMinOneIndex <- rho >= -1
  valueIndex <- as.logical(negativeIndex*greaterThanMinOneIndex)
  result <- rho*0
  result[valueIndex] <- 2*priorRho(rho[valueIndex], kappa)
  return(result)
}

# 2. Likelihood -------------
# These are the functions used for the likelihood
#
aFunction <- function(n, r, rho) {
  hyperTerm <- Re(hypergeo::f15.1.1("A"=1-n/2, "B"=1-n/2, "C"=1/2, "z"=(r*rho)^2))
  result <- exp(((n-1)/2)*log(1-rho^2) + (3/2-n)*log(1-(r*rho)^2))*hyperTerm
  return(result)

  # hyperTerm <- Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=1/2, "z"=(r*rho)^2))
  # result <- (1-rho^2)^((n-1)/2)*(1-(r*rho)^2)^(3/2-n)*hyperTerm
}

bFunction <- function(n, r, rho) {
  hyperTerm <- Re(hypergeo::f15.1.1("A"=(3-n)/2, "B"=(3-n)/2, "C"=3/2, "z"=(r*rho)^2))
  # hyperTerm <- Re(hypergeo::f15.3.3("A"=n/2, "B"=n/2, "C"=3/2, "z"=(r*rho)^2))
  logTerm <- 2*(lgamma(n/2)-lgamma((n-1)/2))+((n-1)/2)*log(1-rho^2)+(3/2-n)*log(1-(r*rho)^2)
  result <- 2*r*rho*exp(logTerm)*hyperTerm
  return(result)
}

hFunction <- function(n, r, rho) {
  result <- aFunction("n"=n, "r"=r, rho) + bFunction("n"=n, "r"=r, rho)
  return(result)
}

hFunctionCombined <- function(nOri, rOri, nRep, rRep, rho) {
  result <- hFunction(n=nOri, r=rOri, rho)*hFunction(n=nRep, r=rRep, rho)
  return(result)
}

hFunctionCombinedTwoSided <- function(nOri, rOri, nRep, rRep, rho) {
  result <- aFunction(n=nOri, r=rOri, rho)*aFunction(n=nRep, r=rRep, rho) +
    bFunction(n=nOri, r=rOri, rho)*bFunction(n=nRep, r=rRep, rho)
  return(result)
}

hJeffreysApprox <- function(n, r, rho) {
  logResult <- ((n-1)/2)*log(1 - rho^2)-(n-3/2)*log(1 - rho*r)
  return(exp(logResult))
  # result <- ((1 - rho^(2))^(0.5*(n - 1)))/((1 - rho*r)^(n - 1 - 0.5))
  # return(result)
}

# 1.1 Explicit marginal likelihood functions
m0MarginalLikelihood <- function(s, t, n) {
  logTerm <- 2*lgamma(0.5*(n-1))
  result <- 1/4*n^(0.5*(1-2*n))*pi^(1-n)*(s*t)^(1-n)*exp(logTerm)
  return(result)
}

m1MarginalLikelihoodNoRho <- function(s, t, n, r, rho) {
  return(m0MarginalLikelihood(s, t, n)*
           (aFunction(n=n, r=r, rho)+bFunction(n=n, r=r, rho)))
}


#' Function returns the p value from correlation.
#'
#' @param r numeric in (-1, 1) representing the sample correlation
#' @param n integer representing sample size
#' @param method character, "pearson", "kendall" or "spearman
#'
#' @return list of three p-values
#' @export
#'
#' @examples
pValueFromCor <- function(n, stat, method=c("pearson","kendall", "spearman")) {
  sidedList <- list("p"=NA)
  result <- list("two.sided"=sidedList,
                 "greater"=sidedList,
                 "less"=sidedList)

  if (n <= 2){
    # result[["two.sided"]] <- sidedList
    # # tau < 0
    # result[["less"]] <- sidedList
    # # tau > 0
    # result[["greater"]] <- sidedList
    return(result)
  }

  if (method[1] == "pearson"){
    # Use t-distribution based on bivariate normal assumption using r to t transformation
    #
    df <- n - 2
    t <- stat*sqrt(df/(1-stat^2))
    result <- pValueFromT("t"=t, "n1"=n-1, "n2"=0, var.equal=TRUE)
  } else if (method[1] == "kendall"){
    if (n > 2 && n < 50) {
      # Exact sampling distribution
      # tau neq 0
      result[["two.sided"]][["p"]] <- 1 - SuppDists::pKendall("q"=abs(stat), N=n) +
        SuppDists::pKendall("q"=-abs(stat), "N"=n)
      # tau < 0
      result[["less"]][["p"]] <- SuppDists::pKendall("q"=stat, "N"=n)
      # tau > 0
      result[["greater"]][["p"]] <- SuppDists::pKendall("q"=stat, "N"=n, "lower.tail" = FALSE)
    } else if (n >= 50){
      # normal approximation
      #
      someSd <- sqrt(2*(2*n+5)/(9*n*(n-1)))

      # tau neq 0
      result[["two.sided"]][["p"]] <- 2 * stats::pnorm(-abs(stat), "sd"=someSd)
      # tau < 0
      result[["less"]][["p"]] <- stats::pnorm(stat, "sd"=someSd)
      # tau > 0
      result[["greater"]][["p"]] <- stats::pnorm(stat, "sd"=someSd, lower.tail = FALSE)
    }
  } else if (method[1] == "spearman"){
    # TODO: Johnny
    # Without code this will print a NULL, if we go through here
  }
  return(result)
}


# TODO(Alexander): Unify ".pValueFromT", and check that this is not used in Alexandra's rewrite or somehwere else
# The important difference is that it now outputs a list
#
#' Function returns the p value from correlation.
#'
#' @param t numeric representing observed t-statistic
#' @param n1 integer representing sample size of the first sample
#' @param n2 integer representing sample size of the second sample, n2=0 implies that it's a one-sample
#' @param var.equal logical, TRUE is equal, cannot have unequal, because we don't have access to s1 and s2
#'
#' @return list of three p-values
#' @export
#'
#' @examples
pValueFromT <- function(t, n1, n2 = 0, var.equal = TRUE) {
  # Function returns the p value from t statistic
  #
  # Args:
  #   t: t value input by user
  #   n1: sample size of group 1
  #   n2: sample size of group 2 (Note the hack by setting n2 = 0)
  #   var.equal: Note: always true: var.equal, we do not have enough info for different
  #              variances. In that case we also need s1 and s2
  #
  # Output:
  #   number in [0, 1] which is the p value
  sidedList <- list("p"=NA)
  result <- list("two.sided"=sidedList,
                 "greater"=sidedList,
                 "less"=sidedList)

  if (n2 > 0) {
    # If n2 > 0, then two-sample
    someDf <- n1 + n2 - 2
  } else {
    # If n2 <= 0, then one-sample
    someDf <- n1 - 1
  }
  # mu \neq 0
  result[["two.sided"]][["p"]] <- 2 * stats::pt(-abs(t), "df" = someDf)
  # mu < 0
  result[["less"]][["p"]] <- stats::pt(t, "df" = someDf)
  # mu > 0
  result[["greater"]][["p"]] <- stats::pt(t, "df" = someDf, "lower.tail" = FALSE)

  return(result)
}

# 3. Bayes factor
# These are the functions used to compute the Bayes factors
#

# 3.1 Two-sided main Bayes factor ----------------------------------------------
## Suit:
#' Title
#'
#'
#'#' Exact bf10, bfPlus0, bfMin0 based on the exact
#'
#'  for Pearson's correlation based on the exact reduced likelihood, see Ly et al. (2018).
#'
#' @param n integer representing sample size
#' @param r numeric in (-1, 1) sample Pearson correlation r
#' @param kappa numeric > 0 sample
#' @param h0 numeric in (-1, 1) which represents the test point of H0
#' @param methodNumber integer either 1, or 2. Methodnumber 1 uses the exact results of Ly et.al (2018),
#' whereas Methodnumber 2 uses the exact integral based on Jeffreys's approximation to the reduced likelihood.
#' @param hyperGeoOverFlowThreshold integer set to 24. Some experiments show that when log(bf10) > 24, that the
#' one-sided Bayes factors become instable [hypergeo v. 1.2-13]. For instance,
#'
#' myN <- 300
#' myR <- 0.42 # 0.415
#'
#' (bf10 <- bf10CorExact(n=myN, stat=myR))
#' (bfPlus0 <- bf10CorExact(n=myN, stat=myR) + mPlusCorExact(n=myN, stat=myR))
#' (bfMin0 <- bf10CorExact(n=myN, stat=myR) + mPlusCorExact(n=myN, stat=-myR))
#'
#' Changing myR <- 0.42. Shows that the result is okay
#'
#' @export
#' @return Returns bf10, bfPlus0, bfMin0 and the stretched beta fit of the posterior.
#'

#' Exact BF10 for Pearson's correlation based on the exact reduced likelihood, see Ly et al. (2018).
#'
#' @inheritParams bcor.test
#'
#' @export
#' @return Returns Bayes factor BF10 in favour of the alternative over the null
#'
bf10CorExact <- function(n, r, kappa=1, h0=0, oneThreshold=1e-3) {
  # Note (Alexander): Input check is done at a higher level: computePearsonBCor
  # Maximum it can take is
  #     r=0.993, which works well up to n = 337, but r=0.992 and n=3 fails
  #     r=0.6, which works up to n=3201

  checkR <- (1 - abs(r) < oneThreshold)

  if (checkR) {
    if (kappa < 2/(n-2)) {
      result <- tryOrFailWithNA(
        exp(lbeta(1/kappa+(n-1)/2, 1/2) - lbeta(1/kappa, 1/2) +
              lgamma(1/kappa) + lgamma(1/kappa -n/2 + 1) - 2*lgamma(1/kappa+1/2))
      )
      return(result)
    } else {
      return(Inf)
    }
  }

  logHyperTerm <- tryOrFailWithNA(
    (1/kappa+1-n/2)*log(1-r^2) +
      log(Re(hypergeo::f15.1.1("A"=(1/kappa+1/2), "B"=(1/kappa+1/2), "C"=(1/kappa+n/2), "z"=r^2)))
  )

  if (is.na(logHyperTerm))
    return(NaN)

  logBetaTerms <- lbeta(1/kappa+(n-1)/2, 1/2)-lbeta(1/kappa, 1/2)

  logResult <- logBetaTerms + logHyperTerm
  realResult <- tryOrFailWithNA(exp(Re(logResult)))

  if (!is.numeric(realResult))
    return(NaN)

  # Failed
  if (realResult < 0)
    return(NaN)

  return(realResult)
}


# 2.2 Two-sided secondairy Bayes factor
#' BF10 for Pearson's correlation based on Jeffreys's approximation to the reduced likelihood,
#' see Jeffreys (1961, pp. 289-292).
#'
#' @inheritParams bcor.test
#'
#' @return Returns approximate Bayes factor BF10 in favour of the alternative over the null
bf10CorJeffreysIntegrate <- function(n, r, kappa=1, h0=0) {
  # Note(Alexander): Input check is done at a higher level: computePearsonBCor

  hyperTerm <- tryOrFailWithNA(Re(hypergeo::f15.1.1(A=3/4+1/kappa, B=1/4+1/kappa, C=(n+2/kappa)/2, z=r^2)))

  # hyperTerm <- tryOrFailWithNA(Re(hypergeo::f15.3.3(A=(2*n-3)/4, B=(2*n-1)/4, C=(n+2/kappa)/2, z=r^2)))

  if (is.na(hyperTerm))
    return(NaN)

  logTerm <- lgamma((n+2/kappa-1)/2) - lgamma((n+2/kappa)/2) - lbeta(1/kappa, 1/kappa) + (1+1/kappa-n/2) * log(1-r^2)
  result <- tryOrFailWithNA(sqrt(pi) * 2^(1-2/kappa) * exp(logTerm) * hyperTerm)

  if (!is.numeric(result))
    return(NaN)

  # Failed
  if (result < 0)
    return(NaN)

  return(result)
}

computeBCorOneSided <- function(bf10, n, r, kappa, methodNumber, betaA, betaB,
                                 hyperGeoOverFlowThreshold=25) {
  failSided <- list("bf"=NA, "tooPeaked"=TRUE)
  result <- list("two.sided"=list("bf"=bf10),
                 "greater"=failSided,
                 "less"=failSided)

  if (methodNumber %in% 1:2) {
    if (log(bf10) < hyperGeoOverFlowThreshold) {
      # Exact method doesn't work well for log(bf10) > 25
      subCounter <- 1L
    } else {
      subCounter <- 2L
    }
  } else {
    subCounter <- 3L
  }

  while (subCounter <= 3) {
    if (subCounter==1L) {
      # (a) Try exact calculations:
      #
      if (methodNumber==1L) {
        bfPlus0 <- tryOrFailWithNA(bf10 + mPlusCorExact(n=n, r=r, kappa))
        bfMin0 <- tryOrFailWithNA(bf10 + mPlusCorExact(n=n, r=-r, kappa))
      } else if (methodNumber==2L) {
        bfPlus0 <- tryOrFailWithNA(bf10 + mPlusCorJeffreysIntegrate(n=n, r=r, kappa=kappa))
        bfMin0 <- tryOrFailWithNA(bf10 + mPlusCorJeffreysIntegrate(n=n, r=-r, kappa=kappa))
      }
    } else if (subCounter==2L) {
      # (b) Compute numerical
      #
      plusSidedIntegrand <- function(x){hJeffreysApprox(n=n, r=r, x)*priorRhoPlus(x, kappa=kappa)}
      minSidedIntegrand <- function(x){hJeffreysApprox(n=n, r=r, x)*priorRhoMin(x, kappa=kappa)}

      bfPlus0 <- tryOrFailWithNA(bf10 + integrate(plusSidedIntegrand, 0, 1)[["value"]])
      bfMin0 <- tryOrFailWithNA(bf10 + integrate(minSidedIntegrand, -1, 0)[["value"]])
    } else if (subCounter==3L) {
      if (!isSomeNA(betaA, betaB)) {
        # Compute bfPlus0 and bfMin0 based on the posterior mass using the found bf10
        #
        tempList <- computeBCorOneSidedSavageDickey("bf10"=bf10, "betaA"=betaA, "betaB"=betaB,
                                                     "h0"=h0, "kappa"=kappa)
        bfPlus0 <- tempList[["bfPlus0"]]
        bfMin0 <- tempList[["bfMin0"]]
      } else {
        return(result)
      }
    }

    if (!isSomeNA(bfPlus0, bfMin0) && isEveryFinite(bfPlus0, bfMin0) &&
        bfPlus0 >= 0 && bfMin0 >= 0) {
      # Result seem good, do consistency check
      if (!(bfPlus0 > 1 && bfMin0 > 1) && bfPlus0 > 0 && bfMin0 > 0 || subCounter==3) {
        tempList <- list("greater"=list("bf"=bfPlus0, "tooPeaked"=FALSE),
                         "less"=list("bf"=bfMin0, "tooPeaked"=FALSE))
        result <- modifyList(result, tempList)
        return(result)
      }
    }
    subCounter <- subCounter + 1
  }
  return(result)
}


computeBCorOneSidedSavageDickey <- function(bf10, betaA, betaB, h0=0, kappa=1) {
  result <- list(bf10=bf10, bfPlus0=NA, bfMin0=NA)

  if (is.finite(bf10)) {
    # bf10 is finite, now calculate one-sided stuff
    #
    rightProportion <- NA
    leftProportion <- tryOrFailWithNA(stats::pbeta(1/2, shape1=betaA, shape2=betaB))

    if (is.na(leftProportion)) {
      return(result)
    }

    if (leftProportion > 0 && leftProportion < 1) {
      bfPlus0 <- 2*bf10*(1-leftProportion)
      bfMin0 <- 2*bf10*leftProportion
    } else if (leftProportion >= 1) {
      bfPlus0 <- 10^(-317)
      bfMin0 <- 2*bf10
    } else {
      rightProportion <- tryOrFailWithNA(stats::pbeta(1/2, "shape1"=betaA, "shape2"=betaB, "lower.tail"=FALSE))

      if (is.na(rightProportion)) {
        return(result)
      } else if (rightProportion > 0 && rightProportion < 1) {
        bfPlus0 <- 2*bf10*rightProportion
        bfMin0 <- 2*bf10*(1-rightProportion)
      } else if (rightProportion >= 1) {
        bfPlus0 <- 2*bf10
        bfMin0 <- 10^(-317)
      }
    }

    # Note(Alexander): Consistency check
    #
    if (bfPlus0 > 1 && bfMin0 > 1) {
      if (leftProportion > 0.5) {
        # Note(Alexander): The problem is machine precision here,
        # because there aren't enough significant figures in leftProportion to have
        # 2 * bf10 (1-leftProportion) < 1
        #
        #  Alternatively, we can change bf10 to the Savage-Dickey approximation
        #
        bfMin0 <- 2*leftProportion*bf10
        bfPlus0 <- 2*(1-leftProportion)
      } else {
        bfMin0 <- 2*leftProportion
        bfPlus0 <- 2*(1-leftProportion)*bf10
      }
    } else if (bfPlus0 < 0 || bfMin0 < 0) {
      return(result)
    }

    #
    result <- list(bf10=bf10, bfPlus0=bfPlus0, bfMin0=bfMin0)
  }
  # return default results
  return(result)
}

# 2.4 The Marsman MH sampler

logCorTarget <- function(rho, n, r, kappa=1) {
  # More correctly, with rho=rhoProp, this is upper term of the acceptance ratio.
  # Recall that in the numerator and denominator of the acceptance ratio are given by
  #
  #   likelihood(rhoProp)*prior(rhoProp)*proposal(rhoCurrent || rhoProp)
  #   likelihood(rhoCurrent)*prior(rhoCurrent)*proposal(rhoProp || rhoCurrent)
  #
  # respectively, see Robert (2015). As the proposal is defined on Fisher transform atanh(rho),
  # we have a Jocobian. The Jacobian in the upper part is (1-rhoCurrent^2)^(-1) and in the lower part
  # is (1-rhoProp^2)^(-1).
  # Note that since prior(rhoProp) \propto (1-rhoProp^2)^(alpha-1), we can drop the one due
  # to the Jacobian term of the denominator ending up in the numerator.
  #
  # For the likelihood we use Jeffreys's approximation
  #   rho \propto (1-rho^2)^((n-1)/2)*(1-rho*r)^((3-2*n)/2)
  #
  # The log prior is (alpha-1)*log(1-rho^2) # where alpha=1/kappa
  #
  return((1/kappa+(n-1)/2)*log(1-rho^2)-(2*n-3)/2*log(1-rho*r))
}


logCorProposal <- function(rho, n, r, z=NULL) {
  # The proposal is defined on Fisher hyperbolic tangent transform of rho,
  # the Jacobian is absorbed in logCorTarget
  if (!is.null(z)) {
    return(-(n-3)/2*(z-atanh(r))^2)
  } else {
    return(-(n-3)/2*(atanh(rho)-atanh(r))^2)
  }
}

#' Fits a beta density to the the posterior samples based on an independent Metropolis sampler, see Tierney (1994).
#'
#' @inheritParams bcor.test
#'
#' @return Returns a list with the beta parameters of a stretched beta distribution on (-1, 1) and the acceptance rate.
#' @export
#'
marsmanCorMHSampler <- function(n, r, kappa=1, nIters=50000L) {
  rhoMetropolisChain <- numeric(nIters)

  if (n <= 3) {
    std <- 1
  } else {
    std <- 1 / sqrt(n - 3)
  }

  zCandidates <- rnorm(nIters, "mean"=atanh(r), "sd"=std)
  rhoCandidates <- tanh(zCandidates)
  logTargetCandidates <- logCorTarget("rho"=rhoCandidates, "n"=n, "r"=r, "kappa"=kappa)
  logPropCandidates <- logCorProposal("z"=zCandidates, "n"=n, "r"=r)
  acceptMechanism <- runif(nIters)
  candidateAcceptance <- numeric(nIters)

  rhoCurrent <- r

  for (iter in 1:nIters) {
    zCurrent <- atanh(rhoCurrent)

    candidateAcceptance[iter] <- logTargetCandidates[iter]+logCorProposal("z"=zCurrent, "n"=n, "r"=r)-
      (logCorTarget("rho"=rhoCurrent, "n"=n, "r"=r, "kappa"=kappa)+logPropCandidates[iter])

    # Accept candidate and update rhoCurrent for next iteration
    if (log(acceptMechanism[iter]) <= candidateAcceptance[iter])
      rhoCurrent <- rhoCandidates[iter]

    rhoMetropolisChain[iter] <- rhoCurrent
  }

  acceptanceRate <- length(unique(rhoMetropolisChain))/nIters

  metropolisVar <- var(rhoMetropolisChain)/2^2
  metropolisMean <- mean((rhoMetropolisChain+1)/2)

  mhFit <- betaParameterEstimates(metropolisMean, metropolisVar)
  mhFit[["acceptanceRate"]] <- acceptanceRate
  return(mhFit)
}


# 3.2 One-sided preparation ----------------------------------------------------

#' Add this to the exact two-sided Bayes factor to get bfPlus0
#'
#' @inherit bf10CorExact
mPlusCorExact <- function(n, r, kappa=1) {
  # Ly et al 2015
  # This is the contribution of one-sided test
  #
  # Note(Alexander): Input check is done at a higher level: computePearsonBCor.
  # In particular the case with n <= 2
  #
  # TODO(Alexander): In hindsight, I'm not so happy with this version, due to instability of 3F2.
  # Try to simplify this expression
  #
  hyperTerm <- Re(hypergeo::genhypergeo(U=c(1, n/2, n/2), L=c(3/2, (2+kappa*(n+1))/(2*kappa)), z=r^2, maxiter=1e5))
  logTerm <- 2*(lgamma(n/2)-lgamma((n-1)/2))-lbeta(1/kappa, 1/kappa)
  result <- 2^((3*kappa-2)/kappa)*kappa*r/(2+(n-1)*kappa)*exp(logTerm)*hyperTerm
  return(result)
}


#' Add this to the semi-exact two-sided Bayes factor to get bfPlus0
#'
#' @inherit bf10CorExact
mPlusCorJeffreysIntegrate <- function(n, r, kappa=1, ...) {
  # Ly et al 2015
  # This is the exact result with symmetric beta prior on rho
  # This is the contribution of one-sided test
  #
  # Note (Alexander): Input check is done at a higher level: computePearsonBCor.
  # In particular the case with n <= 2
  #
  #
  hyperTerm <- Re(
    hypergeo::genhypergeo(U=c(1, (2*n-1)/4, (2*n+1)/4),L=c(3/2, (n+1+2/kappa)/2), z=r^2, ...)
  )
  logTerm <- -lbeta(1/kappa, 1/kappa)
  result <- 2^(1-2/kappa)*r*(2*n-3)/(n+2/kappa-1)*exp(logTerm)*hyperTerm
  return(result)
}

#' Compute Bayes factors bf10, bfPlus0, bfMin0 based on the default stretched beta prior on (-1, 1) for
#' Pearson's correlation coefficient rho
#'
#' @param n numeric > 0, number of samples
#' @param r numeric in (-1, 1), the observed Pearson correlation r in the sample
#' @param h0 numeric in (-1, 1), the hypothesised null value
#' @param kappa numeric > 0, the scale of the beta distribution, i.e., beta(1/kappa, 1/kappa)
#' @param ciValue numeric in (0, 1), the credible value
#' @param hyperGeoOverFlowThreshold numeric > 0, the threshold function for which some computations for which the
#' function genhypergeo in hypergeo [version 1.2-13] used for the one-sided Bayes factors lead to some instablish
#' results.
#' @param methodNumber numeric in reference to the computation method used to calculate Bayes factors: (1) Exact
#' results, see Ly et al. (2018) [when log(bf10) > hyperGeoOverFlowThreshold, the one-sided Bayes factors are
#' calculated using a numerical integrator, or a Savage Dickey method redistribution], (2) Semi-exact,
#' see Wagenmakers et al. (2015), (3) Savage-Dickey: First fit a stretched beta to the posterior based on the exact
#' posterior mean and variance, then compute the ratio of prior to posterior at the restriction point, that is, h0
#' (4) First generate posterior samples using Marsman's IMH sampler, which is the used to fit a stretched beta and
#' again the ratio of prior and posterior is used to compute bf10, (5) use Jeffreys's approximation to bf10 based on
#' kappa=1
#' @param oneThreshold numeric > 0, used to determine when abs(r) is considered indistinguishable from one.
#'
#' @return Returns a list with bf10, bfPlus0, bfMin0, whether the result is plottable, the credible interval.
#' @export
#'
#' @examples
computePearsonBCor <- function(n, r, h0=0, kappa=1, ciValue=0.95, hyperGeoOverFlowThreshold=25,
                                   methodNumber=1L, oneThreshold=1e-3) {
  #
  sidedResult <- list("n"=n, "stat"=r, "bf"=NA, "tooPeaked"=NA,
                      "lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA)

  result <- list("two.sided"=sidedResult,
                 "less"=sidedResult,
                 "greater"=sidedResult,
                 "kappa"=kappa, "ciValue"=ciValue, "acceptanceRate"=1, "h0"=h0,
                 "methodNumber"=methodNumber, "call"=match.call(), "betaA"=NA, "betaB"=NA
  )

  failedSidedResult <- list("n"=n, "stat"=NaN, "bf"=NA, "tooPeaked"=TRUE,
                            "ciValue"=ciValue, "lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA)

  failedResult <- list("two.sided"=failedSidedResult,
                       "less"=failedSidedResult,
                       "greater"=failedSidedResult,
                       "kappa"=kappa, "ciValue"=ciValue, "acceptanceRate"=1,
                       "methodNumber"=6, "call"=match.call(), "betaA"=NA, "betaB"=NA
  )

  # When the prior is trivial (null is alternative) or when the data is predictively matched
  #
  predictiveMatchingList <- list("two.sided"=list("bf"=1, "tooPeaked"=FALSE),
                                 "greater"=list("bf"=1, "tooPeaked"=FALSE),
                                 "less"=list("bf"=1, "tooPeaked"=FALSE),
                                 "methodNumber"=0)

  # Information consistent result
  #
  plusSidedInfList <- list("two.sided"=list("bf"=Inf, "tooPeaked"=TRUE, "posteriorMedian"=r),
                           "less"=list("bf"=0, "tooPeaked"=TRUE, "posteriorMedian"=r),
                           "greater"=list("bf"=Inf, "tooPeaked"=TRUE)
  )

  minSidedInfList <- list("two.sided"=list("bf"=Inf, "tooPeaked"=TRUE, "posteriorMedian"=r),
                          "less"=list("bf"=Inf, "tooPeaked"=TRUE),
                          "greater"=list("bf"=0, "tooPeaked"=TRUE, "posteriorMedian"=r)
  )

  # Note: If log(bf10) then use beta fit for the one sided bfs
  # The 100L is randomlyish chosen based on n=3000, r=-0.3500786
  # hyperGeoOverFlowThreshold <- 100L

  checkTypeInput <- !isEveryNumeric(r, kappa, n)

  if (checkTypeInput) {
    result <- failedResult
    errorMessage <- "Input error: the sample size n, the summary statistic stat, or kappa are not numeric"
    result[["error"]] <- errorMessage
    return(result)
  }

  checkData <- failIfNot(abs(r) <= 1, n >= 0, kappa > 0)

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

  checkR <- (1 - abs(r) < oneThreshold) # check whether 1 - |r| < oneThreshold

  if (n <= 2 || kappa==0) {
    result <- modifyList(result, predictiveMatchingList)
    return(result)
  } else if (checkR && kappa >= 2/(n-2)) {
    if (r > 0) {
      result <- modifyList(result, plusSidedInfList)
    } else if (r <= 0) {
      result <- modifyList(result, minSidedInfList)
    }
    return(result)
  }

  # Note(Alexander): Add stretched beta fitted posterior
  #
  fit <- tryOrFailWithNA(posteriorBetaParameters(n=n, r=r, kappa=kappa))

  # Note(Alexander): Here compute the two sided bf
  #
  while (methodNumber <= 5) {
    # Note: Try all normal methods
    # 1. Exact Ly et al (2018)
    # 2. Exact with Jeffreys's (1961) approximation to the reduced likelihood
    # 3. Savage-Dickey: Via stretched beta fitted to the posterior moments
    # 4. Savage-Dickey: Via Marsman sampler and a stretched beta
    # 5. Jeffreys's approximation to the Bayes factor, only works for h0=0 though
    #     TODO(Alexander): Approximate for h0 \neq 0
    #
    if (methodNumber == 1) {
      bf10 <-  tryOrFailWithNA(bf10CorExact(n=n, r=r, kappa=kappa, h0=h0, oneThreshold=oneThreshold))
    } else if (methodNumber == 2) {
      bf10 <- tryOrFailWithNA(bf10CorJeffreysIntegrate(n=n, r=r, kappa=kappa, h0=h0))
    } else if (methodNumber %in% 3:4) {
      if (methodNumber==4) {
        fit <- marsmanCorMHSampler(n=n, r=r, kappa=kappa)
        result[["acceptanceRate"]] <- fit[["acceptanceRate"]]
      }

      priorAtH0 <- tryOrFailWithNA(priorRho(rho=h0, kappa=kappa))
      posteriorAtH0 <- tryOrFailWithNA(stretchedBeta(rho=h0, betaA=fit[["betaA"]], betaB=fit[["betaB"]]))

      bf10 <- tryOrFailWithNA(priorAtH0/posteriorAtH0)
    } else if (methodNumber == 5) {
      bf10 <- tryOrFailWithNA(
        ((2*n-3)/pi)^(-0.5)*(1-r^2)^((4-n)/2)
      )

      # Note(Alexander): So the methodNumber == 5 is retained
      result[["methodNumber"]] <- 5
      break()
    }

    if (!is.na(bf10) && is.finite(bf10)) {
      result[["two.sided"]][["bf"]] <- bf10
      result[["two.sided"]][["tooPeaked"]] <- FALSE
      result[["methodNumber"]] <- methodNumber
      break()
    }
    methodNumber <- methodNumber+1
  }

  if (is.na(bf10) || bf10 < 0) {
    # Total failure; it's true
    result <- failedResult
    return(result)
  }

  # Note(Alexander): Here compute one-sided bfs
  #
  if (is.infinite(bf10)) {
    if (r >= 0) {
      result <- modifyList(result, plusSidedInfList)
      return(result)
    } else {
      result <- modifyList(result, minSidedInfList)
      return(result)
    }
  } else if (is.finite(bf10)) {
    tempList <- tryOrFailWithNA(
      computeBCorOneSided("bf10"=bf10, "n"=n, "r"=r, "kappa"=kappa,
                           "methodNumber"=methodNumber, "betaA"=fit[["betaA"]], "betaB"=fit[["betaB"]],
                           "hyperGeoOverFlowThreshold" = hyperGeoOverFlowThreshold)
    )

    if (isSomeNA(tempList)) {
      # result <- failedResult
      errorMessage <- "Can't compute one-sided Bayes factors"
      result[["error"]] <- errorMessage
      # return(result)
    }

    result <- modifyList(result, tempList)
  }

  if (!isSomeNA(fit)) {
    tempList <- computePearsonCredibleInterval("betaA"=fit[["betaA"]], "betaB"=fit[["betaB"]],
                                                "ciValue"=ciValue)

    # Add list $two.sided$ci, $plusSided$ci, $less$ci, ciValue
    result <- modifyList(result, tempList)
  } else {
    tempFailedCiResult <- list("tooPeaked"=TRUE)
    failedCiResult <- list("two.sided"=tempFailedCiResult, "greater"=tempFailedCiResult, "less"=tempFailedCiResult)
    result <- modifyList(result, failedCiResult)

    # TODO(Alexander): Consider the case that everything can be computed, except for the cis, then get error.
    # Most of this is covered using Gauss' Vandermonde Identity
    result[["error"]] <- "Can't compute credible intervals"
  }
  return(result)
}


# 4. Posteriors ------------
#


# 4.1 Two-sided
#'
#' @inheritParams computePearsonBCor
#'
#' @return Returns the posterior density
#' @export
#'
#' @examples
posteriorRho <- function(bfObject, rho, alternative="two.sided") {
  sidedObject <- getSidedObject(bfObject, "alternative"=alternative, itemNames=c("kappa"))

  bf <- sidedObject[["bf"]]
  n <- sidedObject[["n"]]
  r <- sidedObject[["stat"]]
  kappa <- sidedObject[["kappa"]]

  if (alternative=="two.sided") {
    return(1/bf*hFunction("n"=n, "r"=r, rho)*priorRho(rho, kappa))
  } else if (alternative=="greater") {
    return(1/bf*hFunction("n"=n, "r"=r, rho)*priorRhoPlus(rho, kappa))
  } else if (alternative=="less") {
    return(1/bf*hFunction("n"=n, "r"=r, rho)*priorRhoMin(rho, kappa))
  }
}

# 4.1 Two-sided
#' ASDF
#'
#' @inheritParams computePearsonBCor
#'
#' @return Returns the posterior density
#' @export
#'
#' @examples
posteriorRhoStat <- function(n, r, rho, kappa=1, alternative="two.sided") {
  if (alternative=="two.sided") {
    if (!is.na(r)) {
      return(1/bf10CorExact("n"=n, "r"=r, "kappa"=kappa)*hFunction("n"=n, "r"=r, rho)*priorRho(rho, kappa))
    } else if (!is.na(r) && r==0) {
      return(1/bf10CorJeffreysIntegrate(n=n, "r"=r, kappa)*hJeffreysApprox(n=n, "r"=r, rho)*priorRho(rho, kappa))
    }
  } else if (alternative=="greater") {
    return(posteriorRhoPlusStat("n"=n, "r"=r, "rho"=rho, "kappa"=kappa))
  } else if (alternative=="less") {
    return(posteriorRhoMinStat("n"=n, "r"=r, "rho"=rho, "kappa"=kappa))
  }
}

posteriorRhoBetaApprox <- function(bfObject, rho, alternative="two.sided") {
  if (alternative == "two.sided") {
    posteriorLine <- stretchedBeta(rho, "betaA"=bfObject[["betaA"]], betaB=bfObject[["betaB"]])
  } else if (alternative == "greater") {
    posteriorLine <- stretchedBeta(rho, "betaA"=bfObject[["betaA"]], betaB=bfObject[["betaB"]]) /
      pbeta("q"=1/2,  "shape1"=bfObject[["betaA"]], "shape2"=bfObject[["betaB"]], "lower.tail"=FALSE)
    posteriorLine[rho < 0] <- 0
  } else if (alternative == "less") {
    posteriorLine <- stretchedBeta(rho, "betaA"=bfObject[["betaA"]], "betaB"=bfObject[["betaB"]]) /
      pbeta("q"=1/2,  "shape1"=bfObject[["betaA"]], "shape2"=bfObject[["betaB"]], "lower.tail"=TRUE)
    posteriorLine[rho > 0] <- 0
  }
  return(posteriorLine)
}

posteriorRhoBetaApproxStat <- function(rho, betaA, betaB, alternative="two.sided") {
  if (alternative == "two.sided") {
    posteriorLine <- stretchedBeta(rho, betaA=betaA, betaB=betaB)
  } else if (alternative == "greater") {
    posteriorLine <- stretchedBeta(rho, betaA=betaA, betaB=betaB) /
      pbeta(q=1/2,  shape1=betaA, shape2=betaB, lower.tail=FALSE)
    posteriorLine[rho < 0] <- 0
  } else if (alternative == "less") {
    posteriorLine <- stretchedBeta(rho, betaA=betaA, betaB=betaB) /
      pbeta(q=1/2,  shape1=betaA, shape2=betaB, lower.tail=TRUE)
    posteriorLine[rho > 0] <- 0
  }
  return(posteriorLine)
}

posteriorRhoPlusStat <- function(n, r, rho, kappa=1) {
  if (!is.na(r)) {
    return(1/computePearsonBCor(n=n, "r"=r, kappa=kappa, methodNumber=1)[["greater"]][["bf"]]
           *hFunction(n=n, "r"=r, rho)*priorRhoPlus(rho, kappa))
  } else if (!is.na(r) && r==0) {
    return(1/computePearsonBCor(n=n, "r"=r, kappa=kappa, methodNumber=2)[["greater"]][["bf"]]
           *hJeffreysApprox(n=n, "r"=r, rho)*priorRhoPlus(rho, kappa))
  }
}

posteriorRhoMinStat <- function(n, r, rho, kappa=1) {
  if (!is.na(r)) {
    return(1/computePearsonBCor("n"=n, "r"=r, "kappa"=kappa, methodNumber=1)[["less"]][["bf"]]
           *hFunction("n"=n, "r"=r, rho)*priorRhoMin(rho, kappa))
  } else if (!is.na(r) && r==0) {
    return(1/computePearsonBCor(n=n, "r"=r, kappa=kappa, methodNumber=2)[["less"]][["bf"]]
           *hJeffreysApprox(n=n, "r"=r, rho)*priorRhoMin(rho, kappa))
  }
}

posteriorRhoFisherApprox <- function(bfObject, rho, alternative="two.sided") {
  sidedObject <- bfObject[[alternative]]
  n <- sidedObject[["n"]]
  r <- sidedObject[["stat"]]

  if (alternative=="two.sided") {
    return(approximatePosteriorRho("rho"=rho, "n"=n, "r"=r))
  } else if (alternative=="greater") {
    return(approximatePosteriorRhoPlus("rho"=rho, "n"=n, "r"=r))
  } else if (alternative=="less") {
    return(approximatePosteriorRhoMin("rho"=rho, "n"=n, "r"=r))
  }
}

approximatePosteriorRho <- function(rho, n, r) {
  if (n <= 3) {
    std <- 1
  } else {
    std <- 1 / sqrt(n - 3)
  }
  return(1/(1-rho^2)*stats::dnorm(atanh(rho), mean=atanh(r), sd=std))
}

approximatePosteriorRhoPlus <- function(rho, n, r) {
  if (n <= 3) {
    std <- 1
  } else {
    std <- 1 / sqrt(n - 3)
  }
  result <- (approximatePosteriorRho(rho, n, r) * (rho > 0)) /
    (stats::pnorm(0, "mean"=atanh(r), "sd"=std, "lower.tail"=FALSE))
  return(result)
}

approximatePosteriorRhoMin <- function(rho, n, r) {
  if (n <= 3) {
    std <- 1
  } else {
    std <- 1 / sqrt(n - 3)
  }
  result <- (approximatePosteriorRho(rho, n, r) * (rho<0)) /
    (stats::pnorm(0, "mean"=atanh(r), "sd"=std))
  return(result)
}


# 4.2
posteriorMean <- function(n, r, kappa=1, old2F1=FALSE, oneThreshold=1e-3) {
  # NEW CODE CAN OFFICIALLY DO posteriorMean(1219, 0.83)
  #
  checkR <- (1 - abs(r) < oneThreshold)

  if (checkR) {
    if (kappa < 2/(n-2)) {
      hypRatio <- exp(
        lgamma(1/kappa+1+n/2)+lgamma(1/kappa+1-n/2) - 2*lgamma(1/kappa+1) -
        (lgamma(1/kappa+n/2) + lgamma(1/kappa-n/2+1) - 2*lgamma(1/kappa+1/2))
      )
    } else {
      # Note(Alexander): Hack due to consistency and posterior mean -> mle, because n grows or even r to 1
      return(r)
    }
  } else {
    logHyperTerm1 <- tryOrFailWithNA(
      # Re(hypergeo::f15.3.3("A"=n/2, "B"=n/2, "C"=(2+(n+2)*kappa)/(2*kappa), "z"=r^2))
      # Re(hypergeo::f15.3.3("A"=n/2, "B"=n/2, "C"=1/kappa+1+n/2, "z"=r^2))
      log(Re(hypergeo::f15.1.1("A"=1+1/kappa, "B"=1+1/kappa, "C"=1/kappa+1+n/2, "z"=r^2)))
    )
    logHyperTerm2 <- tryOrFailWithNA(
      # Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=(2+n*kappa)/(2*kappa), "z"=r^2))
      # Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=1/kappa+n/2, "z"=r^2))
      log(Re(hypergeo::f15.1.1("A"=1/kappa+1/2, "B"=1/kappa+1/2, "C"=1/kappa+n/2, "z"=r^2)))
    )
    hypRatio <- tryOrFailWithNA(exp(logHyperTerm1-logHyperTerm2))
  }

  # Note(Alexander): that this is a bit of a hack here, as but correct due to large samples.
  #
  if (is.na(hypRatio)) {
    return(r)
  }

  logW <-  2*(lgamma(n/2)-lgamma((n-1)/2))
  result <- (2*kappa*r)/(2+n*kappa)*exp(logW)*hypRatio
  return(result)

  # Old code can officially do posteriorMean(1216, 0.83)
  #
  # hyperTerm1 <- tryOrFailWithNA(
  #   Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=c((2+(n+2)*kappa)/(2*kappa)), z=r^2))
  # )
  # hyperTerm2 <- tryOrFailWithNA(
  #   Re(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=c((2+n*kappa)/(2*kappa)), z=r^2))
  # )
}

posteriorSecondMoment <- function(n, r, kappa=1, oneThreshold=1e-3) {
  # New code can do:posteriorSecondMoment(1219, 0.83) n=3 more than old code
  #
  #
  checkR <- (1 - abs(r) < oneThreshold)

  if (checkR) {
    if (kappa < 2/(n-2)) {
      logHypTerm1a <- lgamma(1/kappa+1+n/2) + lgamma(1/kappa - n/2 + 2) - 2*lgamma(1/kappa + 3/2)
      logHypTerm1b <- lgamma(n/2+1/kappa+2) + lgamma(1/kappa - n/2 + 1) - 2*lgamma(1/kappa + 3/2)
      logHypTerm2 <- lgamma(1/kappa+n/2) + lgamma(1/kappa-n/2+1) - 2*lgamma(1/kappa+1/2)

      hypRatioA <- tryOrFailWithNA(exp(logHypTerm1a-logHypTerm2))
      hypRatioB <- tryOrFailWithNA(exp(logHypTerm1b-logHypTerm2))
    } else {
      # Note(Alexander): Quite the hack here. Note that this is the upper bound (Jensen)
      return(r^2)
    }
  } else {
    hyperTerm1a <- tryOrFailWithNA(
      # Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=(2+(n+2)*kappa)/(2*kappa), "z"=r^2))
      # Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=1/kappa+1+n/2, "z"=r^2))
      log(Re(hypergeo::f15.1.1("A"=(1/kappa+3/2), "B"=(1/kappa+3/2), "C"=(n/2+1/kappa+1), "z"=r^2)))
    )
    hyperTerm1b <- tryOrFailWithNA(
      # Re(hypergeo::f15.3.3("A"=(n+1)/2, "B"=(n+1)/2, "C"=(2+(n+2)*kappa)/(2*kappa)+1, "z"=r^2))
      # Re(hypergeo::f15.3.3("A"=(n+1)/2, "B"=(n+1)/2, "C"=(n/2+1/kappa+2), "z"=r^2))
      Re(hypergeo::f15.1.1("A"=(1/kappa+3/2), "B"=(1/kappa+3/2), "C"=n/2+1/kappa+2, "z"=r^2))
    )
    hyperTerm2 <- tryOrFailWithNA(
      # Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=(2+n*kappa)/(2*kappa), "z"=r^2))
      # Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=1/kappa+n/2, "z"=r^2))
      log(Re(hypergeo::f15.1.1("A"=(1/kappa+1/2), "B"=(1/kappa+1/2), "C"=(1/kappa+n/2), "z"=r^2)))
    )

    hypRatioA <- exp(log(1-r^2)+hyperTerm1a-hyperTerm2)
    hypRatioB <- hyperTerm1b/exp(hyperTerm2)
  }

  # TODO(Alexander): Add asymptotic approximation here
  #
  result <- tryOrFailWithNA(
    kappa/(n*kappa+2) * (hypRatioA+ kappa*(n-1)^(2)/(2+(n+2)*kappa)*r^2*hypRatioB)
  )

  if (is.na(result))
    return(r^2)

  return(result)

  # OLD CODE CAN DO:posteriorSecondMoment(1204, 0.83), at 1205 get inf
  #
  # hyperTerm1 <- tryOrFailWithNA(
  #   Re(hypergeo::genhypergeo(U=c(3/2, (n-1)/2, (n-1)/2),
  #                            L=c(1/2, (2+(n+2)*kappa)/(2*kappa)), z=r^2))
  # )
  # hyperTerm2 <- tryOrFailWithNA(
  #   Re(hypergeo::f15.3.3("A"=(n-1)/2, "B"=(n-1)/2, "C"=(2+n*kappa)/(2*kappa), "z"=r^2))
  # )
  #
  # result <- kappa/(n*kappa+2)*hyperTerm1/hyperTerm2
  # return(result)
}

posteriorVariance <- function(n, r, kappa=1, oneThreshold=1e-3) {
  # Posterior mean of the bf10CorExact
  #	That is, (rho+1)/2, thus, on 0,1 scale to estimate a, b in a beta distribution
  #
  # Note(Alexander): Gauss' Vandermonde identity probably catched this case
  #
  #     add safeguard for large n as then hyperTerm1/hyperTerm2 is almost 1
  # 	  and also for logTerm almost being 1
  #
  # 	posteriorVariance(199, 0.8) yields 6808.702
  #
  #
  result <- tryOrFailWithNA(
    posteriorSecondMoment("n"=n, "r"=r, "kappa"=kappa, "oneThreshold"=oneThreshold) -
      (posteriorMean("n"=n, "r"=r, "kappa"=kappa, "oneThreshold"=oneThreshold))^2)

  # Asymptotic approximation Based on Fisher
  #
  # if (is.na(result))
  #   result <- tanh(1/sqrt(n-3))^2

  return(result)
}

betaParameterEstimates <- function(someMean, someVar) {
  # someMean \in (0, 1)
  # Note(Alexander): Gauss' Vandermonde identity covers the case that someMean = 1
  #
  someA <- tryOrFailWithNA(someMean*(someMean*(1-someMean)/someVar-1))
  someB <- tryOrFailWithNA((1-someMean)*(someMean*(1-someMean)/someVar-1))

  result <- list("betaA"=someA, "betaB"=someB)
  return(result)
}

posteriorBetaParameters <- function(n, r, kappa=1, oneThreshold=1e-3) {
  # posteriorBetaParameters
  # Let rho = 2*x - 1 where x \sim beta, thus, x = (rho+1)/2.Hence, someMu.
  # For the variance we have var(rho)/2^2
  #
  someMu <- tryOrFailWithNA((posteriorMean("n"=n, "r"=r, "kappa"=kappa, "oneThreshold"=oneThreshold)+1)/2, )
  someVar <- tryOrFailWithNA(posteriorVariance("n"=n, "r"=r, "kappa"=kappa, "oneThreshold"=oneThreshold)/4)


  if (isSomeNA(someMu, someVar) || isSomeInfinite(someMu, someVar)) {
    # TODO(Alexander): Before doing this try the MH sampler
    return(list(betaA=NA, betaB=NA))
  } else {
    return(betaParameterEstimates(someMu, someVar))
  }
}

computePearsonCredibleInterval <- function(betaA, betaB, ciValue, h0=0) {
  # Compute Pearson's correlation credible interval based on a beta fit
  #
  check <- failIfNot(betaA > 0, betaB > 0, ciValue > 0, ciValue < 1, !isSomeNull(betaA, betaB))

  failedResult <- list("lowerCi"=NA, "upperCi"=NA, "posteriorMedian"=NA, "ciValue"=ciValue)
  result <- list("two.sided"=failedResult, "greater"=failedResult, "less"=failedResult,
                 "ciValue"=ciValue, "betaA"=betaA, "betaB"=betaB)

  if (!is.null(check)) {
    result[["betaA"]] <- NA
    result[["betaB"]] <- NA
    result[["error"]] <- check
    return(result)
  }

  typeOne <- 1-ciValue
  excessLevel <- typeOne/2

  if (isSomeInfinite(betaA, betaB)) {
    result[["betaA"]] <- NA
    result[["betaB"]] <- NA
    result[["error"]] <- "Can't compute credible intervals"
    return(result)
  } else {
    # Note: Zero one refers to the problem on the (0, 1) rather than on (-1, 1)
    lowerCIZeroOne <- tryOrFailWithNA(qbeta(excessLevel, betaA, betaB))
    medianCIZeroOne <-tryOrFailWithNA(qbeta(1/2, betaA, betaB))
    upperCIZeroOne <- tryOrFailWithNA(qbeta(1-excessLevel, betaA, betaB))

    if (isSomeNA(lowerCIZeroOne, medianCIZeroOne, upperCIZeroOne)) {
      return(result)
    } else {
      # Note: This is simply an application of the definition of the stretched beta
      lowerCi <- 2*lowerCIZeroOne-1
      posteriorMedian <- 2*medianCIZeroOne-1
      upperCi <- 2*upperCIZeroOne-1
    }
  }


  tempList <- list("lowerCi"=lowerCi, "upperCi"=upperCi, "posteriorMedian"=posteriorMedian)

  if (abs(upperCi-lowerCi) <= .Machine$double.eps) {
    tempList[["tooPeaked"]] <- TRUE
  }

  result[["two.sided"]] <- modifyList(result[["two.sided"]], tempList)

  # One sided:
  tempCi <- computePearsonMinSidedCredibleInterval("betaA"=betaA, "betaB"=betaB, "ciValue"=ciValue)
  tempList <- list("lowerCi"=tempCi[1], "upperCi"=tempCi[3], "posteriorMedian"=tempCi[2])

  if (abs(tempCi[3]-tempCi[1]) <= .Machine$double.eps) {
    tempList[["tooPeaked"]] <- TRUE
  }
  result[["less"]] <- modifyList(result[["less"]], tempList)

  # The problem is symmetric
  tempCi  <- computePearsonMinSidedCredibleInterval("betaA"=betaB, "betaB"=betaA, "ciValue"=ciValue)
  tempList <- list("lowerCi"=-tempCi[3], "upperCi"=-tempCi[1], "posteriorMedian"=-tempCi[2])

  if (abs(tempCi[3]-tempCi[1]) <= .Machine$double.eps) {
    tempList[["tooPeaked"]] <- TRUE
  }
  result[["greater"]] <- modifyList(result[["greater"]], tempList)

  return(result)
}

computePearsonMinSidedCredibleInterval <- function(betaA, betaB, ciValue) {
  # Compute min sided Pearson's correlation credible interval based on a beta fit
  #
  result <- NA
  typeOne <- 1-ciValue
  excessLevel <- typeOne/2

  if (any(is.na(c(betaA, betaB)), is.infinite(c(betaA, betaB)))) {
    return(result)
  } else {
    leftArea <- pbeta(1/2, betaA, betaB)
    lowerCIZeroOne <- tryOrFailWithNA(qbeta(excessLevel*leftArea, betaA, betaB))
    medianCIZeroOne <- tryOrFailWithNA(qbeta(leftArea/2, betaA, betaB))
    upperCIZeroOne <- tryOrFailWithNA(qbeta((1-excessLevel)*leftArea, betaA, betaB))

    if (isSomeNA(lowerCIZeroOne, medianCIZeroOne, upperCIZeroOne)) {
      return(result)
    } else {
      lowerCI <- 2*lowerCIZeroOne-1
      medianCI <- 2*medianCIZeroOne-1
      upperCI <- 2*upperCIZeroOne-1
    }
  }
  result <- c(lowerCI, medianCI, upperCI)
  return(result)
}

makeKappas <- function(n) {
  someKappas <- sin(seq(1.5*pi, 2*pi, length=n))+1
  someKappas[1] <- someKappas[2]/10
  someKappas[n] <- 1
  someKappas <- 2*someKappas

  return(someKappas)
}

# 5. Replication TODO(Alexander) Needs revising-------
# These are the functions used to compute replication Bayes factors
#
bfCorrieRepJosine <- function(nOri, rOri, nRep, rRep, kappa=1, hyperGeoOverFlowThreshold=25) {
  result <- list(combined=list(bf10=NA, bfPlus0=NA, bfMin0=NA))

  methodNumber <- 1
  while (methodNumber <= 4 && any(is.na(c(result$combined$bf10,
                                          result$combined$bfPlus0,
                                          result$combined$bfMin0)),
                                  is.infinite(result$combined$bf10))) {
    result <- bfCorrieRepJosineKernel(nOri=nOri, rOri=rOri, nRep=nRep, rRep=rRep, kappa=kappa, methodNumber=methodNumber, hyperGeoOverFlowThreshold=hyperGeoOverFlowThreshold)
    methodNumber <- methodNumber+1
  }

  result[["call"]] <-
    paste0("bfCorrieRepJosine(nOri=", nOri, ", rOri=", rOri, ", nRep=", nRep, ", rRep=", rRep, ", kappa=", kappa, ", hyperGeoOverFlowThreshold=", hyperGeoOverFlowThreshold, ")")

  return(result)
}

bfCorrieRepJosineKernel <- function(nOri, rOri, nRep, rRep, kappa=1, methodNumber=1, hyperGeoOverFlowThreshold=25, h0=0) {
  #
  #  Ly, A., Etz, A., Marsman, M., & Wagenmakers, E.--J. (2017) Replication Bayes factors. Manuscript in preparation
  #  Ly, A., Marsman, M., & Wagenmakers, E.-J. (2017) Analytic Posteriors for Pearson’s Correlation Coefficient. Under review
  #  Wagenmakers, E.-J., Verhagen, A. J., & Ly, A. (2016). How to quantify the evidence for the absence of a correlation. Behavior Research Methods, 48, 413-426.
  #
  # Replication BF for the correlation
  #
  # 1:2 are based on the exact reduced likelihood functions
  # 3:4 are based on the beta approximations to the reduced likelihood functions
  #
  #	methodNumber=1: Use exact likelihood Ly, Marsman, Wagenmakers (2017)
  #	methodNumber=2: Use semi-exact result, based on approximation of the likelihood JeffreysExact, see Wagenmakers et al (2015) bathing
  #	methodNumber=3: Savage Dickey beta approximation
  #	methodNumber=4: Marsman's IMH sampler and then Savage Dickey beta approximation
  #
  # Output is a list of the
  #   - original data,
  #   - rep data,
  #   - combined inference
  #   - replication BFs given the original
  #
  #


  # TODO: avoid when pass through object
  oriObj <- computePearsonBCor(n=nOri, r=rOri, h0=h0, method=methodNumber, kappa=kappa)

  # Default is "NA" list
  result <- list(ori=oriObj, rep=list(NULL),
                 combined=list(n=c(nOri, nRep), r=c(rOri, rRep),
                               repMethodNumber=methodNumber,
                               bf10=NA, bfPlus0=NA, bfMin0=NA,
                               betaA=NA, betaB=NA) ,
                 repGivenOri=list(n=c(nOri, nRep), r=c(rOri, rRep),
                                  bf10=NA, bfPlus0=NA, bfMin0=NA),
                 repMethodNumber=methodNumber)

  if (is.infinite(oriObj[["two.sided"]][["bf"]])) {
    # No use, too big too great, it's true
    #
    return(result)
  }

  # Calculate beta fits of the combined likelihood
  # TODO(Alexander): The fact that kappa=1 sohuldn't be necessary
  if (kappa==1) {
    #
    # methods 3 and 4 are highly dependent on the beta fits based on kappa = 1
    if (methodNumber %in% 3:4 && any(is.na(c(oriObj[["betaA"]], oriObj[["betaB"]])))) {
      # Total failure, real sad
      return(result)
    }

    repObj <- computePearsonBCor(n=nRep, r=rRep, method=methodNumber, kappa=kappa)
    result[["rep"]] <- repObj

    if (methodNumber %in% 3:4 && any(is.na(c(repObj[["betaA"]], repObj[["betaB"]])))) {
      # Failed
      return(result)
    }

    result[["combined"]][["betaA"]] <- oriObj[["betaA"]]-1+repObj[["betaA"]]
    result[["combined"]][["betaB"]] <- oriObj[["betaB"]]-1+repObj[["betaB"]]
  } else {
    # kappa \neq 1

    if (methodNumber %in% 1:3) {
      oriLikelihoodFit <- posteriorBetaParameters(n=nOri, r=rOri, kappa=1, expand=FALSE)
      repLikelihoodFit <- posteriorBetaParameters(n=nRep, r=rRep, kappa=1, expand=FALSE)
    }

    if (methodNumber==4) {
      oriLikelihoodFit <- marsmanCorMHSampler(n=nOri, r=rOri, kappa=1)

      if (is.na(oriLikelihoodFit[["betaA"]]) || is.na(oriLikelihoodFit[["betaB"]])) {
        # Total failure, it's sad
        #
        return(result)
      }
      repLikelihoodFit <- marsmanCorMHSampler(n=nRep, r=rRep, kappa=1)
    }

    if (methodNumber %in% 3:4) {
      if (any(is.na(c(oriLikelihoodFit[["betaA"]], oriLikelihoodFit[["betaB"]],
                      repLikelihoodFit[["betaA"]], repLikelihoodFit[["betaB"]])))) {
        # Failure
        return(result)
      }
    }
    # combine here
    result[["combined"]][["betaA"]] <- oriLikelihoodFit[["betaA"]]-1+repLikelihoodFit[["betaA"]]-1+1/kappa
    result[["combined"]][["betaB"]] <- oriLikelihoodFit[["betaB"]]-1+repLikelihoodFit[["betaB"]]-1+1/kappa


    # Here kappa not 1, but still can see what the original default bfs will do for the rep data
    repObj <- computePearsonBCor(n=nRep, r=rRep, method=methodNumber, kappa=kappa)
    result[["rep"]] <- repObj
  }

  if (methodNumber=="exact" || methodNumber==1) {
    twoSidedIntegrand <- function(x){hFunctionCombinedTwoSided(nOri=nOri, rOri=rOri, nRep=nRep, rRep=rRep, x)*priorRho(x, kappa=kappa)}
    plusSidedIntegrand <- function(x){hFunctionCombined(nOri=nOri, rOri=rOri, nRep=nRep, rRep=rRep, x)*priorRhoPlus(x, kappa=kappa)}
    minSidedIntegrand <- function(x){hFunctionCombined(nOri=nOri, rOri=rOri, nRep=nRep, rRep=rRep, x)*priorRhoMin(x, kappa=kappa)}
  } else if (methodNumber=="jeffreysIntegrate" || methodNumber==2) {
    twoSidedIntegrand <- function(x){hJeffreysApprox(nRep, rRep, x)*hJeffreysApprox(nOri, rOri, x)*priorRho(x, kappa=kappa)}
    plusSidedIntegrand <- function(x){hJeffreysApprox(nRep, rRep, x)*hJeffreysApprox(nOri, rOri, x)*priorRhoPlus(x, kappa=kappa)}
    minSidedIntegrand <- function(x){hJeffreysApprox(nRep, rRep, x)*hJeffreysApprox(nOri, rOri, x)*priorRhoMin(x, kappa=kappa)}
  }

  if (methodNumber %in% 1:2) {
    bf10Combined <- tryOrFailWithNA(integrate(twoSidedIntegrand, -1, 1)[["value"]])

    if (is.na(bf10Combined)) {
      # So sad combined bf10 not available
      result[["combined"]][["bf10"]] <- NA
      return(result)
    }

    if (is.infinite(bf10Combined)) {
      # So big, totally infinite
      #
      result$combined$bf10 <- Inf
      result$repGivenOri$bf10 <- Inf

      # TODO(Aelxander): rRep?
      if (r >= 0) {
        result$combined$bfPlus0 <- Inf
        result$combined$bfMin0 <- 0

        result$repGivenOri$bfPlus0 <- Inf
        result$repGivenOri$bfMin0 <- 0
      } else if (r < 0) {
        result$combined$bfPlus0 <- 0
        result$combined$bfMin0 <- Inf

        result$repGivenOri$bfPlus0 <- 0
        result$repGivenOri$bfMin0 <- Inf
      }
      return(result)
    }

    if (is.finite(bf10Combined)) {
      # Total winner, real great, it's the best

      result[["combined"]][["two.sided"]][["bf"]] <- bf10Combined
      result[["repGivenOri"]][["two.sided"]][["bf"]] <- bf10Combined/oriObj[["two.sided"]][["bf"]]

      if (log(bf10Combined) > hyperGeoOverFlowThreshold) {
        # So big like my hands, can't handle it need to adjust
        # TODO(Alexander): recompute one sided Savage-Dickey
        #
        tempList <- computeBCorOneSidedSavageDickey("bf10"=bf10Combined, "betaA"=result$combined$betaA,
                                                    "betaB"=result$combined$betaB, "kappa"=kappa)

        result[["combined"]][["greater"]][["bf"]] <- tempList[["greater"]][["bf"]]
        result[["combined"]][["less"]][["bf"]] <- tempList[["less"]][["bf"]]
      } else {
        # No overflow, thus, try numerically integrate
        #
        bfPlus0Combined <- tryOrFailWithNA(integrate(plusSidedIntegrand, 0, 1)[["value"]])
        bfMin0Combined <- tryOrFailWithNA(integrate(minSidedIntegrand, -1, 0)[["value"]])

        if (isSomeNA(bfPlus0Combined, bfMin0Combined)) {
          # One sided failed
          return(result)
        }

        if ( bfPlus0Combined < 0 || bfMin0Combined < 0) {
          # One sided failed
          return(result)
        }

        if (is.infinite(bfPlus0Combined) || is.infinite(bfMin0Combined) ||
            (bfPlus0Combined > 1 && bfMin0Combined > 1) ||
            (bfPlus0Combined < 1 && bfMin0Combined < 1) ) {
          tempList <- computeBCorOneSidedSavageDickey(bf10=bf10Combined, betaA=result[["combined"]][["betaA"]],
                                                      betaB=result[["combined"]][["betaB"]], kappa=kappa)

          result[["combined"]][["greater"]][["bf"]] <- tempList[["greater"]][["bf"]]
          result[["combined"]][["less"]][["bf"]] <- tempList[["less"]][["bf"]]
        } else {
          # All good, store numerically calculated one-sided bfs

          result[["combined"]][["greater"]][["bf"]] <- bfPlus0Combined
          result[["combined"]][["less"]][["bf"]] <- bfMin0Combined
        }
      }
    }
  }


  if (methodNumber %in% 3:4) {
    # TODO(AlexandeR):
    if (!is.na(result$combined$betaA) && !is.na(result$combined$betaB)) {
      # Use beta fit and Savage-Dickey
      browser()
      tempList <- computeBCorOneSidedSavageDickey(betaA=result$combined$betaA, betaB=result$combined$betaB, kappa=kappa)
      # tempList <- computeCorBf10SavageDickey(betaA=result$combined$betaA, betaB=result$combined$betaB, kappa=kappa,
      #                                         methodNumber=methodNumber)
      result[["combined"]][["two.sided"]][["bf"]] <- tempList$two.sided$bf
      result[["combined"]][["greater"]][["bf"]] <- tempList$greater$bf
      result[["combined"]][["less"]][["bf"]] <- tempList$less$bf
    }
  }

  # TODO(Alexander): checks for bf10Combined, bfPlus0Combined, bfMin0Combined for zeroes and infinities
  # result[["repGivenOri"]][["two.sided"]][["bf"]] <-


  result[["repGivenOri"]][["two.sided"]][["bf"]] <- (result[["combined"]][["two.sided"]][["bf"]] ) / (oriObj[["two.sided"]][["bf"]])
  result[["repGivenOri"]][["greater"]][["bf"]] <- (result[["combined"]][["greater"]][["bf"]] ) / (oriObj[["greater"]][["bf"]])
  result[["repGivenOri"]][["less"]][["bf"]] <- (result[["combined"]][["less"]][["bf"]] ) / (oriObj[["less"]][["bf"]])

  return(result)
}


