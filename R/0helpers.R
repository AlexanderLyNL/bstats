#' Try to evaluate an expression, if not fail with NA (default)
#'
#' @param expr Expression to be evaluated
#' @param value Return value if there is an error, default is NA_real_
#'
#' @return Returns the evaluation of the expression, or value if it doesn't work out
#' @export
#'
#' @examples
#' tryOrFailWithNA(integrate(exp, -Inf, Inf)[["value"]], NA)
#' tryOrFailWithNA(integrate(exp, 0, 3)[["value"]], NA)
tryOrFailWithNA <- function(expr, value=NA_real_) {
  tryCatch(
    error=function(cnd) value,
    expr
  )
}


#' Ensure the Truth of R Expressions and returns TRUE if the expressions are not met.
#'
#' This is basically stopifnot{base}, but instead of stopping it returns TRUE. The following descriptions is
#' adapted from stopifnot{base}: If any of the expressions in ... are not all valid, then instead of stopping a TRUE is
#' returned and an error message is printed indicating the first of the elements of ... which were not true.
#'
#' @param ... any number of (logical) R expressions, which should evaluate to TRUE
#'
#' @return Returns TRUE if the provided expressions are not met
#' @export
#'
#' @examples
#' x <- 0
#' msg <- failIfNot(x > 3)
#' print(msg)
failIfNot <- function (...) {
  # This is equivalent to
  #
  # tryCatch(error=function(cnd){
  #   return(list("failed"=TRUE, "error"=conditionMessage(cnd)))
  # },
  # stopifnot(...)
  # )
  #
  result <- NULL

  ll <- list(...)
  n <- length(ll)

  if (n == 0L) {
    return(result)
  }

  Dparse <- function(call, cutoff = 60L) {
    ch <- deparse(call, width.cutoff = cutoff)
    if (length(ch) > 1L) {
      paste(ch[1L], "....")
    } else {
      ch
    }
  }

  head <- function(x, n = 6L) {
    x[seq_len(if (n < 0L) max(length(x) + n, 0L) else min(n, length(x)))]
  }

  abbrev <- function(ae, n = 3L) {
    paste(c(head(ae, n), if (length(ae) > n) "...."), collapse = "\n  ")
  }

  mc <- match.call()

  for (i in 1L:n) {
    if (!(is.logical(r <- ll[[i]]) && !anyNA(r) && all(r))) {
      cl.i <- mc[[i + 1L]]
      msg <- if (is.call(cl.i) && identical(cl.i[[1]], quote(all.equal)) &&
                 (is.null(ni <- names(cl.i)) || length(cl.i) == 3L ||
                  length(cl.i <- cl.i[!nzchar(ni)]) == 3L)) {
        sprintf(gettext("%s and %s are not equal:\n  %s"),
                Dparse(cl.i[[2]]), Dparse(cl.i[[3]]), abbrev(r))
      } else {
        sprintf(ngettext(length(r), "%s is not TRUE", "%s are not all TRUE"),
                Dparse(cl.i))
      }

      result <- msg
      return(result)
    }
  }
  return(result)
}

# Some/any checks -----------
#'Checks whether some object evaluates to TRUE for a provided criterion function
#'
#' @param ... objects that need testing
#' @param func function used to evaluate the truth
#'
#' @return Returns TRUE if there's some object that evaluates to TRUE according to func, and FALSE if all func
#' evaluations lead to FALSE
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- NA
#'
#' isSome(x, y, z, func=is.numeric)
#' isSome(x, y, z, func=is.na)
#' isSome(x, y, z, func=is.vector)
#' isSome(x, y, func=is.na)
isSome <- function(..., func) {
  # TODO: Make these to return at first find
  obj <- list(...)
  return( purrr::some(obj, func) )
}

#'Checks whether some object is NA
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's an NA, FALSE when all objects are not NA.
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- NA
#'
#' isSomeNA(x, y)
#' isSomeNA(x, y, z)
isSomeNA <- function(...) {
  return(isSome(..., "func"=anyNA, recursive=TRUE))
}

#'Checks whether some object is NULL
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's a NULL, FALSE when all objects are not NULL
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- NULL
#'
#' isSomeNull(x, y)
#' isSomeNull(x, y, z)
isSomeNull <- function(...) {
  return(isSome(..., func=is.null))
}


#'Checks whether some object is infinite
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's infinite, FALSE when all objects are finite
#' @export
#'
#'@examples
#'
#' isSomeInfinite(1:10)
#' isSomeInfinite(c(1e2000, 1:1000000))
#'
#' # Faster than
#' # isSomeInfinite(c(1:1000000, 1e2000))
isSomeInfinite <- function(obj) {
  purrr::some(obj, is.infinite)
}

#'Checks whether some object is TRUE
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's TRUE, FALSE when all objects are not TRUE
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- TRUE
#'
#' isSomeTrue(x, y)
#' isSomeTrue(x, y, z)
isSomeTrue <- function(...) {
  isSome(..., func=isTRUE)
}



#' Checks whether some object is a try error.
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if there's some object that's a try-error, FALSE when all objects are not try-errors
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- try(integrate(exp, -Inf, Inf))
#' isTryError(x, y)
#' isTryError(x, y, z)
isTryError <- function(...) {
  return(isSome(..., func=function(x){inherits(x, "try-error")}))
}


# Every/all checks -----------

#' Checks whether all objects evalutes to TRUE for a provided function criterion
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if any objects are a try error, returns FALSE otherwise
#' @export
#'
#' @examples
#' x <- 1
#' y <- "a"
#' z <- NA
#'
#' isEvery(x, y, z, func=is.numeric)
#' isEvery(x, y, z, func=is.na)
#' isEvery(x, y, z, func=is.vector)
isEvery <- function(..., func) {
  # TODO: Make these to return at first find
  obj <- list(...)
  return(purrr::every(obj, func))
}

#' Checks whether all objects are numeric
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if all objects are numeric, returns FALSE otherwise
#' @export
#'
#' @examples
#' x <- 1
#' y <- 2
#' z <- "3"
#'
#' isEveryNumeric(x, y)
#' isEveryNumeric(x, y, z)
isEveryNumeric <- function(...) {
  # TODO: Make these to return at first find
  return(isEvery(..., func=is.numeric))
}

#' Checks whether all objects are finite
#'
#' @inheritParams isSome
#'
#' @return Returns TRUE if all objects are finite, returns FALSE otherwise
#' @export
#'
#' @examples
#' x <- 1
#' y <- 2
#' z <- 10^(1e10)
#'
#' isEveryFinite(x, y)
#' isEveryFinite(x, y, z)
isEveryFinite <- function(...) {
  isEvery(..., func=is.finite)
}
