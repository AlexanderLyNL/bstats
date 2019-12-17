# 0. Helpers ------------
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


# TODO(Alexander): Ask Don to catch message of stopifnot
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
#'
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




#' Checks every element using the provided function
#'
#' @param ... objects that need testing for try error
#'
#' @return Returns FALSE if there's a single element that does not check out, and TRUE if all elements check out
#' @export
#'
#' @examples
isEvery <- function(..., func) {
  # TODO: Make these to return at first find
  obj <- list(...)
  return(purrr::every(obj, func))
}

isSomeNA <- function(...) {
  return(isSome(..., "func"=anyNA, recursive=TRUE))
}


isSome <- function(..., func) {
  # TODO: Make these to return at first find
  obj <- list(...)
  return( purrr::some(obj, func) )
}

#' Checks for try errors.
#'
#' @param ... objects that need testing for try error
#'
#' @return Returns TRUE whenever there's a single try-error, FALSE otherwise
#' @export
#'
#' @examples
#' kaas <- try(integrate(exp, -Inf, Inf))
#' isTryError(kaas)
#'
isTryError <- function(...) {
  return(isSome(..., func=function(x){inherits(x, "try-error")}))
}


#' Check for any NULL
#'
#' @param ... objects that need testing for try error
#'
#' @return Returns TRUE if there's a single NULL, returns FALSE if no NULL
#' @export
#'
#' @examples
isSomeNull <- function(...) {
  return(isSome(..., func=is.null))
}


#' Check for whether all are numeric
#'
#' @param ... objects that need testing for being numeric
#'
#' @return Returns TRUE if all objects are numeric, returns FALSE otherwise
#' @export
#'
#' @examples
isEveryNumeric <- function(...) {
  # TODO: Make these to return at first find
  return(isEvery(..., func=is.numeric))
}


#' Check for any NA
#'
#' @param ... objects that need testing for NA
#'
#' @return Returns TRUE if any is NA, or its decendants are NA, returns FALSE otherwise
#' @export
#'
#' @examples
isSomeNA <- function(...) {
  # TODO: Make these to return at first find
  return(isSome(..., func=is.na))
}

isEveryFinite <- function(...) {
  isEvery(..., func=is.finite)
}

isSomeInfinite <- function(...) {
  isSome(..., func=is.infinite)
}

isSomeTrue <- function(...) {
  isSome(..., func=isTRUE)
}

