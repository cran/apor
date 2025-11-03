#' @title Argmax Mapping from an Estimated Probability Distribution (EPD) to a Predicted Class
#'
#' @description Deterministically maps each row of an estimated probability distribution (EPD)
#' matrix to a single predicted class by taking the index of the maximum
#' probability. Rows are normalized to sum to one (within tolerance). Ties can
#' be broken by first, last, or at random among maximizers.
#'
#' @param P A numeric matrix of size \eqn{n \times k}, where each row contains
#'   the estimated probabilities \eqn{\hat\pi_{ij}} for subject \eqn{i} and
#'   classes \eqn{j = 1,\ldots,k}. Values must be nonnegative; rows are
#'   normalized to sum to one if needed.
#' @param tie_break Character string indicating how to break ties among
#'   equal maxima. One of \code{"first"} (default), \code{"last"},
#'   or \code{"random"}.
#' @param tol Numeric tolerance used for (i) row-sum checks and (ii) equality
#'   when identifying ties among maximum probabilities. Defaults to \code{1e-12}.
#'
#' @return An integer vector of length \eqn{n} with the predicted class indices
#'   in \eqn{\{1,\ldots,k\}} for each row of \code{P}.
#'
#' @details
#' The function normalizes each row of \code{P} to sum to one (within
#' \code{tol}). Rows with (near) zero total probability trigger an error.
#' \cr
#' If multiple classes achieve the same (within \code{tol}) maximum probability,
#' the returned class depends on \code{tie_break}:
#' \itemize{
#'   \item \code{"first"} — smallest index among maximizers (default).
#'   \item \code{"last"} — largest index among maximizers.
#'   \item \code{"random"} — one index sampled uniformly from the set of
#'         maximizers.
#' }
#'
#' @importFrom stats runif
#' @seealso \code{\link[apor:nopa]{nopa}},
#' \code{\link[apor:ordPredRandom]{ordPredRandom}},
#' \code{\link[apor:opdRef]{opdRef}}
#'
#' @examples
#' P <- rbind(
#'   c(0.05, 0.10, 0.25, 0.60),
#'   c(0.40, 0.40, 0.10, 0.10), # tie between classes 1 and 2
#'   c(NA,   0.20, 0.80, 0.00)  # NA treated as 0
#' )
#'
#' @export
#' @name ordPredArgmax

ordPredArgmax <- function(P, tie_break = c("first", "random", "last"), tol = 1e-12) {
  tie_break <- match.arg(tie_break)
  if (!is.matrix(P)) P <- as.matrix(P)
  if (!is.numeric(P)) stop("P must be a numeric matrix.")
  P[is.na(P)] <- 0

  rs <- rowSums(P)
  if (any(rs <= tol)) stop("Some rows have (near) zero total probability.")
  P <- P / rs

  n <- nrow(P)
  yhat <- integer(n)

  for (i in seq_len(n)) {
    row <- P[i, ]
    m <- max(row)
    idx <- which(abs(row - m) <= tol)

    if (length(idx) == 1L || tie_break == "first") {
      yhat[i] <- idx[1L]
    } else if (tie_break == "last") {
      yhat[i] <- idx[length(idx)]
    } else { # tie_break == "random"
      yhat[i] <- sample(idx, 1L)
    }
  }
  yhat
}

