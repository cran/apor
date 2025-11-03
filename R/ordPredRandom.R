#' @title Randomized Mapping from an Estimated Probability Distribution (EPD) to a Predicted Class
#'
#' @description Stochastically maps each row of an estimated probability distribution (EPD)
#' matrix to a single predicted class by drawing one sample from the row's
#' categorical distribution. Rows are normalized to sum to one (within tolerance),
#' and the cut-points method is used with intervals \eqn{(c_{i,j-1}, c_{i,j}]},
#' ensuring \eqn{z_i=1} maps to class \eqn{k}.
#'
#' @param P A numeric matrix of size \eqn{n \times k}, where each row contains
#'   the estimated probabilities \eqn{\hat\pi_{ij}} for subject \eqn{i} and
#'   classes \eqn{j = 1,\ldots,k}. Values must be nonnegative; rows are
#'   normalized to sum to one if needed.
#' @param z Optional numeric vector of length \eqn{n} with values in \eqn{(0,1]}
#'   providing external uniforms for reproducibility or control. If \code{NULL}
#'   (default), draws are generated internally via \code{runif(n)}.
#' @param tol Numeric tolerance used for row-sum checks and for guarding against
#'   underflow when normalizing. Defaults to \code{1e-12}.
#'
#' @return An integer vector of length \eqn{n} with the predicted class indices
#'   in \eqn{\{1,\ldots,k\}} for each row of \code{P}.
#'
#' @details
#' The mapping follows the cumulative cut-points
#' \eqn{c_{i,0}=0}, \eqn{c_{i,j}=\sum_{\ell=1}^j \hat\pi_{i\ell}} for
#' \eqn{j=1,\ldots,k}, and assigns class \eqn{j} whenever
#' \eqn{c_{i,j-1} < z_i \le c_{i,j}}. When \code{z} is supplied, values are
#' clipped to \eqn{(0,1]} to respect interval boundaries. Rows with (near) zero
#' total probability trigger an error.
#'
#' @seealso \code{\link[apor:nopa]{nopa}},
#' \code{\link[apor:ordPredArgmax]{ordPredArgmax}}
#' \code{\link[apor:opdRef]{opdRef}}
#'
#' @examples
#' set.seed(1)
#' P <- rbind(
#'   c(0.05, 0.10, 0.25, 0.60),
#'   c(0.40, 0.40, 0.10, 0.10),
#'   c(0.00, 0.20, 0.80, 0.00)
#' )
#'
#' # Stochastic draws from each row's EPD
#' ordPredRandom(P)
#'
#' # Reproducible draws using provided uniforms
#' z <- c(0.2, 0.85, 1.0)
#' ordPredRandom(P, z = z)
#'
#' @export
#' @name ordPredRandom

ordPredRandom <- function(P, z = NULL, tol = 1e-12) {
  if (!is.matrix(P)) P <- as.matrix(P)
  if (!is.numeric(P)) stop("P must be a numeric matrix.")
  P[is.na(P)] <- 0

  rs <- rowSums(P)
  if (any(rs <= tol)) stop("Some rows have (near) zero total probability.")
  P <- P / rs

  n <- nrow(P)
  k <- ncol(P)

  if (is.null(z)) {
    z <- runif(n)
  } else {
    if (length(z) != n) stop("Length of z must match number of rows in P.")
  }
  # Clip to (0,1] so that z=1 maps to class k and z>0 never maps to 0
  z <- pmin(pmax(z, .Machine$double.eps), 1.0)

  # Row-wise cumulative sums = cut-points c_{i,j}
  C <- t(apply(P, 1L, cumsum))  # n x k

  # For each i, find j with c_{i,j-1} < z_i <= c_{i,j}
  yhat <- mapply(function(ci, zi) {
    # ci is the cumulative vector for one row; prepend 0 as c_{i,0}
    findInterval(zi, vec = c(0, ci), rightmost.closed = TRUE)
  }, split(C, row(C)), z)

  as.integer(yhat)
}


