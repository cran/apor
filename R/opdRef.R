#' @title OPD Reference Points: Empirical vs Uniform Baselines
#'
#' @description Computes two reference values for the Ordinal Prediction Disagreement (OPD):
#' (i) the expected OPD when the predicted label \eqn{\hat Y} follows the *same*
#' empirical distribution as \eqn{Y}; and (ii) the expected OPD when
#' \eqn{\hat Y} is *uniform* over the \eqn{k} ordered categories while \eqn{Y}
#' retains its empirical distribution. These values are useful as dataset-specific
#' anchors for interpreting raw OPD and for constructing normalized benchmarks.
#'
#' @param p A probability vector of length \eqn{k} giving the empirical
#'   distribution of the observed ordinal outcome \eqn{Y\in\{1,\dots,k\}}.
#'   Each entry must be nonnegative and the entries must sum to 1.
#'
#' @details
#' Let \eqn{p=(p_1,\dots,p_k)} denote the empirical distribution of \eqn{Y}.
#' The function returns two scalars:
#' \itemize{
#'   \item \code{OPDempDist}: \eqn{\mathbb{E}|\,\hat Y-Y\,|} when
#'         \eqn{\hat Y\sim p} independently of \eqn{Y\sim p}.
#'   \item \code{OPDur}: \eqn{\mathbb{E}|\,\hat Y-Y\,|} when
#'         \eqn{\hat Y\sim \mathrm{Unif}\{1,\dots,k\}} independently of
#'         \eqn{Y\sim p}.
#' }
#' Both are computed via the disagreement-level decomposition
#' \deqn{\mathbb{E}|\,\hat Y-Y\,|
#'       = \sum_{d=0}^{k-1} d \;\mathbb{P}(|\hat Y-Y|=d),}
#' where, for the uniform case,
#' \deqn{\mathrm{OPD}_{UR}=\frac{1}{k}\sum_{d=0}^{k-1}
#'       d\Big[\mathbb{P}\{Y\le k-d\}-\mathbb{P}\{Y\le d\} + \mathbb{P}\{Y\ge d+1\}\Big],}
#' which is the discrete-\eqn{\{1,\dots,k\}} version of the expression shown in the manuscript.
#'
#' @return A named numeric vector of length two:
#'   \code{c(OPDempDist = ..., OPDur = ...)}.
#'
#' @examples
#' # Example with k = 5 categories and an empirical distribution p:
#' p <- c(0.10, 0.20, 0.40, 0.20, 0.10)
#' opdRef(p)
#'
#' @seealso \code{\link[apor:nopa]{nopa}},
#' \code{\link[apor:ordPredArgmax]{ordPredArgmax}},
#' \code{\link[apor:ordPredRandom]{ordPredRandom}}
#'
#' @export
#' @name opdRef

opdRef <- function(p) {
  # Produce two reference points for OPD

  # Probability vector must sum 1.
  if (sum(p)!=1)
    stop("Probability vector 'p' must sum 1.")

  k <- length(p)

  # Probability vector must have values between 0 and 1.
  if (sum(p>0 & p<1)!=k)
    stop("Probability vector 'p' must sum 1.")


  # 1.- Reference point for \hat{Y} with the same distribution as
  # the empirical one of observed Y.
  phat <- p
  sum_d <- 0
  for (d in 0:(k-1)) {
    sum_j <- 0
    for (j in 1:(k-d)) {
      sum_j <- sum_j  + (p[d+j]*phat[j] + phat[d+j]*p[j])
    }
    sum_d <- sum_d + d*sum_j
  }
  sum_dSameDist <- sum_d

  # 2.- Reference point for \hat{Y} uniformly distributed
  # and Y obs with its own empirical distribution.
  phat <- rep(1/k,k)
  sum_d <- 0
  for (d in 0:(k-1)) {
    sum_j <- 0
    for (j in 1:(k-d)) {
      sum_j <- sum_j  + (p[d+j]*phat[j] + phat[d+j]*p[j])
    }
    sum_d <- sum_d + d*sum_j
  }
  sum_dYhatRandUnif <- sum_d

  return(list("OPDempDist" = sum_dSameDist,"OPDur" = sum_dYhatRandUnif))
}


