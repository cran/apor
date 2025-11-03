#' @title Normalized Ordinal Prediction Agreement (NOPA)
#'
#' @description Compute the Normalized Ordinal Prediction Agreement (NOPA) metric,
#' a performance measure for models with ordinal-scaled response variables
#' that output estimated probability distributions (EPDs) instead of predicted labels.
#'
#' This function assesses the predictive quality of a model for an ordinal response
#' by aggregating the predicted probability mass as a function of the level of
#' disagreement with respect to the observed category. It provides a normalized
#' and interpretable score between 0 and 1, where 1 indicates perfect agreement
#' and 0 represents the worst possible prediction.
#'
#' NOPA compares the estimated probability distribution produced by a model for
#' each unit of analysis against the observed ordinal response of the same unit. The
#' maximum disagreement is \eqn{k-1}, where \eqn{k} is the number of ordinal categories
#' of the response variable, and the
#' minimum disagreement is 0. Then, aggregates the disagreements of all units of analysis into one single measure.
#'
#' The function internally computes:
#' \itemize{
#'   \item \code{OPD} — Ordinal Prediction Disagreement, the average level of disagreement
#'     between the predicted and observed categories.
#'   \item \code{w} — The worst possible OPD given the dataset, representing the
#'     maximum disagreement achievable.
#'   \item \code{NOPA} — The normalized agreement metric defined as \eqn{1 - OPD / w}.
#'   \item \code{OPDempDist, OPDur, NOPAempDist, NOPAur}: Reference values for
#'     empirical and uniform-random baselines to contextualize model performance assessment
#'     provided by OPD and NOPA.
#' }
#' @param predMat A numeric matrix with \eqn{k} columns and \eqn{n} rows, where \eqn{k} is the number of
#' ordinal categories and \eqn{n} is the number of units of analysis. Each row must be the estimated probability
#' distribution for the unit of analysis to respond each one of the \eqn{k} categories.
#' @param obsVect A numeric or integer vector of observed categories, with values from 1 to \eqn{k},
#'   where \eqn{k} is the number of categories of the ordinal response variable (matching the
#'   number of columns in \code{predMat}).
#' @return A list containing:
#' \describe{
#'   \item{\code{predMat}}{Input matrix of predicted probabilities.}
#'   \item{\code{obsVect}}{Input vector of observed categories.}
#'   \item{\code{disagreementsObs}}{A matrix with \eqn{k} columns (number of ordinal categories
#' of the response variable), and \eqn{n} rows. Each row shows the level of disagreement
#' of each ordinal category with respect to the observed one for the same unit of analysis.}
#'   \item{\code{rearrangedProbObs}}{Matrix of probabilities aggregated by level of disagreement.}
#'   \item{\code{meanDistObs}}{Mean aggregated disagreement profile.}
#'   \item{\code{OPD}}{Observed Ordinal Prediction Disagreement.}
#'   \item{\code{w}}{OPD for the worst prediction possible (maximum disagreement).}
#'   \item{\code{NOPA}}{Normalized Ordinal Prediction Agreement (main metric).}
#'   \item{\code{OPDempDist}}{A version of a reference point for OPD. It considers an
#' ordinal prediction disagreement measure for the case where the estimated probability distribution
#' for the \eqn{k} categories of the ordinal response follows the same distribution as the empirical one.}
#'   \item{\code{OPDur}}{A version of a reference point for OPD. It considers an
#' ordinal prediction disagreement measure for the case where the observed response
#' variable has its own empirical distribution and the estimated probability distribution
#' for the \eqn{k} categories of the ordinal response follows a uniform distribution.}
#'   \item{\code{NOPAempDist}}{A version of a reference point for NOPA. It considers a
#' normalized ordinal prediction agreement measure for the case where the estimated probability distribution
#' for the \eqn{k} categories of the ordinal response follows the same distribution as the empirical one.}
#'   \item{\code{NOPAur}}{A version of a reference point for NOPA. It considers a
#' normalized ordinal prediction agreement measure for the case where the estimated probability distribution
#' for the \eqn{k} categories of the ordinal response follows a uniform distribution.}
#' }
#' @examples
#' EPD <- t(apply(matrix(runif(100),ncol=5),1,function(y) y/sum(y)))
#' sum(rowSums(EPD))==nrow(EPD)
#' ordResponse <- sample(1:5,20, replace=TRUE)
#' nopa(predMat=EPD,obsVect=ordResponse)
#' @seealso \code{\link[apor:ordPredArgmax]{ordPredArgmax}},
#' \code{\link[apor:ordPredRandom]{ordPredRandom}}
#' \code{\link[apor:opdRef]{opdRef}}
#' @references Javier
#' @export
#' @name nopa

nopa <- function(predMat,obsVect) {

  ### Number of categories
  nk <- dim(predMat)[2]

  ### Validations
  # Validate object types
  if (!is.matrix(predMat))
    stop("'predMat' must be a numeric matrix (n x k) of probabilities.")

  if (!is.numeric(predMat))
    stop("'predMat' must contain numerical values only and between 0 and 1.")

  if (!is.numeric(obsVect) && !is.integer(obsVect))
    stop("'obsVect' must be a numeric or integer vector of observed categories
    with values from 1 to ", nk, ".")

  # Validate dimensions
  if (nrow(predMat) != length(obsVect))
    stop("The number of rows in 'predMat' must match the length of 'obsVect'.")

  # Validate category range
  if (any(obsVect < 1 | obsVect > nk))
    stop("Values in 'obsVect' must be between 1 and ", nk, ".")

  # Validate that rows in predMat approximately sum to 1
  row_sums <- rowSums(predMat)
  if (any(abs(row_sums - 1) > 1e-6)) {
    warning("Some rows in 'predMat' do not sum exactly to 1; they will be normalized.")
    predMat <- predMat / row_sums
  }

  # Validate missing or invalid values
  if (any(is.na(predMat)) || any(is.na(obsVect)))
    stop("Missing (NA) values detected in 'predMat' or 'obsVect'.")

  if (any(is.nan(predMat)))
    stop("NaN values detected in 'predMat'.")

  ### Number of observed vectors
  nvect <- dim(predMat)[1]
  ### Number of categories away from the observed one
  disagreementsObs <- abs(matrix(rep(c(1:nk),nvect), nrow=nvect, byrow = TRUE)-
                            matrix(rep(obsVect,nk), ncol=nk, byrow = FALSE))
  ### Rearranged probabilities for each vector according to their level of disagreement
  ### with respect to its corresponding observed category
  rearrangedProbObs <- matrix(0,nrow=nvect,ncol=nk)
  colnames(rearrangedProbObs) <- paste(rep("d",nk),1:nk-1, sep="")
  for (i in 1:nvect) {
    for (j in 1:nk) {
      rearrangedProbObs[i,j] <- sum((disagreementsObs[i,]==j-1)*predMat[i,])
    }
  }
  meanDistObs <- colSums(rearrangedProbObs)/nvect
  ### Ordinal Prediction Disagreement based on the observed category and
  ### the predicted category with the highest probability
  opdObs <- sum(meanDistObs*(1:nk-1))
  ### Ordinal Prediction Disagreement based on the observed category and
  ### the predicted category with the highest disagreement (worst prediction)
  opdWP  <- sum((table(apply(disagreementsObs,1,max))/nvect)*
                  as.numeric(names(table(apply(disagreementsObs,1,max)))))
  NOPA  <- 1-opdObs/opdWP

  p <- prop.table(table(obsVect))
  resul_opdRef <- opdRef(p=p)

  return(list("predMat"=predMat, "obsVect"=obsVect, "disagreementsObs"=disagreementsObs,
              "rearrangedProbObs"=rearrangedProbObs, "meanDistObs"=meanDistObs,
              "OPD"=opdObs, "w"=opdWP, "NOPA"=NOPA,
              "OPDempDist" = resul_opdRef$OPDempDist,
              "OPDur" = resul_opdRef$OPDur,
              "NOPAempDist" = 1-resul_opdRef$OPDempDist/opdWP,
              "NOPAur" = 1-resul_opdRef$OPDur/opdWP))
}


