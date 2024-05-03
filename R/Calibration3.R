#'
#' Calibration method 3 (adjusted cutoff values considering borrowing and non-borrowing cases)
#'
#' @description  This calibration method recalculated an adjusted nominal testing level and split it to borrowing and non-borrowing cases. Two separate
#' cutoff values will be returned in this case. We determine whether to reject the null hypothesis of the primary test based on each case and its
#' corresponding critical value.
#'
#' @usage Calibration3(x1_mean, x2_mean, x1_var, x2_var, cov, delta, p1, alpha_p = 0.05, w = cov/x2_var, W = x1_mean - w*ifelse(abs(x2_mean) < delta, x2_mean, 0))
#'
#' @param x1_mean the sample mean of the treatment effect estimation (mean difference between the treatment arm and control arm)
#' @param x2_mean the sample mean that evaluate the heterogeneity between the control arm and external control arm
#' @param x1_var the estimated variance of x1
#' @param x2_var the estimated variance of x2
#' @param cov the covariance between x1 and x2, which is the variance of the endpoints of the control arm
#' @param delta the equivalence margin, which is the maximum difference that can be treated as clinically acceptable
#' @param p1 the proportion of type I error spitted to borrowing case; p1 can be an array to test multiple split proportion simultaneously
#' @param alpha the significance of the default hypothesis test, which is H0 : μt = μc vs H1 : μt != μc; The default value is 0.05
#' @param w the borrowing weight; we choose the optimal w that minimizes the variance as the default value, which equals cov/x2_var
#' @param W the capital W represents the overall test statistics, which includes the external control data if x2_mean falls within the range between -delta and delta
#'
#' @details \code{Calibration3()} calculate two separate cutoff values based on the conclusion of equivalence test and conclude a decision regarding to reject the
#' null hypothesis or not considering borrowing and non-borrowing cases. The default hypothesis test is a two-sided test with null hypothesis and
#' alternative hypothesis equaling to H0 : μt = μc vs H1 : μt != μc
#'
#' @return \code{Calibration3()} returns (1) cutoffValue1: the cutoff value for the case when the equivalence between the current control arm and external control arm is established (borrowing);
#'  (2) cutoffValue2: the cutoff value for the case when the equivalence between the current control arm and external control arm is not established (no borrowing);
#'  (3) rejectNull: whether we should reject the null hypothesis of the primary test
#'
#' @seealso toBeAdded
#'
#' @export
#'

Calibration3 <- function(x1_mean, x2_mean, x1_var, x2_var, w, delta, p1, alpha = 0.05, cov = w*x2_var, W = x1_mean - w*ifelse(abs(x2_mean) < delta, x2_mean, 0)){

  library("mvtnorm")

  findAlpha = function (cutoffValue, w, x1_var, x2_var, cov, delta) {

    corr.Z1Z2 = cov/sqrt(x1_var)/sqrt(x2_var)

    # Equivalence test H0 : μh − μc ≥ δ or μh − μc ≤ −δ vs H1 : −δ < μh − μc < δ
    # Primary test: H0 : μt - μc >= 0 vs H1 : μt - μc < 0

    # function Z3: the probability of rejecting the null of equivalence test (borrow) and rejecting the null of primary test

    rejectingBorrowing <- pmvnorm(lower = c(cutoffValue, -delta/sqrt(x2_var)), upper = c(Inf, delta/sqrt(x2_var)), corr = matrix(c(1, 0, 0, 1), nrow = 2)) +
      pmvnorm(lower = c(-Inf, -delta/sqrt(x2_var)), upper = c(-cutoffValue, delta/sqrt(x2_var)), corr = matrix(c(1, 0, 0, 1), nrow = 2))

    # function Z1: the probability of accepting the null of equivalence test (do not borrow) and rejecting the null of primary test

    rejectingNonborrowing <- pmvnorm(lower = c(cutoffValue, -Inf), upper = c(Inf, -delta/sqrt(x2_var)), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2)) +
      pmvnorm(lower = c(cutoffValue, delta/sqrt(x2_var)), upper = c(Inf, Inf), mean = c(0, 0), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2)) +
      pmvnorm(lower = c(-Inf, -Inf), upper = c(-cutoffValue, -delta/sqrt(x2_var)), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2)) +
      pmvnorm(lower = c(-Inf, delta/sqrt(x2_var)), upper = c(-cutoffValue, Inf), mean = c(0, 0), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2))

    return(c(rejectingNonborrowing + rejectingBorrowing, rejectingBorrowing, rejectingNonborrowing))
  }

  maxAlpha_borrowing = findAlpha(0, w, x1_var, x2_var, cov, delta)[2]
  maxAlpha_nonborrowing = findAlpha(0, w, x1_var, x2_var, cov, delta)[3]

  typeIerror_borrowing = alpha*p1
  typeIerror_nonborrowing = alpha*(1-p1)

  if (sum(typeIerror_nonborrowing > maxAlpha_nonborrowing)>0) {
    cutoffValue_nonborrowing = sapply(1:length(p1), function(ii) {ifelse((typeIerror_nonborrowing > maxAlpha_nonborrowing)[ii] == TRUE, 0, uniroot(function(cutoffValue){
      findAlpha(cutoffValue, w, x1_var, x2_var, cov, delta)[3] - typeIerror_nonborrowing[ii]
    },lower = 0, upper = 10, tol = 1e-8, maxiter = 1e4)$root)})
    typeIerror_left = sapply(typeIerror_nonborrowing - maxAlpha_nonborrowing, function(ii){ifelse(ii > 0, ii, 0)})

    typeIerror_borrowing = typeIerror_borrowing + typeIerror_left
    cutoffValue_borrowing = sapply(typeIerror_borrowing, function(ii) {
      uniroot(function(cutoffValue){
        findAlpha(cutoffValue, w, x1_var, x2_var, cov, delta)[2] - ii
      },lower = 0, upper = 10, tol = 1e-8, maxiter = 1e4)$root
    })
  }
  else if (sum(typeIerror_borrowing > maxAlpha_borrowing)>0) {
    cutoffValue_borrowing = sapply(1:length(p1), function(ii) {ifelse((typeIerror_borrowing > maxAlpha_borrowing)[ii] == TRUE, 0, uniroot(function(cutoffValue){
      findAlpha(cutoffValue, w, x1_var, x2_var, cov, delta)[2] - typeIerror_borrowing[ii]
    },lower = 0, upper = 10, tol = 1e-8, maxiter = 1e4)$root)})
    typeIerror_left = sapply(typeIerror_borrowing - maxAlpha_borrowing, function(ii){ifelse(ii > 0, ii, 0)})

    typeIerror_nonborrowing = typeIerror_nonborrowing + typeIerror_left
    cutoffValue_nonborrowing = sapply(typeIerror_nonborrowing, function(ii) {
      uniroot(function(cutoffValue){
        findAlpha(cutoffValue, w, x1_var, x2_var, cov, delta)[3] - ii
      },lower = 0, upper = 10, tol = 1e-8, maxiter = 1e4)$root
    })
  }
  else {
    cutoffValue_borrowing = sapply (typeIerror_borrowing, function(ii) {
      uniroot(function(cutoffValue){
        findAlpha(cutoffValue, w, x1_var, x2_var, cov, delta)[2] - ii
      },lower = 0, upper = 10, tol = 1e-8, maxiter = 1e4)$root
    })

    cutoffValue_nonborrowing = sapply (typeIerror_nonborrowing, function(ii) {
      uniroot(function(cutoffValue){
        findAlpha(cutoffValue, w, x1_var, x2_var, cov, delta)[3] - ii
      },lower = 0, upper = 10, tol = 1e-8, maxiter = 1e4)$root
    })
  }

  x3_var = x1_var + w^2*x2_var - 2*w*cov

  if(abs(x2_mean) < delta) {
    rejectNull = abs(W/sqrt(x3_var)) > cutoffValue_borrowing
  }
  else {
    rejectNull = abs(W/sqrt(x1_var)) > cutoffValue_nonborrowing
  }

  return(list(cutoffValue1 = cutoffValue_borrowing,
              cutoffValue2 = cutoffValue_nonborrowing,
              rejectNull = rejectNull))
}
