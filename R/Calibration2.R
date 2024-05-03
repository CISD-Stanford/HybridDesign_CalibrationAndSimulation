#'
#' Calibration method 2 (adjusted nominal testing level)
#'
#' @description  This calibration method recalculated an adjusted nominal testing level adjusted_alpha that makes the actual testing level equals alpha.
#' The adjusted common cutoff value would be (1-adjusted_alpha/2)-th quantile of standard normal distribution
#'
#' @usage Calibration2(x1_mean, x2_mean, x1_var, x2_var, cov, delta, alpha_p = 0.05, w = cov/x2_var, W = x1_mean - w*ifelse(abs(x2_mean) < delta, x2_mean, 0))
#'
#' @param x1_mean the sample mean of the treatment effect estimation (mean difference between the treatment arm and control arm)
#' @param x2_mean the sample mean that evaluate the heterogeneity between the control arm and external control arm
#' @param x1_var the estimated variance of x1
#' @param x2_var the estimated variance of x2
#' @param cov the covariance between x1 and x2, which is the variance of the endpoints of the control arm
#' @param delta the equivalence margin, which is the maximum difference that can be treated as clinically acceptable
#' @param alpha the significance of the default hypothesis test, which is H0 : μt = μc vs H1 : μt != μc; The default value is 0.05
#' @param w the borrowing weight; we choose the optimal w that minimizes the variance as the default value, which equals cov/x2_var
#' @param W the capital W represents the overall test statistics, which includes the external control data if x2_mean falls within the range between -delta and delta
#'
#' @details \code{Calibration2()} calculate the adjusted nominal testing level that makes the actual testing level equals alpha based on several parameters and
#' conclude a decision regarding to reject the null hypothesis or not. The default hypothesis test is a two-sided test with null hypothesis and alternative hypothesis equaling to H0 : μt = μc vs H1 : μt != μc
#'
#' @return \code{Calibration2()} returns (1) cutoffValue: the common cutoff value; (2) rejectNull: whether we should reject the null hypothesis
#'
#' @seealso toBeAdded
#'
#' @export
#'

Calibration2 <- function(x1_mean, x2_mean, x1_var, x2_var, w, delta, alpha = 0.05, cov = w*x2_var, W = x1_mean - w*ifelse(abs(x2_mean) < delta, x2_mean, 0)){

  library("mvtnorm")

  findAlpha = function (cutoffValue, w, x1_var, x2_var, cov, delta) {

    corr.Z1Z2 = cov/sqrt(x1_var)/sqrt(x2_var)

    # Equivalence test H0 : μh − μc ≥ δ or μh − μc ≤ −δ vs H1 : −δ < μh − μc < δ
    # Primary test: H0 : μt - μc >= 0 vs H1 : μt - μc < 0

    # function Z1: the probability of accepting the null of equivalence test (do not borrow) and rejecting the null of primary test

    rejectingNonborrowing <- pmvnorm(lower = c(cutoffValue, -Inf), upper = c(Inf, -delta/sqrt(x2_var)), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2)) +
      pmvnorm(lower = c(cutoffValue, delta/sqrt(x2_var)), upper = c(Inf, Inf), mean = c(0, 0), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2)) +
      pmvnorm(lower = c(-Inf, -Inf), upper = c(-cutoffValue, -delta/sqrt(x2_var)), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2)) +
      pmvnorm(lower = c(-Inf, delta/sqrt(x2_var)), upper = c(-cutoffValue, Inf), mean = c(0, 0), corr = matrix(c(1, corr.Z1Z2, corr.Z1Z2, 1), nrow = 2))


    # function Z3: the probability of rejecting the null of equivalence test (borrow) and rejecting the null of primary test

    rejectingBorrowing <- pmvnorm(lower = c(cutoffValue, -delta/sqrt(x2_var)), upper = c(Inf, delta/sqrt(x2_var)), corr = matrix(c(1, 0, 0, 1), nrow = 2)) +
      pmvnorm(lower = c(-Inf, -delta/sqrt(x2_var)), upper = c(-cutoffValue, delta/sqrt(x2_var)), corr = matrix(c(1, 0, 0, 1), nrow = 2))

    return(c(rejectingNonborrowing + rejectingBorrowing, rejectingNonborrowing, rejectingBorrowing))
  }

  cutoffValue <- uniroot(function(cutoffValue){
    findAlpha(cutoffValue, w, x1_var, x2_var, cov, delta)[1] - alpha
  },lower = 0, upper = 10, tol = 1e-8, maxiter = 1e4)$root

  x3_var = x1_var + w^2*x2_var - 2*w*cov

  rejectNull = abs(W/ifelse(abs(x2_mean) < delta, sqrt(x3_var), sqrt(x1_var))) > cutoffValue

  return(list(cutoffValue = cutoffValue,
              rejectNull = rejectNull))
}
