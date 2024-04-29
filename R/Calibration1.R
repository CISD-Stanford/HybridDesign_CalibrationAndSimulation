#'
#' Calibration method 1 (Normal Approximation)
#'
#' @description  This calibration method estimate the variance of the overall test statistics based on several parameters
#'
#' @usage Calibration1(x1_mean, x2_mean, x1_var, x2_var, cov, delta, alpha_p = 0.05, w = \eqn{\frac{cov}{x2_var}}, \eqn{W = x1_mean - w*ifelse(abs(x2_mean) < delta, x2_mean, 0)})
#'
#' @param x1_mean the sample mean of the treatment effect estimation (mean difference between the treatment arm and control arm)
#' @param x2_mean the sample mean that evaluate the heterogeneity between the control arm and external control arm
#' @param x1_var the estimated variance of x1
#' @param x2_var the estimated variance of x2
#' @param cov the covariance between x1 and x2, which is the variance of the endpoints of the control arm
#' @param delta the equivalence margin, which is the maximum difference that can be treated as clinically acceptable
#' @param alpha the significance of the default hypothesis test, which is H0 : μt = μc vs H1 : μt != μc; The default value is 0.05
#' @param w the borrowing weight; we choose the optimal w that minimizes the variance as the default value, which equals \eqn{\frac{cov}{x2_var}}
#' @param W the capital W represents the overall test statistics, which includes the external control data if x2_mean falls within the range between -delta and delta
#'
#' @details \code{Calibration1()} estimates the variance of the overall test statistics W based on several parameters and conclude a decision regarding to reject the null hypothesis or not.
#'  The default hypothesis test is a two-sided test with null hypothesis and alternative hypothesis equaling to H0 : μt = μc vs H1 : μt != μc
#'
#' @return \code{Calibration1()} returns (1) W_var: the variance of overall test statistics; (2) rejectNull: whether we should reject the null hypothesis
#'
#' @seealso toBeAdded
#'
#' @export
#'

Calibration1 <- function (x1_mean, x2_mean, x1_var, x2_var, cov, delta, alpha = 0.05, w = cov/x2_var, W = x1_mean - w*ifelse(abs(x2_mean) < delta, x2_mean, 0)) {

  force(w)
  force(W)

  # E[I(abs(Z2) < theta) * Z2]

  E1 = integrate(function(z) {z*dnorm(z, x2_mean, sqrt(x2_var))}, -delta, delta)$value

  # E[I(abs(Z2) < theta) * Z2^2]

  E2 = integrate(function(z) {z^2*dnorm(z, x2_mean, sqrt(x2_var))}, -delta, delta)$value

  W_var = x1_var + w^2*(E2 - E1^2) - 2*w*(E2 - x2_mean*E1)*cov/x2_var

  rejectNull = W/sqrt(W_var) > qnorm (1 - alpha/2)

  return(list(W_var = W_var, rejectNull = rejectNull))
}
