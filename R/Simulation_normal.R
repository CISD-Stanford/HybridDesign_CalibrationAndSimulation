#'
#' Run simulation study under normal assumption
#'
#' @description  This function is used to run simulation studies under normal assumption: the endpoints of the treatment arm, control arm and external control arm follow normal distributions 
#' with means equal to mu_T, mu_C, mu_R and standard deviations equal to sigma_T, sigma_C and sigma_R. This function returns a list that contain the simulated data, the proportion of rejecting
#' the null hypothesis of each calibration method.
#'
#' @usage Simulation.normal(mu_T, mu_C, mu_R, sigma_T, sigma_C, sigma_R, alpha_p, delta, alpha_eq, n_T, n_C, n_R, p1, nsim, seed = NULL)
#'
#' @param mu_T the true mean of the treatment arm endpoint
#' @param mu_C the true mean of the control arm endpoint
#' @param mu_R the true mean of the external control arm endpoint
#' @param sigma_T the true standard deviation of the treatment arm endpoint
#' @param sigma_C the true standard deviation of the control arm endpoint
#' @param sigma_R the true standard deviation of the external control arm endpoint
#' @param alpha_p the significance level of the primary hypothesis test, which is H0 : μt = μc vs H1 : μt != μc
#' @param delta the equivalence margin, which is the maximum difference that can be treated as clinically acceptable
#' @param alpha_eq the significance level of the equivalence test, which is H0 : μr - μc >= δeq or μr - μc <= δeq vs H1 : δeq < μr - μc < δeq
#' @param n_T the sample size of the treatment arm
#' @param n_C the sample size of the control arm
#' @param n_R the sample size of the external control arm
#' @param nsim the total number of trials to be simulated
#' @param p1 the proportion of adjusted type I error spitted to the borrowing case, this parameter is only applied in the third calibration method; For details, please refer to the description of function 'calibration3' for details
#' @param seed the random seed for simulation
#' 
#'
#' @details \code{Simulation.normal()} run the simulation study under the normal assumption, which is : the endpoints of the treatment arm, control arm and external control arm follow normal distributions 
#' with means equal to mu_T, mu_C, mu_R and standard deviations equal to sigma_T, sigma_C and sigma_R. This function returns a list that contain the simulated data, the proportion of rejecting
#' the null hypothesis of each calibration method.
#'
#' @return \code{Simulation.normal()} returns a list that contains (1) data: the simulated data; 
#' (2) rejectNull: The proportion of times we reject the null hypothesis when real-world data is not taken into consideration;
#' (3) rejectNull_brw: The proportion of times we reject the null hypothesis when borrowing real-world data only happens when the equivalence test is accepted;
#' (4) rejectNull_c1: The proportion of times we reject the null hypothesis when the first calibration method is applied. Please refer to the description of function 'calibration1' for details;
#' (2) rejectNull_c2: The proportion of times we reject the null hypothesis when the second calibration method is applied. Please refer to the description of function 'calibration2' for details;
#' (2) rejectNull_c3: The proportion of times we reject the null hypothesis when the third calibration method is applied. Please refer to the description of function 'calibration3' for details;
#' (2) rejectNull_pp1: The proportion of times we reject the null hypothesis when borrowing rel-world data by power prior approach;
#'
#' @seealso toBeAdded
#'
#' @export
#'

Simulation = function (mu_T, mu_C, mu_R,
                       sigma_T, sigma_C, sigma_R,
                       alpha_p,
                       delta, alpha_eq,
                       n_T, n_C, n_R, 
                       p1, nsim, seed = NULL) {
  R = c()
  set.seed(seed)
  
  for (i in 1:nsim) {
    
    # Sample treatment, control group and RWD from population
    Treatment  = rnorm(n = n_T, mean = mu_T, sd = sigma_T)
    Control = rnorm(n = n_C, mean = mu_C, sd = sigma_C)
    RWD = rnorm(n = n_R, mean = mu_R, sd = sigma_R)
    
    # Sample mean and standard error
    Y = mean (Treatment)
    Y.se = sd(Treatment)/sqrt(n_T)
    X1 = mean (Control)
    X1.se = sd(Control)/sqrt(n_C)
    X2 = mean (RWD)
    X2.se = sd(RWD)/sqrt(n_R)
    
    Z1 = Y - X1
    Z1_var = Y.se^2 + X1.se^2
    Z2 = X2 - X1
    Z2_var = X2.se^2 + X1.se^2
    
    w = n_R/(n_C+n_R)
    cov = w*sqrt(Z2_var)^2
    theta = delta - qnorm(1-alpha_eq/2)*sqrt(Z2_var)
    borrow = abs (Z2) <= theta
    
    Z3 = Z1 - w*Z2
    Z3_var = Z1_var + w^2*Z2_var - 2*w*cov
    rejectNull   = abs(Z1/sqrt(Z1_var)) > qnorm (1 - alpha_p/2)
    rejectNull_brw = ifelse(borrow, abs(Z3/sqrt(Z3_var)) > qnorm (1 - alpha_p/2), abs(Z1/sqrt(Z1_var)) > qnorm (1 - alpha_p/2))
    
    ############################## Normal approximation (Approach 1) ############################
    
    res1 = Calibration1(Z1, Z2, Z1_var, Z2_var, w, theta, alpha_p)
    rejectNull_c1 = res1$rejectNull
    
    ############################## Common cutoff value (Approach 2) #############################
    
    res2 = Calibration2(Z1, Z2, Z1_var, Z2_var, w, theta, alpha_p)
    rejectNull_c2 = res2$rejectNull
    
    ############################## Split type I error (Approach 3) ##############################
    
    res3 = Calibration3(Z1, Z2, Z1_var, Z2_var, w, theta, p1, alpha_p)
    rejectNull_c3 = res3$rejectNull
    
    ############################## Power prior (Bayesian approach) ##############################
    
    logML <- function(mean1, mean2, sigmaX1, sigmaX2, alpha,sig0,mu0){
      # sig0 is the prior variance and mu0 is the prior of shared mu
      
      # Origin
      signew <- 1/(1/(sigmaX2)^2 * alpha + 1/sig0^2)
      munew <- signew * (1/(sigmaX2)^2 * alpha * mean2 + mu0/sig0^2)
      sigstar <- 1/(1/(sigmaX2)^2 * alpha + 1/sig0^2 + 1/(sigmaX1)^2)
      mustar <- sigstar * (1/(sigmaX2)^2 * alpha * mean2 + mu0/sig0^2 + 1/(sigmaX1)^2 * mean1)
      logllkOrigin <- 0.5 * log(1/signew^2) - 0.5 * munew^2/signew^2 + 0.5 * mustar^2/sigstar^2
      # Deriv
      logllkDeriv <- 1/(1/(sigmaX2)^2 * alpha + 1/sig0^2) -
        (2 * (1/(sigmaX2)^2 * alpha * mean2 + mu0/sig0^2) * mean2 * (1/(sigmaX2)^2 + 1/sig0^2) -
           (1/(sigmaX2)^2 * alpha * mean2 + mu0/sig0^2)^2)/(1/(sigmaX2)^2 * alpha + 1/sig0^2)^2 +
        (2 * (1/(sigmaX2)^2 * alpha * mean2 + mu0/sig0^2 + 1/(sigmaX1)^2 * mean1) * mean2 * (1/(sigmaX2)^2 + 1/sig0^2 + 1/(sigmaX1)^2) -
           (1/(sigmaX2)^2 * alpha * mean2 + mu0/sig0^2 + 1/(sigmaX1)^2 * mean1)^2)/(1/(sigmaX2)^2 * alpha + 1/sig0^2 + 1/(sigmaX1)^2)^2
      logllkDeriv <- logllkDeriv * 1/(sigmaX2)^2
      return(c(logllkOrigin,logllkDeriv))
    }
    
    # sig0 is the prior variance and mu0 is the prior of shared mu
    deriv = logML(X1, X2, X1.se, X2.se, alpha = 1, sig0 = 1e2, mu0 =0)[2]
    
    # Update alpha
    if(deriv > 0){
      alphahat <- 1
    }
    else {
      alphahat <- uniroot(function(alpha){
        logML(X1, X2, X1.se, X2.se, alpha, sig0 = 1e2, mu0 =0)[2]
      },lower = 0,upper = 1,tol = 1e-8, maxiter = 1e4)$root
    }
    
    logML2 = function (alpha, X1, sigma_X1, X2, sigma_X2, mu0, sigma0) {
      numerator = integrate(function(mu) {dnorm(x = mu, mean = X1, sd = sigma_X1)*dnorm(x = mu, mean = X2, sd = sigma_X2)^alpha*dnorm(x = mu, mean = mu0, sd = sigma0)}, -Inf, Inf)$value
      denominator = integrate(function(mu) {dnorm(x = mu, mean = X2, sd = sigma_X2)^alpha*dnorm(x = mu, mean = mu0, sd = sigma0)}, -Inf, Inf)$value
      return(numerator/denominator)
    }
    alphahat2 = optimize(logML2, interval = c(0, 1), X1 = X1, sigma_X1 = X1.se, X2 = X2, sigma_X2 = X2.se,  mu0 = 0, sigma0 = 1e2, maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
    
    # the combined estimator
    w1 = n_R/(n_C + n_R) * alphahat
    w2 = n_R/(n_C + n_R) * alphahat2
    estimator = (1-w1)*X1+w1*X2
    estimator.sd = sqrt((1-w1)^2*X1.se^2 + w1^2*X2.se^2)
    estimator2 = (1-w2)*X1+w2*X2
    estimator.sd2 = sqrt((1-w2)^2*X1.se^2 + w2^2*X2.se^2)
    
    rejectNull_pp1 = abs(Y - estimator)/sqrt(Y.se^2 + estimator.sd^2) > qnorm (1 - alpha_p/2)
    rejectNull_pp2 = abs(Y - estimator2)/sqrt(Y.se^2 + estimator.sd2^2) > qnorm (1 - alpha_p/2)

    R = rbind(R, data.frame(Z1 = Z1,
                            Z2 = Z2,
                            Z3 = Z3,
                            Z1_var = Z1_var,
                            Z2_var = Z2_var,
                            Z3_var = Z3_var,
                            theta = theta,
                            w = w,
                            cov= cov,
                            borrow = borrow,
                            rejectNull = rejectNull,
                            rejectNull_brw = rejectNull_brw,
                            rejectNull_c1 = rejectNull_c1,
                            rejectNull_c2 = rejectNull_c2,
                            rejectNull_c3 = rejectNull_c3,
                            alphahat = alphahat,
                            alphahat2 = alphahat2,
                            rejectNull_pp1 = rejectNull_pp1,
                            rejectNull_pp2 = rejectNull_pp2)
              )
  }
  return (list(data = R,
               rejectNull = mean(R$rejectNull),
               rejectNull_brw = mean(R$rejectNull_brw),
               rejectNull_c1 = mean(R$rejectNull_c1),
               rejectNull_c2 = mean(R$rejectNull_c2),
               rejectNull_c3 = mean(R$rejectNull_c3),
               rejectNull_pp1 = mean(R$rejectNull_pp1),
               rejectNull_pp2 = mean(R$rejectNull_pp2)
               ))
}