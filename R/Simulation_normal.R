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
    Reject   = abs(Z1/sqrt(Z1_var)) > qnorm (1 - alpha_p/2)
    Reject_n = ifelse(borrow, abs(Z3/sqrt(Z3_var)) > qnorm (1 - alpha_p/2), abs(Z1/sqrt(Z1_var)) > qnorm (1 - alpha_p/2))
    
    ############################## Normal approximation (Approach 1) ############################
    
    res1 = Calibration1(Z1, Z2, Z1_var, Z2_var, w, theta, alpha_p)
    Reject_c1 = res1$rejectNull
    
    ############################## Common cutoff value (Approach 2) #############################
    
    res2 = Calibration2(Z1, Z2, Z1_var, Z2_var, w, theta, alpha_p)
    Reject_c2 = res2$rejectNull
    
    ############################## Split type I error (Approach 3) ##############################
    
    res3 = Calibration3(Z1, Z2, Z1_var, Z2_var, w, theta, p1, alpha_p)
    Reject_c3 = res3$rejectNull
    
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
    
    Reject_pp1 = abs(Y - estimator)/sqrt(Y.se^2 + estimator.sd^2) > qnorm (1 - alpha_p/2)
    Reject_pp2 = abs(Y - estimator2)/sqrt(Y.se^2 + estimator.sd2^2) > qnorm (1 - alpha_p/2)

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
                            Reject = Reject,
                            Reject_c1 = Reject_c1,
                            Reject_c2 = Reject_c2,
                            Reject_c3 = Reject_c3,
                            alphahat = alphahat,
                            alphahat2 = alphahat2,
                            Reject_pp1 = Reject_pp1,
                            Reject_pp2 = Reject_pp2)
              )
  }
  return (R)
}