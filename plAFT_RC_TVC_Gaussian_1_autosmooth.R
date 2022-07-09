# based on "plAFT_RC_TVC_Gaussian.R", with different strategy of calculating the knots and sigmas of Gaussian basis functions
# search "comment" and "under construction" for special comments
plAFT_RC_TVC_Gaussian_1_autosmooth <- function(tmat, Xmat, Zmat, beta_initial, gamma_initial, smooth_initial = 1e-3, autosmooth_no = 5, numIntKnt = 7, maxiter = 10000, threshold = 1e-6, knots_min_quantile = 1, knots_max_quantile = 99, sig_coverage = 0.6, Hes_addition = 0, knots_fix_threshold = 1e-4, knots_fix_iter_no = 10000, frequent_knots_update = TRUE, eta = 3e-1, X_option = 2, knots_option = "percentile", beta_Hessian_option = "modified", gamma_Hessian_option = "modified", draw_plot = TRUE) {
  # tmat includes two columns: 1. termination time,
  #   2. censoring indicator (censor_n=1 for event time, censor_n=0 for right-censored time).
  #
  # Xmat includes columns: 1. individual id, 2-. time-fixed covariate(s) X.
  #
  # Zmat includes columns: 1. individual id, 2. time interval index,
  #   time intervals (two columns, 3. interval_L and 4. interval_R),
  #   5-. time-varying covariate(s) Z.
  # browser()
  # handle missing arguments
  if (missing(Xmat)) {stop("Please include the time-fixed covariates (Xmat) in the input arguments.\n")}
  if (missing(tmat)) {stop("Please specify follow up time and censoring indicator (tmat).\n")}
  if (missing(Zmat) & !missing(gamma_initial)) {
    stop("No need to specify initial value for gamma if time-varying covariates (Zmat) are not included.\n")
  }
  if (!missing(Zmat) & missing(gamma_initial)) {
    stop("Please specify intial value for gamma (gamma_initial).\n")
  }
  
  # library(splines2)
  
  # time[is.infinite(time) | time==0] = NA
  
  # information from tmat matrix
  tmat = as.matrix(tmat)
  rownames(tmat) = NULL
  time = tmat[, 1, drop = FALSE] # n by 1 matrix
  censor_n = tmat[, 2, drop = FALSE] # n by 1 matrix
  
  # information from Xmat matrix
  Xmat = as.matrix(Xmat)
  rownames(Xmat) = NULL
  # id_X = Xmat[,1] # n vector
  X = Xmat[, -1, drop = FALSE] # n by p matrix, 20200330 the negative doesn't work with class as data.frame if drop = FALSE
  n = dim(X)[1]
  p = dim(X)[2] # number of time-fixed covariates
  if (X_option==1) {X = X - matrix(1,n,1) %*% colMeans(X,na.rm = TRUE)} # standardisation by subtracting the mean of each column
  if (X_option==2) {X = X} # without standardisation
  n_E = sum(censor_n)
  # n_RC = n - n_E, n_RC denotes number of right-censored time observations
  if (n_E==0) {X_E = matrix(0,1,p)}
  else {X_E = X[censor_n==1L, , drop = FALSE]} # n_E by p matrix, n_E denotes number of event time observations
  
  # information from Zmat matrix
  if (missing(Zmat) | missing(gamma_initial)) {dim_Zmat = c(0,0)} # 20210126 added to accommodate (missing(Zmat) & missing(gamma_initial))
  else {dim_Zmat = dim(Zmat)} # 20210126 added to accommodate (missing(Zmat) & missing(gamma_initial))
  
  # 20201112 added (dim(Zmat)[1]>(dim(Xmat)[1]+1)) condition to avoid computation issue when none or only one individual took Z treatment
  # if (!missing(Zmat) & !missing(gamma_initial) & (dim(Zmat)[1]>(dim(Xmat)[1]+1))) { #20201112
  if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) { # 20210126 modified to accommodate (missing(Zmat) & missing(gamma_initial))
    Zmat = as.matrix(Zmat)
    rownames(Zmat) = NULL
    id_Z = Zmat[,1] # N vector, N is greater than or equal to n
    # interval_index = Zmat[,2] # N vector
    # ni = tapply(interval_index,id_Z,FUN=which.max) # n vector, change FUN from max to which.max on 20220318, robust
    ni = diff(c(0,which(c(diff(id_Z)!=0,TRUE)))) # n vector, 20220318, robust, not needing interval_index any more
    interval_length = Zmat[, 4, drop = FALSE] - Zmat[, 3, drop = FALSE] # N by 1 matrix
    Z = Zmat[, -c(1:4), drop = FALSE] # N by q matrix
    N = dim(Z)[1]
    q = dim(Z)[2] # number of time-varying covariates
    Z_ni = Z[c(diff(id_Z)!=0,TRUE), , drop = FALSE] # n by q matrix
    if (n_E==0) {Z_ni_E = matrix(0,1,q)}
    else {Z_ni_E = Z_ni[censor_n==1L, , drop = FALSE]} # n_E by q matrix
    # for full Hessian of gamma
    # if (gamma_Hessian_option=="full") {
    ni_E = ni[censor_n==1L]
    censor_N = rep(censor_n, times = ni)
    # interval_L_E = Zmat[censor_N==1L, 3, drop = FALSE] # N_E by 1 matrix
    # interval_R_E = Zmat[censor_N==1L, 4, drop = FALSE] # N_E by 1 matrix
    interval_length_E = Zmat[censor_N==1L, 4, drop = FALSE] - Zmat[censor_N==1L, 3, drop = FALSE] # N_E by 1 matrix
    Z_E = Z[censor_N==1, , drop = FALSE] # N_E by 1 matrix
    # }
  }
  
  
  # Gaussian basis part, necessary functions to calculate the derivative and integral of Gaussian basis
  Gaussian_basis <- function(x, mean, sd) {
    if (nrow(x)==1) {
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}))}
    }
  }
  
  Gaussian_basis_integ1 <- function(q, mean, sd) {
    if (nrow(q)==1) {
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd)) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd[i]) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd[i])) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd)) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {(pnorm(q, mean[i], sd[i]) - matrix(1, NROW(q), 1)%*%pnorm(0, mean[i], sd[i])) * sd[i]*sqrt(2*pi)}))}
    }
  }
  
  Gaussian_basis_deriv1 <- function(x, mean, sd) {
    if (nrow(x)==1) {
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd^2 * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd[i]^2 * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd^2 * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {-(x-mean[i])/sd[i]^2 * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}))}
    }
  }
  
  Gaussian_basis_deriv2 <- function(x, mean, sd) {
    if (nrow(x)==1) {
      if (length(sd)==1) {
        return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd^4 - 1/sd^2) * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}), nrow = 1))
      }
      else {return(matrix(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd[i]^4 - 1/sd[i]^2) * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}), nrow = 1))}
    }
    else {
      if (length(sd)==1) {
        return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd^4 - 1/sd^2) * dnorm(x, mean[i], sd) * sd*sqrt(2*pi)}))
      }
      else {return(sapply(X = 1:(numIntKnt+2), FUN = function(i) {((x-mean[i])^2/sd[i]^4 - 1/sd[i]^2) * dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}))}
    }
  }
  
  # outer loop
  smooth_iter = matrix(NA, autosmooth_no, 1)
  smooth = smooth_initial
  df_old = 0
  for (j in 1:autosmooth_no) {
    cat(j,"th outer loop begins. The current smoothing parameter is",smooth,".\n")
    smooth_iter[j] = smooth
    
    # inner loop
    
    # dgrSp = ordSp - 1
    # numSp = numIntKnt + ordSp
    
    # 0. initial value of beta, gamma and theta
    # library(survival)
    # iniFit = survreg(Surv(time = time[,1], time2 = time[,2], type = "interval2") ~ X, dist="weibull")
    # beta_old = matrix(iniFit$coef[2:(p+1)],ncol = 1)
    # beta_old = matrix(0,p,1) # alternative initial value for beta
    beta_old = matrix(beta_initial) # p by 1 matrix
    # theta_old = matrix(1L,numSp,1) # numSp by 1 matrix
    theta_old = matrix(1L,numIntKnt+2,1) # Gaussian basis part, 1 less knots compare to spline basis
    
    # outputs of this function (faster than the alternative assignment)
    grad_beta_iter = matrix(NA,maxiter,p)
    Hes_beta_iter = matrix(NA,maxiter,p^2)
    det_Hes_beta_iter = matrix(NA,maxiter,1)
    eigenvl_Hes_beta_iter = matrix(NA,maxiter,p)
    likelihood_beta_iter_before_ls = matrix(NA,maxiter,1)
    iter_Newton_beta = matrix(NA,maxiter,1)
    likelihood_beta_iter_after_ls = matrix(NA,maxiter,1)
    ts_range_beta_iter = matrix(NA,maxiter,2) # ts_range after beta update
    beta_iter = matrix(NA,maxiter,p)
    
    # grad_theta_iter = matrix(NA,maxiter,numSp)
    grad_theta_iter = matrix(NA,maxiter,numIntKnt+2) # Gaussian basis part, 1 less theta's compare to spline basis
    penlike_theta_iter_before_ls = matrix(NA,maxiter,1)
    iter_MI_theta = matrix(NA,maxiter,1)
    penlike_theta_iter_after_ls = matrix(NA,maxiter,1)
    ts_range_theta_iter = matrix(NA,maxiter,2) # ts_range after theta update
    # theta_iter = matrix(NA,maxiter,numSp)
    theta_iter = matrix(NA,maxiter,numIntKnt+2)# Gaussian basis part, 1 less theta's compare to spline basis
    # Alternative assignment
    # likelihood_beta_iter_after_ls = NULL ## example: likelihood_beta_iter_after_ls = c(likelihood_beta_iter_after_ls,likelihood_beta_iter_1)
    # iter_Newton_beta = NULL
    # penlike_theta_iter_after_ls = NULL
    # iter_MI_theta = NULL
    # grad_beta_iter = NULL
    # ts_range_iter = NULL
    # beta_iter = NULL
    # theta_iter = NULL
    knots_iter = matrix(NA,maxiter,numIntKnt+2)
    
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]<=n+1)) {
      warning("There is NO time-varying covariates in this dataset!\n")
      # gamma_old = matrix(gamma_initial) # comment this line in simulations considering time-varying covariates
      # gamma_new = gamma_old + 1 # comment this line in simulations considering time-varying covariates
    }
    
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      gamma_old = matrix(gamma_initial) # q by 1
      gamma_new = gamma_old + 1
      
      grad_gamma_iter = matrix(NA,maxiter,q)
      Hes_gamma_iter = matrix(NA,maxiter,q^2)
      det_Hes_gamma_iter = matrix(NA,maxiter,1)
      eigenvl_Hes_gamma_iter = matrix(NA,maxiter,q)
      likelihood_gamma_iter_before_ls = matrix(NA,maxiter,1)
      iter_Newton_gamma = matrix(NA,maxiter,1)
      likelihood_gamma_iter_after_ls = matrix(NA,maxiter,1)
      ts_range_gamma_iter = matrix(NA,maxiter,2) # ts_range after gamma update
      gamma_iter = matrix(NA,maxiter,q)
    }
    else {
      gamma_old = 0
      gamma_new = gamma_old + knots_fix_threshold
    }
    
    # 0. accelerated time
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      ts = exp(-X%*%beta_old) * as.matrix(tapply(exp(-Z%*%gamma_old)*interval_length, id_Z, FUN = sum))
    }
    else {ts = exp(-X%*%beta_old) * time}
    # ts_E = as.matrix(ts[censor_n==1])
    ts_E = ts[censor_n==1, , drop = FALSE]
    
    # 0. knots
    # # knots option 1
    # if (knots_option==1) {
    #   bryKnt = c(min(ts), max(ts)+1e-40)
    #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5,95,length.out = numIntKnt)/100, type=1)
    # }
    # # knots option 2
    # if (knots_option==2) {
    #   Knt = quantile(sort(ts), seq(0,1,length.out = numSp-dgrSp+1), type=1)
    #   bryKnt = c(Knt[1],Knt[numSp-dgrSp+1]+1e-40)
    #   IntKnt = Knt[2:(numSp-dgrSp)]
    # }
    # # knots option 3
    # if (knots_option==3) {
    #   ts_min = min(ts)
    #   ts_max = max(ts)
    #   bryKnt = c(ts_min, ts_max+1e-40)
    #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
    # }
    # bryKnt_initial = bryKnt
    # IntKnt_initial = IntKnt
    
    # Gaussian basis part
    if (knots_option=="equal_space") {
      bryKnt = range(ts)
      bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
      gknots = seq(bryKnt[1], bryKnt[2], bin_width)
      sig = (2/3) * bin_width
      # gauss = pnorm(0,gknots,sig)*sqrt(2*pi*sig^2)
      # gauss = pnorm(bryKnt[1],gknots,sig)*sqrt(2*pi*sig^2) # which one to choose?
    }
    else if (knots_option=="percentile") {
      # different strategy of calculating knots and sigmas
      gknots = quantile(sort(ts), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
      dist = sapply(gknots, function(x) {abs(x - ts)})
      sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
    }
    else {stop("knots_option is either percentile or equal_space.\n")}
    
    
    # # Gaussian basis part
    # # bryKnt = range(ts)
    # # basis_coverage = 0.3
    # ts_sorted = sort(ts)
    # basis_coverage = 1 / (numSp - 1)
    # sig = rep(NA,numSp - 1)
    # gknots = quantile(ts_sorted, seq(5,95,length.out = numSp - 1)/100, type=1)
    # sig[1] = (quantile(ts_sorted,basis_coverage)-gknots[1])/1.96
    # # bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
    # # gknots = seq(bryKnt[1], bryKnt[2], bin_width)
    # # gknots = quantile(sort(ts), seq(0,1,length.out = numSp-dgrSp+1), type=1)
    # # sig = (2/3) * bin_width
    # # gauss = pnorm(0,gknots,sig)*sqrt(2*pi*sig^2)
    # 
    # sum(ts>gknots[1]-1.96*sig[1] & ts<gknots[1]+1.96*sig[1])/n
    # 
    # browser()
    
    
    
    
    # 0. spline and cumulative spline
    # if (n_E==0) {phi_E = matrix(0,1,numSp)}
    # else {
    #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
    # }
    # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
    
    # 0. hazard and cumulative hazard
    # if (n_E==0) {h0ts_E = h0tss_E = 1}
    # else {
    #   h0ts_E = phi_E%*%theta_old
    #   h0tss_E = h0ts_E
    #   h0tss_E[h0tss_E<1e-40] = 1e-40
    # }
    # ch0ts = cphi%*%theta_old
    
    
    
    # Gaussian basis part
    h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
    ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_old
    # Jinqing's version
    # h0ts = matrix(0,n,1)
    # ch0ts = matrix(0,n,1)
    # for (i in 1:n) { # will improve using apply function
    #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
    #   ch0ts[i] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)%*%theta_old
    # }
    if (n_E==0) {h0ts_E = h0tss_E = 1}
    else {
      h0ts_E = h0ts[censor_n==1, , drop = FALSE]
      h0tss_E = h0ts_E
      h0tss_E[h0tss_E<1e-40] = 1e-40
    }
    
    # 0. log-likelihood
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      like_beta_old = sum(log(h0tss_E) - X_E%*%beta_old - Z_ni_E%*%gamma_old) - sum(ch0ts) # 20220109 change h0ts_E to h0tss_E inspired by plAFT_RC_TVC_Gaussian_ma.R
    }
    else {like_beta_old = sum(log(h0tss_E) - X_E%*%beta_old) - sum(ch0ts)}
    
    # browser()
    ####################### main iterative loop #########################
    # tryCatch({ #20181129
    for(k in 1L:maxiter) {
      ########################## Newton step for beta ##########################
      # cat(k,"th iteration.\n")
      # 0.5. differential spline for gradient and modified Hessian
      # if (n_E==0) {dphi_E = matrix(0,1,numSp)}
      # else {dphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)}
      # phi = mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
      
      # 0.5. segmentations for gradient and modified Hessian
      # dh0ts_E = dphi_E%*%theta_old
      # h0ts_E = phi_E%*%theta_old
      # h0tss_E = h0ts_E
      # h0tss_E[h0tss_E<1e-40] = 1e-40
      # h0ts = phi%*%theta_old
      
      # Gaussian basis part
      dh0ts = Gaussian_basis_deriv1(x = ts, mean = gknots, sd = sig)%*%theta_old
      h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
      if (beta_Hessian_option=="full") {
        ddh0ts = Gaussian_basis_deriv2(x = ts, mean = gknots, sd = sig)%*%theta_old
      }
      
      # dh0ts = matrix(0,n,1)
      # h0ts = matrix(0,n,1)
      # for (i in 1:n) { # will improve using apply function
      #   dh0ts[i] = t(dnorm_deriv1(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
      #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
      # }
      if (n_E==0) {
        dh0ts_E = h0ts_E = h0tss_E = 1
        if (beta_Hessian_option=="full") {ddh0ts_E = 1}
      }
      else {
        dh0ts_E = dh0ts[censor_n==1, , drop = FALSE]
        h0ts_E = h0ts[censor_n==1, , drop = FALSE]
        h0tss_E = h0ts_E
        h0tss_E[h0tss_E<1e-40] = 1e-40
        if (beta_Hessian_option=="full") {ddh0ts_E = ddh0ts[censor_n==1, , drop = FALSE]}
      }
      
      # 0.5. gradient
      if (n_E==0) {quan_E_grad_beta = matrix(0,1,1)}
      else {quan_E_grad_beta = dh0ts_E*ts_E/h0tss_E + 1}
      quan_grad_beta = h0ts*ts
      grad_beta = t(-X_E)%*%quan_E_grad_beta - t(-X)%*%quan_grad_beta
      grad_beta_iter[k,] = t(grad_beta)
      
      # # 0.5. differential spline for Hessian
      # ddphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=2)
      # dphi_RC = mSpline(ts_RC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # dphi_LC = mSpline(ts_LC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # dphi_IC_L = mSpline(ts_IC_L, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # dphi_IC_R = mSpline(ts_IC_R, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
      # # 0.5. supplement segmentations for Hessian
      # ddh0ts_E = ddphi_E%*%theta_old
      # dh0ts_RC = dphi_RC%*%theta_old
      # dh0ts_LC = dphi_LC%*%theta_old
      # dh0ts_IC_L = dphi_IC_L%*%theta_old
      # dh0ts_IC_R = dphi_IC_R%*%theta_old
      # 0.5. Hessian
      # quan_E_Hes = (ddh0ts_E*(ts_E)^2*h0ts_E-(dh0ts_E*ts_E)^2)/(h0tss_E)^2 - dh0ts_E*(ts_E)^2 - h0ts_E*ts_E
      # quan_RC_Hes = dh0ts_RC*(ts_RC)^2 - h0ts_RC*ts_RC
      # quan_LC_Hes = exp(-ch0ts_LC)*ts_LC*(-h0ts_LC^2*ts_LC + dh0ts_LC*ts_LC + h0ts_LC)/(1-exp(-ch0ts_LC)) - (exp(-ch0ts_LC)*h0ts_LC*ts_LC/(1-exp(-ch0ts_LC)))^2
      # quan_IC_Hes = (exp(-ch0ts_IC_L)*ts_IC_L*(h0ts_IC_L^2*ts_IC_L - dh0ts_IC_L*ts_IC_L - h0ts_IC_L) + exp(-ch0ts_IC_R)*ts_IC_R*(-h0ts_IC_R^2*ts_IC_R+dh0ts_IC_R*ts_IC_R+h0ts_IC_R))/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)) - ((-exp(-ch0ts_IC_L)*h0ts_IC_L*ts_IC_L+exp(-ch0ts_IC_R)*h0ts_IC_R*ts_IC_R)/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)))^2
      # 
      # Hes = t(X_E)%*%diag(x = as.vector(quan_E_Hes), nrow = n_E, ncol = n_E)%*%X_E + t(X_RC)%*%diag(x = as.vector(quan_RC_Hes), nrow = n_RC, ncol = n_RC)%*%X_RC + t(X_LC)%*%diag(x= as.vector(quan_LC_Hes), nrow = n_LC, ncol = n_LC)%*%X_LC + t(X_IC)%*%diag(x = as.vector(quan_IC_Hes), nrow = n_IC, ncol = n_IC)%*%X_IC
      
      if (beta_Hessian_option=="full") {
        # 0.5. full Hessian
        if (n_E==0) {quan_E_Hes_beta = matrix(0,1,1)}
        else {quan_E_Hes_beta = (ddh0ts_E / h0tss_E - dh0ts_E^2 / h0tss_E^2) * ts_E^2 + dh0ts_E / h0tss_E * ts_E}
        quan_Hes_beta = dh0ts * ts^2 + h0ts * ts
        # Hes_beta = t(-X_E)%*%(quan_E_Hes_beta%*%matrix(1,1,p)*(-X_E)) - t(-X)%*%(quan_Hes_beta%*%matrix(1,1,p)*(-X))
        Hes_beta = t(-X_E)%*%(as.numeric(quan_E_Hes_beta)*(-X_E)) - t(-X)%*%(as.numeric(quan_Hes_beta)*(-X)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
      }
      else if (beta_Hessian_option=="modified") {
        # 0.5. modified Hessian
        if (n_E==0) {quan_E_Hes_beta = matrix(0,1,1)}
        else {quan_E_Hes_beta = (dh0ts_E*ts_E/h0tss_E)^2} # only pick the negative definite terms from the full Hessian and make them positive
        quan_Hes_beta = h0ts*ts # only pick the negative definite terms from the full Hessian and make them positive
        # Hes_beta = t(X_E)%*%diag(x = as.vector(quan_E_Hes_beta), nrow = dim(X_E)[1], ncol = dim(X_E)[1])%*%X_E + t(X)%*%diag(x = as.vector(h0ts*ts), nrow = dim(X)[1], ncol = dim(X)[1])%*%X
        # Hes_beta = t(-X_E)%*%(quan_E_Hes_beta%*%matrix(1,1,p)*(-X_E)) + t(-X)%*%(quan_Hes_beta%*%matrix(1,1,p)*(-X)) # faster alternative, but cause tiny difference in some elements
        Hes_beta = t(-X_E)%*%(as.numeric(quan_E_Hes_beta)*(-X_E)) + t(-X)%*%(as.numeric(quan_Hes_beta)*(-X)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
      }
      else {stop("beta_Hessian_option is either full or modified.\n")}
      
      Hes_beta_iter[k,] = c(t(Hes_beta)) # convert the Hes_beta into vector by row, and save
      det_Hes_beta_iter[k] = det(Hes_beta)
      eigenvl_Hes_beta_iter[k,] = eigen(Hes_beta)$values
      Hes_beta = Hes_beta + diag(Hes_addition,p) # 20190812 to avoid singular in solve()
      # eigen(Hes_beta)
      # if (k==27) {browser()}
      # 0.5. update beta temporarily
      beta_inc = solve(Hes_beta, grad_beta)
      beta_new = beta_old + beta_inc
      
      # 0.5. accelerated time
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        ts = exp(-X%*%beta_new) * as.matrix(tapply(exp(-Z%*%gamma_old)*interval_length, id_Z, FUN = sum))
      }
      else {ts = exp(-X%*%beta_new) * time}
      # ts_E = as.matrix(ts[censor_n==1])
      ts_E = ts[censor_n==1, , drop = FALSE]
      
      # 0.5. knots
      # knots option 1
      # if (knots_option==1) {
      #   bryKnt = c(min(ts), max(ts)+1e-40)
      #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5,95,length.out = numIntKnt)/100, type=1)
      # }
      # # knots option 2
      # if (knots_option==2) {
      #   Knt = quantile(sort(ts), seq(0,1,length.out = numSp-dgrSp+1), type=1)
      #   bryKnt = c(Knt[1],Knt[numSp-dgrSp+1]+1e-40)
      #   IntKnt = Knt[2:(numSp-dgrSp)]
      # }
      # # knots option 3
      # if (knots_option==3) {
      #   ts_min = min(ts)
      #   ts_max = max(ts)
      #   bryKnt = c(ts_min, ts_max+1e-40)
      #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
      # }
      
      # if (all(abs(beta_new-beta_old)<1e-3)) {browser()}
      # Gaussian basis part
      # if (all(abs(beta_new-beta_old)<1e-4) & abs(gamma_new-gamma_old)<1e-4) {browser()}
      if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
        if (knots_option=="equal_space") {
          bryKnt = range(ts)
          bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
          gknots = seq(bryKnt[1], bryKnt[2], bin_width)
          sig = (2/3) * bin_width
        }
        else if (knots_option=="percentile") {
          # different strategy of calculating knots and sigmas
          gknots = quantile(sort(ts), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
          dist = sapply(gknots, function(x) {abs(x - ts)})
          sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
        }
      }
      
      # gauss = pnorm(0,gknots,sig)*sqrt(2*pi*sig^2)
      # gauss = pnorm(bryKnt[1],gknots,sig)*sqrt(2*pi*sig^2) # which one to choose?
      
      # 0.5. spline and cumulative spline
      # if (n_E==0) {phi_E = matrix(0,1,numSp)}
      # else {
      #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
      # }
      # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
      
      # 0.5. hazard and cumulative hazard
      # if (n_E==0) {h0ts_E = h0tss_E = 1}
      # else {
      #   h0ts_E = phi_E%*%theta_old
      #   h0tss_E = h0ts_E
      #   h0tss_E[h0tss_E<1e-40] = 1e-40
      # }
      # ch0ts = cphi%*%theta_old
      
      # Gaussian basis part
      h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
      ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_old
      # h0ts = matrix(0,n,1)
      # ch0ts = matrix(0,n,1)
      # for (i in 1:n) { # will improve using apply function
      #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
      #   ch0ts[i] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)%*%theta_old
      # }
      if (n_E==0) {h0ts_E = h0tss_E = 1}
      else {
        h0ts_E = h0ts[censor_n==1, , drop = FALSE]
        h0tss_E = h0ts_E
        h0tss_E[h0tss_E<1e-40] = 1e-40
      }
      # if(k==27) {browser()}
      # 0.5. log-likelihood
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_old) - sum(ch0ts)
      }
      else {like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new) - sum(ch0ts)}
      
      likelihood_beta_iter_before_ls[k] = like_beta_new
      # 1. backtracking line search for beta
      alpha_N_beta = 1
      iteration_beta = 0L
      # message("for loop No. is ",k)
      while (like_beta_new <= like_beta_old) {
        iteration_beta = iteration_beta + 1L
        # message("like in while loop equals ",like)
        if (alpha_N_beta >= 1e-2) {alpha_N_beta = alpha_N_beta*0.6} # original step size multiplier 0.6, alternative 0.1
        else if (alpha_N_beta < 1e-2 & alpha_N_beta >= 1e-5) {alpha_N_beta = alpha_N_beta*5e-2} # original step size multiplier 5e-2, alternative 1e-3
        else if (alpha_N_beta < 1e-5 & alpha_N_beta >= 1e-20) {alpha_N_beta = alpha_N_beta*1e-5} # original step size multiplier 1e-5, alternative 1e-8
        # else if (alpha_N_beta<1e-20 & alpha_N_beta>=1e-100) {alpha_N_beta=alpha_N_beta*1e-10}
        # if (any(abs(grad_beta)>0.1)) {beta_new = beta_old + 0.01}
        else {break}
        # else if (all(abs(grad_beta)<0.01)) {break}
        
        # 1. update beta with search direction
        beta_new = beta_old + alpha_N_beta * beta_inc # bug of version 0.1.0 fixed
        
        # 1. accelerated time
        if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
          ts = exp(-X%*%beta_new) * as.matrix(tapply(exp(-Z%*%gamma_old)*interval_length, id_Z, FUN = sum))
        }
        else {ts = exp(-X%*%beta_new) * time}
        # ts_E = as.matrix(ts[censor_n==1])
        ts_E = ts[censor_n==1, , drop = FALSE]
        
        # 1. knots
        # # knots option 1
        # if (knots_option==1) {
        #   bryKnt = c(min(ts), max(ts)+1e-40)
        #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5,95,length.out = numIntKnt)/100, type=1)
        # }
        # # knots option 2
        # if (knots_option==2) {
        #   Knt = quantile(sort(ts), seq(0,1,length.out = numSp-dgrSp+1), type=1)
        #   bryKnt = c(Knt[1],Knt[numSp-dgrSp+1]+1e-40)
        #   IntKnt = Knt[2:(numSp-dgrSp)]
        # }
        # # knots option 3
        # if (knots_option==3) {
        #   ts_min = min(ts)
        #   ts_max = max(ts)
        #   bryKnt = c(ts_min, ts_max+1e-40)
        #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
        # }
        
        # Gaussian basis part
        if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
          if (knots_option=="equal_space") {
            bryKnt = range(ts)
            bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
            gknots = seq(bryKnt[1], bryKnt[2], bin_width)
            sig = (2/3) * bin_width
          }
          else if (knots_option=="percentile") {
            # different strategy of calculating knots and sigmas
            gknots = quantile(sort(ts), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
            dist = sapply(gknots, function(x) {abs(x - ts)})
            sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
          }
        }
        
        # gauss = pnorm(0,gknots,sig)*sqrt(2*pi*sig^2)
        # gauss = pnorm(bryKnt[1],gknots,sig)*sqrt(2*pi*sig^2) # which one to choose?
        
        # 1. spline and cumulative spline
        # if (n_E==0) {phi_E = matrix(0,1,numSp)}
        # else {
        #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        # }
        # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        
        # 1. hazard and cumulative hazard
        # if (n_E==0) {h0ts_E = h0tss_E = 1}
        # else {
        #   h0ts_E = phi_E%*%theta_old
        #   h0tss_E = h0ts_E
        #   h0tss_E[h0tss_E<1e-40] = 1e-40
        # }
        # ch0ts = cphi%*%theta_old
        
        # Gaussian basis part
        h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
        ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_old
        # h0ts = matrix(0,n,1)
        # ch0ts = matrix(0,n,1)
        # for (i in 1:n) { # will improve using apply function
        #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
        #   ch0ts[i] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)%*%theta_old
        # }
        if (n_E==0) {h0ts_E = h0tss_E = 1}
        else {
          h0ts_E = h0ts[censor_n==1, , drop = FALSE]
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
        }
        
        # 1. log-likelihood
        if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
          like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_old) - sum(ch0ts)
        }
        else {like_beta_new = sum(log(h0tss_E) - X_E%*%beta_new) - sum(ch0ts)}
      }
      # # one way to escape from local maximum 20210608
      # if (any(abs(grad_beta)>=0.01) & iteration_beta == 7) {beta_new = beta_old - 1e-1*abs(beta_old)*sign(grad_beta)}
      
      ts_range_beta_iter[k,] = range(ts)
      likelihood_beta_iter_after_ls[k] = like_beta_new
      iter_Newton_beta[k] = iteration_beta
      beta_iter[k,] = beta_new
      
      ########################## Newton step for gamma ##########################
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        like_gamma_old = like_beta_new
        # 0.5. differential spline for gradient and modified Hessian
        # if (n_E==0) {dphi_E = matrix(0,1,numSp)}
        # else {dphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)}
        # phi = mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        
        # 0.5. segmentations for gradient and modified Hessian
        # dh0ts_E = dphi_E%*%theta_old
        # h0ts_E = phi_E%*%theta_old
        # h0tss_E = h0ts_E
        # h0tss_E[h0tss_E<1e-40] = 1e-40 # to avoid 0 denominator
        # h0ts = phi%*%theta_old
        
        # Gaussian basis part
        dh0ts = Gaussian_basis_deriv1(x = ts, mean = gknots, sd = sig)%*%theta_old
        h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
        # dh0ts = matrix(0,n,1)
        # h0ts = matrix(0,n,1)
        # for (i in 1:n) { # will improve using apply function
        #   dh0ts[i] = t(dnorm_deriv1(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
        #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
        # }
        # for full Hessian of gamma
        if (gamma_Hessian_option=="full") {ddh0ts = Gaussian_basis_deriv2(x = ts, mean = gknots, sd = sig)%*%theta_old}
        
        if (n_E==0) {
          dh0ts_E = h0ts_E = h0tss_E = 1
          # for full Hessian of gamma
          if (gamma_Hessian_option=="full") {ddh0ts_E = 1}
        }
        else {
          dh0ts_E = dh0ts[censor_n==1, , drop = FALSE]
          h0ts_E = h0ts[censor_n==1, , drop = FALSE]
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
          if (gamma_Hessian_option=="full") {ddh0ts_E = ddh0ts[censor_n==1, , drop = FALSE]} # for full Hessian of gamma}
        }
        
        # 0.5. first derivative of t_star w.r.t. gamma
        if (q==1) {
          # dts_gamma = (exp(-X%*%beta_new)%*%matrix(1,1,q)) * as.matrix(tapply((-Z)*((exp(-Z%*%gamma_old)*interval_length)%*%matrix(1,1,q)), id_Z, FUN = sum)) # n by 1
          dts_gamma = as.numeric(exp(-X%*%beta_new)) * as.matrix(tapply((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length), id_Z, FUN = sum)) # 20220119 faster alternative of above line, use R recycling rule, may cause trouble
          # dts_gamma = as.numeric(exp(-X%*%beta_new)) * apply(as.matrix((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length)), 2, function(x) tapply(x, id_Z, sum)) # 20220119 slower alternative of above line, but works for q > 0
        }
        else {
          # dts_gamma = (exp(-X%*%beta_new)%*%matrix(1,1,q)) * as.matrix(tapply((-Z)*((exp(-Z%*%gamma_old)*interval_length)%*%matrix(1,1,q)), id_Z, FUN = colSums)) # n by q
          dts_gamma = as.numeric(exp(-X%*%beta_new)) * apply((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length), 2, function(x) tapply(x, id_Z, sum))
        }
        dts_E_gamma = dts_gamma[censor_n==1, , drop = FALSE] # n_E by q
        
        # 0.5. gradient
        if (n_E==0) {quan_E_grad_gamma = matrix(0,1,1)}
        # else {quan_E_grad_gamma = (dh0ts_E/h0tss_E)%*%matrix(1,1,q)} # n_E by q
        else {quan_E_grad_gamma = dh0ts_E/h0tss_E} # n_E by q
        # grad_gamma = matrix(colSums(quan_E_grad_gamma*dts_E_gamma - Z_ni_E) - colSums(h0ts%*%matrix(1,1,q)*dts_gamma)) # q by 1
        grad_gamma = matrix(t(dts_E_gamma)%*%quan_E_grad_gamma - colSums(Z_ni_E) - t(dts_gamma)%*%h0ts) # q by 1
        grad_gamma_iter[k,] = t(grad_gamma)
        
        # # 0.5. differential spline for Hessian
        # ddphi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=2)
        # dphi_RC = mSpline(ts_RC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # dphi_LC = mSpline(ts_LC, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # dphi_IC_L = mSpline(ts_IC_L, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # dphi_IC_R = mSpline(ts_IC_R, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=1)
        # # 0.5. supplement segmentations for Hessian
        # ddh0ts_E = ddphi_E%*%theta_old
        # dh0ts_RC = dphi_RC%*%theta_old
        # dh0ts_LC = dphi_LC%*%theta_old
        # dh0ts_IC_L = dphi_IC_L%*%theta_old
        # dh0ts_IC_R = dphi_IC_R%*%theta_old
        # 0.5. Hessian
        # quan_E_Hes = (ddh0ts_E*(ts_E)^2*h0ts_E-(dh0ts_E*ts_E)^2)/(h0tss_E)^2 - dh0ts_E*(ts_E)^2 - h0ts_E*ts_E
        # quan_RC_Hes = dh0ts_RC*(ts_RC)^2 - h0ts_RC*ts_RC
        # quan_LC_Hes = exp(-ch0ts_LC)*ts_LC*(-h0ts_LC^2*ts_LC + dh0ts_LC*ts_LC + h0ts_LC)/(1-exp(-ch0ts_LC)) - (exp(-ch0ts_LC)*h0ts_LC*ts_LC/(1-exp(-ch0ts_LC)))^2
        # quan_IC_Hes = (exp(-ch0ts_IC_L)*ts_IC_L*(h0ts_IC_L^2*ts_IC_L - dh0ts_IC_L*ts_IC_L - h0ts_IC_L) + exp(-ch0ts_IC_R)*ts_IC_R*(-h0ts_IC_R^2*ts_IC_R+dh0ts_IC_R*ts_IC_R+h0ts_IC_R))/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)) - ((-exp(-ch0ts_IC_L)*h0ts_IC_L*ts_IC_L+exp(-ch0ts_IC_R)*h0ts_IC_R*ts_IC_R)/(exp(-ch0ts_IC_L)-exp(-ch0ts_IC_R)))^2
        # 
        # Hes = t(X_E)%*%diag(x = as.vector(quan_E_Hes), nrow = n_E, ncol = n_E)%*%X_E + t(X_RC)%*%diag(x = as.vector(quan_RC_Hes), nrow = n_RC, ncol = n_RC)%*%X_RC + t(X_LC)%*%diag(x= as.vector(quan_LC_Hes), nrow = n_LC, ncol = n_LC)%*%X_LC + t(X_IC)%*%diag(x = as.vector(quan_IC_Hes), nrow = n_IC, ncol = n_IC)%*%X_IC
        
        if (gamma_Hessian_option=="full") {
          # 0.5. full Hessian
          if (n_E==0) {
            quan_E_Hes_gamma_1 = quan_E_Hes_gamma_2 = matrix(0,1,1)
          }
          else {
            quan_E_Hes_gamma_1 = ddh0ts_E / h0tss_E - dh0ts_E^2 / h0tss_E^2 # n_E by 1
            # quan_E_Hes_gamma_2 = rep(dh0ts_E*exp(-X_E%*%beta_new),ni_E)*(exp(-Z_E%*%gamma_old)*interval_length_E)
            quan_E_Hes_gamma_2 = rep(dh0ts_E / h0tss_E * exp(-X_E%*%beta_new),ni_E)*(exp(-Z_E%*%gamma_old)*interval_length_E)
          }
          quan_Hes_gamma_1 = dh0ts
          quan_Hes_gamma_2 = rep(h0ts*exp(-X%*%beta_new),ni)*(exp(-Z%*%gamma_old)*interval_length) # N by 1, rep(c(),c()) suggested by Ali Shariati 20190804
          # Hes_gamma = t(dts_E_gamma)%*%diag(x = as.vector(quan_E_Hes_gamma), nrow = dim(X_E)[1], ncol = dim(X_E)[1])%*%dts_E_gamma + t(Z)%*%diag(x = as.vector(quan_Hes_gamma), nrow = dim(Z)[1], ncol = dim(Z)[1])%*%(Z)
          # Hes_gamma = t(dts_E_gamma)%*%(quan_E_Hes_gamma_1%*%matrix(1,1,q)*dts_E_gamma) +
          #             t(Z_E)%*%(quan_E_Hes_gamma_2%*%matrix(1,1,q)*Z_E) -
          #             t(dts_gamma)%*%(quan_Hes_gamma_1%*%matrix(1,1,q)*dts_gamma) -
          #             t(Z)%*%(quan_Hes_gamma_2%*%matrix(1,1,q)*Z) # faster alternative
          Hes_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hes_gamma_1)*dts_E_gamma) +
                      t(Z_E)%*%(as.numeric(quan_E_Hes_gamma_2)*Z_E) -
                      t(dts_gamma)%*%(as.numeric(quan_Hes_gamma_1)*dts_gamma) -
                      t(Z)%*%(as.numeric(quan_Hes_gamma_2)*Z) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
        }
        else if (gamma_Hessian_option=="modified") {
          # 0.5. modified Hessian
          if (n_E==0) {quan_E_Hes_gamma = matrix(0,1,1)}
          else {quan_E_Hes_gamma = (dh0ts_E/h0tss_E)^2} # n_E by 1
          quan_Hes_gamma = rep(h0ts*exp(-X%*%beta_new),ni)*(exp(-Z%*%gamma_old)*interval_length) # N by 1, rep(c(),c()) suggested by Ali Shariati 20190804
          # Hes_gamma = t(dts_E_gamma)%*%diag(x = as.vector(quan_E_Hes_gamma), nrow = dim(X_E)[1], ncol = dim(X_E)[1])%*%dts_E_gamma + t(Z)%*%diag(x = as.vector(quan_Hes_gamma), nrow = dim(Z)[1], ncol = dim(Z)[1])%*%(Z)
          # Hes_gamma = t(dts_E_gamma)%*%(quan_E_Hes_gamma%*%matrix(1,1,q)*dts_E_gamma) +
          #             t(-Z)%*%(quan_Hes_gamma%*%matrix(1,1,q)*(-Z)) # faster alternative
          Hes_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hes_gamma)*dts_E_gamma) +
                      t(-Z)%*%(as.numeric(quan_Hes_gamma)*(-Z)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
        }
        else {stop("gamma_Hessian_option is either full or modified.\n")}
        Hes_gamma_iter[k,] = c(t(Hes_gamma)) # convert the Hes_gamma into vector by row, and save
        det_Hes_gamma_iter[k] = det(Hes_gamma)
        eigenvl_Hes_gamma_iter[k,] = eigen(Hes_gamma)$values
        Hes_gamma = Hes_gamma + diag(Hes_addition,q) # to avoid singular in solve()
        # eigen(Hes_gamma)
        
        # 0.5. update gamma temporarily
        gamma_inc = solve(Hes_gamma, grad_gamma)
        gamma_new = gamma_old + gamma_inc
        
        # 0.5. accelerated time
        ts = exp(-X%*%beta_new) * as.matrix(tapply(exp(-Z%*%gamma_new)*interval_length, id_Z, FUN = sum))
        # ts_E = as.matrix(ts[censor_n==1])
        ts_E = ts[censor_n==1, , drop = FALSE]
        
        # 0.5. knots
        # # knots option 1
        # if (knots_option==1) {
        #   bryKnt = c(min(ts), max(ts)+1e-40)
        #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5,95,length.out = numIntKnt)/100, type=1)
        # }
        # # knots option 2
        # if (knots_option==2) {
        #   Knt = quantile(sort(ts), seq(0,1,length.out = numSp-dgrSp+1), type=1)
        #   bryKnt = c(Knt[1],Knt[numSp-dgrSp+1]+1e-40)
        #   IntKnt = Knt[2:(numSp-dgrSp)]
        # }
        # # knots option 3
        # if (knots_option==3) {
        #   ts_min = min(ts)
        #   ts_max = max(ts)
        #   bryKnt = c(ts_min, ts_max+1e-40)
        #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
        # }
        
        # Gaussian basis part
        if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
          if (knots_option=="equal_space") {
            bryKnt = range(ts)
            bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
            gknots = seq(bryKnt[1], bryKnt[2], bin_width)
            sig = (2/3) * bin_width
          }
          else if (knots_option=="percentile") {
            # different strategy of calculating knots and sigmas
            gknots = quantile(sort(ts), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
            dist = sapply(gknots, function(x) {abs(x - ts)})
            sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
          }
        }
        
        # gauss = pnorm(0,gknots,sig)*sqrt(2*pi*sig^2)
        # gauss = pnorm(bryKnt[1],gknots,sig)*sqrt(2*pi*sig^2) # which one to choose?
        
        # 0.5. spline and cumulative spline
        # if (n_E==0) {phi_E = matrix(0,1,numSp)}
        # else {
        #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        # }
        # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
        
        # 0.5. hazard and cumulative hazard
        # if (n_E==0) {h0ts_E = h0tss_E = 1}
        # else {
        #   h0ts_E = phi_E%*%theta_old
        #   h0tss_E = h0ts_E
        #   h0tss_E[h0tss_E<1e-40] = 1e-40
        # }
        # ch0ts = cphi%*%theta_old
        
        # Gaussian basis part
        h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
        ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_old
        # h0ts = matrix(0,n,1)
        # ch0ts = matrix(0,n,1)
        # for (i in 1:n) { # will improve using apply function
        #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
        #   ch0ts[i] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)%*%theta_old
        # }
        if (n_E==0) {h0ts_E = h0tss_E = 1}
        else {
          h0ts_E = h0ts[censor_n==1, , drop = FALSE]
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
        }
        
        # 0.5. log-likelihood
        like_gamma_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new) - sum(ch0ts)
        
        likelihood_gamma_iter_before_ls[k] = like_gamma_new
        # 1. backtracking line search for gamma
        alpha_N_gamma = 1
        iteration_gamma = 0L
        # message("for loop No. is ",k)
        while (like_gamma_new <= like_gamma_old) {
          iteration_gamma = iteration_gamma + 1L
          # if(iteration_gamma==7) {browser()}
          # message("like in while loop equals ",like)
          if (alpha_N_gamma >= 1e-2) {alpha_N_gamma = alpha_N_gamma*0.6} # original step size multiplier 0.6
          else if (alpha_N_gamma < 1e-2 & alpha_N_gamma >= 1e-5) {alpha_N_gamma = alpha_N_gamma*5e-2} # original step size multiplier 5e-2
          else if (alpha_N_gamma < 1e-5 & alpha_N_gamma >= 1e-20) {alpha_N_gamma = alpha_N_gamma*1e-5} # original step size multiplier 1e-5
          # else if (alpha_N_gamma<1e-20 & alpha_N_gamma>=1e-100) {alpha_N_beta=alpha_N_gamma*1e-10}
          else {break}
          
          # iteration_gamma==1, alpha_N_gamma == 1e-1
          # iteration_gamma==2, alpha_N_gamma == 1e-2
          # iteration_gamma==3, alpha_N_gamma == 1e-3
          # iteration_gamma==4, alpha_N_gamma == 1e-6
          # iteration_gamma==5, alpha_N_gamma == 1e-14
          # iteration_gamma==6, alpha_N_gamma == 1e-22
          
          # 1. update gamma with search direction
          gamma_new = gamma_old + alpha_N_gamma * gamma_inc # bug of version 0.1.0 fixed
          
          # 1. accelerated time
          ts = exp(-X%*%beta_new) * as.matrix(tapply(exp(-Z%*%gamma_new)*interval_length, id_Z, FUN = sum))
          # ts_E = as.matrix(ts[censor_n==1])
          ts_E = ts[censor_n==1, , drop = FALSE]
          
          # 1. knots
          # # knots option 1
          # if (knots_option==1) {
          #   bryKnt = c(min(ts), max(ts)+1e-40)
          #   IntKnt = quantile(sort(ts)[2:(n-1)], seq(5,95,length.out = numIntKnt)/100, type=1)
          # }
          # # knots option 2
          # if (knots_option==2) {
          #   Knt = quantile(sort(ts), seq(0,1,length.out = numSp-dgrSp+1), type=1)
          #   bryKnt = c(Knt[1],Knt[numSp-dgrSp+1]+1e-40)
          #   IntKnt = Knt[2:(numSp-dgrSp)]
          # }
          # # knots option 3
          # if (knots_option==3) {
          #   ts_min = min(ts)
          #   ts_max = max(ts)
          #   bryKnt = c(ts_min, ts_max+1e-40)
          #   IntKnt = seq(ts_min, ts_max, length.out = numSp-dgrSp+1)[2:(numSp-dgrSp)]
          # }
          
          # Gaussian basis part
          if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==TRUE) { # 20210630 # 20220208 add frequent_knots_update check
            if (knots_option=="equal_space") {
              bryKnt = range(ts)
              bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
              gknots = seq(bryKnt[1], bryKnt[2], bin_width)
              sig = (2/3) * bin_width
            }
            else if (knots_option=="percentile") {
              # different strategy of calculating knots and sigmas
              gknots = quantile(sort(ts), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
              dist = sapply(gknots, function(x) {abs(x - ts)})
              sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
            }
          }
          
          # gauss = pnorm(0,gknots,sig)*sqrt(2*pi*sig^2)
          # gauss = pnorm(bryKnt[1],gknots,sig)*sqrt(2*pi*sig^2) # which one to choose?
          
          # 1. spline and cumulative spline
          # if (n_E==0) {phi_E = matrix(0,1,numSp)}
          # else {
          #   phi_E = mSpline(ts_E, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
          # }
          # cphi = iSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt)
          
          # 1. hazard and cumulative hazard
          # if (n_E==0) {h0ts_E = h0tss_E = 1}
          # else {
          #   h0ts_E = phi_E%*%theta_old
          #   h0tss_E = h0ts_E
          #   h0tss_E[h0tss_E<1e-40] = 1e-40
          # }
          # ch0ts = cphi%*%theta_old
          
          # Gaussian basis part
          h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_old
          ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_old
          # h0ts = matrix(0,n,1)
          # ch0ts = matrix(0,n,1)
          # for (i in 1:n) { # will improve using apply function
          #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_old
          #   ch0ts[i] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)%*%theta_old
          # }
          if (n_E==0) {h0ts_E = h0tss_E = 1}
          else {
            h0ts_E = h0ts[censor_n==1, , drop = FALSE]
            h0tss_E = h0ts_E
            h0tss_E[h0tss_E<1e-40] = 1e-40
          }
          
          # 1. log-likelihood
          like_gamma_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new) - sum(ch0ts)
        }
        
        # # one way to escape from local maximum 20210608
        # if (any(abs(grad_gamma)>=0.01) & iteration_gamma == 7) {gamma_new = gamma_old - 1e-1*abs(gamma_old)*sign(grad_gamma)}
        
        ts_range_gamma_iter[k,] = range(ts)
        likelihood_gamma_iter_after_ls[k] = like_gamma_new
        iter_Newton_gamma[k] = iteration_gamma
        gamma_iter[k,] = gamma_new
        
        like_theta_old = like_gamma_new
      }
      else {like_theta_old = like_beta_new}
      
      
      ################ roughness penalty matrix ##################
      # 1. penalty equals theta_transpose %*% R %*% theta, R equals to phi_sd %*% phi_sd_transpose
      # R=matrix(0, nrow=numSp, ncol=numSp) 
      # xknots = c(rep(bryKnt[1], ordSp), IntKnt, rep(bryKnt[2], ordSp))
      # # browser()
      # for (ii in 1:numSp){
      #   for (jj in ii:numSp){
      #     if (jj - ii<ordSp){
      #       kntset = xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]]
      #       kntsum = 0
      #       for (kk in 1:(length(kntset)-1)){
      #         kntsum = kntsum + mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, Boundary.knots=bryKnt, 
      #                                   derivs=dgrSp)[ii]*mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, 
      #                                                             Boundary.knots=bryKnt,derivs=dgrSp)[jj]*(kntset[kk+1]-kntset[kk])
      #       }
      #       R[ii, jj] = kntsum
      #     }
      #   }
      # }
      # R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
      
      
      # try alternative
      
      # R=matrix(0, nrow=numSp, ncol=numSp) 
      # xknots = c(rep(bryKnt[1], ordSp), IntKnt, rep(bryKnt[2], ordSp))
      # for (ii in 1:numSp){
      #   for (jj in ii:numSp){
      #     if (jj - ii<ordSp){
      #       kntset = xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]]
      #       ddphi_kntset = mSpline(kntset, knots=IntKnt, degree=dgrSp, intercept=T, Boundary.knots=bryKnt,derivs=dgrSp)
      #       kntsum = 0
      #       for (kk in 1:(length(kntset)-1)){
      #         kntsum = kntsum + ddphi_kntset[kk,ii]*ddphi_kntset[kk,jj]*(kntset[kk+1]-kntset[kk])
      #       }
      #       R[ii, jj] = kntsum
      #     }
      #   }
      # }
      # R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
      
      
      if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no) {
        # Ding's optimised version
        num_dots = 5*n
        ts_range = range(ts)
        bin_edges = seq(ts_range[1],ts_range[2],length.out=num_dots+1)
        bin_middle_points = as.matrix((bin_edges[2:(num_dots+1)] + bin_edges[1:num_dots])/2)
        ddphi_bin_middle_points = Gaussian_basis_deriv2(x = bin_middle_points, mean = gknots, sd = sig)
        # R = matrix(0,numIntKnt+2,numIntKnt+2)
        # for (u in 1:(numIntKnt+2)) {
        #   for (r in u:(numIntKnt+2)) {
        #     R[u,r] = (ts_range[2]-ts_range[1])/num_dots * sum(ddphi_bin_middle_points[,u]*ddphi_bin_middle_points[,r])
        #   }
        # }
        # R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
        R = t(ddphi_bin_middle_points)%*%ddphi_bin_middle_points * (ts_range[2]-ts_range[1])/num_dots # 20211102
      }
      
      # ddphi_bind = mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt,derivs = 2)
      # int_ddphi = matrix(NA,length(ts)-1,1)
      # for(interval in 1:(length(ts)-1)) {
      #   int_ddphi[interval] = ddphi_bind[interval,]%*%t(ddphi_bind[interval,])*(ts[interval+1]-ts[interval])
      # }
      # browser()
      # penalty = t(theta_old)%*%(h0ts_bind)%*%theta_old
      
      ########################## MI step for theta ##########################
      # 1. penalised likelihood
      # if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {like_theta_old = like_gamma_new}
      # else {like_theta_old = like_beta_new}
      plike_old = like_theta_old - smooth*t(theta_old)%*%R%*%theta_old
      
      # 1.5. terms of theta update algorithm
      # if (n_E==0) {quan_E_nume = matrix(0,1,numSp)}
      # else {quan_E_nume = phi_E/(h0tss_E%*%matrix(1,1,numSp))}
      # quan_deno = cphi
      
      # Gaussian basis part
      if (n_E==0) {quan_E_nume = matrix(0,1,numIntKnt+2)}
      else {
        # phi = sapply(X = gknots, FUN = dnorm, x = ts, sd = sig)*sqrt(2*pi*sig^2)
        # phi = matrix(0,n,numIntKnt+2)
        # for (i in 1:n) {
        #   phi[i,] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))
        # }
        # phi_E = phi[censor_n==1, , drop = FALSE]
        phi_E = Gaussian_basis(x = ts_E, mean = gknots, sd = sig)
        # quan_E_nume = phi_E/(h0tss_E%*%matrix(1,1,numIntKnt+2))
        quan_E_nume = phi_E / as.numeric(h0tss_E) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
      }
      # cphi = sapply(X = gknots, FUN = pnorm, q = ts, sd = sig)*sqrt(2*pi*sig^2) - matrix(1,n,1)%*%pnorm(0,gknots,sig)*sqrt(2*pi*sig^2)
      # cphi = matrix(0,n,numIntKnt+2)
      # for (i in 1:n) {
      #   cphi[i,] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)
      # }
      cphi = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)
      quan_deno = cphi
      
      nume = colSums(quan_E_nume) - 2*smooth*R%*%theta_old*(R%*%theta_old<=0)+eta
      deno = colSums(quan_deno) + 2*smooth*R%*%theta_old*(R%*%theta_old>0)+eta
      
      # 1.5. update theta temporarily
      theta_new = theta_old*(nume/deno)
      theta_inc = theta_new-theta_old
      
      grad_theta = nume - deno
      grad_theta_iter[k,] = t(grad_theta)
      
      # 1.5. hazard and cumulative hazard
      # if (n_E==0) {h0ts_E = h0tss_E = 1}
      # else {
      #   h0ts_E = phi_E%*%theta_new
      #   h0tss_E = h0ts_E
      #   h0tss_E[h0tss_E<1e-40] = 1e-40
      # }
      # ch0ts = cphi%*%theta_new
      
      # Gaussian basis part
      h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_new
      ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_new
      # h0ts = matrix(0,n,1)
      # ch0ts = matrix(0,n,1)
      # for (i in 1:n) { # will improve using apply function
      #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_new
      #   ch0ts[i] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)%*%theta_new
      # }
      if (n_E==0) {h0ts_E = h0tss_E = 1}
      else {
        h0ts_E = h0ts[censor_n==1, , drop = FALSE]
        h0tss_E = h0ts_E
        h0tss_E[h0tss_E<1e-40] = 1e-40
      }
      
      # 1.5. log-likelihood and penalised log-likelihood
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new) - sum(ch0ts)
      }
      else {like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new) - sum(ch0ts)}
      plike_new = like_theta_new - smooth*t(theta_new)%*%R%*%theta_new
      
      penlike_theta_iter_before_ls[k] = plike_new
      # 2. backtracking line search for theta
      alpha_MI = 1
      iteration_theta = 0L
      while (plike_new <= plike_old) {
        iteration_theta = iteration_theta + 1L
        if (alpha_MI >= 1e-2) {alpha_MI = alpha_MI*0.6}
        else if (alpha_MI < 1e-2 & alpha_MI >= 1e-5) {alpha_MI = alpha_MI*5e-2}
        else if (alpha_MI < 1e-5 & alpha_MI >= 1e-20) {alpha_MI = alpha_MI*1e-5}
        # else if (alpha_MI<1e-20 & alpha_MI>=1e-100) {alpha_MI=alpha_MI*1e-10}
        else {break}
        
        # 2. update theta with search direction
        theta_new = theta_old + alpha_MI * theta_inc
        
        # 2. hazard and cumulative hazard
        # if (n_E==0) {h0ts_E = h0tss_E = 1}
        # else {
        #   h0ts_E = phi_E%*%theta_new
        #   h0tss_E = h0ts_E
        #   h0tss_E[h0tss_E<1e-40] = 1e-40
        # }
        # ch0ts = cphi%*%theta_new
        
        # Gaussian basis part
        h0ts = Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_new
        ch0ts = Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_new
        # h0ts = matrix(0,n,1)
        # ch0ts = matrix(0,n,1)
        # for (i in 1:n) { # will improve using apply function
        #   h0ts[i] = t(dnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2))%*%theta_new
        #   ch0ts[i] = t(pnorm(ts[i],gknots,sig)*sqrt(2*pi*sig^2) - gauss)%*%theta_new
        # }
        if (n_E==0) {h0ts_E = h0tss_E = 1}
        else {
          h0ts_E = h0ts[censor_n==1, , drop = FALSE]
          h0tss_E = h0ts_E
          h0tss_E[h0tss_E<1e-40] = 1e-40
        }
        
        # 2. log-likelihood and penalised log-likelihood
        if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
          like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new - Z_ni_E%*%gamma_new) - sum(ch0ts)
        }
        else {like_theta_new = sum(log(h0tss_E) - X_E%*%beta_new) - sum(ch0ts)}
        plike_new = like_theta_new - smooth*t(theta_new)%*%R%*%theta_new
      }
      
      penlike_theta_iter_after_ls[k] = plike_new
      iter_MI_theta[k] = iteration_theta
      theta_iter[k,] = theta_new
      ts_range_theta_iter[k,] = range(ts)
      
      knots_iter[k,] = gknots
      
      iteration_no = k
      
      if ((all(abs(beta_new-beta_old)>knots_fix_threshold) | all(abs(gamma_new-gamma_old)>knots_fix_threshold)) & k<=knots_fix_iter_no & frequent_knots_update==FALSE) { # 20220208 add frequent_knots_update check
        if (knots_option=="equal_space") {
          bryKnt = range(ts)
          bin_width = (bryKnt[2] - bryKnt[1]) / (numIntKnt+1)
          gknots = seq(bryKnt[1], bryKnt[2], bin_width)
          sig = (2/3) * bin_width
        }
        else if (knots_option=="percentile") {
          # different strategy of calculating knots and sigmas
          gknots = quantile(sort(ts), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
          dist = sapply(gknots, function(x) {abs(x - ts)})
          sig = apply(dist, 2, function(x) {quantile(x, sig_coverage) / 2})
        }
      }
      ########################## stopping criteria ##########################
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
        beta_hat = beta_new
        gamma_hat = gamma_new
        theta_hat = theta_new
        if (all(abs(beta_new-beta_old)<threshold) & all(abs(gamma_new-gamma_old)<threshold) & all(abs(theta_new-theta_old)<threshold)) {
          cat("The stopping criteria are met and the results are obtained after",iteration_no,"iterations.\n")
          break
        }
        else if (iteration_no==maxiter) {
          warning("The stopping criteria are NOT met after ",iteration_no," itreations!\n")
        }
        else {
          beta_old = beta_new
          gamma_old = gamma_new
          if (k==knots_fix_iter_no) {theta_old = matrix(1, numIntKnt+2, 1)} # 20220208 inspired by Jun's modification
          else {theta_old = theta_new}
          like_beta_old = like_theta_new
          # print(paste("The penalised likelihood after iteration No.",k," is ",plike_new))
        }
      }
      else {
        beta_hat = beta_new
        theta_hat = theta_new
        if (all(abs(beta_new-beta_old)<threshold) & all(abs(theta_new-theta_old)<threshold)) {
          cat("The stopping criteria are met and the results are obtained after",iteration_no,"iterations.\n")
          break
        }
        else if (iteration_no==maxiter) {
          warning("The stopping criteria are NOT met after ",iteration_no," itreations!\n")
        }
        else {
          beta_old = beta_new
          if (k==knots_fix_iter_no) {theta_old = matrix(1, numIntKnt+2, 1)} # 20220208 inspired by Jun's modification
          else {theta_old = theta_new}
          like_beta_old = like_theta_new
          # print(paste("The penalised likelihood after iteration No.",k," is ",plike_new))
        }
      }
      # test if nonparametric estimates of AFT model are noisy with only 5 iterations
      # if (k==5){break}
    }
    
    # }, error=function(e){ #20181129
    #   theta_new=matrix(NA,numSp,1)
    #   beta_new=matrix(NA,p,1)
    # })
    
    ########################## estimation results ##########################
    time_range = rbind(range(time),range(ts))
    rownames(time_range) <- c("original time","accelerated time")
    colnames(time_range) <- c("min","max")
    ts_range_beta_iter = ts_range_beta_iter[1:k,]
    ts_range_theta_iter = ts_range_theta_iter[1:k,]
    knots_iter = knots_iter[1:k,]
    if (any(abs(diff(knots_iter[,1]))<=1e-12 & abs(diff(knots_iter[,ncol(knots_iter)]))<=1e-12)) {
      cat("The knots location stops changing since the",max(which(abs(diff(knots_iter[,1]))>1e-12 | abs(diff(knots_iter[,ncol(knots_iter)]))>1e-12))+1,"th iteration.\n")
      # cat("The knots location stops changing since the",max(which(abs(diff(knots_iter[,1]))>.Machine$double.eps | abs(diff(knots_iter[,ncol(knots_iter)]))>.Machine$double.eps))+1,"th iteration.\n")
    }
    multiple_ts_time = diff(range(ts))/diff(range(time))
    if (multiple_ts_time>1) {
      warning("This is a decelerated failure time case! The range of the decelerated time is ",multiple_ts_time," times that of the original time!\n")
    }
    
    grad_beta_iter = grad_beta_iter[1:k,]
    Hes_beta_iter = Hes_beta_iter[1:k,]
    det_Hes_beta_iter = det_Hes_beta_iter[1:k]
    eigenvl_Hes_beta_iter = eigenvl_Hes_beta_iter[1:k,]
    beta_iter = beta_iter[1:k,]
    
    grad_theta_iter = grad_theta_iter[1:k,]
    theta_iter = theta_iter[1:k,]
    
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      cvg = cbind(likelihood_beta_iter_before_ls,iter_Newton_beta,likelihood_beta_iter_after_ls,likelihood_gamma_iter_before_ls,iter_Newton_gamma,likelihood_gamma_iter_after_ls,penlike_theta_iter_before_ls,iter_MI_theta,penlike_theta_iter_after_ls)[1:k,]
      colnames(cvg) = c("like_beta_iter_before_ls","iters_ls_beta","like_beta_iter_after_ls","like_gamma_iter_before_ls","iters_ls_gamma","like_gamma_iter_after_ls","penlike_theta_iter_before_ls","iters_ls_theta","penlike_theta_iter_after_ls")
      grad_gamma_iter = grad_gamma_iter[1:k,]
      Hes_gamma_iter = Hes_gamma_iter[1:k,]
      det_Hes_gamma_iter = det_Hes_gamma_iter[1:k]
      eigenvl_Hes_gamma_iter = eigenvl_Hes_gamma_iter[1:k,]
      ts_range_gamma_iter = ts_range_gamma_iter[1:k,]
      gamma_iter = gamma_iter[1:k,]
    }
    else {
      cvg = cbind(likelihood_beta_iter_before_ls,iter_Newton_beta,likelihood_beta_iter_after_ls,penlike_theta_iter_before_ls,iter_MI_theta,penlike_theta_iter_after_ls)[1:k,]
      colnames(cvg) = c("like_beta_iter_before_ls","iters_ls_beta","like_beta_iter_after_ls","penlike_theta_iter_before_ls","iters_ls_theta","penlike_theta_iter_after_ls")
    }
    
    
    # browser()
    ########################## Asymptotic Normality ##########################
    # Hessian calculation, 3 Full Hessians, 3 Cross Hessians
    dts_beta = -X * as.numeric(ts) # n by p
    dts_E_beta = dts_beta[censor_n==1, , drop = FALSE] # n_E by p
    
    phi = Gaussian_basis(x = ts, mean = gknots, sd = sig)
    phi_E = phi[censor_n==1, , drop = FALSE]
    h0ts = phi%*%theta_hat
    h0ts_E = h0ts[censor_n==1, , drop = FALSE]
    h0tss_E = h0ts_E
    h0tss_E[h0tss_E<1e-40] = 1e-40
    dh0ts = Gaussian_basis_deriv1(x = ts, mean = gknots, sd = sig)%*%theta_hat
    dh0ts_E = dh0ts[censor_n==1, , drop = FALSE]
    dphi_E = Gaussian_basis_deriv1(x = ts_E, mean = gknots, sd = sig)
    ddh0ts_E = Gaussian_basis_deriv2(x = ts_E, mean = gknots, sd = sig)%*%theta_hat
    
    # full Hessian beta
    if (n_E==0) {quan_E_Hessian_beta = matrix(0,1,1)}
    else {quan_E_Hessian_beta = (ddh0ts_E / h0tss_E - dh0ts_E^2 / h0tss_E^2) * ts_E^2 + dh0ts_E / h0tss_E * ts_E}
    quan_Hessian_beta = dh0ts * ts^2 + h0ts * ts
    
    # Hessian_beta = t(-X_E)%*%(quan_E_Hessian_beta%*%matrix(1,1,p)*(-X_E)) - t(-X)%*%(quan_Hessian_beta%*%matrix(1,1,p)*(-X)) # p by p
    Hessian_beta = t(-X_E)%*%(as.numeric(quan_E_Hessian_beta)*(-X_E)) - t(-X)%*%(as.numeric(quan_Hessian_beta)*(-X)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
    
    # full Hessian beta theta
    # quan_E_Hessian_beta_theta = dphi_E / (h0tss_E%*%matrix(1,1,numIntKnt+2)) - phi_E * (dh0ts_E / h0tss_E^2)%*%matrix(1,1,numIntKnt+2) # n_E by (numIntKnt+2)
    quan_E_Hessian_beta_theta = dphi_E / as.numeric(h0tss_E) - phi_E * as.numeric(dh0ts_E / h0tss_E^2) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
    quan_Hessian_beta_theta = phi # n by (numIntKnt+2)
    
    Hessian_beta_theta = t(dts_E_beta)%*%quan_E_Hessian_beta_theta - t(dts_beta)%*%quan_Hessian_beta_theta # p by (numIntKnt+2)
    
    # full Hessian theta
    if (n_E==0) {Hessian_theta_no_penalty = matrix(0,1,numIntKnt+2)}
    else {
      # quan_E_Hessian_theta = t(phi_E / (h0tss_E%*%matrix(1,1,numIntKnt+2)))%*%(phi_E / (h0tss_E%*%matrix(1,1,numIntKnt+2)))
      Hessian_theta_no_penalty = -t(phi_E / as.numeric(h0tss_E))%*%(phi_E / as.numeric(h0tss_E)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
    }
    
    # Hessian_theta = Hessian_theta_no_penalty - 2*smooth*R # (numIntKnt+2) by (numIntKnt+2)
    
    # browser()
    # full Hessian gamma
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      if (q==1) {
        # dts_gamma = (exp(-X%*%beta_new)%*%matrix(1,1,q)) * as.matrix(tapply((-Z)*((exp(-Z%*%gamma_old)*interval_length)%*%matrix(1,1,q)), id_Z, FUN = sum)) # n by 1
        dts_gamma = as.numeric(exp(-X%*%beta_new)) * as.matrix(tapply((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length), id_Z, FUN = sum)) # 20220119 faster alternative of above line, use R recycling rule, may cause trouble
        # dts_gamma = as.numeric(exp(-X%*%beta_new)) * apply(as.matrix((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length)), 2, function(x) tapply(x, id_Z, sum)) # 20220119 slower alternative of above line, but works for q > 0
      }
      else {
        # dts_gamma = (exp(-X%*%beta_new)%*%matrix(1,1,q)) * as.matrix(tapply((-Z)*((exp(-Z%*%gamma_old)*interval_length)%*%matrix(1,1,q)), id_Z, FUN = colSums)) # n by q
        dts_gamma = as.numeric(exp(-X%*%beta_new)) * apply((-Z)*as.numeric(exp(-Z%*%gamma_old)*interval_length), 2, function(x) tapply(x, id_Z, sum))
      }
      dts_E_gamma = dts_gamma[censor_n==1, , drop = FALSE] # n_E by q
      
      if (n_E==0) {quan_E_Hessian_gamma_1 = quan_E_Hessian_gamma_2 = matrix(0,1,1)}
      else {
        quan_E_Hessian_gamma_1 = ddh0ts_E / h0tss_E - dh0ts_E^2 / h0tss_E^2 # n_E by 1
        quan_E_Hessian_gamma_2 = rep(dh0ts_E / h0tss_E * exp(-X_E%*%beta_hat),ni_E) * (exp(-Z_E%*%gamma_hat)*interval_length_E)
      }
      quan_Hessian_gamma_1 = dh0ts
      quan_Hessian_gamma_2 = rep(h0ts * exp(-X%*%beta_hat),ni) * (exp(-Z%*%gamma_hat)*interval_length) # N by 1, rep(c(),c()) suggested by Ali Shariati 20190804
      
      # Hessian_gamma = t(dts_E_gamma)%*%(quan_E_Hessian_gamma_1%*%matrix(1,1,q)*dts_E_gamma) +
      #                 t(-Z_E)%*%(quan_E_Hessian_gamma_2%*%matrix(1,1,q)*(-Z_E)) -
      #                 t(dts_gamma)%*%(quan_Hessian_gamma_1%*%matrix(1,1,q)*dts_gamma) -
      #                 t(-Z)%*%(quan_Hessian_gamma_2%*%matrix(1,1,q)*(-Z)) # q by q
      Hessian_gamma = t(dts_E_gamma)%*%(as.numeric(quan_E_Hessian_gamma_1)*dts_E_gamma) +
        t(-Z_E)%*%(as.numeric(quan_E_Hessian_gamma_2)*(-Z_E)) -
        t(dts_gamma)%*%(as.numeric(quan_Hessian_gamma_1)*dts_gamma) -
        t(-Z)%*%(as.numeric(quan_Hessian_gamma_2)*(-Z)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
      # browser()
      # full Hessian beta gamma
      if (n_E==0) {quan_E_Hessian_beta_gamma_1 = quan_E_Hessian_beta_gamma_2 = matrix(0,1,1)}
      else {
        quan_E_Hessian_beta_gamma_1 = ddh0ts_E / h0tss_E - dh0ts_E^2 / h0tss_E^2 # n_E by 1
        quan_E_Hessian_beta_gamma_2 = rep(dh0ts_E / h0tss_E * exp(-X_E%*%beta_hat),ni_E)*(exp(-Z_E%*%gamma_hat)*interval_length_E) # N_E by 1
        # quan_E_Hessian_beta_gamma = (ddh0ts_E / h0tss_E - dh0ts_E^2 / h0tss_E^2) * ts_E + dh0ts_E / h0tss_E
      }
      quan_Hessian_gamma_1 = dh0ts # n by 1
      quan_Hessian_gamma_2 = rep(h0ts * exp(-X%*%beta_hat),ni)*(exp(-Z%*%gamma_hat)*interval_length) # N by 1
      # quan_Hessian_beta_gamma = dh0ts * ts + h0ts
      
      # use data[rep(seq_len(nrow(data)), a_vector), ] to create repetitions of rows of data according to a_vector # https://statisticsglobe.com/repeat-rows-of-data-frame-n-times-in-r for reference
      # Hessian_beta_gamma = t(dts_E_beta)%*%(quan_E_Hessian_gamma_1%*%matrix(1,1,q)*dts_E_gamma) +
      #                      t(-X_E[rep(seq_len(n_E),ni_E), ])%*%(quan_E_Hessian_gamma_2%*%matrix(1,1,q)*(-Z_E)) - # under construction, think twice of this line
      #                      t(dts_beta)%*%(quan_Hessian_gamma_1%*%matrix(1,1,q)*dts_gamma) -
      #                      t(-X[rep(seq_len(n),ni), ])%*%(quan_Hessian_gamma_2%*%matrix(1,1,q)*(-Z)) # p by q
      
      if (n_E==0) {X_E_N = matrix(0, 1, p)}
      else {X_E_N = X_E[rep(seq_len(n_E), ni_E), , drop = FALSE]}
      X_N = X[rep(seq_len(n), ni), , drop = FALSE]
      
      Hessian_beta_gamma = t(dts_E_beta)%*%(as.numeric(quan_E_Hessian_gamma_1)*dts_E_gamma) +
        t(-X_E_N)%*%(as.numeric(quan_E_Hessian_gamma_2)*(-Z_E)) - # under construction, think twice of this line
        t(dts_beta)%*%(as.numeric(quan_Hessian_gamma_1)*dts_gamma) -
        t(-X_N)%*%(as.numeric(quan_Hessian_gamma_2)*(-Z)) # 20220117 faster alternative of above line, use R recycling rule, may cause trouble
      # Hessian_beta_gamma_1 = t(-X_E)%*%(as.numeric(quan_E_Hessian_beta_gamma)*dts_E_gamma) -
      # t(-X)%*%(as.numeric(quan_Hessian_beta_gamma)*dts_gamma) # under construction, inconsistent with Hessian_beta_gamma
      # browser()
      # term_Hessian_beta_gamma_E_Ding = t(dts_E_beta)%*%(quan_E_Hessian_gamma_1%*%matrix(1,1,q)*dts_E_gamma) +
      #   t(-X_E[rep(seq_len(n_E),ni_E), ])%*%(quan_E_Hessian_gamma_2%*%matrix(1,1,q)*(-Z_E))
      # term_Hessian_beta_gamma_E_Jun = t(-X_E)%*%(dts_E_gamma*((ddh0ts_E*ts_E+dh0ts_E) / h0tss_E - dh0ts_E^2*ts_E / h0tss_E^2))
      
      # full Hessian gamma theta
      Hessian_gamma_theta = t(dts_E_gamma)%*%quan_E_Hessian_beta_theta - t(dts_gamma)%*%quan_Hessian_beta_theta # q by (numIntKnt+2)
      
      # full Hessian
      F_mat = -rbind(cbind(Hessian_beta,Hessian_beta_gamma,Hessian_beta_theta),
                     cbind(t(Hessian_beta_gamma),Hessian_gamma,Hessian_gamma_theta),
                     cbind(t(Hessian_beta_theta),t(Hessian_gamma_theta),Hessian_theta_no_penalty)) # (p+q+numIntKnt+2) by (p+q+numIntKnt+2)
      Q_mat = matrix(0, nrow = p+q+numIntKnt+2, ncol = p+q+numIntKnt+2)
      Q_mat[(p+q+1):(p+q+numIntKnt+2),(p+q+1):(p+q+numIntKnt+2)] = 2*smooth*R
      F_mat_penalised = F_mat + Q_mat
      # F_mat_penalised = -rbind(cbind(Hessian_beta,Hessian_beta_gamma,Hessian_beta_theta),
      #                          cbind(t(Hessian_beta_gamma),Hessian_gamma,Hessian_gamma_theta),
      #                          cbind(t(Hessian_beta_theta),t(Hessian_gamma_theta),Hessian_theta)) # (p+q+numIntKnt+2) by (p+q+numIntKnt+2)
      
      # if (any(grad_theta<(-1e-1))) { # if else to deal with 0 active constraints
      if (any(grad_theta<0 & theta_hat<1e-2)) { # 20220213
        # U_index_active_constraints = which(grad_theta<(-1e-1))+p+q
        U_index_active_constraints = which(grad_theta<0 & theta_hat<1e-2)+p+q # 20220213
        # U_index_active_constraints = c(4,5)
        U_mat = diag(1L, nrow = p+q+numIntKnt+2)[, -U_index_active_constraints]
      }
      else {U_mat = diag(1L, nrow = p+q+numIntKnt+2)}
    }
    else {
      # full Hessian
      F_mat = -rbind(cbind(Hessian_beta,Hessian_beta_theta),
                     cbind(t(Hessian_beta_theta),Hessian_theta_no_penalty)) # (p+numIntKnt+2) by (p+numIntKnt+2)
      Q_mat = matrix(0, nrow = p+numIntKnt+2, ncol = p+numIntKnt+2)
      Q_mat[(p+1):(p+numIntKnt+2),(p+1):(p+numIntKnt+2)] = 2*smooth*R
      F_mat_penalised = F_mat + Q_mat
      # F_mat_penalised = -rbind(cbind(Hessian_beta,Hessian_beta_theta),
      #                          cbind(t(Hessian_beta_theta),Hessian_theta)) # (p+numIntKnt+2) by (p+numIntKnt+2)
      # if (any(grad_theta<(-1e-1))) {
      if (any(grad_theta<0 & theta_hat<1e-2)) { # 20220213
        # U_index_active_constraints = which(grad_theta<(-1e-1))+p
        U_index_active_constraints = which(grad_theta<0 & theta_hat<1e-2)+p # 20220213
        U_mat = diag(1L, nrow = p+numIntKnt+2)[, -U_index_active_constraints]
      }
      else {U_mat = diag(1L, nrow = p+numIntKnt+2)}
    }
    
    # variance-covariance matrix calculation
    B_mat_inv = U_mat%*%solve(t(U_mat)%*%F_mat_penalised%*%U_mat)%*%t(U_mat)
    cov_mat = t(B_mat_inv)%*%F_mat%*%B_mat_inv
    if (any(diag(cov_mat)<0)) {cov_mat = B_mat_inv} # 20220208
    # table of estimates
    var_diag = diag(cov_mat)
    asym_var_beta = as.matrix(var_diag[1:p])
    asym_std_beta = sqrt(asym_var_beta)
    Z_score_beta = beta_hat / sqrt(asym_var_beta)
    p_value_beta = 2*pnorm(abs(Z_score_beta), lower.tail = FALSE)
    # significance_beta = rep(TRUE,p)
    # significance_beta[p_value_beta>=0.05] = FALSE
    if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
      asym_var_gamma = as.matrix(var_diag[(p+1):(p+q)])
      asym_std_gamma = sqrt(asym_var_gamma)
      Z_score_gamma = gamma_hat / sqrt(asym_var_gamma)
      p_value_gamma = 2*pnorm(abs(Z_score_gamma), lower.tail = FALSE)
      # significance_gamma = rep(TRUE,1)
      # significance_gamma[p_value_gamma>=0.05] = FALSE
      asym_var_theta = as.matrix(var_diag[(p+q+1):(p+q+numIntKnt+2)])
      asym_std_theta = sqrt(asym_var_theta)
      Z_score_theta = theta_hat / sqrt(asym_var_theta)
      p_value_theta = 2*pnorm(abs(Z_score_theta), lower.tail = FALSE)
      # significance_theta = rep(TRUE,numIntKnt+2)
      # significance_theta[p_value_theta>=0.05] = FALSE
      table_coef = rbind(cbind(beta_hat, asym_std_beta, Z_score_beta, p_value_beta),
                        cbind(gamma_hat, asym_std_gamma, Z_score_gamma, p_value_gamma),
                        cbind(theta_hat, asym_std_theta, Z_score_theta, p_value_theta))
      rownames(table_coef)[(p+q+1):(p+q+numIntKnt+2)] <- c(sprintf("basis_function_%d", 1:(numIntKnt+2)))
    }
    else {
      asym_var_theta = as.matrix(var_diag[(p+1):(p+numIntKnt+2)])
      asym_std_theta = sqrt(asym_var_theta)
      Z_score_theta = theta_hat / sqrt(asym_var_theta)
      p_value_theta = 2*pnorm(abs(Z_score_theta), lower.tail = FALSE)
      # significance_theta = rep(TRUE,numIntKnt+2)
      # significance_theta[p_value_theta>=0.05] = FALSE
      table_coef = rbind(cbind(beta_hat, asym_std_beta, Z_score_beta, p_value_beta),
                        cbind(theta_hat, asym_std_theta, Z_score_theta, p_value_theta))
      rownames(table_coef)[(p+1):(p+numIntKnt+2)] <- c(sprintf("basis_function_%d", 1:(numIntKnt+2)))
    }
    colnames(table_coef) <- c("coef", "se(coef)", "z", "p(>|z|)")
    # colnames(table_coef) <- c("est", "asym_std", "Z", "p", "significant? (<0.05)")
    
    # browser()
    df = numIntKnt+2 - sum(diag(solve(t(U_mat)%*%F_mat_penalised%*%U_mat)%*%(t(U_mat)%*%Q_mat%*%U_mat)))
    if (df < 0) {df = 1}
    smooth_new = as.numeric(df / (2 * t(theta_hat)%*%R%*%theta_hat))
    
    # # Asymptotics
    # dbeta=grad
    # dtheta=nume-deno
    # 
    # vv=as.numeric(exp(colMeans(Xmat)%*%beta_hat))
    # if (any(abs(theta_hat)<1e-5)) {
    #   id0theta=intersect(which(abs(theta_hat)<1e-5),which(dbeta<1e-5))
    #   # id0theta=which(abs(theta_hat)<1e-5)
    # }
    # else{id0theta = c()}
    # idtheta=setdiff(1:numSp, id0theta)
    # 
    # 
    # ddphi= mSpline(ts, knots=IntKnt, degree=dgrSp, intercept=TRUE, Boundary.knots=bryKnt, derivs=2)
    # ddh0ts=as.matrix(ddphi%*%theta_hat)
    # # 2nd derivative beta, p by p matrix
    # pl_beta=t(X)%*% ((( censor_n*((ddh0ts*ts^2+dh0ts*ts)/h0tss-(dh0ts*ts/h0tss)^2) - (ddh0ts*ts^2+dh0ts*ts) )%*%matrix(1,1,p))*X)
    # # 2nd derivative theta, m by m matrix
    # l_theta=-t(h0tsM*phi)%*%((1/h0tss)%*%matrix(1,1,numSp)*phi)
    # l_theta_s=l_theta[idtheta,idtheta]
    # pl_theta=l_theta-smooth*R
    # pl_theta_s=pl_theta[idtheta,idtheta]
    # # 2nd derivative beta and theta, p by m matrix
    # pl_beta_theta= -t(X)%*% ( (censor_n*ts/h0tss)%*%matrix(1,1,numSp)*dphi - (censor_n*ts*dh0ts/h0tss^2)%*%matrix(1,1,numSp)*phi  -  (ts%*%matrix(1,1,numSp))*phi)
    # pl_beta_theta_s = pl_beta_theta[,idtheta]
    # 
    # # F-1G
    # # d2PLs=rbind(cbind(pl_beta,pl_beta_theta_s),cbind(t(pl_beta_theta_s),pl_theta_s))
    # # dfHs=solve(d2PLs, rbind(cbind(pl_beta,pl_beta_theta_s),cbind(t(pl_beta_theta_s),l_theta_s)))
    # # Asym_V_s = solve(d2PLs,dfHs);
    # # F_matrix_inv=solve(-rbind(cbind(pl_beta,pl_beta_theta),cbind(t(pl_beta_theta),pl_theta)) )
    # # G_matrix=-rbind(cbind(pl_beta,pl_beta_theta),cbind(t(pl_beta_theta),l_theta))
    # # Asym_V=F_matrix_inv%*%G_matrix%*%t(F_matrix_inv)
    # F_matrix_inv=solve(rbind(cbind(pl_beta,pl_beta_theta_s),cbind(t(pl_beta_theta_s),pl_theta_s)) )
    # G_matrix=-rbind(cbind(pl_beta,pl_beta_theta_s),cbind(t(pl_beta_theta_s),l_theta_s))
    # Asym_V_s=F_matrix_inv%*%G_matrix%*%t(F_matrix_inv)
    # # browser()
    # B11=Asym_V_s[1:p,1:p]
    # B12=matrix(0,p,numSp)
    # B12[,idtheta]=Asym_V_s[1:p,(p+1):(p+length(idtheta))]
    # B22=matrix(0,numSp,numSp)
    # B22[idtheta,idtheta]=Asym_V_s[(p+1):(p+length(idtheta)),(p+1):(p+length(idtheta))]
    # 
    # Asym_V=rbind(cbind(B11,B12),cbind(t(B12),B22))
    # Asym_V_diag=diag(Asym_V)
    # 
    # Asym_V_beta=Asym_V_diag[1:p]
    # Asym_V_theta=Asym_V_diag[-(1:p)]
    # # Asym_v_h0=vv^2*diag(phi%*%diag(Asym_V_theta)%*%t(phi))
    # 
    # ASD_beta=sqrt(Asym_V_beta)
    # ASD_theta=sqrt(Asym_V_theta)
    # 
    # 
    # # S0=exp(-ch0ts*vv)
    # 
    # time_axis=seq(min(time),max(time),length.out=500)
    # phi_axis=mSpline(time_axis, knots = IntKnt, degree = dgrSp,intercept=TRUE)
    # Asym_v_h0=vv^2*diag(phi_axis%*%diag(Asym_V_theta)%*%t(phi_axis))
    # # Asym_v_h0=diag(phi_axis%*%diag(Asym_V_theta)%*%t(phi_axis))
    # 
    # ASD_h0=sqrt(Asym_v_h0)
    # h0_hat_plot=vv*phi_axis%*%theta_hat
    # h0_hat_plot=phi_axis%*%theta_hat
    # h0_hat_plot_a_upper=h0_hat_plot+1.96*ASD_h0
    # h0_hat_plot_a_lower=h0_hat_plot-1.96*ASD_h0
    # 
    # plot(time_axis,h0_hat_plot,type="l",main=bquote(paste("MPL approach with ",lambda==.(smooth))),xlab="Time",ylab="Estimated baseline hazard")
    # lines(time_axis,h0_hat_plot_a_upper,type="l",col="blue")
    # lines(time_axis,h0_hat_plot_a_lower,type="l",col="red")
    
    
    smooth_hat = smooth_new
    if (abs(df - df_old)<1) {
      break
    }
    else {
      smooth = smooth_new
      df_old = df
      beta_initial = beta_hat
      if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {gamma_initial = gamma_hat}
    }
  }
  smooth_iter = smooth_iter[1:j,]
  
  
  # browser()
  ################################# plots #################################
  if (draw_plot==TRUE) {
    # hazard plot
    # plots of final accelerated time
    # plot_hazard = cbind(ts, Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_hat)
    # plot_hazard = plot_hazard[order(plot_hazard[,1]), ] # sort according to ascending order of ts
    # plot(plot_hazard[,1], plot_hazard[,2], xlim = c(0,max(plot_hazard[,1])), ylim = c(0,max(plot_hazard[,2])), main = "hazard plot", xlab = "accelerated/decelerated time", ylab = "hazard")
    plot(ts, Gaussian_basis(x = ts, mean = gknots, sd = sig)%*%theta_hat, xlim = c(0,max(ts)), main = "hazard plot", xlab = "accelerated/decelerated time", ylab = "hazard")
    # smooth curve connecting data points
    lines(seq(0,max(ts),length.out=max(1000,10*n)), Gaussian_basis(x = as.matrix(seq(0,max(ts),length.out=max(1000,10*n))), mean = gknots, sd = sig)%*%theta_hat)
    
    # survival plot
    # plots of final accelerated time
    # plot_survival = cbind(ts, exp(-Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_hat))
    # plot_survival = plot_survival[order(plot_survival[,1]), ] # sort according to ascending order of ts
    # plot(plot_survival[,1], plot_survival[,2], xlim = c(0,max(plot_hazard[,1])), ylim = c(0,1), main = "Survival plot", xlab = "accelerated/decelerated time", ylab = "Survival")
    plot(ts, exp(-Gaussian_basis_integ1(q = ts, mean = gknots, sd = sig)%*%theta_hat), xlim = c(0,max(ts)), ylim = c(0,1), main = "Survival plot", xlab = "accelerated/decelerated time", ylab = "Survival")
    # smooth curve connecting data points
    lines(seq(0,max(ts),length.out=max(1000,10*n)), exp(-Gaussian_basis_integ1(q = as.matrix(seq(0,max(ts),length.out=max(1000,10*n))), mean = gknots, sd = sig)%*%theta_hat))
    
    # below plots for original time may be incorrect
    if (knots_option=="equal_space") {
      bryKnt_0 = range(time)
      bin_width_0 = (bryKnt_0[2] - bryKnt_0[1]) / (numIntKnt+1)
      gknots_0 = seq(bryKnt_0[1], bryKnt_0[2], bin_width_0)
      sig_0 = (2/3) * bin_width_0
    }
    else if (knots_option=="percentile") {
      gknots_0 = quantile(sort(time), seq(knots_min_quantile, knots_max_quantile, length.out = numIntKnt+2)/100, type=1)
      dist_0 = sapply(gknots_0, function(x) {abs(x - time)})
      sig_0 = apply(dist_0, 2, function(x) {quantile(x, sig_coverage) / 2})
    }
    # baseline hazard plot
    # plots of original time
    # plot_hazard_0 = cbind(time,Gaussian_basis(x = time, mean = gknots_0, sd = sig_0)%*%theta_hat)
    # plot_hazard_0 = plot_hazard_0[order(plot_hazard_0[,1]), ] # sort according to ascending order of time
    # plot(plot_hazard_0[,1], plot_hazard_0[,2], xlim = c(0,max(plot_hazard_0[,1])), ylim = c(0,max(plot_hazard_0[,2])), main = "baseline hazard plot", xlab = "time", ylab = "baseline hazard")
    plot(time, Gaussian_basis(x = time, mean = gknots_0, sd = sig_0)%*%theta_hat, xlim = c(0,max(time)), main = "baseline hazard plot", xlab = "time", ylab = "baseline hazard")
    # smooth curve connecting data points
    lines(seq(0,max(time),length.out=max(1000,10*n)), Gaussian_basis(x = as.matrix(seq(0,max(time),length.out=max(1000,10*n))), mean = gknots_0, sd = sig_0)%*%theta_hat)
    
    # baseline survival plot
    # plots of original time
    # plot_survival_0 = cbind(time,exp(-Gaussian_basis_integ1(q = time, mean = gknots_0, sd = sig_0)%*%theta_hat))
    # plot_survival_0 = plot_survival_0[order(plot_survival_0[,1]), ] # sort according to ascending order of time
    # plot(plot_survival_0[,1], plot_survival_0[,2], xlim = c(0,max(plot_hazard_0[,1])), ylim = c(0,1), main = "baseline Survival plot", xlab = "time", ylab = "baseline Survival")
    plot(time, exp(-Gaussian_basis_integ1(q = time, mean = gknots_0, sd = sig_0)%*%theta_hat), xlim = c(0,max(time)), ylim = c(0,1), main = "baseline Survival plot", xlab = "time", ylab = "baseline Survival")
    # smooth curve connecting data points
    lines(seq(0,max(time),length.out=max(1000,10*n)), exp(-Gaussian_basis_integ1(q = as.matrix(seq(0,max(time),length.out=max(1000,10*n))), mean = gknots_0, sd = sig_0)%*%theta_hat))
  }
  
  # browser()
  ############################### Output ###############################
  if (!missing(Zmat) & !missing(gamma_initial) & (dim_Zmat[1]>n+1)) {
    return(list(beta_hat=beta_hat,gamma_hat=gamma_hat,theta_hat=theta_hat,
                asym_var_beta=asym_var_beta,asym_var_gamma=asym_var_gamma,asym_var_theta=asym_var_theta,
                asym_std_beta=asym_std_beta,asym_std_gamma=asym_std_gamma,asym_std_theta=asym_std_theta,
                Z_score_beta=Z_score_beta,Z_score_gamma=Z_score_gamma,Z_score_theta=Z_score_theta,
                p_value_beta=p_value_beta,p_value_gamma=p_value_gamma,p_value_theta=p_value_theta,
                table_coef=table_coef,
                smooth_hat=smooth_hat, smooth_iter=smooth_iter,
                cvg=cvg,time_range=time_range,multiple_ts_time=multiple_ts_time,knots_iter=knots_iter,
                ts_range_beta_iter=ts_range_beta_iter,ts_range_gamma_iter=ts_range_gamma_iter,ts_range_theta_iter=ts_range_theta_iter,
                iteration_no=iteration_no,
                grad_beta_iter=grad_beta_iter,grad_gamma_iter=grad_gamma_iter,grad_theta_iter=grad_theta_iter,
                Hes_beta_iter=Hes_beta_iter,det_Hes_beta_iter=det_Hes_beta_iter,eigenvl_Hes_beta_iter=eigenvl_Hes_beta_iter,
                Hes_gamma_iter=Hes_gamma_iter,det_Hes_gamma_iter=det_Hes_gamma_iter,eigenvl_Hes_gamma_iter=eigenvl_Hes_gamma_iter,
                Hes_beta_final=Hes_beta,Hes_gamma_final=Hes_gamma,
                beta_iter=beta_iter,gamma_iter=gamma_iter,theta_iter=theta_iter)) #20181129
  }
  else {
    return(list(beta_hat=beta_hat,theta_hat=theta_hat,
                asym_var_beta=asym_var_beta,asym_var_theta=asym_var_theta,
                asym_std_beta=asym_std_beta,asym_std_theta=asym_std_theta,
                Z_score_beta=Z_score_beta,Z_score_theta=Z_score_theta,
                p_value_beta=p_value_beta,p_value_theta=p_value_theta,
                table_coef=table_coef,
                smooth_hat=smooth_hat, smooth_iter=smooth_iter,
                cvg=cvg,time_range=time_range,multiple_ts_time=multiple_ts_time,knots_iter=knots_iter,
                ts_range_beta_iter=ts_range_beta_iter,ts_range_theta_iter=ts_range_theta_iter,
                iteration_no=iteration_no,
                grad_beta_iter=grad_beta_iter,grad_theta_iter=grad_theta_iter,
                Hes_beta_iter=Hes_beta_iter,det_Hes_beta_iter=det_Hes_beta_iter,eigenvl_Hes_beta_iter=eigenvl_Hes_beta_iter,
                Hes_beta_final=Hes_beta,
                beta_iter=beta_iter,theta_iter=theta_iter)) #20181129
  }
  # return(list(beta=beta_hat,theta=theta_hat,beta_ASD=ASD_beta,h0_ASD=ASD_h0,h0=h0_hat_plot,theta_ASD=ASD_theta,h0_lower=h0_hat_plot_a_lower,h0_upper=h0_hat_plot_a_upper,beta_gradient=dbeta,theta_gradient=dtheta,penlike=penlike_theta_iter_after_ls))  #,"tstar"=t))
  # message("No. of iterations to convergence is ",k)
}
