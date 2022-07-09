# a fix of datagen_RC_TVC_discrete.R 20210915
datagen_RC_TVC_discrete_1 <- function(n = 100, beta = c(1,-1), gamma = -0.1, t_treat_option = 1, t_treat_min = 0, t_treat_max, c_min = 0, c_max, datagen_option = "RC_TVC") {
  # 20201110 updated to incorporate a separate regulator for t_treat, i.e. t_treat_min != c_min and t_treat_max != c_max
  # if(missing(c_max)) {stop("Please specify a maximum value (c_max) in order to generate right-censored time")}
  
  #################### generate time-fixed covariate matrix Xmat ####################
  id_X = 1:n
  X_1 = rbinom(n,1,0.5)
  # X_2 = rlnorm(n)
  X_3 = runif(n,min=0,max=3)
  # X = cbind(X_1,X_2,X_3)
  X = cbind(X_1,X_3)
  Xmat = cbind(id_X,X)
  
  ############################ generate time matrix tmat ############################
  S = runif(n,min = 0,max = 1) # survival
  
  if (datagen_option=="RC_TVC") {
    if (t_treat_option == 1) {t_treat = runif(n, min = t_treat_min, max = t_treat_max)} # treatment time (uniform distribution)
    else if (t_treat_option == 2) {t_treat = rep(t_treat_max,n)} # treatment time (all same and fixed)
    # # slower version of calculating t using for loop
    # # Z_i(t) == 0 if t < t_treat[i], Z_i(t) == 1 if t >= t_treat[i]
    # t = rep(NA, times = n)
    # for (i in 1:n) {
    #   if (-log(S[i]) < (exp(-X[i,]%*%beta) * t_treat[i])^3) {t[i] = exp(X[i,]%*%beta) * (-log(S[i]))^(1/3)}
    #   else {t[i] = exp(X[i,]%*%beta + gamma) * (-log(S[i]) - (exp(-X[i,]%*%beta) * t_treat[i])^3 + (exp(-X[i,]%*%beta - gamma) * t_treat[i])^3)^(1/3)}
    # }
    # faster alternative of calculating t using matrix multiplication
    less_or_more = (-log(S) < (exp(-X%*%beta) * t_treat)^3)
  }
  
  t = exp(X%*%beta) * (-log(S))^(1/3)
  
  if (datagen_option=="RC_TVC") {
    # # Jun's version
    # t_alternative = exp(X%*%beta + matrix(1,n,1)%*%gamma) * (-log(S) - (exp(-X%*%beta) * t_treat)^3 + (exp(-X%*%beta - matrix(1,n,1)%*%gamma) * t_treat)^3)^(1/3)
    # alternative version
    t_alternative = exp(X%*%beta + matrix(1,n,1)%*%gamma) * (-log(S))^(1/3) - exp(matrix(1,n,1)%*%gamma) * t_treat + t_treat
    t[less_or_more==FALSE] = t_alternative[less_or_more==FALSE]
  }
  t = as.numeric(t)
  
  c = runif(n, min = c_min, max = c_max)
  y = pmin(t,c) # y stands for the last observed time (the minimum between the event time and right-censored time)
  indicator = (t < c)*1 # indicator=1 for event time, indicator=0 for right-censored time
  tmat = cbind(y,indicator)
  censor_prop = 1 - sum(indicator)/n
  range_t = range(t)
  range_y = range(y)
  
  ################## generate time-varying covariates matrix Zmat ##################
  if (datagen_option=="RC_TVC") {
    id_Z = c(rbind(id_X,id_X))
    interval_index = rep(1:2, times = n)
    Z = rep(0:1, times = n)
    
    Zmat_deletion_mark_interval_1 = numeric(n)
    Zmat_deletion_mark_interval_2 = numeric(n)
    
    # if 0 < y < t_treat, need to mark the rows of second intervals to delete later
    Zmat_deletion_mark_interval_2[y>0 & y<t_treat] = 1
    # if 0 < y < t_treat, need to replace the relevant t_treat's by y
    t_treat[y>0 & y<t_treat] = y[y>0 & y<t_treat]
    
    Zmat_deletion_interval_indicator = c(rbind(Zmat_deletion_mark_interval_1,Zmat_deletion_mark_interval_2))
    
    interval_L = c(rbind(numeric(n),t_treat))
    interval_R = c(rbind(t_treat,y))
    
    Zmat = cbind(id_Z,interval_index,interval_L,interval_R,Z)
    
    Zmat = Zmat[Zmat_deletion_interval_indicator==0,]
    
    treat_prop = sum(Zmat[,2]==2)/n
    
    indicator_N = rep(indicator, times = tapply(Zmat[,2],Zmat[,1],FUN = max))
    event_treated_prop = sum(indicator_N==1 & Zmat[,5]==1) / n
    
    if (event_treated_prop==0) {warning("The generated dataset doesn't contain any event cases with time varying covariates! Please rerun the code and generate a different dataset!\n")}
    return(list(tmat=tmat, Xmat=Xmat, Zmat=Zmat, range_t=range_t, range_y=range_y, censor_prop=censor_prop, treat_prop=treat_prop, event_treated_prop=event_treated_prop))
  }
  else {
    return(list(tmat=tmat, Xmat=Xmat, range_t=range_t, range_y=range_y, censor_prop=censor_prop))
  }
}

