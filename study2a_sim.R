library(dplyr)
library(splines2)

rJM <- function(hazard, censoring, x, r, 
                subdivisions = 1000, tmin = 0, tmax, file = NULL, ...){
  ## compute hazard for every person i at time
  nsub <- nrow(x)
  
  time <- rep(NA, nsub) 
  
  Hazard <- function(hazard, time, x, r) { 
    integrate(hazard, 0, time, x = x, r = r,
              subdivisions = subdivisions)$value 
  } 
  
  InvHazard <- function(Hazard, hazard, x, r, tmin, tmax) { 
    negLogU <- -log(runif(1, 0, 1)) 
    # check if first Lambda value is smaller than sample
    rootfct <- function(time) { 
      negLogU - Hazard(hazard, time, x, r) 
    } 
    if(rootfct(tmin)<0){
      return(0)
    } else {
      root <- try(uniroot(rootfct, interval = c(0, tmax))$root, silent=TRUE)
      root <- if(inherits(root, "try-error")) {
        # if root not within [0, tmax] --> error --> set it to tmax + 0.01 (will be censored)
        tmax + 0.01
      } else {root}
    }
    return(root)
  }
  
  # Finding Survival Times
  cumhaz <- rep(NA, nsub)
  survprob <- rep(NA, nsub)
  for(i in 1:nsub) { 
    time[i] <- InvHazard(Hazard, hazard, x[i,], r[i,], tmin, tmax)
    cumhaz[i] <- Hazard(hazard, time[i], x[i,], r[i,])
    survprob[i] <- exp((-1)*cumhaz[i])
  } 
  
  time_event <- censoring(time, tmax, length(time)) 
  
  # Make data (long format)
  data_short <- data.frame(t_L = time_event[, 2], t_R = time_event[, 3], 
                           event = time_event[, 4], interval = time_event[, 5],
                           right = time_event[, 6], left = time_event[, 7],
                           x, r, cumhaz = cumhaz, true = time)
  names(data_short) <- gsub(".", "", names(data_short), fixed = TRUE) 
  
  return(data_short)
}

hazard_function_poly2 = function(event_time, x, r){
  
  kappa = 4
  alpha = 1
  
  #h0_t = 3*event_time^2
  h0_t = ((kappa/alpha) * (event_time/alpha)^(kappa-1))/(1 + (event_time/alpha)^kappa)
  
  w_yi = (0.5 + r$g0) + (0.5 + r$g1)*event_time - 
    0.1*event_time^2 + 0.2*event_time^3 - x[3]*0.2 - 0.5*x[3]*event_time - 
    0.2*x[3]*event_time^2
  
  true_beta1 = 0.5
  true_beta2 = -1
  true_gamma = -0.5
  
  ht = h0_t * (exp(c(true_beta1*x[1] + true_beta2*x[2] + true_gamma*w_yi)))
  return(ht)
  
}



censoring <- function(time, tmax, n){
  #set alphas to control width of interval censoring intervals
  a1 = 0.9
  a2 = 1.2
  
  #uniform variables
  U_E = runif(n)
  U_L = runif(n,0,2)
  U_R = runif(n,U_L,2)
  pi_E = 0
  events = as.numeric(U_E < pi_E & time < tmax)
  interval = as.numeric(a1*U_L <= time & time <= a2*U_R & U_E >= pi_E & time < tmax)
  right = as.numeric((a2*U_R < time & U_E >= pi_E) | time >= tmax)
  left = as.numeric(time < a1*U_L & U_E >= pi_E & time < tmax)
  
  t_L = t_R = rep(0, n)
  
  t_L[which(events==1)] = time[which(events==1)]
  t_R[which(events==1)] = time[which(events==1)]
  
  t_L[which(interval==1)] = (a1*U_L)[which(interval==1)]
  t_R[which(interval==1)] = (a2*U_R)[which(interval==1)]
  
  t_L[which(right==1)] = (a2*U_R)[which(right==1)]
  t_R[which(right==1)] = rep(Inf, sum(right))
  
  t_R[which(left==1)] = (a1*U_L)[which(left==1)]
  
  
  return(cbind(time, t_L, t_R, events, interval, right, left, U_E))
}



gen_interval_censoring_jointmodel = function(n, nobs){
  
  g0 = rnorm(n, 0, 0.3)
  g1 = rnorm(n, 0, 0.3)
  
  r = data.frame(g0, g1)
  
  fixed_x1 = rnorm(n, 0, 1)
  #fixed_x1 = fixed_x1 - mean(fixed_x1)
  fixed_x2 = rbinom(n, 1, 0.5)
  #fixed_x2 = fixed_x2 - mean(fixed_x2)
  fixed_x3 = rbinom(n, 1, 0.5)
  #fixed_x3 = fixed_x3 - mean(fixed_x3)
  fixed_x = as.matrix(cbind(fixed_x1, fixed_x2, fixed_x3), ncol = 3)
  
  dat.baseline = rJM(hazard_function_poly2, censoring, x = fixed_x, r = r, tmin = 0, tmax = 3)
  
  
  dat.baseline$id = c(1:n)
  dat.baseline$latest_time = dat.baseline$t_R
  dat.baseline$latest_time[which(dat.baseline$latest_time == Inf)] = dat.baseline$t_L[which(dat.baseline$latest_time == Inf)]
  
  dat.baseline$mid_point = 0
  dat.baseline$mid_point[which(dat.baseline$right == 1)] = dat.baseline$t_L[which(dat.baseline$right == 1)]
  dat.baseline$mid_point[which(dat.baseline$event == 1)] = dat.baseline$t_L[which(dat.baseline$event == 1)]
  dat.baseline$mid_point[which(dat.baseline$interval == 1)] = dat.baseline$t_L[which(dat.baseline$interval == 1)] + 
    (dat.baseline$t_R[which(dat.baseline$interval == 1)] - dat.baseline$t_L[which(dat.baseline$interval == 1)])/2
  dat.baseline$mid_point[which(dat.baseline$left == 1)] = dat.baseline$t_R[which(dat.baseline$left == 1)]/2
  
  
  #for each i, generate a random sequence of observation times for w(t)
  #num.obs = rpois(n, 8) + 1
  if(!is.null(nobs)){
    num.obs = rep(nobs, n)
  }else{
    num.obs = rpois(n, 9) + 1
  }
  #sample random obs times from tv and pick up "true" w(t) values for each i,
  #truncated at event/censoring time
  id_long = NULL
  wt_samp = NULL
  wt_i_samp = NULL
  obs_time_samp = NULL
  e_i_t = NULL
  visit_num = NULL
  W_long = NULL
  before_midpoint = NULL
  
  for(i in 1:n){
    obs.times = c(cumsum(runif(100000, 0, min(dat.baseline$latest_time)/3)))
    obs.times = obs.times[which(obs.times <= dat.baseline$latest_time[i])]
    if(length(obs.times) > num.obs[i]){
      obs.incl = sample(1:length(obs.times), num.obs[i])
      obs.times = c(0, obs.times[obs.incl[order(obs.incl)]])
    }else{
      obs.times = c(0, obs.times)
    }
    
    #time_diff = runif((length(obs.times)-2), 0.01, 0.04)
    #obs.times[3:length(obs.times)] = obs.times[3:length(obs.times)] - time_diff
    obs.times = obs.times[which(obs.times <= dat.baseline$latest_time[i])]
    
    id_long = c(id_long, rep(i, length(obs.times)))
    obs_time_samp = c(obs_time_samp, obs.times)
    
    e_i = rnorm(length(obs.times), 0, 0.05)
    e_i_t = c(e_i_t, e_i)
    
    wt_tia = (0.5 + r$g0[i]) +
      (0.5 + r$g1[i])*obs.times -
      0.1*obs.times^2 + 0.2*obs.times^3 - 
      0.2*fixed_x3[i] - 0.5*fixed_x3[i]*obs.times - 
      0.2*fixed_x3[i]*obs.times^2 + e_i
    
    #(0.5+ r$g1[i])*obs.times - (1)*obs.times^2 + (0.5)*obs.times^3 +
    #(fixed_x3[i] > 0)*(0.5+ r$g2[i])*obs.times + (fixed_x3[i] > 0)*(0.5)*obs.times^2 +
    #e_i #measurement error
    
    wt_samp = c(wt_samp, wt_tia)
    wt_i_samp = c(wt_i_samp, wt_tia)
    visit_num = c(visit_num, c(1:length(wt_tia)))
    W_long = c(W_long, rep(fixed_x3[i],length(wt_tia)))
    before_midpoint = c(before_midpoint, as.numeric(obs.times < dat.baseline$mid_point[i]))
    
  }
  
  
  dat.long = data.frame(cbind(id_long, visit_num, wt_samp, wt_i_samp, obs_time_samp, e_i_t, W_long, before_midpoint))
  
  dat.baseline$n.obs = (dat.long %>% group_by(id_long) %>% tally())$n
  
  out = list(dat.baseline = dat.baseline, dat.long = dat.long)
  
  return(out)
  
}


dat = gen_interval_censoring_jointmodel(n=200, nobs = NULL)
