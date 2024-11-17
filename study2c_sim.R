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
  #h0_t = ((kappa/alpha) * (event_time/alpha)^(kappa-1))/(1 + (event_time/alpha)^kappa)
  
  h0_t = ((1/(sqrt(2*pi)*0.5*event_time+1e-4)) * exp(-0.5 * ((log(event_time+1e-4)+0.3)/0.5)^2))/(1 - pnorm((log(event_time+1e-4)+0.3)/0.5))
  
  w_yi1 =  -0.6 + (1-(0.75)*exp(4*event_time)/(1 + exp(4*event_time)) ) + r$g0 + +r$g5*event_time
  #w_yi1 = (r$g0 + 0.2) + (0.5 + r$g5)* event_time
  
  w_yi2 = (r$g6 + 0.1) + (0.2 + r$g7)* event_time - 0.2 * event_time^2 + x[1]*0.2
  
  true_beta1 = 0.2
  true_beta2 = -0.5
  true_gamma1 = 1
  true_gamma2 = 0
  
  #ht = h0_t * (exp(c(true_beta1*x[1] + true_beta2*x[2] + true_gamma1*w_yi1 + true_gamma2*w_yi2 )))
  ht = h0_t * (exp(c(true_beta1*x[1] + true_beta2*x[2] + true_gamma1*w_yi1 + + true_gamma2*w_yi2)))
  return(ht)
  
}



censoring <- function(time, tmax, n){
  #set alphas to control width of interval censoring intervals
  a1 = 0.7
  a2 = 1.1
  
  #uniform variables
  U_E = runif(n)
  U_E = runif(n)
  U_L = runif(n,0,1.5)
  #U_L = rweibull(n, 2, 0.5)
  U_R = runif(n,U_L,1.5)
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
  g0 = rnorm(n, 0, 0.05)
  #g0 = rnorm(n, 0, 0.1)
  g1 = rnorm(n, 0, 0.15)
  g2 = rnorm(n, 0, 0.1)
  g3 = rnorm(n, 0, 0.1)
  g4 = rnorm(n, 0, 0.05)
  g5 = rnorm(n, 0, 0.1)
  #g5 = rnorm(n, 0, 0.05)
  
  g6 = rnorm(n, 0, 0.2)
  g7 = rnorm(n, 0, 0.1)
  
  #r = data.frame(g0, g1, g2, g3, g4, g5, g6, g7)
  
  
  
  
  fixed_x1 = rnorm(n, 0, 1)
  #fixed_x1 = fixed_x1 - mean(fixed_x1)
  fixed_x2 = rbinom(n, 1, 0.5)
  #fixed_x2 = fixed_x2 - mean(fixed_x2)
  fixed_x3 = rep(1, n)
  #fixed_x3 = fixed_x3 - mean(fixed_x3)
  fixed_x = as.matrix(cbind(fixed_x1, fixed_x2, fixed_x3), ncol = 3)
  
  
  a0 = 1.5 + g0
  a1 = 1 + g1
  a2 = 0.5 + g2
  a3 = 0.2 + g3
  a4 = 0.2 + g4
  a5 = 0.1 + g5
  
  r = data.frame(g0, g1, g2, g3, g4, g5, g6, g7, a0, a1, a2, a3, a4, a5)
  
  
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
    num.obs = rpois(n, 5) + 1
  }
  #sample random obs times from tv and pick up "true" w(t) values for each i,
  #truncated at event/censoring time
  id_long = NULL
  wt_samp = NULL
  wt_i_samp2 = NULL
  obs_time_samp = NULL
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
    
    e_i1 = rnorm(length(obs.times), 0, 0.02)
    #wt_tia1 = (r$g0[i] + 0.2) + (0.5 + r$g5[i])* obs.times +  e_i1
    wt_tia1 =  -0.6 + (1-(0.75)*exp(4*obs.times)/(1 + exp(4*obs.times)) ) + r$g0[i] + r$g5[i]*obs.times +  e_i1
    wt_samp = c(wt_samp, wt_tia1)
  
    
    e_i2 = rnorm(length(obs.times), 0, 0.05)
    wt_tia2 = (r$g6[i] + 0.1) + (0.2 + r$g7[i])* obs.times - 0.2 * obs.times^2 + fixed_x3[i]*0.2 + e_i2
    wt_i_samp2 = c(wt_i_samp2, wt_tia2)
    
    
    visit_num = c(visit_num, c(1:length(wt_tia1)))
    W_long = c(W_long, rep(fixed_x3[i],length(wt_tia1)))
    before_midpoint = c(before_midpoint, as.numeric(obs.times < dat.baseline$mid_point[i]))
    
  }
  
  
  dat.long = data.frame(cbind(id_long, visit_num, wt_samp, wt_i_samp2, obs_time_samp, W_long, before_midpoint))
  
  dat.baseline$n.obs = (dat.long %>% group_by(id_long) %>% tally())$n
  
  out = list(dat.baseline = dat.baseline, dat.long = dat.long)
  
  return(out)
  
}

