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

hazard_function_compare = function(event_time, x, r){
  
  #kappa = 4
  #alpha = 1
  
  #h0_t = 3*event_time^2
  #h0_t = ((kappa/alpha) * (event_time/alpha)^(kappa-1))/(1 + (event_time/alpha)^kappa)
  h0_t = 0.5 * exp(2 * event_time)
  
  w_yi = (-0.1 + r$g0) - (0.1 + r$g1)*event_time - x[3]*0.3 
  
  true_beta1 = -1
  true_gamma = -0.3
  
  ht = h0_t * (exp(c(true_beta1*x[1] + true_gamma*w_yi)))
  return(ht)
  
}

interval_censoring <- function(time, tmax, n){
  #observation_times = cumsum(runif(10, 0.1, 0.8))
  observation_times = seq(0.25, 3, 0.25)
  #observation_times = seq(0.05, 3, 0.05)
  right_cens_i = runif(n, 0.5, 2.5)
  
  events = right = left = interval = rep(0, n)
  t_L = t_R = rep(0, n)
  
  for(i in 1:n){
    #observation_times = cumsum(runif(10, 0.2, 0.5))
    observation_times_i = observation_times[which(observation_times < right_cens_i[i])]
    #observation_times_i = c(observation_times_i, right_cens_i[i])
    if(time[i] < min(observation_times_i)){
      #left censoring
      left[i] = 1
      t_R[i] = min(observation_times_i)
    }else if(time[i] > max(observation_times_i)){
      #right censoring
      right[i] = 1
      t_L[i] = max(observation_times_i)
      t_R[i] = Inf
    }else{
      interval[i] = 1
      t_L[i] = max(observation_times_i[which(observation_times_i < time[i])])
      t_R[i] = min(observation_times_i[which(observation_times_i > time[i])])
    }
  }
  
  return(cbind(time, t_L, t_R, events, interval, right, left, right_cens_i))
}

gen_interval_censoring_jointmodel = function(n){
  
  g0 = rnorm(n, 0, 0.2)
  g1 = rnorm(n, 0, 0.4)
  
  r = data.frame(g0, g1)
  
  fixed_x1 = rnorm(n, 0, 1)
  #fixed_x1 = fixed_x1 - mean(fixed_x1)
  fixed_x2 = rnorm(n, 0, 1)
  #fixed_x2 = fixed_x2 - mean(fixed_x2)
  fixed_x3 = rbinom(n, 1, 0.5)
  #fixed_x3 = fixed_x3 - mean(fixed_x3)
  fixed_x = as.matrix(cbind(fixed_x1, fixed_x2, fixed_x3), ncol = 3)
  
  dat.baseline = rJM(hazard_function_compare, interval_censoring, x = fixed_x, r = r, tmin = 0, tmax = 3)
  
  #for each i, generate a random sequence of observation times for w(t)
  #num.obs = rpois(n, 19) + 1
  
  #sample random obs times from tv and pick up "true" w(t) values for each i,
  #truncated at event/censoring time
  id_long = NULL
  wt_samp = NULL
  wt_i_samp = NULL
  obs_time_samp = NULL
  e_i_t = NULL
  visit_num = NULL
  W_long = NULL
  
  for(i in 1:n){
    max.obs.time = ifelse(dat.baseline$right[i] == 1, dat.baseline$t_L[i], dat.baseline$t_R[i])
    obs.times = seq(0, max.obs.time, by = 0.25)
    
    id_long = c(id_long, rep(i, length(obs.times)))
    obs_time_samp = c(obs_time_samp, obs.times)
    
    e_i = rnorm(length(obs.times), 0, 0.1)
    e_i_t = c(e_i_t, e_i)
    
    wt_tia = (-0.1 + r$g0[i]) -
      (0.1 + r$g1[i])*obs.times -
      0.3*fixed_x3[i] + e_i
    
    wt_samp = c(wt_samp, wt_tia)
    wt_i_samp = c(wt_i_samp, wt_tia)
    visit_num = c(visit_num, c(1:length(wt_tia)))
    W_long = c(W_long, rep(fixed_x3[i],length(wt_tia)))
    
    
  }
  
  
  dat.long = data.frame(cbind(visit.ID = id_long, visit.num = visit_num, visit.time = obs_time_samp, W = wt_i_samp, Z.1 = W_long))
  
  dat.baseline$n.obs = (dat.long %>% group_by(visit.ID) %>% tally())$n
  
  out = list(dat.baseline = dat.baseline, dat.long = dat.long)
  
  return(out)
  
}


dat = gen_interval_censoring_jointmodel(n=200)
