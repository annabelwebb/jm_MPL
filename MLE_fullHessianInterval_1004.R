jm_fitting_interval_MLE_new = function(dat.baseline, dat.long, bmm, 
                                        max.iter = c(10, 1000), 
                                        n.knots.h0 = 2, step_size = 1, mod_mat_long_f,
                                        re_ind, range, init_re, init_sigma_eps, init_sigma_re){
  
  
  kappa = 1/0.5
  
  event.time = dat.baseline$latest_time
  obs.time = dat.long$obs_time_samp
  id = dat.baseline$id
  id_long = dat.long$id_long
  cont = dat.long$wt_samp
  fixed = as.matrix(bmm)
  W_long = (dat.long$W_long)
  
  fixed_e = (fixed[which(dat.baseline$event == 1),])
  fixed_r = (fixed[which(dat.baseline$right == 1),])
  fixed_l = (fixed[which(dat.baseline$left == 1),])
  fixed_i = (fixed[which(dat.baseline$interval == 1),])
  
  #initial set up
  
  #m-splines for baseline hazard function
  int.knots.event = quantile(c(dat.baseline$mid_point)[(c(dat.baseline$mid_point) > 0)], seq(range[1], range[2], length.out=n.knots.h0))
  bound.knots.event = c(- 1e-4, max(c(obs.time, event.time)[(c(obs.time, event.time) > 0)]) + 1e-4)
  event.kn = list(int = int.knots.event, bound = bound.knots.event)
  
  psi_e = mSpline(dat.baseline$t_L, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[dat.baseline$event == 1,]
  psi_r = mSpline(dat.baseline$t_L[dat.baseline$right == 1], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)
  psi_l = mSpline(dat.baseline$t_L, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[dat.baseline$left == 1,]
  psi_i1 = mSpline(dat.baseline$t_L, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[dat.baseline$interval == 1]
  psi_i2 = mSpline(dat.baseline$t_R[dat.baseline$interval == 1], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)
  
  tm = max(c(obs.time, event.time)[(c(obs.time, event.time) > 0)])
  #b-splines for TVC function
  phi.mu = mod_mat_long_f(obs.time, tmax = tm, W_long = W_long)
  
  #b-splines for TVC function
  phi = mod_mat_long_f(event.time, tmax = tm, W_long =  dat.baseline$fixed_x3)
  phi_e = phi[dat.baseline$event == 1,]
  phi_r = phi[dat.baseline$right == 1,]
  phi_l = phi[dat.baseline$left == 1,]
  phi_i2 = phi[dat.baseline$interval == 1,]
  
  phi_i1 = mod_mat_long_f(dat.baseline$t_L[which(dat.baseline$interval == 1)], 
                          tmax = tm, W_long =  dat.baseline$fixed_x3[which(dat.baseline$interval == 1)])
  
  
  #set up dimensions
  p = ncol(fixed)
  q = 1
  m = length(int.knots.event) + 3
  r = ncol(phi.mu)
  #v = 1
  w = length(re_ind)
  n = length(id)
  n_long = length(id_long)
  
  #initalise parameter vectors
  beta = matrix(rep(0, p), ncol = 1)
  gamma = matrix(0)
  theta = matrix(rep(1,m), ncol = 1)
  alpha = matrix(rep(0.1,r), ncol = 1)
  a_re_long = init_re
  a_re_cols = matrix(a_re_long, ncol = w)
  a_re_cols_pad = matrix(rep(0,r*n), ncol = r)
  a_re_cols_pad_quad = matrix(rep(0,r*n*15), ncol = r)
  a_re_cols_pad[,re_ind] = as.matrix(a_re_cols)
  a_re_cols_pad_quad[,re_ind] = matrix(sapply(a_re_cols, rep, 15), ncol = 2)
  
  i1_ind_quad = c(sapply(dat.baseline$interval, rep, 15))
  
  #variance components
  sigma2_Et = init_sigma_eps
  sigma2_re = init_sigma_re
  
  df.theta = n
  df.alpha = n
  df.epsilon = n_long
  df.a_re = rep(n, w)
  
  #penalty matrices and smoothing parameters
  theta.G = theta_penalty_f(3, int.knots.event, bound.knots.event)
  #alpha.G = zeta_penalty_f(3, int.knots.mu, bound.knots.mu)
  alpha.G = matrix(0, r, r)
  
  theta.lambda = 0
  alpha.lambda = 0
  
  #quadrature rules etc.
  quad.lambda = (event.time - 0)/2
  quad.mu = (event.time + 0)/2
  quad.y = t(as.matrix(quad.lambda) %*% rules[[15]]$x + quad.mu) #one row per time, one column per i
  quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event))
  quad.phi.event = mod_mat_long_f(c(quad.y), tmax = tm,W_long =  c(sapply(dat.baseline$fixed_x3, rep, 15)))
  quad.w = rules[[15]]$w
  
  quad.lambda_i1 = (dat.baseline$t_L[dat.baseline$interval == 1] - 0)/2
  quad.mu_i1 = (dat.baseline$t_L[dat.baseline$interval == 1] + 0)/2
  quad.y_i1 = t(as.matrix(quad.lambda_i1) %*% rules[[15]]$x + quad.mu_i1) #one row per time, one column per i
  quad.psi.y_i1 = t(sapply(quad.y_i1, mSpline, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event))
  quad.phi.y_i1 = mod_mat_long_f(c(quad.y_i1), tmax = tm,W_long =  c(sapply(dat.baseline$fixed_x3[dat.baseline$interval == 1], rep, 15)))
  i1_ind_quad = c(sapply(dat.baseline$interval, rep, 15))
  
  quad.lambda_i2 = (dat.baseline$t_R[dat.baseline$interval == 1] - dat.baseline$t_L[dat.baseline$interval == 1])/2
  quad.mu_i2 = (dat.baseline$t_L[dat.baseline$interval == 1] + dat.baseline$t_R[dat.baseline$interval == 1])/2
  quad.y_i2 = t(as.matrix(quad.lambda_i2) %*% rules[[15]]$x + quad.mu_i2) #one row per time, one column per i
  quad.psi.y_i2 = t(sapply(quad.y_i2, mSpline, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event))
  quad.phi.y_i2 = mod_mat_long_f(c(quad.y_i2), tmax = tm,W_long =  c(sapply(dat.baseline$fixed_x3[dat.baseline$interval == 1], rep, 15)))
  
  
  
  for(it in 1:max.iter[1]){
    for(iter in 1:max.iter[2]){
      
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      H0_t_l = H0_t[dat.baseline$left == 1]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      pl_l = sum(log(1 - S_t_l))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      pl_i = sum(log(S_t_i1 - S_t_i2 ))
      
      #least squares
      z_t_ia = NULL
      for(i in 1:n){
        alpha_re = alpha
        alpha_re[re_ind] = alpha_re[re_ind] + a_re_cols[i,re_ind]
        ind.long = which(id_long == i)
        z_t_ia = c(z_t_ia, phi.mu[ind.long,] %*% (alpha_re))
      }
      
      pl_ls = (1/(2*sigma2_Et))*sum((cont - z_t_ia)^2)
      
      #log likelihood
      log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #alpha
      B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event # h0(t) * exp(zTg) * phi(t) for each i at each quad t, (n*15) x r
      df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* B_t_quad) #multiply above with quad.w (weight) for each t, then each row is labelled with i, (n*15) x (r+1)
      B_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1]) # aggregate() sums down each column for each i (removing i label column), then multiply each i with its own quad.lamba
      
      Bphi_t_e = B_t[dat.baseline$event == 1,]
      Bphi_t_r = B_t[dat.baseline$right == 1,]
      Bphi_t_l = B_t[dat.baseline$left == 1,]
      Bphi_t_i2 = B_t[dat.baseline$interval == 1,]
      
      B_t_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * quad.phi.y_i1 #same as above but only for tL for interval censored indiviudals
      df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* B_t_quad_i1)
      Bphi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
      
      alpha_score_e = apply(as.numeric(gamma) * phi_e - as.numeric(gamma) * as.numeric(eXtB_e) * Bphi_t_e, 2, sum)
      alpha_score_r = apply(- as.numeric(gamma) * as.numeric(eXtB_r) * Bphi_t_r, 2, sum)
      alpha_score_l = apply(as.numeric(gamma) * as.numeric(eXtB_l * S_t_l/(1-S_t_l)) * Bphi_t_l, 2, sum)
      alpha_score_i = apply(- as.numeric(gamma) * as.numeric(eXtB_i * S_t_i1 /(S_t_i1-S_t_i2)) *  Bphi_t_i1 +  
                              as.numeric(gamma) * as.numeric(eXtB_i * S_t_i2 /(S_t_i1-S_t_i2)) *  Bphi_t_i2, 2, sum)
      
      alpha_score_ls = as.numeric( c(1/sigma2_Et) * t(as.matrix((cont - z_t_ia))) %*% phi.mu)
      
      alpha_score = alpha_score_e + alpha_score_r + alpha_score_l + alpha_score_i + alpha_score_ls - 2*alpha.lambda * alpha.G %*% alpha
      
      #B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
      B_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, 15)) * B_t_quad
      B_t_quad2_e = B_t_quad2[c(sapply(dat.baseline$event == 1, rep, 15)),]
      B_t_quad2_r = B_t_quad2[c(sapply(dat.baseline$right == 1, rep, 15)),]
      B_t_quad2_l = B_t_quad2[c(sapply(dat.baseline$left == 1, rep, 15)),]
      B_t_quad2_i2 = B_t_quad2[c(sapply(dat.baseline$interval == 1, rep, 15)),]
      B_t_quad2_i1 = c(rep(quad.w, sum(dat.baseline$interval == 1))) * 
        c(sapply(quad.lambda_i1, rep, 15)) * B_t_quad_i1 
      
      if(sum(dat.baseline$event)>0){
        alpha_hess_e = t(-c(sapply(c(gamma^2) * c(eXtB_e), rep, 15)) * B_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, 15)),]
      }else{
        alpha_hess_e = 0
      }
      
      alpha_hess_r = t(-c(sapply(c(gamma^2) * c(eXtB_r), rep, 15)) * B_t_quad2_r) %*% quad.phi.event[c(sapply(dat.baseline$right == 1, rep, 15)),]
      alpha_hess_l = t(c(sapply(c(gamma^2) * c(eXtB_l*(S_t_l)/(1-S_t_l)), rep, 15)) * B_t_quad2_l) %*% quad.phi.event[c(sapply(dat.baseline$left == 1, rep, 15)),] -
        t(c(c(gamma^2) * (eXtB_l^2 * S_t_l/(1-S_t_l)^2)) * Bphi_t_l) %*% Bphi_t_l
      alpha_hess_i = t(-c(sapply(c(gamma^2) * c(eXtB_i*S_t_i1 /(S_t_i1-S_t_i2)), rep, 15)) * B_t_quad2_i1) %*% quad.phi.y_i1 + 
        t(c(sapply(c(gamma^2) * c(eXtB_i*S_t_i2 /(S_t_i1-S_t_i2)), rep, 15)) * B_t_quad2_i2) %*% quad.phi.event[c(sapply(dat.baseline$interval == 1, rep, 15)),]-
        t(c(c(gamma^2) * (eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i1 + 
        t(c(c(gamma^2) * (eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i2
      
      alpha_hess_partial_new = alpha_hess_e + alpha_hess_r + alpha_hess_l + alpha_hess_i
      
      alpha_hess = alpha_hess_partial_new  - c(1/sigma2_Et) * t(phi.mu) %*% phi.mu - 2*alpha.lambda * alpha.G
      alpha_hess_neg = -alpha_hess
      
      alpha_old = alpha
      alpha = alpha_old + solve(alpha_hess_neg)%*%alpha_score
      
      #step size for alpha
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      H0_t_l = H0_t[dat.baseline$left == 1]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      pl_l = sum(log(1 - S_t_l))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      pl_i = sum(log(S_t_i1 - S_t_i2))
      
      #least squares
      z_t_ia = NULL
      for(i in 1:n){
        alpha_re = alpha
        alpha_re[re_ind] = alpha_re[re_ind] + a_re_cols[i,re_ind]
        ind.long = which(id_long == i)
        z_t_ia = c(z_t_ia, phi.mu[ind.long,] %*% (alpha_re))
      }
      
      pl_ls = (1/(2*sigma2_Et))*sum((cont - z_t_ia)^2)
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #print(c("alpha", log_lik_old, log_lik))
      
      if(((log_lik < log_lik_old) & step_size == 1)){
        ii = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          alpha = alpha_old + omega1 * solve(alpha_hess_neg)%*%alpha_score
          
          #update log-likelihood
          #likelihood
          eXtB = exp(fixed %*% beta)
          
          h0_t_quad = quad.psi.event %*% theta
          exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          H0_t_l = H0_t[dat.baseline$left == 1]
          #H0_t_i2 = H0_t[dat.baseline$interval == 1]
          
          h0_t_quad_i1 = quad.psi.y_i1 %*% theta
          exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
          h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
          h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
          
          h0_t_quad_i2 = quad.psi.y_i2 %*% theta
          exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
          h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
          h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #left
          eXtB_l = exp(fixed_l %*% beta)
          S_t_l = exp(-eXtB_l * H0_t_l)
          pl_l = sum(log(1 - S_t_l ))
          
          #interval
          eXtB_i = exp(fixed_i %*% beta)
          S_t_i1 = exp(-eXtB_i * H0_t_i1)
          S_t_i2 = exp(-eXtB_i * H0_t_i2)
          pl_i = sum(log(S_t_i1 - S_t_i2 ))
          
          #least squares
          z_t_ia = NULL
          for(i in 1:n){
            alpha_re = alpha
            alpha_re[re_ind] = alpha_re[re_ind] + a_re_cols[i,re_ind]
            ind.long = which(id_long == i)
            z_t_ia = c(z_t_ia, phi.mu[ind.long,] %*% (alpha_re) ) 
          }
          
          pl_ls = (1/(2*sigma2_Et))*sum((cont - z_t_ia)^2)
          
          #log likelihood
          log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
            sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
            theta.lambda* t(theta)%*%theta.G%*%theta - 
            alpha.lambda* t(alpha)%*%alpha.G%*%alpha
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          ii = ii+1
          if(ii>500){break}
          
          
        }
      }
      
      #gamma
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      A_t_quad = h0_t_quad * exp_zTg_quad * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
      A_t_quad = matrix(A_t_quad, ncol = n, nrow = 15, byrow = FALSE)
      A_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  A_t_quad, 2, sum)
      
      A_t_e = A_t[dat.baseline$event == 1]
      A_t_r = A_t[dat.baseline$right == 1]
      A_t_l = A_t[dat.baseline$left == 1]
      A_t_i2 = A_t[dat.baseline$interval == 1]
      
      exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
      A_t_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1 * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)
      A_t_quad_i1 = matrix(A_t_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      A_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  A_t_quad_i1, 2, sum)
      
      gamma_score_e = sum(zt_e - eXtB_e * A_t_e)
      gamma_score_r = sum(- eXtB_r * A_t_r)
      gamma_score_l = sum(S_t_l * eXtB_l * A_t_l / (1 - S_t_l))
      gamma_score_i = sum(- eXtB_i * (S_t_i1 * A_t_i1 - S_t_i2 * A_t_i2) / (S_t_i1 - S_t_i2 ))
      
      gamma_score = gamma_score_e + gamma_score_r + gamma_score_l + gamma_score_i
      
      #second derivative of gamma
      z_t_quad = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
      A2_t_quad_long = h0_t_quad * exp_zTg_quad * (z_t_quad)^2
      A_t_quad2 = matrix(A2_t_quad_long, ncol = n, nrow = 15, byrow = FALSE)
      A2_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  A_t_quad2, 2, sum)
      
      A_t_quad2_e = A2_t[dat.baseline$event == 1]
      A_t_quad2_r = A2_t[dat.baseline$right == 1]
      A_t_quad2_l = A2_t[dat.baseline$left == 1]
      A_t_quad2_i2 = A2_t[dat.baseline$interval == 1]
      
      z_t_quad_i1 = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)
      A_t_quad_i1_long = h0_t_quad_i1 * exp_zTg_quad_i1 * (z_t_quad_i1)^2
      A_t_quad2_i1 = matrix(A_t_quad_i1_long, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      A_t_quad2_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  A_t_quad2_i1, 2, sum)
      
      gamma_hess_e = t(-eXtB_e) %*% A_t_quad2_e
      gamma_hess_r = t(-eXtB_r) %*% A_t_quad2_r
      #gamma_hess_l = t(S_t_l * eXtB_l/(1-S_t_l)) %*% A_t_quad2_l - 
      #  t(c( S_t_l * eXtB_l^2/(1-S_t_l)^2) * A_t_l) %*% A_t_l 
      #gamma_hess_i = t(-eXtB_i *S_t_i1 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i1 +
      #  t(eXtB_i *S_t_i2 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i2 -
      #  t(c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2 ) * (A_t_i1 - A_t_i2)) %*% (A_t_i1 - A_t_i2)
      
      
      gamma_hess_l = - t(c( S_t_l * eXtB_l^2/(1-S_t_l)^2) * A_t_l) %*% A_t_l 
      gamma_hess_i = - t(eXtB_i *S_t_i1 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i1 -
        t(c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2 ) * (A_t_i1 - A_t_i2)) %*% (A_t_i1 - A_t_i2)
      
      
      # t(eXtB_i *S_t_i2 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i2
      
      gamma_hess = gamma_hess_e + gamma_hess_r + gamma_hess_l + gamma_hess_i
      gamma_hess_neg = -gamma_hess
      
      gamma_old = gamma
      gamma = gamma_old + solve(gamma_hess_neg)%*%gamma_score
      
      #step size for gamma
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      H0_t_l = H0_t[dat.baseline$left == 1]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      pl_l = sum(log(1 - S_t_l ))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      pl_i = sum(log(S_t_i1 - S_t_i2 ))
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      step_size_gamma = 1
      if(((log_lik < log_lik_old) & step_size_gamma == 1)){
        ii = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          gamma = gamma_old +  omega1 * solve(gamma_hess_neg)%*%gamma_score
          
          #likelihood
          eXtB = exp(fixed %*% beta)
          
          h0_t_quad = quad.psi.event %*% theta
          exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          H0_t_l = H0_t[dat.baseline$left == 1]
          #H0_t_i2 = H0_t[dat.baseline$interval == 1]
          
          h0_t_quad_i1 = quad.psi.y_i1 %*% theta
          exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
          h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
          h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
          
          h0_t_quad_i2 = quad.psi.y_i2 %*% theta
          exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
          h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
          h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #left
          eXtB_l = exp(fixed_l %*% beta)
          S_t_l = exp(-eXtB_l * H0_t_l)
          pl_l = sum(log(1 - S_t_l ))
          
          #interval
          eXtB_i = exp(fixed_i %*% beta)
          S_t_i1 = exp(-eXtB_i * H0_t_i1)
          S_t_i2 = exp(-eXtB_i * H0_t_i2)
          pl_i = sum(log(S_t_i1 - S_t_i2 ))
          
          #log likelihood
          log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
            sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
            theta.lambda* t(theta)%*%theta.G%*%theta - 
            alpha.lambda* t(alpha)%*%alpha.G%*%alpha
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          ii = ii+1
          if(ii>1000){break}
          
          
        }
      }
      
      
      #theta
      psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
      df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* psi_t_star_quad)
      Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
      
      Psi_t_e = Psi_t[dat.baseline$event == 1,]
      Psi_t_r = Psi_t[dat.baseline$right == 1,]
      Psi_t_l = Psi_t[dat.baseline$left == 1,]
      #Psi_t_i2 = Psi_t[dat.baseline$interval == 1,]
      
      psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
      df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* psi_t_star_quad_i1)
      Psi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
      
      psi_t_star_quad_i2 = c(exp_zTg_quad_i2) * quad.psi.y_i2
      df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* psi_t_star_quad_i2)
      Psi_t_i2 = Psi_t_i1 + c(quad.lambda_i2) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
      
      
      TwoLRtheta = as.numeric(2*theta.lambda)*(theta.G%*%theta)
      
      theta_score_neg = apply(as.numeric(eXtB_e)* Psi_t_e, 2, sum) + 
        apply(as.numeric(eXtB_r)* Psi_t_r, 2, sum) +
        apply(as.numeric(eXtB_i * S_t_i1/(S_t_i1-S_t_i2))* Psi_t_i1, 2, sum) +
        TwoLRtheta*(TwoLRtheta>0) + 1e-3
      
      theta_score_pos = apply(as.numeric(1/h0_t_e)* psi_e, 2, sum) +
        apply(as.numeric(eXtB_i * S_t_i2/(S_t_i1-S_t_i2))* Psi_t_i2, 2, sum) +
        apply(as.numeric(eXtB_l * S_t_l/(1-S_t_l))* Psi_t_l, 2, sum) - 
        TwoLRtheta*(TwoLRtheta<0) + 1e-3
      
      theta_score = (theta_score_pos - theta_score_neg)
      
      D_matrix = theta/(theta_score_neg)
      
      theta_old = theta
      theta = theta_old + D_matrix*(theta_score)
      
      
      #step size for theta
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      H0_t_l = H0_t[dat.baseline$left == 1]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(quad.w *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      pl_l = sum(log(1 - S_t_l ))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      pl_i = sum(log(S_t_i1 - S_t_i2 ))
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      
      step_size_theta = 1
      
      if(((log_lik < log_lik_old) & step_size_theta == 1)){
        i = 0
        omega1 = 0.5
        
        
        while(log_lik<log_lik_old){
          theta = theta_old +  omega1 * D_matrix*(theta_score)
          
          #likelihood
          eXtB = exp(fixed %*% beta)
          
          h0_t_quad = quad.psi.event %*% theta
          exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          H0_t_l = H0_t[dat.baseline$left == 1]
          #H0_t_i2 = H0_t[dat.baseline$interval == 1]
          
          h0_t_quad_i1 = quad.psi.y_i1 %*% theta
          exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
          h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
          h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i1 = quad.lambda_i1 * apply(quad.w *  h0_t_star_quad_i1, 2, sum)
          
          h0_t_quad_i2 = quad.psi.y_i2 %*% theta
          exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
          h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
          h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #left
          eXtB_l = exp(fixed_l %*% beta)
          S_t_l = exp(-eXtB_l * H0_t_l)
          pl_l = sum(log(1 - S_t_l))
          
          #interval
          eXtB_i = exp(fixed_i %*% beta)
          S_t_i1 = exp(-eXtB_i * H0_t_i1)
          S_t_i2 = exp(-eXtB_i * H0_t_i2)
          pl_i = sum(log(S_t_i1 - S_t_i2 ))
          
          #log likelihood
          log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
            sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
            theta.lambda* t(theta)%*%theta.G%*%theta - 
            alpha.lambda* t(alpha)%*%alpha.G%*%alpha
          
          #print(c(log_lik_old, log_lik, omega1))
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>500){break}
          
          
        }
      }
      
      
      
      #beta
      
      if(sum(dat.baseline$event) > 0){
        beta_score_e = t(fixed_e) %*% c(1- eXtB_e * H0_t_e)
      }else{
        beta_score_e = 0
      }
      beta_score_r = - t(fixed_r) %*% c(eXtB_r * H0_t_r)
      beta_score_l = t(fixed_l) %*% c(S_t_l * eXtB_l * H0_t_l / (1 - S_t_l))
      beta_score_i = - t(fixed_i) %*% c((S_t_i1 * eXtB_i * H0_t_i1 - S_t_i2 * eXtB_i * H0_t_i2) / (S_t_i1 - S_t_i2))
      
      beta_score = beta_score_e + beta_score_r + beta_score_l + beta_score_i
      
      beta_hess_e = t(fixed_e) %*% diag(c(-eXtB_e * H0_t_e)) %*% fixed_e
      beta_hess_r = t(fixed_r) %*% diag(c(-eXtB_r * H0_t_r)) %*% fixed_r
      beta_hess_l = t(fixed_l) %*% diag(c((eXtB_l*H0_t_l*S_t_l - ((eXtB_l*H0_t_l)^2)*S_t_l) /(1 - S_t_l))) %*% fixed_l -
        t(fixed_l) %*% diag(c(((eXtB_l*H0_t_l*S_t_l)^2)/(1 - S_t_l)^2))%*% fixed_l
      beta_hess_i = -t(fixed_i) %*% diag(c((S_t_i1 * eXtB_i * H0_t_i1 ) / (S_t_i1 - S_t_i2))) %*% fixed_i + 
        t(fixed_i) %*% diag(c((S_t_i2 * eXtB_i * H0_t_i2 ) / (S_t_i1 - S_t_i2))) %*% fixed_i -
        t(fixed_i) %*% diag(c(S_t_i1 * S_t_i2 * (eXtB_i * H0_t_i1 - eXtB_i *  H0_t_i2)^2/ (S_t_i1 - S_t_i2)^2)) %*% fixed_i
      
      
      beta_hess = beta_hess_e + beta_hess_r + beta_hess_l + beta_hess_i
      beta_neg_hess = -beta_hess
      
      beta_old = beta
      
      beta = beta_old + solve(beta_neg_hess)%*%beta_score
      
      #step size for beta
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      pl_l = sum(log(1 - S_t_l))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      pl_i = sum(log(S_t_i1 - S_t_i2))
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #print(c("beta", log_lik_old, log_lik))
      step_size_beta =1
      if(((log_lik < log_lik_old) & step_size_beta == 1)){
        i = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          beta = beta_old +  omega1 * solve(beta_neg_hess)%*%beta_score
          
          #likelihood
          eXtB = exp(fixed %*% beta)
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #left
          eXtB_l = exp(fixed_l %*% beta)
          S_t_l = exp(-eXtB_l * H0_t_l)
          pl_l = sum(log(1 - S_t_l))
          
          #interval
          eXtB_i = exp(fixed_i %*% beta)
          S_t_i1 = exp(-eXtB_i * H0_t_i1)
          S_t_i2 = exp(-eXtB_i * H0_t_i2)
          pl_i = sum(log(S_t_i1 - S_t_i2))
          
          #log likelihood
          log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
            sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
            theta.lambda* t(theta)%*%theta.G%*%theta - 
            alpha.lambda* t(alpha)%*%alpha.G%*%alpha
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>500){break}
          
          
        }
        #print(c(omega1))
        
      }
      
      a_re_cols_old = a_re_cols
      
      #random effects
      
      B_t_re_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event[,1:w]
      df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* B_t_re_quad)
      B_t_re = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
      
      Bphi_t_re_e = B_t_re[dat.baseline$event == 1,]
      Bphi_t_re_r = B_t_re[dat.baseline$right == 1,]
      Bphi_t_re_l = B_t_re[dat.baseline$left == 1,]
      Bphi_t_re_i2 = B_t_re[dat.baseline$interval == 1,]
      
      B_t_re_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * quad.phi.y_i1[,1:w]
      df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), 
                      c(rep(quad.w, sum(dat.baseline$interval == 1)))* B_t_re_quad_i1)
      B_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
      
      a1i_score_e = as.numeric(gamma) * phi_e[,1:w] - as.numeric(gamma) * as.numeric(eXtB_e) * (Bphi_t_re_e)
      a1i_score_r = - as.numeric(gamma) * as.numeric(eXtB_r) * (Bphi_t_re_r)
      a1i_score_l = as.numeric(gamma) * as.numeric(eXtB_l * S_t_l/(1-S_t_l)) * (Bphi_t_re_l)
      a1i_score_i = - as.numeric(gamma) * as.numeric(eXtB_i * S_t_i1/(S_t_i1-S_t_i2)) * (B_t_re_i1) + 
        as.numeric(gamma) * as.numeric(eXtB_i * S_t_i2/(S_t_i1-S_t_i2)) * (Bphi_t_re_i2)
      
      id_censor_order = c(id[dat.baseline$event == 1], id[dat.baseline$right == 1], id[dat.baseline$left == 1], id[dat.baseline$interval == 1])
      a1i_score_s = rbind(as.matrix(a1i_score_e), as.matrix(a1i_score_r), as.matrix(a1i_score_l), as.matrix(a1i_score_i))
      a1i_score_s = data.frame(id_censor_order, a1i_score_s)[order(id_censor_order),-1]
      
      phi.mu.1 = phi.mu[,1:w]
      a1i_score_ls = data.frame(id_long, (cont-z_t_ia)*phi.mu.1)
      a1i_score_ls = as.matrix(aggregate( a1i_score_ls[,2:(w+1)], list(a1i_score_ls[,1]), FUN = sum )[,-1])
      
      a1i_score = a1i_score_s + (1/sigma2_Et)*a1i_score_ls - a_re_cols/matrix(rep(sigma2_re, n), ncol = w, byrow = T)
      
      a_cols_update = matrix(0, nrow = n, ncol = w)
      
      S_t = exp(-H0_t*eXtB)
      
      S_t_i1_long = rep(0, length(S_t))
      S_t_i1_long[which(dat.baseline$interval == 1)] = S_t_i1
      
      B_t_re_quad_i1_long = matrix(0, nrow = nrow(as.matrix(B_t_re_quad)), ncol = ncol(as.matrix(B_t_re_quad)))
      B_t_re_quad_i1_long[i1_ind_quad,] = B_t_re_quad_i1
      
      B_t_re_i1_long = matrix(0, nrow = nrow(B_t_re), ncol = ncol(B_t_re))
      B_t_re_i1_long[which(dat.baseline$interval == 1),] = B_t_re_i1
      
      quad.phi.y_i1_long = matrix(0, nrow = nrow(quad.phi.event), ncol = ncol(quad.phi.event))
      quad.phi.y_i1_long[i1_ind_quad,] = quad.phi.y_i1
      
      phi.mu.1 = as.matrix(phi.mu.1)
      #quad.phi.event = as.matrix(quad.phi.event)
      B_t_re_quad = as.matrix(B_t_re_quad)
      quad.lambda_i1_long = rep(0, n)
      quad.lambda_i1_long[which(dat.baseline$interval == 1)] = quad.lambda_i1
      
      
      for(i in 1:n){
      
        ind.long = which(id_long == i)
        a_hess_i = t(c(dat.baseline$event[i] + dat.baseline$right[i]) * rep(- as.numeric(gamma^2) * as.numeric(eXtB[i] ), 15) * B_t_re_quad[(15*(i-1)+1):(15*(i-1)+15),]) %*%  as.matrix(quad.lambda[i] * quad.w * quad.phi.event[(15*(i-1)+1):(15*(i-1)+15),1:w]) +
          t(c(dat.baseline$left[i]) *  rep(as.numeric(gamma^2) * as.numeric(eXtB[i] * (S_t[i]/(1-S_t[i])) ), 15) * B_t_re_quad[(15*(i-1)+1):(15*(i-1)+15),]) %*%  (quad.lambda[i] * quad.w * quad.phi.event[(15*(i-1)+1):(15*(i-1)+15),1:w]) - 
          ((c(dat.baseline$left[i]) *  as.numeric(gamma^2) * as.numeric((eXtB[i]^2) * (S_t[i]/(1-S_t[i])^2) ))* B_t_re[i,])  %*%  t(B_t_re[i,]) - 
          t(c(dat.baseline$interval[i]) *  rep(as.numeric(gamma^2) * as.numeric(eXtB[i] * (S_t_i1_long[i]/(S_t_i1_long[i]-S_t[i])) ), 15) * B_t_re_quad_i1_long[(15*(i-1)+1):(15*(i-1)+15),]) %*%  (quad.lambda_i1_long[i] * quad.w * quad.phi.y_i1_long[(15*(i-1)+1):(15*(i-1)+15),1:w]) + 
          t(c(dat.baseline$interval[i]) *  rep(as.numeric(gamma^2) * as.numeric(eXtB[i] * (S_t[i]/(S_t_i1_long[i]-S_t[i])) ), 15) * B_t_re_quad[(15*(i-1)+1):(15*(i-1)+15),]) %*%  (quad.lambda[i] * quad.w * quad.phi.event[(15*(i-1)+1):(15*(i-1)+15),1:w]) -
          (c(dat.baseline$interval[i]) *  as.numeric(gamma^2) * as.numeric((eXtB[i]^2) * (S_t[i]*S_t_i1_long[i]/(S_t_i1_long[i]-S_t[i])^2) )* (B_t_re_i1_long[i,] - B_t_re[i,]))  %*%  t(B_t_re_i1_long[i,] - B_t_re[i,]) -
          t((1/sigma2_Et) * phi.mu.1[ind.long,]) %*% (phi.mu.1[ind.long,]) - 
         diag(c(1/sigma2_re))
        
        
        a_cols_update[i,] = solve(-a_hess_i)%*%t(as.matrix(a1i_score[i,]))
    #    a_cols_update[i,] = solve(-a_hess_i)%*%t(as.matrix(a1i_score[i,]))
        
      }
      
      a_re_cols = a_re_cols_old + a_cols_update
      #a_re_cols = a_re_cols_old + as.matrix(a1i_score)
      
      
      a_re_long_old = a_re_long
      a_re_cols_pad_old = a_re_cols_pad
      a_re_cols_pad_quad_old = a_re_cols_pad_quad
      
      a_re_long = unlist(c(a_re_cols))
      a_re_cols_pad[,re_ind] = as.matrix(a_re_cols)
      a_re_cols_pad_quad[,re_ind] = matrix(sapply(a_re_cols, rep, 15), ncol = 2)
      
      
      #step size for random effects
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      H0_t_l = H0_t[dat.baseline$left == 1]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(quad.w *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      pl_l = sum(log(1 - S_t_l ))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      pl_i = sum(log(S_t_i1 - S_t_i2 ))
      
      #least squares
      z_t_ia = NULL
      for(i in 1:n){
        alpha_re = alpha
        alpha_re[re_ind] = alpha_re[re_ind] + a_re_cols[i,re_ind]
        ind.long = which(id_long == i)
        z_t_ia = c(z_t_ia, phi.mu[ind.long,] %*% (alpha_re) ) 
      }
      
      pl_ls = (1/(2*sigma2_Et))*sum((cont - z_t_ia)^2)
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #print(c("re", log_lik_old, log_lik))
      step_size_a1 = 1
      if(((log_lik < log_lik_old) & step_size_a1 == 1)){
        ii = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          
          a_re_cols = a_re_cols_old + omega1 * a_cols_update
          
          a_re_long = unlist(c(a_re_cols))
          a_re_cols_pad[,re_ind] = as.matrix(a_re_cols)
          a_re_cols_pad_quad[,re_ind] = matrix(sapply(a_re_cols, rep, 15), ncol = 2)
          
          #step size for random effects
          #likelihood
          eXtB = exp(fixed %*% beta)
          
          h0_t_quad = quad.psi.event %*% theta
          exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          H0_t_l = H0_t[dat.baseline$left == 1]
          #H0_t_i2 = H0_t[dat.baseline$interval == 1]
          
          h0_t_quad_i1 = quad.psi.y_i1 %*% theta
          exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
          h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
          h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i1 = quad.lambda_i1 * apply(quad.w *  h0_t_star_quad_i1, 2, sum)
          
          h0_t_quad_i2 = quad.psi.y_i2 %*% theta
          exp_zTg_quad_i2 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i2, 1, sum))
          h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
          h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
          H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #left
          eXtB_l = exp(fixed_l %*% beta)
          S_t_l = exp(-eXtB_l * H0_t_l)
          pl_l = sum(log(1 - S_t_l))
          
          #interval
          eXtB_i = exp(fixed_i %*% beta)
          S_t_i1 = exp(-eXtB_i * H0_t_i1)
          S_t_i2 = exp(-eXtB_i * H0_t_i2)
          pl_i = sum(log(S_t_i1 - S_t_i2 ))
          
          #least squares
          z_t_ia = NULL
          for(i in 1:n){
            alpha_re = alpha
            alpha_re[re_ind] = alpha_re[re_ind] + a_re_cols[i,re_ind]
            ind.long = which(id_long == i)
            z_t_ia = c(z_t_ia, phi.mu[ind.long,] %*% (alpha_re)) 
          }
          
          pl_ls = (1/(2*sigma2_Et))*sum((cont - z_t_ia)^2)
          
          #log likelihood
          log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
            sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
            theta.lambda* t(theta)%*%theta.G%*%theta - 
            alpha.lambda* t(alpha)%*%alpha.G%*%alpha
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          ii = ii+1
          if(ii>500){break}
          
          
        }
        #print(c(omega1))
        
      }
      
      #hist(a_re_cols[,1])
      #hist(a_re_cols[,2])
      
      #pos_temp = rep(1, m)
      for(u in 1:m){
        if((theta[u]< (1e-3) & theta_score[u] > 0)){
          #pos_temp[u] = 0
          theta[u] = 0.5
        }
      }
      
      print(c(iter, beta, gamma, alpha, log_lik))
      
      #print(c(iter, theta_score_pos, theta_score_neg, theta_score_pos/theta_score_neg, log_lik))
      
      
      if(all(abs(c(beta - beta_old, gamma -gamma_old, alpha - alpha_old, theta - theta_old )) < 1e-5) & 
         all(abs(c(a_re_cols - a_re_cols_old)) < 1e-5)){
        break
      }
      
      #if(it == 1 & iter == 100){
      #  break
      #}
      
      
      
    } #end inner loop
    
    
    print(c(iter, beta, gamma, log_lik))
    
    
    
    Hessian = matrix(0, nrow = (p + q + m + r + w*n), ncol = (p + q + m + r + w*n))
    #Hessian = matrix(0, nrow = (p + q + m + r + n), ncol = (p + q + m + r + n))
    
    phi.mu.1 = phi.mu[,re_ind]
    phi.mu.summed  = data.frame(id_long, phi.mu.1)
    phi.mu.summed = as.matrix(aggregate( phi.mu.summed[,2:(w+1)], list(phi.mu.summed[,1]), FUN = sum )[,-1])
    phi.mu.summed_alpha  = data.frame(id_long, phi.mu)
    phi.mu.summed_alpha = as.matrix(aggregate( phi.mu.summed_alpha[,2:(r+1)], list(phi.mu.summed_alpha[,1]), FUN = sum )[,-1])
    
    
    #least squares
    zt_i_samp_est = NULL
    for(i in 1:n){
      alpha_re = alpha
      alpha_re[re_ind] = alpha_re[re_ind] + unlist(c(a_re_cols[i,re_ind]))
      ind.long = which(id_long == i)
      
      zt_i_samp_est = c(zt_i_samp_est, phi.mu[ind.long,] %*% (alpha_re)) 
    }
    
    eXtB_e = exp(fixed_e %*% beta)
    eXtB_r = exp(fixed_r %*% beta)
    eXtB_l = exp(fixed_l %*% beta)
    eXtB_i = exp(fixed_i %*% beta)
    
    h0_t_e = psi_e %*% theta
    
    h0_t_quad = quad.psi.event %*% theta
    exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
    h0_t_star_quad = h0_t_quad * exp_zTg_quad
    h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
    H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
    
    H0_t_e = H0_t[dat.baseline$event == 1]
    H0_t_r = H0_t[dat.baseline$right == 1]
    H0_t_l = H0_t[dat.baseline$left == 1]
    #H0_t_i2 = H0_t[dat.baseline$interval == 1]
    
    h0_t_quad_i1 = quad.psi.y_i1 %*% theta
    exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
    h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
    h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
    H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
    
    H0_t_i1_long = rep(0, n)
    H0_t_i1_long[which(dat.baseline$interval == 1)] = H0_t_i1
    
    S_t = exp(-eXtB * H0_t)
    S_t_i1_long = rep(0, n)
    
    S_t_l = exp(-eXtB_l * H0_t_l)
    S_t_i1 = exp(-eXtB_i * H0_t_i1)
    S_t_i2 = exp(-eXtB_i * H0_t_i2)
    
    S_t = exp(-eXtB * H0_t)
    S_t_i1_long = rep(0, n)
    S_t_i1_long[which(dat.baseline$interval == 1)] = S_t_i1
    
    #first derivative of gamma
    exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
    A_t_quad = h0_t_quad * exp_zTg_quad * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
    A_t_quad = matrix(A_t_quad, ncol = n, nrow = 15, byrow = FALSE)
    A_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  A_t_quad, 2, sum)
    
    A_t_e = A_t[dat.baseline$event == 1]
    A_t_r = A_t[dat.baseline$right == 1]
    A_t_l = A_t[dat.baseline$left == 1]
    A_t_i2 = A_t[dat.baseline$interval == 1]
    
    exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
    A_t_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1 * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)
    A_t_quad_i1 = matrix(A_t_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
    A_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  A_t_quad_i1, 2, sum)
    
    A_t_il_long = rep(0, n)
    A_t_il_long[which(dat.baseline$interval == 1)] = A_t_i1
    
    #second derivative of gamma
    z_t_quad = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
    A2_t_quad_long = h0_t_quad * exp_zTg_quad * (z_t_quad)^2
    A_t_quad2 = matrix(A2_t_quad_long, ncol = n, nrow = 15, byrow = FALSE)
    A2_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  A_t_quad2, 2, sum)
    
    A_t_quad2_e = A2_t[dat.baseline$event == 1]
    A_t_quad2_r = A2_t[dat.baseline$right == 1]
    A_t_quad2_l = A2_t[dat.baseline$left == 1]
    A_t_quad2_i2 = A2_t[dat.baseline$interval == 1]
    
    z_t_quad_i1 = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)
    A_t_quad_i1_long = h0_t_quad_i1 * exp_zTg_quad_i1 * (z_t_quad_i1)^2
    A_t_quad2_i1 = matrix(A_t_quad_i1_long, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
    A_t_quad2_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  A_t_quad2_i1, 2, sum)
    
    A2_t_i1_long = rep(0, length(A2_t))
    A2_t_i1_long[which(dat.baseline$interval == 1)] = A_t_quad2_i1
    
    #d RE d2 gamma
    A_dRE_t_quad2_long = c(c(gamma) * h0_t_quad * exp_zTg_quad * (z_t_quad)^2 + 2 * h0_t_quad * exp_zTg_quad * z_t_quad) * quad.phi.event[,re_ind]
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* A_dRE_t_quad2_long)
    A2_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
    
    A2_dRE_t_e = A2_dRE[dat.baseline$event == 1,]
    A2_dRE_t_r = A2_dRE[dat.baseline$right == 1,]
    A2_dRE_t_l = A2_dRE[dat.baseline$left == 1,]
    A2_dRE_t_i2 = A2_dRE[dat.baseline$interval == 1,]
    
    A_dRE_t_quad_long = c(h0_t_quad * exp_zTg_quad + h0_t_quad * exp_zTg_quad * z_t_quad * c(gamma)) * quad.phi.event[,re_ind]
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* A_dRE_t_quad_long)
    A_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
    
    #first derivative of theta
    psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* psi_t_star_quad)
    Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
    
    Psi_t_e = Psi_t[dat.baseline$event == 1,]
    Psi_t_r = Psi_t[dat.baseline$right == 1,]
    Psi_t_l = Psi_t[dat.baseline$left == 1,]
    Psi_t_i2 = Psi_t[dat.baseline$interval == 1,]
    
    psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
    df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* psi_t_star_quad_i1)
    Psi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
    
    Psi_t_i1_long = matrix(0, nrow = nrow(Psi_t), ncol = ncol(Psi_t))
    Psi_t_i1_long[which(dat.baseline$interval==1),] = Psi_t_i1
    
    #first derivative of alpha
    B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* B_t_quad)
    B_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    Bphi_t_e = B_t[dat.baseline$event == 1,]
    Bphi_t_r = B_t[dat.baseline$right == 1,]
    Bphi_t_l = B_t[dat.baseline$left == 1,]
    Bphi_t_i2 = B_t[dat.baseline$interval == 1,]
    
    B_t_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * quad.phi.y_i1
    df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* B_t_quad_i1)
    Bphi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    B_t_i1_long = matrix(0, nrow = nrow(B_t), ncol = ncol(B_t))
    B_t_i1_long[which(dat.baseline$interval == 1),] = Bphi_t_i1
    
    #second derivative of alpha
    B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
    B_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, 15)) * B_t_quad 
    B_t_quad2_e = B_t_quad2[c(sapply(dat.baseline$event == 1, rep, 15)),]
    B_t_quad2_r = B_t_quad2[c(sapply(dat.baseline$right == 1, rep, 15)),]
    B_t_quad2_l = B_t_quad2[c(sapply(dat.baseline$left == 1, rep, 15)),]
    B_t_quad2_i2 = B_t_quad2[c(sapply(dat.baseline$interval == 1, rep, 15)),]
    B_t_quad_i1_2 = c(rep(quad.w, sum(dat.baseline$interval == 1))) * 
      c(sapply(quad.lambda_i1, rep, 15)) * B_t_quad_i1 
    
    # gamma alpha hessian
    z_est_quad = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
    E_t_quad = c(gamma) *  c(h0_t_quad * exp_zTg_quad * z_est_quad ) * quad.phi.event + 
      c(h0_t_quad * exp_zTg_quad ) * quad.phi.event
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* E_t_quad)
    E_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    Ephi_t_e = E_t[dat.baseline$event == 1,]
    Ephi_t_r = E_t[dat.baseline$right == 1,]
    Ephi_t_l = E_t[dat.baseline$left == 1,]
    Ephi_t_i2 = E_t[dat.baseline$interval == 1,]
    
    E_t_quad_i1 = c(gamma) *  c(h0_t_quad_i1 * exp_zTg_quad_i1 *apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)) * quad.phi.y_i1 + 
      c(h0_t_quad_i1 * exp_zTg_quad_i1 ) * quad.phi.y_i1
    df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* E_t_quad_i1)
    Ephi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    Ephi_t_i1_long = matrix(0, nrow = nrow(E_t), ncol = ncol(E_t))
    Ephi_t_i1_long[which(dat.baseline$interval == 1),] = Ephi_t_i1
    
    # theta gamma hessian (need to multiply this with z_quad)
    psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
    psi_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, 15)) * psi_t_star_quad 
    psi_t_quad2_e = psi_t_quad2[c(sapply(dat.baseline$event == 1, rep, 15)),]
    psi_t_quad2_r = psi_t_quad2[c(sapply(dat.baseline$right == 1, rep, 15)),]
    psi_t_quad2_l = psi_t_quad2[c(sapply(dat.baseline$left == 1, rep, 15)),]
    psi_t_quad2_i2 = psi_t_quad2[c(sapply(dat.baseline$interval == 1, rep, 15)),]
    
    psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
    psi_t_quad2_i1 = c(rep(quad.w, sum(dat.baseline$interval == 1))) * 
      c(sapply(quad.lambda_i1, rep, 15)) * psi_t_star_quad_i1 
    
    psi_t_i1_quad2 = matrix(0, nrow = nrow(psi_t_quad2), ncol = ncol(psi_t_quad2))
    psi_t_i1_quad2[i1_ind_quad,] = psi_t_quad2_i1
    
    #theta alpha hessian
    ##HERE NEED TO MULTIPLY TOGETHER FIRST DERIVATIVE OF THETA WITH PHI MATRIX IN MULTIPLICATION PART
    
    
    
    #beta
    #beta beta
    beta_hess_e = t(fixed_e) %*% diag(c(-eXtB_e * H0_t_e)) %*% fixed_e
    beta_hess_r =  t(fixed_r) %*% diag(c(-eXtB_r * H0_t_r)) %*% fixed_r
    beta_hess_l = t(fixed_l) %*% diag(c(((eXtB_l*H0_t_l)*S_t_l - (eXtB_l*H0_t_l)^2*S_t_l - eXtB_l*H0_t_l*S_t_l^2)  / (1 - S_t_l)^2)) %*% fixed_l
    beta_hess_i = t(fixed_i) %*% diag(c(-eXtB_i * (S_t_i1 * H0_t_i1 ) / (S_t_i1 - S_t_i2))) %*% fixed_i + 
      t(fixed_i) %*% diag(c(eXtB_i * (S_t_i2 * H0_t_i2 ) / (S_t_i1 - S_t_i2))) %*% fixed_i -
      t(fixed_i) %*% diag(c(S_t_i1 * S_t_i2 * (eXtB_i * H0_t_i1 - eXtB_i *  H0_t_i2)^2/ (S_t_i1 - S_t_i2)^2)) %*% fixed_i
    beta_hess = beta_hess_e + beta_hess_r + beta_hess_l + beta_hess_i
    Hessian[1:p, 1:p] = beta_hess
    
    #beta gamma
    betagamma_hess_e = t(fixed_e) %*% (c(-eXtB_e * A_t_e))
    betagamma_hess_r = t(fixed_r) %*% (c(-eXtB_r * A_t_r))
    betagamma_hess_l = t(fixed_l) %*% c(eXtB_l * (S_t_l * A_t_l)/ (1 - S_t_l)) -
      t(fixed_l) %*% c(eXtB_l^2 * S_t_l * A_t_l * H0_t_l/ (1 - S_t_l)^2)
    betagamma_hess_i = t(fixed_i) %*% c(-eXtB_i * S_t_i1 * A_t_i1 /(S_t_i1 - S_t_i2)) +
      t(fixed_i) %*% c(eXtB_i * S_t_i2 * A_t_i2 /(S_t_i1 - S_t_i2)) -
      t(fixed_i) %*% c(eXtB_i * S_t_i1 * S_t_i2 * (A_t_i1 - A_t_i2) * (H0_t_i1 * eXtB_i - H0_t_i2 * eXtB_i)/(S_t_i1 - S_t_i2)^2)
    betagamma_hess = betagamma_hess_e + betagamma_hess_r + betagamma_hess_l + betagamma_hess_i
    Hessian[1:p, (p+1):(p+q)] = betagamma_hess
    Hessian[(p+1):(p+q), 1:p] = t(Hessian[1:p, (p+1):(p+q)])
    
    #beta theta
    betatheta_hess_e = t(c(-eXtB_e) * fixed_e) %*% Psi_t_e
    betatheta_hess_r = t(c(-eXtB_r) * fixed_r) %*% Psi_t_r
    betatheta_hess_l = t(fixed_l * c(eXtB_l*(S_t_l - eXtB_l*H0_t_l*S_t_l - S_t_l^2)  / (1 - S_t_l)^2)) %*% Psi_t_l
    betatheta_hess_i = t(-fixed_i * c(eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2))) %*% Psi_t_i1 + 
      t(fixed_i * c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2))) %*% Psi_t_i2 -
      t(fixed_i * c(eXtB_i^2 * S_t_i1 * S_t_i2 * (H0_t_i1 - H0_t_i2)/(S_t_i1 - S_t_i2)^2)) %*% (Psi_t_i1 - Psi_t_i2) 
    betatheta_hess = betatheta_hess_e + betatheta_hess_r + betatheta_hess_l + betatheta_hess_i
    Hessian[1:p, (p+q+1):(p+q+m)] = betatheta_hess
    Hessian[(p+q+1):(p+q+m), 1:p] = t(Hessian[1:p, (p+q+1):(p+q+m)])
    
    #beta alpha
    betaalpha_hess_e = t(c(gamma) * c(-eXtB_e) * fixed_e) %*% Bphi_t_e
    betaalpha_hess_r = t(c(gamma) * c(-eXtB_r) * fixed_r) %*% Bphi_t_r
    betaalpha_hess_l = t(c(gamma) * fixed_l * c(eXtB_l * (S_t_l - eXtB_l*H0_t_l*S_t_l - S_t_l^2)  / (1 - S_t_l)^2) ) %*% Bphi_t_l
    betaalpha_hess_i = t(fixed_i) %*% (c(gamma) * c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)) * Bphi_t_i1) + 
      t(fixed_i) %*% (c(gamma) * c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)) * Bphi_t_i2) -
      t(c(gamma) * fixed_i * c(eXtB_i * S_t_i1 * S_t_i2 * (eXtB_i*H0_t_i1 - eXtB_i*H0_t_i2)/(S_t_i1 - S_t_i2)^2) ) %*% (Bphi_t_i1 - Bphi_t_i2)
    betaalpha_hess = betaalpha_hess_e + betaalpha_hess_r + betaalpha_hess_l + betaalpha_hess_i
    Hessian[1:p, (p+q+m+1):(p+q+m+r)] = betaalpha_hess
    Hessian[(p+q+m+1):(p+q+m+r), 1:p] = t(Hessian[1:p, (p+q+m+1):(p+q+m+r)])
    
    #beta kappa
    betaa1_hess = do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
      c(c(dat.baseline$event + dat.baseline$right) * c(gamma) * c(eXtB) * B_t[,re_ind]) +
      do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
      c(c(dat.baseline$left) *  c(gamma) * c(eXtB * S_t/(1-S_t)) * B_t[,re_ind]) -
      do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
      c(c(dat.baseline$left) *  c(gamma) * c(eXtB * S_t * (H0_t*eXtB)/(1-S_t)^2) * B_t[,re_ind]) -
      do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
      c(c(dat.baseline$interval) *  c(gamma) * c(eXtB * S_t_i1_long/(S_t_i1_long-S_t)) * B_t_i1_long[,re_ind]) + 
      do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
      c(c(dat.baseline$interval) *  c(gamma) * c(eXtB * S_t/(S_t_i1_long-S_t)) * B_t[,re_ind]) - 
      do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
      c(c(dat.baseline$interval) *  c(gamma) * c(eXtB * S_t_i1_long * S_t * (H0_t_i1_long*eXtB - H0_t*eXtB)/(S_t_i1_long-S_t)^2) * (B_t_i1_long[,re_ind] - B_t[,re_ind]))
    Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)] = betaa1_hess
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), 1:p] = t(Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)])
    
    #gamma
    #gamma gamma
    gamma_hess_e = t(-eXtB_e) %*% A_t_quad2_e
    gamma_hess_r = t(-eXtB_r) %*% A_t_quad2_r
    gamma_hess_l = t(S_t_l * eXtB_l/(1-S_t_l)) %*% A_t_quad2_l - 
      t(c( S_t_l* A_t_l * eXtB_l^2/(1-S_t_l)^2)) %*% A_t_l
    gamma_hess_i = t(-eXtB_i * S_t_i1 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i1 +
      t(eXtB_i* S_t_i2 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i2 -
      t(c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2 ) * (A_t_i1 - A_t_i2)) %*% (A_t_i1 - A_t_i2)
    gamma_hess = gamma_hess_e + gamma_hess_r + gamma_hess_l + gamma_hess_i
    Hessian[(p+1):(p+q), (p+1):(p+q)] = gamma_hess
    
    #gamma theta
    if(sum(dat.baseline$event) > 0){
      gammatheta_hess_e = t(c(sapply(c(-eXtB_e), rep, 15)) * psi_t_quad2_e)%*% (c(sapply(quad.lambda, rep, 15)) * rep(quad.w, n) * z_t_quad)[c(sapply(dat.baseline$event == 1, rep, 15))]
    }else{
      gammatheta_hess_e = 0
    }
    gammatheta_hess_r = t(c(sapply(c(-eXtB_r), rep, 15)) * psi_t_quad2_r) %*% (c(sapply(quad.lambda, rep, 15)) * rep(quad.w, n) * z_t_quad)[c(sapply(dat.baseline$right == 1, rep, 15))]
    gammatheta_hess_l = t(c(sapply(c(eXtB_l * (S_t_l)/(1-S_t_l)), rep, 15)) *psi_t_quad2_l) %*% (c(sapply(quad.lambda, rep, 15)) * rep(quad.w, n) *  z_t_quad)[dat.baseline$left == 1] -
      t(Psi_t_l) %*% (c(eXtB_l^2 *(S_t_l)/(1-S_t_l)^2)* A_t_l) 
    gammatheta_hess_i = t(c(sapply(c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)), rep, 15)) *psi_t_quad2_i1) %*% (c(sapply(quad.lambda_i1, rep, 15)) * rep(quad.w, sum(dat.baseline$interval)) * z_t_quad_i1) + 
      t(c(sapply(c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)), rep, 15)) *psi_t_quad2_i2) %*% (c(sapply(quad.lambda, rep, 15)) * rep(quad.w, n) * z_t_quad)[dat.baseline$interval == 1] -
      t(Psi_t_i1 - Psi_t_i2) %*% (c(eXtB_i^2 * (A_t_i1 - A_t_i2) * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2) ) 
    gammatheta_hess = gammatheta_hess_e + gammatheta_hess_r + gammatheta_hess_l + gammatheta_hess_i
    Hessian[(p+1):(p+q), (p+q+1):(p+q+m)] = gammatheta_hess
    Hessian[(p+q+1):(p+q+m), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+1):(p+q+m)])
    
    #gamma alpha
    if(sum(dat.baseline$event)>0){
      gammaalpha_hess_e = apply(phi_e, 2, sum) - t(eXtB_e) %*% (Ephi_t_e)
    }else{
      gammaalpha_hess_e = 0
    }
    gammaalpha_hess_r = t(-eXtB_r) %*% (Ephi_t_r)
    gammaalpha_hess_l = t(eXtB_l * (S_t_l )/((1-S_t_l))) %*% (Ephi_t_l) -
      t(as.numeric(gamma) * as.numeric(eXtB_l^2 * A_t_l * S_t_l/((1-S_t_l)^2))) %*% Bphi_t_l
    gammaalpha_hess_i = t(-eXtB_i * S_t_i1/(S_t_i1 - S_t_i2)) %*% (Ephi_t_i1) +
      t(eXtB_i * S_t_i2/(S_t_i1 - S_t_i2)) %*% (Ephi_t_i2) -
      t(as.numeric(gamma) * as.numeric(eXtB_i^2 * S_t_i1 * S_t_i2 * (A_t_i1 - A_t_i2)/(S_t_i1 - S_t_i2)^2)) %*% (Bphi_t_i1-Bphi_t_i2)
    gammaalpha_hess = gammaalpha_hess_e + gammaalpha_hess_r + gammaalpha_hess_l + gammaalpha_hess_i
    Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)] = gammaalpha_hess
    Hessian[(p+q+m+1):(p+q+m+r), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)])
    
    #gamma kappa
    gammaa1_hess = c(dat.baseline$event * phi[,re_ind] - c(dat.baseline$event + dat.baseline$right) * (c(eXtB) * E_t[,re_ind])) +
      c(c(dat.baseline$left) * (c(eXtB * S_t/(1-S_t)) * E_t[,re_ind])) -
      c(c(dat.baseline$left) * c(gamma) * (c((eXtB)^2 * A_t * S_t/(1-S_t)^2) * B_t_re[,re_ind])) - 
      c(c(dat.baseline$interval) * (c(eXtB * S_t_i1_long/(S_t_i1_long-S_t)) * Ephi_t_i1_long[,re_ind])) + 
      c(c(dat.baseline$interval) * (c(eXtB * S_t/(S_t_i1_long-S_t)) * E_t[,re_ind])) - 
      c(c(dat.baseline$interval) * c(gamma) * (c((eXtB)^2 * (A_t_il_long - A_t) * S_t_i1_long * S_t/(S_t_i1_long-S_t)^2) * (B_t_i1_long[,re_ind] - B_t_re[,re_ind])))
    Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)] = gammaa1_hess
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)])
    
    #theta
    theta_hess_e = - t(as.numeric(1/(h0_t_e^2)) * psi_e) %*% psi_e
    theta_hess_l = - t(as.numeric(eXtB_l^2 * S_t_l/(1 - S_t_l)^2) * Psi_t_l) %*% Psi_t_l
    theta_hess_i = - t(as.numeric(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2) * (Psi_t_i1 - Psi_t_i2)) %*% (Psi_t_i1 - Psi_t_i2)
    theta_hess = theta_hess_e + theta_hess_l + theta_hess_i
    Hessian[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = theta_hess
    
    #theta alpha
    if(sum(dat.baseline$event)>0){
      thetaalpha_hess_e = t(c(sapply(c(gamma) * c(-eXtB_e), rep, 15)) *psi_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, 15)),]
    }else{
      thetaalpha_hess_e = 0
    }
    
    thetaalpha_hess_r = t(c(sapply(c(gamma) * c(-eXtB_r), rep, 15)) * psi_t_quad2_r) %*% (c(rep(quad.w, sum(dat.baseline$right))) * c(sapply(quad.lambda[dat.baseline$right == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$right == 1, rep, 15)),])
    thetaalpha_hess_l = t(c(sapply(c(gamma) * c(eXtB_l * (S_t_l)/((1-S_t_l))), rep, 15)) * psi_t_quad2_l) %*% (c(rep(quad.w, sum(dat.baseline$left))) * c(sapply(quad.lambda[dat.baseline$left == 1], rep, 15)) * quad.phi.event[dat.baseline$left == 1,]) -
      t(c(gamma) * c(eXtB_l^2 * S_t_l/((1-S_t_l)^2)) * Psi_t_l) %*% Bphi_t_l
    thetaalpha_hess_i = t(c(sapply(c(gamma) * c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)), rep, 15)) * psi_t_quad2_i1) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda_i1, rep, 15)) * quad.phi.y_i1) +
      t(c(sapply(c(gamma) * c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)), rep, 15)) * psi_t_quad2_i2) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda[dat.baseline$interval == 1], rep, 15)) * quad.phi.event[dat.baseline$interval == 1,]) -
      t(c(c(gamma) * eXtB_i^2 * S_t_i1 * S_t_i2/((S_t_i1-S_t_i2)^2)) * (Psi_t_i1 - Psi_t_i2)) %*% (Bphi_t_i1 - Bphi_t_i2)
    Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)] = thetaalpha_hess_e + thetaalpha_hess_r + thetaalpha_hess_l + thetaalpha_hess_i
    Hessian[(p+q+m+1):(p+q+m+r), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)])
    
    #theta a1
    temp = c(quad.phi.event[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) psi_t_quad2))
    df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, 15)), w),
                    re_no = c(sapply(re_ind, rep, 15*n)), c(rep(quad.w, n)) * temp)
    F_t = c(quad.lambda) * as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    temp = c(quad.phi.y_i1[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) psi_t_quad2_i1))
    df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)),
                    re_no = c(sapply(re_ind, rep, 15*sum(dat.baseline$interval))), c(rep(quad.w, sum(dat.baseline$interval == 1))) * temp)
    F_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    F_t_i1_long = matrix(0, nrow = nrow(F_t), ncol = ncol(F_t))
    F_t_i1_long[rep(which(dat.baseline$interval == 1), w),] = F_t_i1
    
    thetaa1_hess = rep((- c(dat.baseline$event + dat.baseline$right) * eXtB * c(gamma)), w) * F_t + 
      rep((c(dat.baseline$left) * eXtB * S_t/(1-S_t) * c(gamma)), w) * F_t -
      rep((c(dat.baseline$left) * (eXtB^2) * S_t/((1-S_t)^2) * c(gamma)), w) * c(B_t_re) * do.call(rbind, lapply(seq_len(w), function(X) Psi_t)) -
      rep((c(dat.baseline$interval) * eXtB * S_t_i1_long/(S_t_i1_long-S_t) * c(gamma)), w) * F_t_i1_long + 
      rep((c(dat.baseline$interval) * eXtB * S_t/(S_t_i1_long-S_t) * c(gamma)), w) * F_t -
      rep((c(dat.baseline$interval) * (eXtB^2) * S_t * S_t_i1_long/((S_t_i1_long-S_t)^2) * c(gamma)), w) * c(B_t_re_i1_long - B_t_re) * do.call(rbind, lapply(seq_len(w), function(X) (Psi_t_i1_long - Psi_t)))
      
    Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)] = t(thetaa1_hess)
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)])
    
    
    #alpha
    if(sum(dat.baseline$event)>0){
      alpha_hess_e = t(-c(sapply(c(gamma^2) * c(eXtB_e), rep, 15)) * B_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, 15)),]
    }else{
      alpha_hess_e = 0
    }
    alpha_hess_r = t(-c(sapply(c(gamma^2) * c(eXtB_r), rep, 15)) * B_t_quad2_r) %*% (c(rep(quad.w, sum(dat.baseline$right))) * c(sapply(quad.lambda[dat.baseline$right == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$right == 1, rep, 15)),])
    alpha_hess_l = t(c(sapply(c(gamma^2) * c(eXtB_l*(S_t_l)/(1-S_t_l)), rep, 15)) * B_t_quad2_l) %*% (c(rep(quad.w, sum(dat.baseline$left))) * c(sapply(quad.lambda[dat.baseline$left == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$left == 1, rep, 15)),]) -
      t(c(c(gamma^2) * (eXtB_l^2 * S_t_l/(1-S_t_l)^2)) * Bphi_t_l) %*% Bphi_t_l
    alpha_hess_i = t(-c(sapply(c(gamma^2) * c(eXtB_i*S_t_i1 /(S_t_i1-S_t_i2)), rep, 15)) * B_t_quad_i1_2) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda_i1, rep, 15)) * quad.phi.y_i1)  + 
      t(c(sapply(c(gamma^2) * c(eXtB_i*S_t_i2 /(S_t_i1-S_t_i2)), rep, 15)) * B_t_quad2_i2) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda[dat.baseline$interval == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$interval == 1, rep, 15)),]) -
      t(c(c(gamma^2) * (eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i1 + 
      t(c(c(gamma^2) * (eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i2
    Hessian[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = alpha_hess_e + alpha_hess_r + alpha_hess_l + alpha_hess_i
    
    #alpha a1
    temp = c(quad.phi.event[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) B_t_quad2))
    df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, 15)), w),
                    re_no = c(sapply(re_ind, rep, 15*n)), c(rep(quad.w, n)) * temp)
    D_t_re = c(quad.lambda) * as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    temp = c(quad.phi.y_i1[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) B_t_quad2_i1))
    df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)),
                    re_no = c(sapply(re_ind, rep, 15*sum(dat.baseline$interval))), c(rep(quad.w, sum(dat.baseline$interval == 1))) * temp)
    D_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    D_t_re_i1_long = matrix(0, nrow = nrow(D_t_re), ncol = ncol(D_t_re))
    D_t_re_i1_long[rep(which(dat.baseline$interval == 1), w),] = D_t_re_i1
    
    alphaa1_hess = - rep((c(dat.baseline$event + dat.baseline$right) * eXtB * c(gamma^2)), w) * D_t_re + 
      rep((c(dat.baseline$left) * eXtB * (S_t/(1-S_t)) * c(gamma^2)), w) * D_t_re - 
      c(rep((c(dat.baseline$left) * (eXtB^2) * (S_t/(1-S_t)^2) * c(gamma^2)), w) * B_t_re) * do.call(rbind, lapply(seq_len(w), function(X) B_t)) -
      rep((c(dat.baseline$interval) * eXtB * (S_t_i1_long/(S_t_i1_long-S_t)) * c(gamma^2)), w) * D_t_re_i1_long + 
      rep((c(dat.baseline$interval) * eXtB * (S_t/(S_t_i1_long-S_t)) * c(gamma^2)), w) * D_t_re - 
      c(rep((c(dat.baseline$interval) * (eXtB^2) * (S_t_i1_long*S_t/(S_t_i1_long-S_t)^2) * c(gamma^2)), w) * (B_t_re_i1_long - B_t_re)) * do.call(rbind, lapply(seq_len(w), function(X) (B_t_i1_long - B_t)))
      
    Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)] = t(alphaa1_hess)
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+1):(p+q+m+r)] = t(Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)])
    
    #a1 a1
    #temp = (c(quad.phi.event[,re_ind]) * 
    #          do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma^2), rep, 15)) * B_t_quad2[,re_ind])))
    #df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, 15)), w),
    #                re_no = c(sapply(re_ind, rep, 15*n)), temp)
    a1a1_hess = alphaa1_hess[,re_ind]
    
    for(re in 1:w){
      diag(Hessian[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = as.matrix(a1a1_hess)[((re-1)*n + 1):((re-1)*n + n),re]
      
      if(re > 1){
        diag(Hessian[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess[((re-1)*n + 1):((re-1)*n + n), (re-1)]
        diag(Hessian[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n)]) = diag(Hessian[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)])
      }
      
    }
    
    #add penalty part
    Q_theta = Q_alpha = Q_a1 = H_epsilon = matrix(0, nrow = (p+q+m+r+w*n), ncol = (p+q+m+r+w*n))
    
    #least squares part
    H_epsilon[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = - c(1/sigma2_Et) * t(phi.mu) %*% phi.mu #alpha alpha
    
    temp = c(phi.mu[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) phi.mu))
    df = data.frame(id_long,
                    re_no = c(sapply(re_ind, rep, nrow(dat.long))), temp)
    alphaa1_hess_ls = - c(1/sigma2_Et) *  as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    H_epsilon[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)] = t(alphaa1_hess_ls) #alpha re
    H_epsilon[ (p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+1):(p+q+m+r)] = t(H_epsilon[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)])
    
    temp = c(phi.mu[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) phi.mu[,re_ind]))
    df = data.frame(id_long,
                    re_no = c(sapply(re_ind, rep, nrow(dat.long))), temp)
    a1a1_hess_ls = - c(1/sigma2_Et) *  as.matrix(aggregate( df[,3:(w+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    for(re in 1:w){
      diag(H_epsilon[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n),re]
      
      if(re > 1){
        diag(H_epsilon[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n), (re-1)]
        diag(H_epsilon[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n)]) = diag(H_epsilon[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)])
      }
      
    }
    
    #theta part
    Q_theta[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = - theta.lambda*theta.G
    
    #alpha part
    Q_alpha[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = - alpha.lambda*alpha.G
    
    #re part
    diag(Q_a1[(p+q+m+r+1):(p+q+m+r+w*n),(p+q+m+r+1):(p+q+m+r+w*n)]) = -c(1/c(sapply(sigma2_re, rep, n)))
    
    H_full = Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1
    
    #Hessian_trunc = H_full[1:(p+q+m+r), 1:(p+q+m+r)]
    #Hessian_RE = H_full[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+r+1):(p+q+m+r+w*n)]
    
    #approximate Hessian for fixed parameters
    
    
    
    #new iteration - a^(k-1) using eta^(k-1), a^(k-1)
    
    #d_Ai_matrix = matrix(c(a_re_cols) - c(a_re_cols_old))
    #d_Omega_matrix = matrix(c(beta, gamma, theta, alpha) - c(beta_old, gamma_old, theta_old, alpha_old))
    #d_Omega_matrix = d_Omega_matrix + 1e-5
    
    #d_Ai_d_Omega = d_Ai_matrix %*% t(1/d_Omega_matrix) 
    
    #Hessian_est = Hessian_trunc + t(d_Ai_d_Omega) %*% Hessian_RE %*% d_Ai_d_Omega
    
    #update smoothing parameters
    H_full_inv = solve(-(Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1))
    
    df.theta.old = df.theta
    df.alpha.old = df.alpha
    df.epsilon.old = df.epsilon
    df.a_re.old = df.a_re
    
    df.theta = sum(diag( -H_full_inv %*% (-Q_theta)))
    df.theta = m - df.theta
    theta.sigma2 = t(theta) %*% theta.G %*% theta / df.theta
    theta.lambda = c(1/(2*theta.sigma2))
    #theta.lambda = 0
    
    df.alpha = sum(diag(H_full_inv %*% (-Q_alpha)))
    df.alpha = r - df.alpha
    alpha.sigma2 = t(alpha) %*% alpha.G %*% alpha / df.alpha
    alpha.lambda = c(1/(2*alpha.sigma2))
    alpha.lambda = 0
    
    df.epsilon = sum(diag(H_full_inv %*% (-H_epsilon)))
    df.epsilon = n_long - df.epsilon
    sigma2_Et = as.numeric(sum((cont - zt_i_samp_est)^2)/ df.epsilon)
    #sigma2_Et = 0.0025
    
    sigma2_re_old = sigma2_re
    sigma2_re = NULL
    for(re in 1:w){
      sigma2_temp = sum((a_re_cols[,re])^2)/(n)
      sigma2_re = c(sigma2_re, sigma2_temp)
    }
    
    #sigma2_re = NULL
    #for(re in 1:w){
    #  Q_a_temp = Q_a1
    #  pos_a1 = rep(0, p+q+m+r+w*n)
    #  pos_a1[(p+q+m+r+(n*(re-1))+1):(p+q+m+r+(n*(re-1))+n)] = 1
    #  diag(Q_a_temp) = diag(Q_a_temp) * pos_a1
    
    #  df.temp = sum(diag(H_full_inv %*% (-Q_a_temp)))
    #  df.temp = n - df.temp
    #  sigma2_temp = as.numeric(t(a_re_cols[,re]) %*% a_re_cols[,re] / df.temp)
    #  df.a_re[re] = c(1/(2*sigma2_temp))
    #  sigma2_re = c(sigma2_re, sigma2_temp)
    #}
    
    #if(any(c(alpha.sigma2, theta.sigma2, sigma2_a0, sigma2_a1, sigma2_Et) < 0)){
    #  print(c(alpha.sigma2, theta.sigma2, sigma2_a0, sigma2_a1, sigma2_Et))
    #  print(c("Variance estimate is negative"))
    #  break
    #}
    
    print(c(it, theta.lambda, sigma2_re, sigma2_Et))
    
    if(all(abs(c(df.theta.old - df.theta, df.alpha.old - df.alpha )) < 1) & 
       abs(df.epsilon.old - df.epsilon) < 1 & all(abs(sigma2_re - sigma2_re_old) < 0.1) ){
      break
    }
    
    
    
  } #end outer loop
  
  #full Hessian matrix for beta, gamma, theta, alpha, kappa
  Hessian = matrix(0, nrow = (p + q + m + r + w*n), ncol = (p + q + m + r + w*n))
  #Hessian = matrix(0, nrow = (p + q + m + r + n), ncol = (p + q + m + r + n))
  
  phi.mu.1 = phi.mu[,re_ind]
  phi.mu.summed  = data.frame(id_long, phi.mu.1)
  phi.mu.summed = as.matrix(aggregate( phi.mu.summed[,2:(w+1)], list(phi.mu.summed[,1]), FUN = sum )[,-1])
  phi.mu.summed_alpha  = data.frame(id_long, phi.mu)
  phi.mu.summed_alpha = as.matrix(aggregate( phi.mu.summed_alpha[,2:(r+1)], list(phi.mu.summed_alpha[,1]), FUN = sum )[,-1])
  
  
  #least squares
  zt_i_samp_est = NULL
  for(i in 1:n){
    alpha_re = alpha
    alpha_re[re_ind] = alpha_re[re_ind] + unlist(c(a_re_cols[i,re_ind]))
    ind.long = which(id_long == i)
    
    zt_i_samp_est = c(zt_i_samp_est, phi.mu[ind.long,] %*% (alpha_re)) 
  }
  
  eXtB_e = exp(fixed_e %*% beta)
  eXtB_r = exp(fixed_r %*% beta)
  eXtB_l = exp(fixed_l %*% beta)
  eXtB_i = exp(fixed_i %*% beta)
  
  h0_t_e = psi_e %*% theta
  
  h0_t_quad = quad.psi.event %*% theta
  exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  h0_t_star_quad = h0_t_quad * exp_zTg_quad
  h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = 15, byrow = FALSE)
  H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
  
  H0_t_e = H0_t[dat.baseline$event == 1]
  H0_t_r = H0_t[dat.baseline$right == 1]
  H0_t_l = H0_t[dat.baseline$left == 1]
  H0_t_i2 = H0_t[dat.baseline$interval == 1]
  
  h0_t_quad_i1 = quad.psi.y_i1 %*% theta
  exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
  h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
  h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
  H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
  
  H0_t_i1_long = rep(0, n)
  H0_t_i1_long[which(dat.baseline$interval == 1)] = H0_t_i1
  
  S_t = exp(-eXtB * H0_t)
  S_t_i1_long = rep(0, n)
  
  S_t_l = exp(-eXtB_l * H0_t_l)
  S_t_i1 = exp(-eXtB_i * H0_t_i1)
  S_t_i2 = exp(-eXtB_i * H0_t_i2)
  
  S_t = exp(-eXtB * H0_t)
  S_t_i1_long = rep(0, n)
  S_t_i1_long[which(dat.baseline$interval == 1)] = S_t_i1
  
  #first derivative of gamma
  exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  A_t_quad = h0_t_quad * exp_zTg_quad * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
  A_t_quad = matrix(A_t_quad, ncol = n, nrow = 15, byrow = FALSE)
  A_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE)  *  A_t_quad, 2, sum)
  
  A_t_e = A_t[dat.baseline$event == 1]
  A_t_r = A_t[dat.baseline$right == 1]
  A_t_l = A_t[dat.baseline$left == 1]
  A_t_i2 = A_t[dat.baseline$interval == 1]
  
  exp_zTg_quad_i1 = exp(c(gamma) * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum))
  A_t_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1 * apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)
  A_t_quad_i1 = matrix(A_t_quad_i1, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
  A_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  A_t_quad_i1, 2, sum)
  
  A_t_il_long = rep(0, n)
  A_t_il_long[which(dat.baseline$interval == 1)] = A_t_i1
  
  #second derivative of gamma
  z_t_quad = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
  A2_t_quad_long = h0_t_quad * exp_zTg_quad * (z_t_quad)^2
  A_t_quad2 = matrix(A2_t_quad_long, ncol = n, nrow = 15, byrow = FALSE)
  A2_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = 15, byrow = FALSE) *  A_t_quad2, 2, sum)
  
  A_t_quad2_e = A2_t[dat.baseline$event == 1]
  A_t_quad2_r = A2_t[dat.baseline$right == 1]
  A_t_quad2_l = A2_t[dat.baseline$left == 1]
  A_t_quad2_i2 = A2_t[dat.baseline$interval == 1]
  
  z_t_quad_i1 = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)
  A_t_quad_i1_long = h0_t_quad_i1 * exp_zTg_quad_i1 * (z_t_quad_i1)^2
  A_t_quad2_i1 = matrix(A_t_quad_i1_long, ncol = sum(dat.baseline$interval), nrow = 15, byrow = FALSE)
  A_t_quad2_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(dat.baseline$interval)), nrow = 15, byrow = FALSE) *  A_t_quad2_i1, 2, sum)
  
  A2_t_i1_long = rep(0, length(A2_t))
  A2_t_i1_long[which(dat.baseline$interval == 1)] = A_t_quad2_i1
  
  #d RE d2 gamma
  A_dRE_t_quad2_long = c(c(gamma) * h0_t_quad * exp_zTg_quad * (z_t_quad)^2 + 2 * h0_t_quad * exp_zTg_quad * z_t_quad) * quad.phi.event[,re_ind]
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* A_dRE_t_quad2_long)
  A2_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
  
  A2_dRE_t_e = A2_dRE[dat.baseline$event == 1,]
  A2_dRE_t_r = A2_dRE[dat.baseline$right == 1,]
  A2_dRE_t_l = A2_dRE[dat.baseline$left == 1,]
  A2_dRE_t_i2 = A2_dRE[dat.baseline$interval == 1,]
  
  A_dRE_t_quad_long = c(h0_t_quad * exp_zTg_quad + h0_t_quad * exp_zTg_quad * z_t_quad * c(gamma)) * quad.phi.event[,re_ind]
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* A_dRE_t_quad_long)
  A_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
  
  #first derivative of theta
  psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* psi_t_star_quad)
  Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
  
  Psi_t_e = Psi_t[dat.baseline$event == 1,]
  Psi_t_r = Psi_t[dat.baseline$right == 1,]
  Psi_t_l = Psi_t[dat.baseline$left == 1,]
  Psi_t_i2 = Psi_t[dat.baseline$interval == 1,]
  
  psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
  df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* psi_t_star_quad_i1)
  Psi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
  
  Psi_t_i1_long = matrix(0, nrow = nrow(Psi_t), ncol = ncol(Psi_t))
  Psi_t_i1_long[which(dat.baseline$interval==1),] = Psi_t_i1
  
  #first derivative of alpha
  B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* B_t_quad)
  B_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
  
  Bphi_t_e = B_t[dat.baseline$event == 1,]
  Bphi_t_r = B_t[dat.baseline$right == 1,]
  Bphi_t_l = B_t[dat.baseline$left == 1,]
  Bphi_t_i2 = B_t[dat.baseline$interval == 1,]
  
  B_t_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * quad.phi.y_i1
  df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* B_t_quad_i1)
  Bphi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
  
  B_t_i1_long = matrix(0, nrow = nrow(B_t), ncol = ncol(B_t))
  B_t_i1_long[which(dat.baseline$interval == 1),] = Bphi_t_i1
  
  #second derivative of alpha
  B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
  B_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, 15)) * B_t_quad 
  B_t_quad2_e = B_t_quad2[c(sapply(dat.baseline$event == 1, rep, 15)),]
  B_t_quad2_r = B_t_quad2[c(sapply(dat.baseline$right == 1, rep, 15)),]
  B_t_quad2_l = B_t_quad2[c(sapply(dat.baseline$left == 1, rep, 15)),]
  B_t_quad2_i2 = B_t_quad2[c(sapply(dat.baseline$interval == 1, rep, 15)),]
  B_t_quad_i1_2 = c(rep(quad.w, sum(dat.baseline$interval == 1))) * 
    c(sapply(quad.lambda_i1, rep, 15)) * B_t_quad_i1 
  
  # gamma alpha hessian
  z_est_quad = apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
  E_t_quad = c(gamma) *  c(h0_t_quad * exp_zTg_quad * z_est_quad ) * quad.phi.event + 
    c(h0_t_quad * exp_zTg_quad ) * quad.phi.event
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, 15)), c(rep(quad.w, n))* E_t_quad)
  E_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
  
  Ephi_t_e = E_t[dat.baseline$event == 1,]
  Ephi_t_r = E_t[dat.baseline$right == 1,]
  Ephi_t_l = E_t[dat.baseline$left == 1,]
  Ephi_t_i2 = E_t[dat.baseline$interval == 1,]
  
  E_t_quad_i1 = c(gamma) *  c(h0_t_quad_i1 * exp_zTg_quad_i1 *apply((matrix(rep(alpha, 15*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad)[i1_ind_quad ==1,] * quad.phi.y_i1, 1, sum)) * quad.phi.y_i1 + 
    c(h0_t_quad_i1 * exp_zTg_quad_i1 ) * quad.phi.y_i1
  df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* E_t_quad_i1)
  Ephi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
  
  Ephi_t_i1_long = matrix(0, nrow = nrow(E_t), ncol = ncol(E_t))
  Ephi_t_i1_long[which(dat.baseline$interval == 1),] = Ephi_t_i1
  
  # theta gamma hessian (need to multiply this with z_quad)
  psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
  psi_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, 15)) * psi_t_star_quad 
  psi_t_quad2_e = psi_t_quad2[c(sapply(dat.baseline$event == 1, rep, 15)),]
  psi_t_quad2_r = psi_t_quad2[c(sapply(dat.baseline$right == 1, rep, 15)),]
  psi_t_quad2_l = psi_t_quad2[c(sapply(dat.baseline$left == 1, rep, 15)),]
  psi_t_quad2_i2 = psi_t_quad2[c(sapply(dat.baseline$interval == 1, rep, 15)),]
  
  psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
  psi_t_quad2_i1 = c(rep(quad.w, sum(dat.baseline$interval == 1))) * 
    c(sapply(quad.lambda_i1, rep, 15)) * psi_t_star_quad_i1 
  
  psi_t_i1_quad2 = matrix(0, nrow = nrow(psi_t_quad2), ncol = ncol(psi_t_quad2))
  psi_t_i1_quad2[i1_ind_quad,] = psi_t_quad2_i1
  
  #theta alpha hessian
  ##HERE NEED TO MULTIPLY TOGETHER FIRST DERIVATIVE OF THETA WITH PHI MATRIX IN MULTIPLICATION PART
  
  
  
  
  #beta
  #beta beta
  beta_hess_e = t(fixed_e) %*% diag(c(-eXtB_e * H0_t_e)) %*% fixed_e
  beta_hess_r =  t(fixed_r) %*% diag(c(-eXtB_r * H0_t_r)) %*% fixed_r
  beta_hess_l = t(fixed_l) %*% diag(c(((eXtB_l*H0_t_l)*S_t_l - (eXtB_l*H0_t_l)^2*S_t_l - eXtB_l*H0_t_l*S_t_l^2)  / (1 - S_t_l)^2)) %*% fixed_l
  beta_hess_i = t(fixed_i) %*% diag(c(-eXtB_i * (S_t_i1 * H0_t_i1 ) / (S_t_i1 - S_t_i2))) %*% fixed_i + 
    t(fixed_i) %*% diag(c(eXtB_i * (S_t_i2 * H0_t_i2 ) / (S_t_i1 - S_t_i2))) %*% fixed_i -
    t(fixed_i) %*% diag(c(S_t_i1 * S_t_i2 * (eXtB_i * H0_t_i1 - eXtB_i *  H0_t_i2)^2/ (S_t_i1 - S_t_i2)^2)) %*% fixed_i
  beta_hess = beta_hess_e + beta_hess_r + beta_hess_l + beta_hess_i
  Hessian[1:p, 1:p] = beta_hess
  
  #beta gamma
  betagamma_hess_e = t(fixed_e) %*% (c(-eXtB_e * A_t_e))
  betagamma_hess_r = t(fixed_r) %*% (c(-eXtB_r * A_t_r))
  betagamma_hess_l = t(fixed_l) %*% c(eXtB_l * (S_t_l * A_t_l)/ (1 - S_t_l)) -
    t(fixed_l) %*% c(eXtB_l^2 * S_t_l * A_t_l * H0_t_l/ (1 - S_t_l)^2)
  betagamma_hess_i = t(fixed_i) %*% c(-eXtB_i * S_t_i1 * A_t_i1 /(S_t_i1 - S_t_i2)) +
    t(fixed_i) %*% c(eXtB_i * S_t_i2 * A_t_i2 /(S_t_i1 - S_t_i2)) -
    t(fixed_i) %*% c(eXtB_i * S_t_i1 * S_t_i2 * (A_t_i1 - A_t_i2) * (H0_t_i1 * eXtB_i - H0_t_i2 * eXtB_i)/(S_t_i1 - S_t_i2)^2)
  betagamma_hess = betagamma_hess_e + betagamma_hess_r + betagamma_hess_l + betagamma_hess_i
  Hessian[1:p, (p+1):(p+q)] = betagamma_hess
  Hessian[(p+1):(p+q), 1:p] = t(Hessian[1:p, (p+1):(p+q)])
  
  #beta theta
  betatheta_hess_e = t(c(-eXtB_e) * fixed_e) %*% Psi_t_e
  betatheta_hess_r = t(c(-eXtB_r) * fixed_r) %*% Psi_t_r
  betatheta_hess_l = t(fixed_l * c(eXtB_l*(S_t_l - eXtB_l*H0_t_l*S_t_l - S_t_l^2)  / (1 - S_t_l)^2)) %*% Psi_t_l
  betatheta_hess_i = t(-fixed_i * c(eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2))) %*% Psi_t_i1 + 
    t(fixed_i * c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2))) %*% Psi_t_i2 -
    t(fixed_i * c(eXtB_i^2 * S_t_i1 * S_t_i2 * (H0_t_i1 - H0_t_i2)/(S_t_i1 - S_t_i2)^2)) %*% (Psi_t_i1 - Psi_t_i2) 
  betatheta_hess = betatheta_hess_e + betatheta_hess_r + betatheta_hess_l + betatheta_hess_i
  Hessian[1:p, (p+q+1):(p+q+m)] = betatheta_hess
  Hessian[(p+q+1):(p+q+m), 1:p] = t(Hessian[1:p, (p+q+1):(p+q+m)])
  
  #beta alpha
  betaalpha_hess_e = t(c(gamma) * c(-eXtB_e) * fixed_e) %*% Bphi_t_e
  betaalpha_hess_r = t(c(gamma) * c(-eXtB_r) * fixed_r) %*% Bphi_t_r
  betaalpha_hess_l = t(c(gamma) * fixed_l * c(eXtB_l * (S_t_l - eXtB_l*H0_t_l*S_t_l - S_t_l^2)  / (1 - S_t_l)^2) ) %*% Bphi_t_l
  betaalpha_hess_i = t(fixed_i) %*% (c(gamma) * c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)) * Bphi_t_i1) + 
    t(fixed_i) %*% (c(gamma) * c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)) * Bphi_t_i2) -
    t(c(gamma) * fixed_i * c(eXtB_i * S_t_i1 * S_t_i2 * (eXtB_i*H0_t_i1 - eXtB_i*H0_t_i2)/(S_t_i1 - S_t_i2)^2) ) %*% (Bphi_t_i1 - Bphi_t_i2)
  betaalpha_hess = betaalpha_hess_e + betaalpha_hess_r + betaalpha_hess_l + betaalpha_hess_i
  Hessian[1:p, (p+q+m+1):(p+q+m+r)] = betaalpha_hess
  Hessian[(p+q+m+1):(p+q+m+r), 1:p] = t(Hessian[1:p, (p+q+m+1):(p+q+m+r)])
  
  #beta kappa
  betaa1_hess = do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
    c(c(dat.baseline$event + dat.baseline$right) * c(gamma) * c(eXtB) * B_t[,re_ind]) +
    do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
    c(c(dat.baseline$left) *  c(gamma) * c(eXtB * S_t/(1-S_t)) * B_t[,re_ind]) -
    do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
    c(c(dat.baseline$left) *  c(gamma) * c(eXtB * S_t * (H0_t*eXtB)/(1-S_t)^2) * B_t[,re_ind]) -
    do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
    c(c(dat.baseline$interval) *  c(gamma) * c(eXtB * S_t_i1_long/(S_t_i1_long-S_t)) * B_t_i1_long[,re_ind]) + 
    do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
    c(c(dat.baseline$interval) *  c(gamma) * c(eXtB * S_t/(S_t_i1_long-S_t)) * B_t[,re_ind]) - 
    do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
    c(c(dat.baseline$interval) *  c(gamma) * c(eXtB * S_t_i1_long * S_t * (H0_t_i1_long*eXtB - H0_t*eXtB)/(S_t_i1_long-S_t)^2) * (B_t_i1_long[,re_ind] - B_t[,re_ind]))
  Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)] = betaa1_hess
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), 1:p] = t(Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)])
  
  #gamma
  #gamma gamma
  gamma_hess_e = t(-eXtB_e) %*% A_t_quad2_e
  gamma_hess_r = t(-eXtB_r) %*% A_t_quad2_r
  gamma_hess_l = t((S_t_l) * eXtB_l/(1-S_t_l)) %*% A_t_quad2_l - 
    t(c( S_t_l* A_t_l * eXtB_l^2/(1-S_t_l)^2)) %*% A_t_l
  gamma_hess_i = t(-eXtB_i * S_t_i1 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i1 +
    t(eXtB_i* S_t_i2 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i2 -
    t(c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2 ) * (A_t_i1 - A_t_i2)) %*% (A_t_i1 - A_t_i2)
  gamma_hess = gamma_hess_e + gamma_hess_r + gamma_hess_l + gamma_hess_i
  Hessian[(p+1):(p+q), (p+1):(p+q)] = gamma_hess
  
  #gamma theta
  if(sum(dat.baseline$event) > 0){
    gammatheta_hess_e = t(c(sapply(c(-eXtB_e), rep, 15)) * psi_t_quad2_e)%*% (z_t_quad)[c(sapply(dat.baseline$event == 1, rep, 15))]
  }else{
    gammatheta_hess_e = 0
  }
  gammatheta_hess_r = t(c(sapply(c(-eXtB_r), rep, 15)) * psi_t_quad2_r) %*% (z_t_quad)[c(sapply(dat.baseline$right == 1, rep, 15))]
  gammatheta_hess_l = t(c(sapply(c(eXtB_l * (S_t_l)/(1-S_t_l)), rep, 15)) *psi_t_quad2_l) %*% ( z_t_quad)[dat.baseline$left == 1] -
    t(Psi_t_l) %*% (c(eXtB_l^2 *(S_t_l)/(1-S_t_l)^2)* A_t_l) 
  gammatheta_hess_i = t(c(sapply(c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)), rep, 15)) *psi_t_quad2_i1) %*% (z_t_quad)[dat.baseline$interval == 1] + 
    t(c(sapply(c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)), rep, 15)) *psi_t_quad2_i2) %*% (z_t_quad)[dat.baseline$interval == 1] -
    t(Psi_t_i1 - Psi_t_i2) %*% (c(eXtB_i^2 * (A_t_i1 - A_t_i2) * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2) ) 
  gammatheta_hess = gammatheta_hess_e + gammatheta_hess_r + gammatheta_hess_l + gammatheta_hess_i
  Hessian[(p+1):(p+q), (p+q+1):(p+q+m)] = gammatheta_hess
  Hessian[(p+q+1):(p+q+m), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+1):(p+q+m)])
  
  #gamma alpha
  if(sum(dat.baseline$event)>0){
    gammaalpha_hess_e = apply(phi_e, 2, sum) - t(eXtB_e) %*% (Ephi_t_e)
  }else{
    gammaalpha_hess_e = 0
  }
  gammaalpha_hess_r = t(-eXtB_r) %*% (Ephi_t_r)
  gammaalpha_hess_l = t(eXtB_l * (S_t_l )/((1-S_t_l))) %*% (Ephi_t_l) -
    t(as.numeric(gamma) * as.numeric(eXtB_l^2 * A_t_l * S_t_l/((1-S_t_l)^2))) %*% Bphi_t_l
  gammaalpha_hess_i = t(-eXtB_i * S_t_i1/(S_t_i1 - S_t_i2)) %*% (Ephi_t_i1) +
    t(eXtB_i * S_t_i2/(S_t_i1 - S_t_i2)) %*% (Ephi_t_i2) -
    t(as.numeric(gamma) * as.numeric(eXtB_i^2 * S_t_i1 * S_t_i2 * (A_t_i1 - A_t_i2)/(S_t_i1 - S_t_i2)^2)) %*% (Bphi_t_i1-Bphi_t_i2)
  gammaalpha_hess = gammaalpha_hess_e + gammaalpha_hess_r + gammaalpha_hess_l + gammaalpha_hess_i
  Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)] = gammaalpha_hess
  Hessian[(p+q+m+1):(p+q+m+r), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)])
  
  #gamma kappa
  gammaa1_hess = c(dat.baseline$event * phi[,re_ind] - c(dat.baseline$event + dat.baseline$right) * (c(eXtB) * E_t[,re_ind])) +
    c(c(dat.baseline$left) * (c(eXtB * S_t/(1-S_t)) * E_t[,re_ind])) -
    c(c(dat.baseline$left) * c(gamma) * (c((eXtB)^2 * A_t * S_t/(1-S_t)^2) * B_t_re[,re_ind])) - 
    c(c(dat.baseline$interval) * (c(eXtB * S_t_i1_long/(S_t_i1_long-S_t)) * Ephi_t_i1_long[,re_ind])) + 
    c(c(dat.baseline$interval) * (c(eXtB * S_t/(S_t_i1_long-S_t)) * E_t[,re_ind])) - 
    c(c(dat.baseline$interval) * c(gamma) * (c((eXtB)^2 * (A_t_il_long - A_t) * S_t_i1_long * S_t/(S_t_i1_long-S_t)^2) * (B_t_i1_long[,re_ind] - B_t_re[,re_ind])))
  Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)] = gammaa1_hess
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)])
  
  #theta
  theta_hess_e = - t(as.numeric(1/(h0_t_e^2)) * psi_e) %*% psi_e
  theta_hess_l = - t(as.numeric(eXtB_l^2 * S_t_l/(1 - S_t_l)^2) * Psi_t_l) %*% Psi_t_l
  theta_hess_i = - t(as.numeric(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2) * (Psi_t_i1 - Psi_t_i2)) %*% (Psi_t_i1 - Psi_t_i2)
  theta_hess = theta_hess_e + theta_hess_l + theta_hess_i
  Hessian[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = theta_hess
  
  #theta alpha
  if(sum(dat.baseline$event)>0){
    thetaalpha_hess_e = t(c(sapply(c(gamma) * c(-eXtB_e), rep, 15)) *psi_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, 15)),]
  }else{
    thetaalpha_hess_e = 0
  }
  
  thetaalpha_hess_r = t(c(sapply(c(gamma) * c(-eXtB_r), rep, 15)) * psi_t_quad2_r) %*% (c(rep(quad.w, sum(dat.baseline$right))) * c(sapply(quad.lambda[dat.baseline$right == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$right == 1, rep, 15)),])
  thetaalpha_hess_l = t(c(sapply(c(gamma) * c(eXtB_l * (S_t_l)/((1-S_t_l))), rep, 15)) * psi_t_quad2_l) %*% (c(rep(quad.w, sum(dat.baseline$left))) * c(sapply(quad.lambda[dat.baseline$left == 1], rep, 15)) * quad.phi.event[dat.baseline$left == 1,]) -
    t(c(gamma) * c(eXtB_l^2 * S_t_l/((1-S_t_l)^2)) * Psi_t_l) %*% Bphi_t_l
  thetaalpha_hess_i = t(c(sapply(c(gamma) * c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)), rep, 15)) * psi_t_quad2_i1) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda_i1, rep, 15)) * quad.phi.y_i1) +
    t(c(sapply(c(gamma) * c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)), rep, 15)) * psi_t_quad2_i2) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda[dat.baseline$interval == 1], rep, 15)) * quad.phi.event[dat.baseline$interval == 1,]) -
    t(c(c(gamma) * eXtB_i^2 * S_t_i1 * S_t_i2/((S_t_i1-S_t_i2)^2)) * (Psi_t_i1 - Psi_t_i2)) %*% (Bphi_t_i1 - Bphi_t_i2)
  Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)] = thetaalpha_hess_e + thetaalpha_hess_r + thetaalpha_hess_l + thetaalpha_hess_i
  Hessian[(p+q+m+1):(p+q+m+r), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)])
  
  #theta a1
  temp = c(quad.phi.event[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) psi_t_quad2))
  df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, 15)), w),
                  re_no = c(sapply(re_ind, rep, 15*n)), c(rep(quad.w, n)) * temp)
  F_t = c(quad.lambda) * as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  temp = c(quad.phi.y_i1[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) psi_t_quad2_i1))
  df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)),
                  re_no = c(sapply(re_ind, rep, 15*sum(dat.baseline$interval))), c(rep(quad.w, sum(dat.baseline$interval == 1))) * temp)
  F_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  F_t_i1_long = matrix(0, nrow = nrow(F_t), ncol = ncol(F_t))
  F_t_i1_long[rep(which(dat.baseline$interval == 1), w),] = F_t_i1
  
  thetaa1_hess = rep((- c(dat.baseline$event + dat.baseline$right) * eXtB * c(gamma)), w) * F_t + 
    rep((c(dat.baseline$left) * eXtB * S_t/(1-S_t) * c(gamma)), w) * F_t -
    rep((c(dat.baseline$left) * (eXtB^2) * S_t/((1-S_t)^2) * c(gamma)), w) * c(B_t_re) * do.call(rbind, lapply(seq_len(w), function(X) Psi_t)) -
    rep((c(dat.baseline$interval) * eXtB * S_t_i1_long/(S_t_i1_long-S_t) * c(gamma)), w) * F_t_i1_long + 
    rep((c(dat.baseline$interval) * eXtB * S_t/(S_t_i1_long-S_t) * c(gamma)), w) * F_t -
    rep((c(dat.baseline$interval) * (eXtB^2) * S_t * S_t_i1_long/((S_t_i1_long-S_t)^2) * c(gamma)), w) * c(B_t_re_i1_long - B_t_re) * do.call(rbind, lapply(seq_len(w), function(X) (Psi_t_i1_long - Psi_t)))
  
  Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)] = t(thetaa1_hess)
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)])
  
  
  #alpha
  if(sum(dat.baseline$event)>0){
    alpha_hess_e = t(-c(sapply(c(gamma^2) * c(eXtB_e), rep, 15)) * B_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, 15)),]
  }else{
    alpha_hess_e = 0
  }
  alpha_hess_r = t(-c(sapply(c(gamma^2) * c(eXtB_r), rep, 15)) * B_t_quad2_r) %*% (c(rep(quad.w, sum(dat.baseline$right))) * c(sapply(quad.lambda[dat.baseline$right == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$right == 1, rep, 15)),])
  alpha_hess_l = t(c(sapply(c(gamma^2) * c(eXtB_l*(S_t_l)/(1-S_t_l)), rep, 15)) * B_t_quad2_l) %*% (c(rep(quad.w, sum(dat.baseline$left))) * c(sapply(quad.lambda[dat.baseline$left == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$left == 1, rep, 15)),]) -
    t(c(c(gamma^2) * (eXtB_l^2 * S_t_l/(1-S_t_l)^2)) * Bphi_t_l) %*% Bphi_t_l
  alpha_hess_i = t(-c(sapply(c(gamma^2) * c(eXtB_i*S_t_i1 /(S_t_i1-S_t_i2)), rep, 15)) * B_t_quad_i1_2) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda_i1, rep, 15)) * quad.phi.y_i1)  + 
    t(c(sapply(c(gamma^2) * c(eXtB_i*S_t_i2 /(S_t_i1-S_t_i2)), rep, 15)) * B_t_quad2_i2) %*% (c(rep(quad.w, sum(dat.baseline$interval))) * c(sapply(quad.lambda[dat.baseline$interval == 1], rep, 15)) * quad.phi.event[c(sapply(dat.baseline$interval == 1, rep, 15)),]) -
    t(c(c(gamma^2) * (eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i1 + 
    t(c(c(gamma^2) * (eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i2
  Hessian[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = alpha_hess_e + alpha_hess_r + alpha_hess_l + alpha_hess_i
  
  #alpha a1
  temp = c(quad.phi.event[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) B_t_quad2))
  df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, 15)), w),
                  re_no = c(sapply(re_ind, rep, 15*n)), c(rep(quad.w, n)) * temp)
  D_t_re = c(quad.lambda) * as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  temp = c(quad.phi.y_i1[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) B_t_quad2_i1))
  df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, 15)),
                  re_no = c(sapply(re_ind, rep, 15*sum(dat.baseline$interval))), c(rep(quad.w, sum(dat.baseline$interval == 1))) * temp)
  D_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  D_t_re_i1_long = matrix(0, nrow = nrow(D_t_re), ncol = ncol(D_t_re))
  D_t_re_i1_long[rep(which(dat.baseline$interval == 1), w),] = D_t_re_i1
  
  alphaa1_hess = - rep((c(dat.baseline$event + dat.baseline$right) * eXtB * c(gamma^2)), w) * D_t_re + 
    rep((c(dat.baseline$left) * eXtB * (S_t/(1-S_t)) * c(gamma^2)), w) * D_t_re - 
    c(rep((c(dat.baseline$left) * (eXtB^2) * (S_t/(1-S_t)^2) * c(gamma^2)), w) * B_t_re) * do.call(rbind, lapply(seq_len(w), function(X) B_t)) -
    rep((c(dat.baseline$interval) * eXtB * (S_t_i1_long/(S_t_i1_long-S_t)) * c(gamma^2)), w) * D_t_re_i1_long + 
    rep((c(dat.baseline$interval) * eXtB * (S_t/(S_t_i1_long-S_t)) * c(gamma^2)), w) * D_t_re - 
    c(rep((c(dat.baseline$interval) * (eXtB^2) * (S_t_i1_long*S_t/(S_t_i1_long-S_t)^2) * c(gamma^2)), w) * (B_t_re_i1_long - B_t_re)) * do.call(rbind, lapply(seq_len(w), function(X) (B_t_i1_long - B_t)))
  
  Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)] = t(alphaa1_hess)
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+1):(p+q+m+r)] = t(Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)])
  
  #a1 a1
  #temp = (c(quad.phi.event[,re_ind]) * 
  #          do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma^2), rep, 15)) * B_t_quad2[,re_ind])))
  #df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, 15)), w),
  #                re_no = c(sapply(re_ind, rep, 15*n)), temp)
  a1a1_hess = alphaa1_hess[,re_ind]
  
  for(re in 1:w){
    diag(Hessian[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = as.matrix(a1a1_hess)[((re-1)*n + 1):((re-1)*n + n),re]
    
    if(re > 1){
      diag(Hessian[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess[((re-1)*n + 1):((re-1)*n + n), (re-1)]
      diag(Hessian[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n)]) = diag(Hessian[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)])
    }
    
  }
  
  #add penalty part
  Q_theta = Q_alpha = Q_a1 = H_epsilon = matrix(0, nrow = (p+q+m+r+w*n), ncol = (p+q+m+r+w*n))
  
  #least squares part
  H_epsilon[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = - c(1/sigma2_Et) * t(phi.mu) %*% phi.mu #alpha alpha
  
  temp = c(phi.mu[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) phi.mu))
  df = data.frame(id_long,
                  re_no = c(sapply(re_ind, rep, nrow(dat.long))), temp)
  alphaa1_hess_ls = - c(1/sigma2_Et) *  as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  H_epsilon[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)] = t(alphaa1_hess_ls) #alpha re
  H_epsilon[ (p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+1):(p+q+m+r)] = t(H_epsilon[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)])
  
  temp = c(phi.mu[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) phi.mu[,re_ind]))
  df = data.frame(id_long,
                  re_no = c(sapply(re_ind, rep, nrow(dat.long))), temp)
  a1a1_hess_ls = - c(1/sigma2_Et) *  as.matrix(aggregate( df[,3:(w+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  for(re in 1:w){
    diag(H_epsilon[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n),re]
    
    if(re > 1){
      diag(H_epsilon[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n), (re-1)]
      diag(H_epsilon[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n)]) = diag(H_epsilon[(p+q+m+r+(re-2)*n + 1):(p+q+m+r+(re-2)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)])
    }
    
  }
  
  #theta part
  Q_theta[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = - theta.lambda*theta.G
  
  #alpha part
  Q_alpha[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = - alpha.lambda*alpha.G
  
  #re part
  diag(Q_a1[(p+q+m+r+1):(p+q+m+r+w*n),(p+q+m+r+1):(p+q+m+r+w*n)]) = -c(1/c(sapply(sigma2_re, rep, n)))
  
  H_full = Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1
  
  
  pos = rep(1, nrow(H_full))
  for(u in 1:m){
    if((theta[u]< (1e-3) & theta_score[u] < (0))){
      pos[p+q+u] = 0
    }
  }
  
  
  #calculate asymptotic variance
  #M2_trunc = (-(Hessian[1:(p+q+m+r),1:(p+q+m+r)] + H_epsilon[1:(p+q+m+r),1:(p+q+m+r)] + 
  #                Q_theta[1:(p+q+m+r),1:(p+q+m+r)] + Q_alpha[1:(p+q+m+r),1:(p+q+m+r)]))
  #
  
  corr = M2_inv = matrix(0, nrow(H_full), nrow(H_full))
  diag(corr) = pos
  #corr[!pos,]=0
  
  M2_inv[(pos==1), (pos==1)] = solve(-H_full[(pos==1), (pos==1)])
  
  A_omega = corr %*% M2_inv %*% t(corr)
  
  cov_H = A_omega %*% (-((Hessian + H_epsilon))) %*% A_omega
  cov_H_RE = A_omega %*% (-((Hessian + H_epsilon + Q_a1))) %*% A_omega
  
  cov_H_RE_2 = solve((-((Hessian + H_epsilon + Q_theta + Q_alpha)))[1:(p+q+m+r), 1:(p+q+m+r)] - 
                       (-((Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1)))[1:(p+q+m+r), (p+q+m+r+1):(p+q+m+r+n*w)] %*%
                       solve((-((Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1)))[(p+q+m+r+1):(p+q+m+r+n*w), (p+q+m+r+1):(p+q+m+r+n*w)]) %*%
                       (-((Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1)))[(p+q+m+r+1):(p+q+m+r+n*w), 1:(p+q+m+r)])
  
  
  
  #param = c(beta, gamma, theta, alpha)
  
  #output
  pars = list(beta = beta, gamma = gamma, theta = theta, alpha = c(alpha), 
              a_re_cols = a_re_cols)
  score = list(beta_score = beta_score, gamma_score = gamma_score, theta_score = theta_score,
               alpha_score = alpha_score)
  knots = list(theta.knots = event.kn)
  variance = list(sigma2_re = sigma2_re, sigma2_Et = sigma2_Et, 
                  theta.lambda = theta.lambda, alpha.lambda = alpha.lambda)
  penalty = list(R = theta.G)
  out = list(pars = pars, knots = knots, variance = variance, penalty = penalty, 
             score = score, cov_H = cov_H, M2_inv = M2_inv, cov_H_RE = cov_H_RE, H_full = H_full, cov_H_RE_2 = cov_H_RE_2)
  return(out)
  
}



