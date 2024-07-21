jm_fitting_right_profile_new = function(dat.baseline, dat.long, bmm, 
                                        max.iter = c(10, 1000), 
                                        n.knots.h0 = 2, step_size = 1, mod_mat_long_f,
                                        re_ind, range, init_re, init_sigma_eps, init_sigma_re,
                                        gauss_rules){
  
  
  kappa = 1/0.5
  
  event.time = dat.baseline$latest_time
  obs.time = dat.long$obs_time_samp
  id = dat.baseline$id
  id_long = dat.long$id_long
  cont = dat.long$wt_samp
  fixed = bmm
  W_long = (dat.long$W_long)
  
  fixed_e = (fixed[which(dat.baseline$event == 1),])
  fixed_r = (fixed[which(dat.baseline$right == 1),])
  
  #initial set up
  
  #m-splines for baseline hazard function
  int.knots.event = quantile(c(obs.time, event.time)[(c(obs.time, event.time) > 0)], seq(range[1], range[2], length.out=n.knots.h0))
  bound.knots.event = c(- 1e-4, max(c(obs.time, event.time)[(c(obs.time, event.time) > 0)]) + 1e-4)
  event.kn = list(int = int.knots.event, bound = bound.knots.event)
  
  psi_e = mSpline(dat.baseline$t_L, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[dat.baseline$event == 1,]
  psi_r = mSpline(dat.baseline$t_L[dat.baseline$right == 1], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)
  
  tm = max(c(obs.time, event.time)[(c(obs.time, event.time) > 0)])
  #b-splines for TVC function
  phi.mu = mod_mat_long_f(obs.time, tmax = tm, W_long = W_long)
  
  #b-splines for TVC function
  phi = mod_mat_long_f(event.time, tmax = tm, W_long =  dat.baseline$fixed_x3)
  phi_e = phi[which(dat.baseline$event == 1),]
  phi_r = phi[dat.baseline$right == 1,]
  
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
  alpha = matrix(rep(0,r), ncol = 1)
  #alpha_v = matrix(rep(0,v), ncol = 1)
  
  a_re_cols = init_re
  a_re_long = unlist(c(init_re))
  
  a_re_cols_pad = matrix(rep(0,r*n), ncol = r)
  a_re_cols_pad_quad = matrix(rep(0,r*n*gauss_rules), ncol = r)
  a_re_cols_pad[,re_ind] = a_re_cols
  a_re_cols_pad_quad[,re_ind] = matrix(sapply(a_re_cols, rep, gauss_rules))
  
  i1_ind_quad = c(sapply(dat.baseline$interval, rep, gauss_rules))
  
  #variance components
  sigma2_Et = init_sigma_eps
  sigma2_re = c(init_sigma_re)
  
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
  quad.y = t(as.matrix(quad.lambda) %*% rules[[gauss_rules]]$x + quad.mu) #one row per time, one column per i
  quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event))
  quad.phi.event = mod_mat_long_f(c(quad.y), tmax = tm,W_long =  c(sapply(dat.baseline$fixed_x3, rep, gauss_rules)))
  quad.w = rules[[gauss_rules]]$w
  
  for(it in 1:max.iter[1]){
    for(iter in 1:max.iter[2]){
      
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      
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
      log_lik = pl_e + pl_r - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #beta
      
      if(sum(dat.baseline$event) > 0){
        beta_score_e = t(fixed_e) %*% c(1- eXtB_e * H0_t_e)
      }else{
        beta_score_e = 0
      }
      beta_score_r = - t(fixed_r) %*% c(eXtB_r * H0_t_r)
      
      beta_score = beta_score_e + beta_score_r
      
      beta_hess_e = t(fixed_e) %*% diag(c(-eXtB_e * H0_t_e)) %*% fixed_e
      beta_hess_r = t(fixed_r) %*% diag(c(-eXtB_r * H0_t_r)) %*% fixed_r
      
      beta_hess = beta_hess_e + beta_hess_r
      beta_neg_hess = -beta_hess
      
      beta_old = beta
      
      beta = beta_old + solve(beta_neg_hess)%*%beta_score
      
      #step size for beta
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1,] %*% alpha_v))
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1,] %*% alpha_v))
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #print(c("beta", log_lik_old, log_lik))
      step_size_beta = 1
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
          #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1,] %*% alpha_v))
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1,] %*% alpha_v))
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #log likelihood
          log_lik = pl_e + pl_r - pl_ls - 
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
          if(i>1000){break}
          
          
        }
        #print(c(omega1))
        
      }
      
      
      #alpha
      B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event # h0(t) * exp(zTg) * phi(t) for each i at each quad t, (n*gauss_rules) x r
      df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* B_t_quad) #multiply above with quad.w (weight) for each t, then each row is labelled with i, (n*gauss_rules) x (r+1)
      B_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1]) # aggregate() sums down each column for each i (removing i label column), then multiply each i with its own quad.lamba
      
      Bphi_t_e = B_t[dat.baseline$event == 1,]
      Bphi_t_r = B_t[dat.baseline$right == 1,]
      
      alpha_score_e = apply(as.numeric(gamma) * phi_e - as.numeric(gamma) * as.numeric(eXtB_e) * Bphi_t_e, 2, sum)
      alpha_score_r = apply(- as.numeric(gamma) * as.numeric(eXtB_r) * Bphi_t_r, 2, sum)
      alpha_score_ls = as.numeric( t(as.matrix((cont - z_t_ia))) %*% phi.mu * (1/sigma2_Et))
      
      alpha_score = alpha_score_e + alpha_score_r + alpha_score_ls - 2*alpha.lambda * alpha.G %*% alpha
      
      #B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
      B_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, gauss_rules)) * B_t_quad
      B_t_quad2_e = B_t_quad2[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
      B_t_quad2_r = B_t_quad2[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
      
      if(sum(dat.baseline$event)>0){
        alpha_hess_e = t(-c(sapply(c(gamma^2) * c(eXtB_e), rep, gauss_rules)) * B_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
      }else{
        alpha_hess_e = 0
      }
      
      alpha_hess_r = t(-c(sapply(c(gamma^2) * c(eXtB_r), rep, gauss_rules)) * B_t_quad2_r) %*% quad.phi.event[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
      alpha_hess_partial_new = alpha_hess_e + alpha_hess_r
      
      alpha_hess = alpha_hess_partial_new  - (1/sigma2_Et) * t(phi.mu) %*% phi.mu - 2*alpha.lambda * alpha.G
      alpha_hess_neg = -alpha_hess
      
      alpha_old = alpha
      alpha = alpha_old + solve(alpha_hess_neg)%*%alpha_score
      
      #step size for alpha
      #likelihood
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      
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
      log_lik = pl_e + pl_r - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #print(c("alpha", log_lik_old, log_lik))
      
      if(((log_lik < log_lik_old) & step_size == 1)){
        i = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          alpha = alpha_old + omega1 * solve(alpha_hess_neg)%*%alpha_score
          
          #update log-likelihood
          #likelihood
          exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          
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
          log_lik = pl_e + pl_r - pl_ls - 
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
          if(i>1000){break}
          
          
        }
      }
      
      #gamma
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      A_t_quad = h0_t_quad * exp_zTg_quad * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
      A_t_quad = matrix(A_t_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
      A_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE) *  A_t_quad, 2, sum)
      
      A_t_e = A_t[dat.baseline$event == 1]
      A_t_r = A_t[dat.baseline$right == 1]
      
      gamma_score_e = sum(zt_e - eXtB_e * A_t_e)
      gamma_score_r = sum(- eXtB_r * A_t_r)
      
      gamma_score = gamma_score_e + gamma_score_r 
      
      
      #second derivative of gamma
      z_t_quad = apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
      A2_t_quad_long = h0_t_quad * exp_zTg_quad * (z_t_quad)^2
      A_t_quad2 = matrix(A2_t_quad_long, ncol = n, nrow = gauss_rules, byrow = FALSE)
      A2_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE) *  A_t_quad2, 2, sum)
      
      A_t_quad2_e = A2_t[dat.baseline$event == 1]
      A_t_quad2_r = A2_t[dat.baseline$right == 1]
      
      gamma_hess_e = t(-eXtB_e) %*% A_t_quad2_e
      gamma_hess_r = t(-eXtB_r) %*% A_t_quad2_r
      
      gamma_hess = gamma_hess_e + gamma_hess_r
      gamma_hess_neg = -gamma_hess
      
      gamma_old = gamma
      gamma = gamma_old + solve(gamma_hess_neg)%*%gamma_score
      
      #step size for gamma
      #likelihood
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1,] %*% alpha_v))
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1,] %*% alpha_v))
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #print(c("gamma", log_lik_old, log_lik))
      step_size_gamma = 1
      if(((log_lik < log_lik_old) & step_size_gamma == 1)){
        i = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          gamma = gamma_old +  omega1 * solve(gamma_hess_neg)%*%gamma_score
          
          #likelihood
          exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1,] %*% alpha_v))
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum) 
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1,] %*% alpha_v))
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #log likelihood
          log_lik = pl_e + pl_r - pl_ls - 
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
          if(i>1000){break}
          
          
        }
      }
      
      
      #theta
      
      #switch Psi1 and Psi2 calculation to use 0 -> tL and tL -> tR like cum haz
      psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
      df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* psi_t_star_quad)
      Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
      
      Psi_t_e = Psi_t[dat.baseline$event == 1,]
      Psi_t_r = Psi_t[dat.baseline$right == 1,]
      
      TwoLRtheta = as.numeric(2*theta.lambda)*(theta.G%*%theta)
      
      theta_score_neg = apply(as.numeric(eXtB_e) * Psi_t_e, 2, sum) + 
        apply(as.numeric(eXtB_r)* Psi_t_r, 2, sum) +
        TwoLRtheta*(TwoLRtheta>0) + 1e-3
      
      theta_score_pos = apply(as.numeric(1/h0_t_e)* psi_e, 2, sum) - 
        TwoLRtheta*(TwoLRtheta<0) + 1e-3
      
      theta_score = (theta_score_pos - theta_score_neg)
      
      D_matrix = theta/(theta_score_neg)
      
      theta_old = theta
      theta = theta_old + D_matrix*(theta_score)
      
      
      #step size for theta
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1,] %*% alpha_v))
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1,] %*% alpha_v))
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_e + pl_r - pl_ls - 
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
          #eXtB = exp(fixed %*% beta)
          
          h0_t_quad = quad.psi.event %*% theta
          #exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1,] %*% alpha_v))
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1,] %*% alpha_v))
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #log likelihood
          log_lik = pl_e + pl_r - pl_ls - 
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
          if(i>1000){break}
          
          
        }
      }
      
      
      a_re_cols_old = a_re_cols
      
      #random effects
      
      B_t_re_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event[,1:w]
      df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* B_t_re_quad)
      B_t_re = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
      
      Bphi_t_re_e = B_t_re[dat.baseline$event == 1,]
      Bphi_t_re_r = B_t_re[dat.baseline$right == 1,]
      
      a1i_score_e = as.numeric(gamma) * phi_e[,1:w] - as.numeric(gamma) * as.numeric(eXtB_e) * (Bphi_t_re_e)
      a1i_score_r = - as.numeric(gamma) * as.numeric(eXtB_r) * (Bphi_t_re_r)
      
      id_censor_order = c(id[dat.baseline$event == 1], id[dat.baseline$right == 1])
      a1i_score_s = rbind(a1i_score_e, a1i_score_r)
      a1i_score_s = data.frame(id_censor_order, a1i_score_s)[order(id_censor_order),-1]
      
      phi.mu.1 = phi.mu[,1:w]
      a1i_score_ls = data.frame(id_long, (cont-z_t_ia)*phi.mu.1)
      a1i_score_ls = as.matrix(aggregate( a1i_score_ls[,2:(w+1)], list(a1i_score_ls[,1]), FUN = sum )[,-1])
      
      a1i_score = a1i_score_s + (1/sigma2_Et)*a1i_score_ls - a_re_cols/matrix(rep(sigma2_re, n), ncol = w, byrow = T)
      
      a_cols_update = matrix(0, nrow = n, ncol = w)
      
      for(i in 1:n){
        
        ind.long = which(id_long == i)
        a_hess_i = t(rep(- as.numeric(gamma^2) * as.numeric(eXtB[i] ), gauss_rules) * B_t_re_quad[(gauss_rules*(i-1)+1):(gauss_rules*(i-1)+gauss_rules),]) %*%  quad.phi.event[(gauss_rules*(i-1)+1):(gauss_rules*(i-1)+gauss_rules),1:w] - 
          t((1/sigma2_Et) * phi.mu.1[ind.long,]) %*% (phi.mu.1[ind.long,]) - 
          diag(1/sigma2_re)
        
        
        a_cols_update[i,] = solve(-a_hess_i)%*%t(as.matrix(a1i_score[i,]))
        
      }
      
      a_re_cols = a_re_cols_old + a_cols_update
      
      a_re_long_old = a_re_long
      a_re_cols_pad_old = a_re_cols_pad
      a_re_cols_pad_quad_old = a_re_cols_pad_quad
      
      a_re_long = unlist(c(a_re_cols))
      a_re_cols_pad[,re_ind] = as.matrix(a_re_cols)
      a_re_cols_pad_quad[,re_ind] = matrix(sapply(a_re_cols, rep, gauss_rules), ncol = 2)
      
      
      #step size for random effects
      #likelihood
      exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[dat.baseline$event == 1]
      H0_t_r = H0_t[dat.baseline$right == 1]
      
      #event
      zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
      zt_e = zt_e[which(dat.baseline$event == 1)]
      zt_g_e = zt_e %*% gamma
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
      
      #right
      pl_r = - sum(eXtB_r * H0_t_r)
      
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
      log_lik = pl_e + pl_r - pl_ls - 
        sum((1/(2*sigma2_re)) * apply(a_re_cols * a_re_cols, 2, sum)) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha)%*%alpha.G%*%alpha
      
      #print(c("re", log_lik_old, log_lik))
      step_size_a1 = 1
      if(((log_lik < log_lik_old) & step_size_a1 == 1)){
        i = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          
          a_re_cols = a_re_cols_old + omega1 * a_cols_update
          
          a_re_long = unlist(c(a_re_cols))
          a_re_cols_pad[,re_ind] = as.matrix(a_re_cols)
          a_re_cols_pad_quad[,re_ind] = matrix(sapply(a_re_cols, rep, gauss_rules), ncol = 2)
          
          #step size for random effects
          #likelihood
          exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[dat.baseline$event == 1]
          H0_t_r = H0_t[dat.baseline$right == 1]
          
          #event
          zt_e = apply((matrix(rep(alpha, n), ncol = r, byrow = TRUE) + a_re_cols_pad) * phi, 1, sum)
          zt_e = zt_e[which(dat.baseline$event == 1)]
          zt_g_e = zt_e %*% gamma
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + zt_g_e - eXtB_e*H0_t_e)
          
          #right
          pl_r = - sum(eXtB_r * H0_t_r)
          
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
          log_lik = pl_e + pl_r - pl_ls - 
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
          if(i>100){break}
          
          
        }
        #print(c(omega1))
        
      }
      
      #hist(a_re_cols[,1])
      #hist(a_re_cols[,2])
      
      pos_temp = rep(1, m)
      for(u in 1:m){
        if((theta[u]< (1e-2) & theta_score[u] < (-1e-2))){
          pos_temp[u] = 0
        }
      }
      
      #print(c(iter, beta, gamma, alpha, log_lik))
      
      print(c(iter, theta_score_pos, theta_score_neg, theta_score_pos/theta_score_neg, log_lik))
      
      if(all(abs(c(beta - beta_old, gamma -gamma_old, alpha - alpha_old, theta - theta_old )) < 1e-6) & 
         all(abs(c(a_re_cols - a_re_cols_old)) < 1e-5)){
        break
      }
      
      #if(it == 1 & iter == 2){
      #  break
      #}
      
      
      
    } #end inner loop
    
    
    print(c(iter, beta, gamma, log_lik))
    
    
    
    Hessian = matrix(0, nrow = (p + q + m + r + w*n), ncol = (p + q + m + r + w*n))
    #Hessian = matrix(0, nrow = (p + q + m + r + n), ncol = (p + q + m + r + n))
    
    phi.mu.summed  = data.frame(id_long, phi.mu.1)
    phi.mu.summed = as.matrix(aggregate( phi.mu.summed[,2:(w+1)], list(phi.mu.summed[,1]), FUN = sum )[,-1])
    phi.mu.summed_alpha  = data.frame(id_long, phi.mu)
    phi.mu.summed_alpha = as.matrix(aggregate( phi.mu.summed_alpha[,2:(r+1)], list(phi.mu.summed_alpha[,1]), FUN = sum )[,-1])
    
    
    #least squares
    zt_i_samp_est = NULL
    for(i in 1:n){
      alpha_re = alpha
      alpha_re[re_ind] = alpha_re[re_ind] + a_re_cols[i,re_ind]
      ind.long = which(id_long == i)
      
      zt_i_samp_est = c(zt_i_samp_est, phi.mu[ind.long,] %*% (alpha_re)) 
    }
    
    eXtB_e = exp(fixed_e %*% beta)
    eXtB_r = exp(fixed_r %*% beta)
    
    #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
    #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
    
    
    h0_t_e = psi_e %*% theta
    
    h0_t_quad = quad.psi.event %*% theta
    exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
    h0_t_star_quad = h0_t_quad * exp_zTg_quad
    h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
    H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
    
    H0_t_e = H0_t[dat.baseline$event == 1]
    H0_t_r = H0_t[dat.baseline$right == 1]
    
    #first derivative of gamma
    exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
    A_t_quad = h0_t_quad * exp_zTg_quad * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
    A_t_quad = matrix(A_t_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
    A_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  A_t_quad, 2, sum)
    
    A_t_e = A_t[dat.baseline$event == 1]
    A_t_r = A_t[dat.baseline$right == 1]
    
    #second derivative of gamma
    z_t_quad = apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
    A2_t_quad_long = h0_t_quad * exp_zTg_quad * (z_t_quad)^2
    A_t_quad2 = matrix(A2_t_quad_long, ncol = n, nrow = gauss_rules, byrow = FALSE)
    A2_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE) *  A_t_quad2, 2, sum)
    
    A_t_quad2_e = A2_t[dat.baseline$event == 1]
    A_t_quad2_r = A2_t[dat.baseline$right == 1]
    
    #first derivative of theta
    psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* psi_t_star_quad)
    Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
    
    Psi_t_e = Psi_t[dat.baseline$event == 1,]
    Psi_t_r = Psi_t[dat.baseline$right == 1,]
    
    #first derivative of alpha
    B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* B_t_quad)
    B_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    Bphi_t_e = B_t[dat.baseline$event == 1,]
    Bphi_t_r = B_t[dat.baseline$right == 1,]
    
    #second derivative of alpha
    B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
    B_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, gauss_rules)) * B_t_quad 
    B_t_quad2_e = B_t_quad2[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
    B_t_quad2_r = B_t_quad2[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
    
    # gamma alpha hessian
    E_t_quad = c(gamma) *  c(h0_t_quad * exp_zTg_quad *apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)) * quad.phi.event + 
      c(h0_t_quad * exp_zTg_quad ) * quad.phi.event
    df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* E_t_quad)
    E_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    Ephi_t_e = E_t[dat.baseline$event == 1,]
    Ephi_t_r = E_t[dat.baseline$right == 1,]
    
    # theta gamma hessian (need to multiply this with z_quad)
    psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
    psi_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, gauss_rules)) * psi_t_star_quad 
    psi_t_quad2_e = psi_t_quad2[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
    psi_t_quad2_r = psi_t_quad2[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
    
    #theta alpha hessian
    ##HERE NEED TO MULTIPLY TOGETHER FIRST DERIVATIVE OF THETA WITH PHI MATRIX IN MULTIPLICATION PART
    
    #beta
    #beta beta
    beta_hess_e = t(fixed_e) %*% diag(c(-eXtB_e *  H0_t_e)) %*% fixed_e
    beta_hess_r =  t(fixed_r) %*% diag(c(-eXtB_r * H0_t_r)) %*% fixed_r
    beta_hess = beta_hess_e + beta_hess_r
    Hessian[1:p, 1:p] = beta_hess
    
    #beta gamma
    betagamma_hess_e = t(fixed_e) %*% (c(-eXtB_e *  A_t_e))
    betagamma_hess_r = t(fixed_r) %*% (c(-eXtB_r * A_t_r))
    betagamma_hess = betagamma_hess_e + betagamma_hess_r
    Hessian[1:p, (p+1):(p+q)] = betagamma_hess
    Hessian[(p+1):(p+q), 1:p] = t(Hessian[1:p, (p+1):(p+q)])
    
    #beta theta
    betatheta_hess_e = t(c(-eXtB_e ) * fixed_e) %*% Psi_t_e
    betatheta_hess_r = t(c(-eXtB_r ) * fixed_r) %*% Psi_t_r
    betatheta_hess = betatheta_hess_e + betatheta_hess_r
    Hessian[1:p, (p+q+1):(p+q+m)] = betatheta_hess
    Hessian[(p+q+1):(p+q+m), 1:p] = t(Hessian[1:p, (p+q+1):(p+q+m)])
    
    #beta alpha
    betaalpha_hess_e = t(c(gamma) * c(-eXtB_e ) * fixed_e) %*% Bphi_t_e
    betaalpha_hess_r = t(c(gamma) * c(-eXtB_r ) * fixed_r) %*% Bphi_t_r
    betaalpha_hess = betaalpha_hess_e + betaalpha_hess_r
    Hessian[1:p, (p+q+m+1):(p+q+m+r)] = betaalpha_hess
    Hessian[(p+q+m+1):(p+q+m+r), 1:p] = t(Hessian[1:p, (p+q+m+1):(p+q+m+r)])
    
    #beta kappa
    betaa1_hess = do.call(rbind, lapply(seq_len(w), function(X) fixed))*
      c(- c(gamma) * c(eXtB) * B_t[,re_ind])
    Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)] = betaa1_hess
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), 1:p] = t(Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)])
    
    #gamma
    #gamma gamma
    gamma_hess_e = t(-eXtB_e) %*% A_t_quad2_e
    gamma_hess_r = t(-eXtB_r) %*% A_t_quad2_r
    gamma_hess = gamma_hess_e + gamma_hess_r
    Hessian[(p+1):(p+q), (p+1):(p+q)] = gamma_hess
    
    #gamma theta
    if(sum(dat.baseline$event) > 0){
      gammatheta_hess_e = t(c(sapply(c(-eXtB_e), rep, gauss_rules)) *psi_t_quad2_e)%*% ( z_t_quad)[c(sapply(dat.baseline$event == 1, rep, gauss_rules))]
    }else{
      gammatheta_hess_e = 0
    }
    gammatheta_hess_r = t(c(sapply(c(-eXtB_r), rep, gauss_rules)) *psi_t_quad2_r) %*% ( z_t_quad)[c(sapply(dat.baseline$right == 1, rep, gauss_rules))]
    gammatheta_hess = gammatheta_hess_e + gammatheta_hess_r
    Hessian[(p+1):(p+q), (p+q+1):(p+q+m)] = gammatheta_hess
    Hessian[(p+q+1):(p+q+m), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+1):(p+q+m)])
    
    #gamma alpha
    if(sum(dat.baseline$event)>0){
      gammaalpha_hess_e = t(matrix(1, nrow(eXtB_e))) %*% phi_e - t(eXtB_e) %*% (Ephi_t_e)
    }else{
      gammaalpha_hess_e = 0
    }
    gammaalpha_hess_r = t(-eXtB_r) %*% (Ephi_t_r)
    gammaalpha_hess = gammaalpha_hess_e + gammaalpha_hess_r
    Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)] = gammaalpha_hess
    Hessian[(p+q+m+1):(p+q+m+r), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)])
    
    #gamma kappa
    gammaa1_hess = c(dat.baseline$event * phi[,re_ind] - (c(eXtB) * E_t[,re_ind]))
    Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)] = gammaa1_hess
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)])
    
    #theta
    theta_hess_e = - t(as.numeric(1/(h0_t_e^2)) * psi_e) %*% psi_e
    theta_hess = theta_hess_e
    Hessian[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = theta_hess
    
    #theta alpha
    if(sum(dat.baseline$event)>0){
      thetaalpha_hess_e = t(c(sapply(c(gamma) * c(-eXtB_e), rep, gauss_rules)) *psi_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
    }else{
      thetaalpha_hess_e = 0
    }
    thetaalpha_hess_r = t(c(sapply(c(gamma) * c(-eXtB_r), rep, gauss_rules)) *psi_t_quad2_r) %*% quad.phi.event[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
    Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)] = thetaalpha_hess_e + thetaalpha_hess_r
    Hessian[(p+q+m+1):(p+q+m+r), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)])
    
    #theta a1
    temp = c(quad.phi.event[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma), rep, gauss_rules)) * psi_t_quad2))
    df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, gauss_rules)), w),
                    re_no = c(sapply(re_ind, rep, gauss_rules*n)), temp)
    thetaa1_hess = as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)] = t(thetaa1_hess)
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)])
    
    #alpha
    if(sum(dat.baseline$event)>0){
      alpha_hess_e = t(-c(sapply(c(gamma^2) * c(eXtB_e), rep, gauss_rules)) * B_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
    }else{
      alpha_hess_e = 0
    }
    alpha_hess_r = t(-c(sapply(c(gamma^2) * c(eXtB_r), rep, gauss_rules)) * B_t_quad2_r) %*% quad.phi.event[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
    Hessian[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = alpha_hess_e + alpha_hess_r
    
    
    #alpha a1
    temp = c(quad.phi.event[,re_ind]) * 
      do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma^2), rep, gauss_rules)) * B_t_quad2))
    df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, gauss_rules)), w),
                    re_no = c(sapply(re_ind, rep, gauss_rules*n)), temp)
    alphaa1_hess = as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)] = t(alphaa1_hess)
    Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+1):(p+q+m+r)] = t(Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)])
    
    #a1 a1
    temp = (c(quad.phi.event[,re_ind]) * 
              do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma^2), rep, gauss_rules)) * B_t_quad2[,re_ind])))
    df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, gauss_rules)), w),
                    re_no = c(sapply(re_ind, rep, gauss_rules*n)), temp)
    a1a1_hess = as.matrix(aggregate( df[,3:(w+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    for(re in 1:w){
      diag(Hessian[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess[((re-1)*n + 1):((re-1)*n + n),re]
      
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
  
  h0_t_e = psi_e %*% theta
  
  h0_t_quad = quad.psi.event %*% theta
  exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  h0_t_star_quad = h0_t_quad * exp_zTg_quad
  h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
  H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  h0_t_star_quad, 2, sum)
  
  H0_t_e = H0_t[dat.baseline$event == 1]
  H0_t_r = H0_t[dat.baseline$right == 1]
  
  #first derivative of gamma
  exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  A_t_quad = h0_t_quad * exp_zTg_quad * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
  A_t_quad = matrix(A_t_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
  A_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE)  *  A_t_quad, 2, sum)
  
  A_t_e = A_t[dat.baseline$event == 1]
  A_t_r = A_t[dat.baseline$right == 1]
  
  #second derivative of gamma
  z_t_quad = apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)
  A2_t_quad_long = h0_t_quad * exp_zTg_quad * (z_t_quad)^2
  A_t_quad2 = matrix(A2_t_quad_long, ncol = n, nrow = gauss_rules, byrow = FALSE)
  A2_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = gauss_rules, byrow = FALSE) *  A_t_quad2, 2, sum)
  
  A_t_quad2_e = A2_t[dat.baseline$event == 1]
  A_t_quad2_r = A2_t[dat.baseline$right == 1]
  
  #d RE d2 gamma
  A_dRE_t_quad2_long = c(c(gamma) * h0_t_quad * exp_zTg_quad * (z_t_quad)^2 + 2 * h0_t_quad * exp_zTg_quad * z_t_quad) * quad.phi.event[,re_ind]
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* A_dRE_t_quad2_long)
  A2_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
  
  A2_dRE_t_e = A2_dRE[dat.baseline$event == 1,]
  A2_dRE_t_r = A2_dRE[dat.baseline$right == 1,]
  
  A_dRE_t_quad_long = c(h0_t_quad * exp_zTg_quad + h0_t_quad * exp_zTg_quad * z_t_quad * c(gamma)) * quad.phi.event[,re_ind]
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* A_dRE_t_quad_long)
  A_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
  
  #first derivative of theta
  psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* psi_t_star_quad)
  Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
  
  Psi_t_e = Psi_t[dat.baseline$event == 1,]
  Psi_t_r = Psi_t[dat.baseline$right == 1,]
  
  #first derivative of alpha
  B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* B_t_quad)
  B_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
  
  Bphi_t_e = B_t[dat.baseline$event == 1,]
  Bphi_t_r = B_t[dat.baseline$right == 1,]
  
  #second derivative of alpha
  B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
  B_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, gauss_rules)) * B_t_quad 
  B_t_quad2_e = B_t_quad2[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
  B_t_quad2_r = B_t_quad2[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
  
  #second derivative of alpha - version 2
  B_t_quad2_v2 = c(rep(quad.w, n)) * B_t_quad 
  B_t_quad2_e_v2 = B_t_quad2_v2[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
  B_t_quad2_r_v2 = B_t_quad2_v2[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
  
  # gamma alpha hessian
  E_t_quad = c(gamma) *  c(h0_t_quad * exp_zTg_quad *apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)) * quad.phi.event + 
    c(h0_t_quad * exp_zTg_quad ) * quad.phi.event
  df = data.frame(id_quad = c(sapply(dat.baseline$id, rep, gauss_rules)), c(rep(quad.w, n))* E_t_quad)
  E_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
  
  Ephi_t_e = E_t[dat.baseline$event == 1,]
  Ephi_t_r = E_t[dat.baseline$right == 1,]
  
  # theta gamma hessian (need to multiply this with z_quad)
  psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
  psi_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, gauss_rules)) * psi_t_star_quad 
  psi_t_quad2_e = psi_t_quad2[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
  psi_t_quad2_r = psi_t_quad2[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
  
  #theta alpha hessian
  ##HERE NEED TO MULTIPLY TOGETHER FIRST DERIVATIVE OF THETA WITH PHI MATRIX IN MULTIPLICATION PART
  
  
  #beta
  #beta beta
  beta_hess_e = t(fixed_e) %*% diag(c(-eXtB_e * H0_t_e)) %*% fixed_e
  beta_hess_r =  t(fixed_r) %*% diag(c(-eXtB_r * H0_t_r)) %*% fixed_r
  beta_hess = beta_hess_e + beta_hess_r
  Hessian[1:p, 1:p] = beta_hess
  
  #beta gamma
  betagamma_hess_e = t(fixed_e) %*% (c(-eXtB_e * A_t_e))
  betagamma_hess_r = t(fixed_r) %*% (c(-eXtB_r * A_t_r))
  betagamma_hess = betagamma_hess_e + betagamma_hess_r 
  Hessian[1:p, (p+1):(p+q)] = betagamma_hess
  Hessian[(p+1):(p+q), 1:p] = t(Hessian[1:p, (p+1):(p+q)])
  
  #beta theta
  betatheta_hess_e = t(c(-eXtB_e) * fixed_e) %*% Psi_t_e
  betatheta_hess_r = t(c(-eXtB_r) * fixed_r) %*% Psi_t_r
  betatheta_hess = betatheta_hess_e + betatheta_hess_r
  Hessian[1:p, (p+q+1):(p+q+m)] = betatheta_hess
  Hessian[(p+q+1):(p+q+m), 1:p] = t(Hessian[1:p, (p+q+1):(p+q+m)])
  
  #beta alpha
  betaalpha_hess_e = t(c(gamma) * c(-eXtB_e) * fixed_e) %*% Bphi_t_e
  betaalpha_hess_r = t(c(gamma) * c(-eXtB_r) * fixed_r) %*% Bphi_t_r
  betaalpha_hess = betaalpha_hess_e + betaalpha_hess_r
  Hessian[1:p, (p+q+m+1):(p+q+m+r)] = betaalpha_hess
  Hessian[(p+q+m+1):(p+q+m+r), 1:p] = t(Hessian[1:p, (p+q+m+1):(p+q+m+r)])
  
  #beta kappa
  betaa1_hess = do.call(rbind, lapply(seq_len(w), function(X) fixed)) * 
    c(c(dat.baseline$event + dat.baseline$right) * c(gamma) * c(eXtB) * B_t[,re_ind])
  Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)] = betaa1_hess
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), 1:p] = t(Hessian[1:p, (p+q+m+r+1):(p+q+m+r+w*n)])
  
  #gamma
  #gamma gamma
  gamma_hess_e = t(-eXtB_e) %*% A_t_quad2_e
  gamma_hess_r = t(-eXtB_r) %*% A_t_quad2_r
  gamma_hess = gamma_hess_e + gamma_hess_r
  Hessian[(p+1):(p+q), (p+1):(p+q)] = gamma_hess
  
  #gamma theta
  if(sum(dat.baseline$event) > 0){
    gammatheta_hess_e = t(c(sapply(c(-eXtB_e), rep, gauss_rules)) *psi_t_quad2_e)%*% (c(exp_zTg_quad) * z_t_quad)[c(sapply(dat.baseline$event == 1, rep, gauss_rules))]
  }else{
    gammatheta_hess_e = 0
  }
  gammatheta_hess_r = t(c(sapply(c(-eXtB_r), rep, gauss_rules)) *psi_t_quad2_r) %*% (c(exp_zTg_quad) * z_t_quad)[c(sapply(dat.baseline$right == 1, rep, gauss_rules))]
  gammatheta_hess = gammatheta_hess_e + gammatheta_hess_r 
  Hessian[(p+1):(p+q), (p+q+1):(p+q+m)] = gammatheta_hess
  Hessian[(p+q+1):(p+q+m), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+1):(p+q+m)])
  
  #gamma alpha
  if(sum(dat.baseline$event)>0){
    gammaalpha_hess_e = t(matrix(1, nrow(eXtB_e))) %*% phi_e - t(eXtB_e) %*% (Ephi_t_e)
  }else{
    gammaalpha_hess_e = 0
  }
  gammaalpha_hess_r = t(-eXtB_r) %*% (Ephi_t_r)
  gammaalpha_hess = gammaalpha_hess_e + gammaalpha_hess_r
  Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)] = gammaalpha_hess
  Hessian[(p+q+m+1):(p+q+m+r), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r)])
  
  #gamma kappa
  gammaa1_hess = c(dat.baseline$event * phi[,re_ind] - (c(eXtB) * E_t[,re_ind]))
  Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)] = gammaa1_hess
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w*n)])
  
  #theta
  theta_hess_e = - t(as.numeric(1/(h0_t_e^2)) * psi_e) %*% psi_e
  theta_hess = theta_hess_e 
  Hessian[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = theta_hess
  
  #theta alpha
  if(sum(dat.baseline$event)>0){
    thetaalpha_hess_e = t(c(sapply(c(gamma) * c(-eXtB_e), rep, gauss_rules)) * psi_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
  }else{
    thetaalpha_hess_e = 0
  }
  thetaalpha_hess_r = t(c(sapply(c(gamma) * c(-eXtB_r), rep, gauss_rules)) * psi_t_quad2_r) %*% quad.phi.event[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
  Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)] = thetaalpha_hess_e + thetaalpha_hess_r
  Hessian[(p+q+m+1):(p+q+m+r), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r)])
  
  #theta a1
  temp = c(quad.phi.event[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma), rep, gauss_rules)) * psi_t_quad2))
  df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, gauss_rules)), w),
                  re_no = c(sapply(re_ind, rep, gauss_rules*n)), temp)
  thetaa1_hess = as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)] = t(thetaa1_hess)
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w*n)])
  
  #alpha
  if(sum(dat.baseline$event)>0){
    alpha_hess_e = t(-c(sapply(c(gamma^2) * c(eXtB_e), rep, gauss_rules)) * B_t_quad2_e) %*% quad.phi.event[c(sapply(dat.baseline$event == 1, rep, gauss_rules)),]
  }else{
    alpha_hess_e = 0
  }
  alpha_hess_r = t(-c(sapply(c(gamma^2) * c(eXtB_r), rep, gauss_rules)) * B_t_quad2_r) %*% quad.phi.event[c(sapply(dat.baseline$right == 1, rep, gauss_rules)),]
  Hessian[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = alpha_hess_e + alpha_hess_r
  
  #alpha a1
  temp = c(quad.phi.event[,re_ind]) * 
    do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma^2), rep, gauss_rules)) * B_t_quad2))
  df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, gauss_rules)), w),
                  re_no = c(sapply(re_ind, rep, gauss_rules*n)), temp)
  alphaa1_hess = as.matrix(aggregate( df[,3:(r+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)] = t(alphaa1_hess)
  Hessian[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+1):(p+q+m+r)] = t(Hessian[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w*n)])
  
  #a1 a1
  temp = (c(quad.phi.event[,re_ind]) * 
            do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma^2), rep, gauss_rules)) * B_t_quad2[,re_ind])))
  df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, gauss_rules)), w),
                  re_no = c(sapply(re_ind, rep, gauss_rules*n)), temp)
  a1a1_hess = as.matrix(aggregate( df[,3:(w+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  for(re in 1:w){
    diag(Hessian[(p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n), (p+q+m+r+(re-1)*n + 1):(p+q+m+r+(re-1)*n + n)]) = a1a1_hess[((re-1)*n + 1):((re-1)*n + n),re]
    
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
  
  #H_epsilon[(p+q+m+r+1):(p+q+m+r+w*n), (p+q+m+r+1):(p+q+m+r+w*n)] = -diag(as.numeric(c(phi.mu.summed^2)/sigma2_Et)) #re re
  
  #theta part
  Q_theta[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = - theta.lambda*theta.G
  
  #alpha part
  Q_alpha[(p+q+m+1):(p+q+m+r),(p+q+m+1):(p+q+m+r)] = - alpha.lambda*alpha.G
  
  #re part
  diag(Q_a1[(p+q+m+r+1):(p+q+m+r+w*n),(p+q+m+r+1):(p+q+m+r+w*n)]) = -c(1/c(sapply(sigma2_re, rep, n)))
  
  H_full = Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1
  
  #H_full_inv = solve(-(Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1)[1:(p+q+m+r),1:(p+q+m+r)])
  
  
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
  
  param = c(beta, gamma, theta, alpha)

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
             score = score, cov_H = cov_H, M2_inv = M2_inv, cov_H_RE = cov_H_RE, H_full = H_full)
  return(out)
  
}



