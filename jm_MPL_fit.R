
# control function ---- 
jm_MPL_control = function(kappa = 1/0.6, max.iter = c(1,1000), n.knots.h0 = 5, range = c(0.1,0.9),
                          par_initial = c(0, 0, 0), re_ind = NULL, step_size = 1, gs.rules = NULL, 
                          theta.lambda.init = 0, alpha.lambda.init = 0){
  
  #check omega's will be between 0 and 1
  if(1/kappa < 0 | 1/kappa > 1){
    stop("kappa value is mis-specified")
  }else{
    kappa = kappa
  }
  
  ## should put creation of gaussian quad rules in here too!!

  if(is.null(gs.rules)){
    gs.rules = 30
  }
  
  out = list(kappa = kappa, max.iter = max.iter, n.knots.h0 = n.knots.h0, range = range, 
             par_initial = par_initial, re_ind = re_ind, step_size = step_size, gs.rules = gs.rules, 
             theta.lambda.init = theta.lambda.init, alpha.lambda.init = alpha.lambda.init)
  return(out)
  
}


jm_MPL = function(time.formula, long.formula, data, id, time_func, control){
  
  mc = match.call(expand.dots = FALSE)
  m=match(c("time.formula","long.formula","data", "id"), names(mc),0)
  mc = mc[c(1,m)]
  if(m[1]==0){stop("A Cox regression formula argument is required.")}
  data.name = if(m[3]!=0){deparse(match.call()[[4]])}else{"-"}      
  form=lapply(list(mc$time.formula,mc$long.formula),as.formula)
  time_mf=model.frame(form[[1]], data)
  long_mf = list()
  for(cov in 1:length(form[[2]])){
    long_mf[[cov]] = model.frame(form[[2]][[cov]], data)
    
  }
  
  surv = model.extract(time_mf,"response")
  type=attr(surv,"class")
  if(!inherits(surv, "Surv")) {stop("Response must be a survival object.")}
  if(attr(surv,which="type")=="right"){
    left=surv[,1]
    right=rep(NA,nrow(surv))
    icase=which(surv[,2]==1)
    right[icase]=surv[icase,1]
    surv=Surv(left,right,type="interval2")
  }
  
  n_long = nrow(surv)
  n = length(unique(id))
  
  #set up baseline model matrix
  n.obs = (as.data.frame(id) %>% group_by(id) %>% tally())$n
  last_record = rep(0,n_long)
  
  for(ii in unique(id)){
    ind.long = which(id == ii)
    last_record[ind.long[n.obs[ii]]] = 1
    
  }
  
  fixed = as.matrix(time_mf[which(last_record == 1),-1])
  
  #set up censoring indicators
  censor_long = matrix(FALSE,nrow=n_long,ncol=4)
  for(ii in 1:n_long){
    censor_long[ii,c(surv[ii,3] + 1)]=TRUE
  }
  
  surv_short = surv[which(last_record == 1),]
  censor = matrix(FALSE,nrow=n,ncol=4)
  for(ii in 1:n){
    censor[ii,(surv_short[ii,3]+1)]=TRUE
  }
  
  fixed_e = (fixed[which(censor[,2] == 1),])
  fixed_r = (fixed[which(censor[,1] == 1),])
  fixed_l = (fixed[which(censor[,3] == 1),])
  fixed_i = (fixed[which(censor[,4] == 1),])
  
  #initial set up
  
  obs.time = long_mf[[1]][,2]
  cont = data.frame(lapply(long_mf, "[", , 1))
  ll_save = NULL
  
  #m-splines for baseline hazard function
  int.knots.event = quantile(c(surv_short[which(censor[,2]), 1],
                               surv_short[which(censor[,1]),1],
                               surv_short[which(censor[,3]),1]/2,
                               (surv_short[which(censor[,4]),1] + (surv_short[which(censor[,4]),2] - surv_short[which(censor[,4]),1])/2)), 
                             seq(control$range[1], control$range[2], length.out=control$n.knots.h0))
  bound.knots.event = c(- 1e-4, max(c(surv_short[,1], surv_short[which(censor[,4]),2], obs.time)) + 1e-4)
  event.kn = list(int = int.knots.event, bound = bound.knots.event)
  
  psi_e = mSpline(surv_short[, 1], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[which(censor[,2]),]
  psi_r = mSpline(surv_short[, 1], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[which(censor[,1]),]
  psi_l = mSpline(surv_short[, 1], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[which(censor[,3]),]
  psi_i1 = mSpline(surv_short[, 1], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[which(censor[,4]),]
  psi_i2 = mSpline(surv_short[, 2], degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event, intercept = F)[which(censor[,4]),]
  
  #model matrix for longitudinal variable - TVC function
  latest_time = surv_short[,1]
  latest_time[which(censor[,4])] = surv_short[which(censor[,4]),2]
  
  phi.mu = phi.mu.re = list()
  phi_s1 = phi_s1.re = list()
  phi_s2 = phi_s2.re = list()
  phi = phi.re = list()
  
  for(cov in 1:ncol(cont)){
    tf = time_func[[cov]]
    
    phi.mu[[cov]] = tf(long_mf[[cov]][,-1])[[1]]
    phi.mu.re[[cov]] = tf(long_mf[[cov]][,-1])[[2]]
    
    phi_s1[[cov]] = tf(cbind(surv_short[,1], long_mf[[cov]][which(last_record == 1),-c(1:2)]))[[1]]
    phi_s1.re[[cov]] = tf(cbind(surv_short[,1], long_mf[[cov]][which(last_record == 1),-c(1:2)]))[[2]]
    
    phi_s2[[cov]] = tf(cbind(surv_short[,2], long_mf[[cov]][which(last_record == 1),-c(1:2)]))[[1]]
    phi_s2.re[[cov]] = tf(cbind(surv_short[,2], long_mf[[cov]][which(last_record == 1),-c(1:2)]))[[2]]
    
    phi[[cov]] = tf(cbind(latest_time, long_mf[[cov]][which(last_record == 1),-c(1:2)]))[[1]]
    phi.re[[cov]] = tf(cbind(latest_time, long_mf[[cov]][which(last_record == 1),-c(1:2)]))[[2]]
    
    
  }
  
  
  phi_e = lapply(phi_s1, "[", which(censor[,2]),)
  phi_r = lapply(phi_s1, "[", which(censor[,1]),)
  phi_l = lapply(phi_s1, "[", which(censor[,3]),)
  phi_i1 = lapply(phi_s1, "[", which(censor[,4]),)
  phi_i2 = lapply(phi_s2, "[", which(censor[,4]),)
  
  phi_e.re = lapply(phi_s1.re, "[", which(censor[,2]),)
  phi_r.re =  lapply(phi_s1.re, "[", which(censor[,1]),)
  phi_l.re =  lapply(phi_s1.re, "[", which(censor[,3]),)
  phi_i1.re =  lapply(phi_s1.re, "[", which(censor[,4]),)
  phi_i2.re =  lapply(phi_s2.re, "[", which(censor[,4]),)
  
  #set up dimensions
  p = ncol(fixed)
  q = length(long_mf)
  m = length(int.knots.event) + 3
  r = lapply(phi.mu, ncol)
  w = lapply(phi.mu.re, ncol)
  
  #initalise parameter vectors
  beta = matrix(rep(0, p), ncol = 1)
  gamma = matrix(rep(0, q))
  theta = matrix(rep(0.1,m), ncol = 1)
  alpha = list()
  a_re_long = a_re_cols = a_re_cols_quad = sigma2_Et = sigma2_re = df.a_re = list()
  
  for(cov in 1:q){
    
    
    alpha[[cov]] = matrix(rep(0.1,r[[cov]]), ncol = 1)
    
    if(length(control$par_initial[[cov]][[1]]) == 1){
      a_re_long[[cov]] = rep(control$par_initial[[cov]][[1]], w[[cov]]*n)
      a_re_cols[[cov]] = matrix(control$par_initial[[cov]][[1]], ncol = w[[cov]], nrow = n)
      #a_re_cols_pad = matrix(rep(control$par_initial[[1]],r*n), ncol = r)
      a_re_cols_quad[[cov]] = matrix(rep(control$par_initial[[1]],w[[cov]]*n*control$gs.rules), ncol = w[[cov]])
    }
    else{
      a_re_long[[cov]] = c(control$par_initial[[cov]][[1]])
      a_re_cols[[cov]] = matrix(a_re_long[[cov]], ncol = w[[cov]])
      #a_re_cols_pad = matrix(rep(0,r*n), ncol = r)
      a_re_cols_quad[[cov]] = matrix(sapply(a_re_cols[[cov]], rep, control$gs.rules), ncol = w[[cov]])
      
    }
    
    
    #variance components
    sigma2_Et[[cov]] = control$par_initial[[cov]][[2]]
    
    if(length(control$par_initial[[cov]][[3]]) < w[[cov]]){
      sigma2_re[[cov]] = rep(control$par_initial[[cov]][[3]], w[[cov]])
    }else{
      sigma2_re[[cov]] = control$par_initial[[cov]][[3]]
    }
    
    
    df.a_re[[cov]] = rep(n, w[[cov]])
    
    
  }
  
  df.theta = n
  df.alpha = n
  df.epsilon = rep(n_long,q)
  
  #penalty matrices and smoothing parameters
  theta.G = theta_penalty_f(3, int.knots.event, bound.knots.event)
  #theta.G = NULL
  
  
  ### FIX THISSSS
  #C = NULL
  #for(uu in 1:(r[[1]]-2)){
  #  row = rep(0, r[[1]])
  #  if(uu >4){
  #    row[uu:(uu+2)] = c(1, -2, 1)
  #  }
  #  C = rbind(C, row)
  #}
  
  #theta.G = t(C) %*% C
  
  #alpha.G = zeta_penalty_f(3, int.knots.mu, bound.knots.mu)
  alpha.G = diag(c(rep(0, r[[1]])))
  #alpha.G = t(C) %*% C
  
  
  theta.lambda = control$theta.lambda.init
  alpha.lambda = control$alpha.lambda.init
  ### FIX THISSS
  
  
  i1_ind_quad = c(sapply(as.numeric(censor[,4]), rep, control$gs.rules))
  
  rules = gaussquad::legendre.quadrature.rules(control$gs.rules)
  #quadrature rules etc.
  
  quad.lambda = (latest_time - 0)/2
  quad.mu = (latest_time + 0)/2
  quad.y = t(as.matrix(quad.lambda) %*% rules[[control$gs.rules]]$x + quad.mu) #one row per time, one column per i
  quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event))
  quad.w = rules[[control$gs.rules]]$w
  
  quad.lambda_i1 = (surv_short[which(censor[,4]),1] - 0)/2
  quad.mu_i1 = (surv_short[which(censor[,4]),1] + 0)/2
  quad.y_i1 = t(as.matrix(quad.lambda_i1) %*% rules[[control$gs.rules]]$x + quad.mu_i1) #one row per time, one column per i
  quad.psi.y_i1 = t(sapply(quad.y_i1, mSpline, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event))
  i1_ind_quad = c(sapply(as.numeric(censor[,4]), rep, control$gs.rules))
  
  quad.lambda_i2 = (surv_short[which(censor[,4]),2] - surv_short[which(censor[,4]),1])/2
  quad.mu_i2 = (surv_short[which(censor[,4]),2] + surv_short[which(censor[,4]),1])/2
  quad.y_i2 = t(as.matrix(quad.lambda_i2) %*% rules[[control$gs.rules]]$x + quad.mu_i2) #one row per time, one column per i
  quad.psi.y_i2 = t(sapply(quad.y_i2, mSpline, degree = 3, knots = int.knots.event, Boundary.knots = bound.knots.event))
  
  
  quad.phi.event = quad.phi.event.re = quad.phi.y_i1 = quad.phi.y_i1.re = quad.phi.y_i2 = quad.phi.y_i2.re = list()
  
  
  for(cov in 1:q){
    tf = time_func[[cov]]
    
    quad.phi.event[[cov]] = tf(cbind(c(quad.y), c(sapply(long_mf[[cov]][which(last_record == 1),-c(1:2)], rep, control$gs.rules))))[[1]]
    quad.phi.event.re[[cov]] = tf(cbind(c(quad.y), c(sapply(long_mf[[cov]][which(last_record == 1),-c(1:2)], rep, control$gs.rules))))[[2]]
    
    quad.phi.y_i1[[cov]] = tf(cbind(c(quad.y_i1), c(sapply(long_mf[[cov]][which(last_record == 1 & censor_long[,4] == 1),-c(1:2)], rep, control$gs.rules))))[[1]]
    quad.phi.y_i1.re[[cov]] = tf(cbind(c(quad.y_i1), c(sapply(long_mf[[cov]][which(last_record == 1 & censor_long[,4] == 1),-c(1:2)], rep, control$gs.rules))))[[2]]
    
    quad.phi.y_i2[[cov]] = tf(cbind(c(quad.y_i2), c(sapply(long_mf[[cov]][which(last_record == 1 & censor_long[,4] == 1),-c(1:2)], rep, control$gs.rules))))[[1]]
    quad.phi.y_i2.re[[cov]] = tf(cbind(c(quad.y_i2), c(sapply(long_mf[[cov]][which(last_record == 1 & censor_long[,4] == 1),-c(1:2)], rep, control$gs.rules))))[[2]]
    
    
  }
  
  
  
  for(it in 1:control$max.iter[1]){
    for(iter in 1:control$max.iter[2]){
      
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
      nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
      z_t_quad = mu_t_quad + nu_t_quad
      exp_zTg_quad = exp(z_t_quad %*% gamma)
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[which(censor[,2])]
      H0_t_r = H0_t[which(censor[,1])]
      H0_t_l = H0_t[which(censor[,3])]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
      nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
      z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
      exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
      nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
      z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
      exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
      zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
      #zt_e = zt_e[which(censor[,2])]
      zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      diff_l = 1 - S_t_l
      diff_l[which(diff_l < 1e-5)] = 1e-5
      pl_l = sum(log(diff_l))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      diff_i = S_t_i1 - S_t_i2
      diff_i[which(diff_i < 1e-5)] = 1e-5
      pl_i = sum(log(diff_i))
      
      #least squares
      mu_t_ia = sapply(seq_along(phi.mu), function(x) phi.mu[[x]] %*% alpha[[x]])
      
      nu_t_ia = sapply(seq_along(phi.mu.re), function(x) apply(phi.mu.re[[x]] * sapply(data.frame(a_re_cols[[x]]), rep, n.obs), 1, sum))
      
      z_t_ia = mu_t_ia + nu_t_ia
      
      pl_ls = sum((sapply(seq_along(cont), function(x) (1/(2*sigma2_Et[[x]]))*sum((cont[[x]] - z_t_ia[,x])^2))))
      
      #log likelihood
      ll_surv =  pl_e + pl_r + pl_l + pl_i
      log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
        sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum)))) - 
        theta.lambda* t(theta)%*%theta.G%*%theta - 
        alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
      
      
      
      
      #alpha
      alpha_old = alpha
      for(cov in 1:q){
        #print(c("update alpha"))
        B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event[[cov]]
        df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* B_t_quad)
        B_t = gamma[cov] * c(quad.lambda) * as.matrix(aggregate( df[,2:(r[[cov]]+1)], list(df[,1]), FUN = sum )[,-1])
        
        
        Bphi_t_e = B_t[which(censor[,2]),]
        Bphi_t_r = B_t[which(censor[,1]),]
        Bphi_t_l = B_t[which(censor[,3]),]
        Bphi_t_i2 = B_t[which(censor[,4]),]
        
        B_t_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * quad.phi.y_i1[[cov]]
        df = data.frame(id_quad = c(sapply(unique(id[which(censor_long[,4])]), rep, control$gs.rules)), c(rep(quad.w, sum((censor[,4]))))* B_t_quad_i1) 
        Bphi_t_i1 = gamma[cov] *  c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r[[cov]]+1)], list(df[,1]), FUN = sum )[,-1])
        
        alpha_score_e = apply((phi_e[[cov]]) - as.numeric(eXtB_e) * Bphi_t_e, 2, sum)
        alpha_score_r = apply(-  as.numeric(eXtB_r) * Bphi_t_r, 2, sum)
        alpha_score_l = apply(as.numeric(eXtB_l * S_t_l/(1-S_t_l)) * Bphi_t_l, 2, sum)
        alpha_score_i = apply(- as.numeric(eXtB_i * S_t_i1 /(S_t_i1-S_t_i2)) *  Bphi_t_i1 +  
                                as.numeric(eXtB_i * S_t_i2 /(S_t_i1-S_t_i2)) *  Bphi_t_i2, 2, sum)
        
        alpha_score_ls = as.numeric( c(1/sigma2_Et[[cov]]) * t(as.matrix((cont[[cov]] - z_t_ia[,cov]))) %*% phi.mu[[cov]])
         
        alpha_score = alpha_score_e + alpha_score_r + alpha_score_l + alpha_score_i + alpha_score_ls #- 2*alpha.lambda * alpha.G %*% alpha
        
        #B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
        B_t_quad2 = c(gamma[cov]) * c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, control$gs.rules)) * B_t_quad
        B_t_quad2_e = B_t_quad2[c(sapply(censor[,2], rep, control$gs.rules)),]
        B_t_quad2_r = B_t_quad2[c(sapply(censor[,1], rep, control$gs.rules)),]
        B_t_quad2_l = B_t_quad2[c(sapply(censor[,3], rep, control$gs.rules)),]
        B_t_quad2_i2 = B_t_quad2[c(sapply(censor[,4], rep, control$gs.rules)),]
        B_t_quad2_i1 = c(gamma[cov]) * c(rep(quad.w, sum(censor[,4]))) * 
          c(sapply(quad.lambda_i1, rep, control$gs.rules)) * B_t_quad_i1
        
        if(sum(censor[,2])>0){
          alpha_hess_e = t(-c(sapply(c(eXtB_e), rep, control$gs.rules)) * B_t_quad2_e) %*% as.matrix(gamma[[cov]] * quad.phi.event[[cov]])[c(sapply(censor[,2], rep, control$gs.rules)),]
        }else{
          alpha_hess_e = 0
        }
        
        alpha_hess_r = t(-c(sapply(c(eXtB_r), rep, control$gs.rules)) * B_t_quad2_r) %*% as.matrix(gamma[[cov]] * quad.phi.event[[cov]])[c(sapply(censor[,1], rep, control$gs.rules)),]
        alpha_hess_l = -t(c((eXtB_l^2 * S_t_l/(1-S_t_l)^2)) * as.matrix(Bphi_t_l)) %*% as.matrix(Bphi_t_l)
        alpha_hess_i = t(-c(sapply(c(eXtB_i*S_t_i1 /(S_t_i1-S_t_i2)), rep, control$gs.rules)) * data.frame(B_t_quad2_i1)) %*% as.matrix(gamma[[cov]] * quad.phi.y_i1[[cov]]) -
          t(c((eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * as.matrix(Bphi_t_i1 - Bphi_t_i2)) %*% as.matrix(Bphi_t_i1 - Bphi_t_i2)
        
        alpha_hess_partial_new = alpha_hess_e + alpha_hess_r + alpha_hess_l + alpha_hess_i
        
        alpha_hess_ls = c(1/sigma2_Et[[cov]]) * t(phi.mu[[cov]]) %*% phi.mu[[cov]]
        
        alpha_hess = alpha_hess_partial_new  - alpha_hess_ls # - 2*alpha.lambda * alpha.G
        alpha_hess_neg = -alpha_hess
        
        alpha_update_save = solve(alpha_hess_neg)%*%alpha_score
        alpha[[cov]] = alpha_old[[cov]] + alpha_update_save
        
        
        #step size for alpha
        #likelihood
        mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
        nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
        z_t_quad = mu_t_quad + nu_t_quad
        exp_zTg_quad = exp(z_t_quad %*% gamma)
        h0_t_star_quad = h0_t_quad * exp_zTg_quad
        h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
        H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
        
        H0_t_e = H0_t[which(censor[,2])]
        H0_t_r = H0_t[which(censor[,1])]
        H0_t_l = H0_t[which(censor[,3])]
        #H0_t_i2 = H0_t[dat.baseline$interval == 1]
        
        mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
        nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
        z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
        exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
        h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
        h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
        H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
        
        mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
        nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
        z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
        exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
        h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
        h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
        H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
        
        #event
        eXtB_e = exp(fixed_e %*% beta)
        #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
        zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
        #zt_e = zt_e[which(censor[,2])]
        zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
        pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
        
        #right
        #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
        pl_r = - sum(eXtB_r * H0_t_r)
        
        #left
        S_t_l = exp(-eXtB_l * H0_t_l)
        diff_l = 1 - S_t_l
        diff_l[which(diff_l < 1e-5)] = 1e-5
        pl_l = sum(log(diff_l))
        
        #interval
        S_t_i1 = exp(-eXtB_i * H0_t_i1)
        S_t_i2 = exp(-eXtB_i * H0_t_i2)
        diff_i = S_t_i1 - S_t_i2
        diff_i[which(diff_i < 1e-5)] = 1e-5
        pl_i = sum(log(diff_i))
        
        #least squares
        
        mu_t_ia = sapply(seq_along(phi.mu), function(x) phi.mu[[x]] %*% alpha[[x]])
        nu_t_ia = sapply(seq_along(phi.mu.re), function(x) apply(phi.mu.re[[x]] * sapply(data.frame(a_re_cols[[x]]), rep, n.obs), 1, sum))
        z_t_ia = mu_t_ia + nu_t_ia
        pl_ls = sum((sapply(seq_along(cont), function(x) (1/(2*sigma2_Et[[x]]))*sum((cont[[x]] - z_t_ia[,x])^2))))
        
        log_lik_old = log_lik
        ll_old_surv = ll_surv
        #log likelihood
        ll_surv = pl_e + pl_r + pl_l + pl_i
        log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
          sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum)))) - 
          theta.lambda* t(theta)%*%theta.G%*%theta -
          alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
        
        
        #print(c("alpha", log_lik_old, log_lik))
        step_size_a = 1
        if(((log_lik- log_lik_old) < (-1e-6) & step_size_a == 1)){
          ii = 0
          omega1 = 1/control$kappa
          
          
          while((log_lik- log_lik_old) < (-1e-6)){
            alpha[[cov]] = alpha_old[[cov]] + omega1 * alpha_update_save
            
            #print(c("alpha ss", cov, log_lik- log_lik_old))
            
            #likelihood
            mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
            nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
            z_t_quad = mu_t_quad + nu_t_quad
            exp_zTg_quad = exp(z_t_quad %*% gamma)
            h0_t_star_quad = h0_t_quad * exp_zTg_quad
            h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
            H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
            
            H0_t_e = H0_t[which(censor[,2])]
            H0_t_r = H0_t[which(censor[,1])]
            H0_t_l = H0_t[which(censor[,3])]
            #H0_t_i2 = H0_t[dat.baseline$interval == 1]
            
            mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
            nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
            z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
            exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
            h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
            h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
            H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
            
            mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
            nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
            z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
            exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
            h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
            h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
            H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
            
            
            #event
            #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
            zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
            #zt_e = zt_e[which(censor[,2])]
            zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
            pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
            
            #right
            #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
            pl_r = - sum(eXtB_r * H0_t_r)
            
            #left
            S_t_l = exp(-eXtB_l * H0_t_l)
            diff_l = 1 - S_t_l
            diff_l[which(diff_l < 1e-5)] = 1e-5
            pl_l = sum(log(diff_l))
            
            #interval
            S_t_i1 = exp(-eXtB_i * H0_t_i1)
            S_t_i2 = exp(-eXtB_i * H0_t_i2)
            diff_i = S_t_i1 - S_t_i2
            diff_i[which(diff_i < 1e-5)] = 1e-5
            pl_i = sum(log(diff_i))
            
            #least squares
            
            mu_t_ia = sapply(seq_along(phi.mu), function(x) phi.mu[[x]] %*% alpha[[x]])
            
            #nu_t_ia = sapply(seq_along(phi.mu.re), function(x) apply(phi.mu.re[[x]] * sapply(data.frame(a_re_cols[[x]]), rep, n.obs), 1, sum))
            
            z_t_ia = mu_t_ia + nu_t_ia
            pl_ls = sum((sapply(seq_along(cont), function(x) (1/(2*sigma2_Et[[x]]))*sum((cont[[x]] - z_t_ia[,x])^2))))
            
            
            #log_lik_old = log_lik
            #log likelihood
            ll_surv = pl_e + pl_r + pl_l + pl_i
            log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
              sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum)))) - 
              theta.lambda* t(theta)%*%theta.G%*%theta -
            alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
            
            #update value of omega1
            if(omega1>=1e-2){
              omega1 = omega1/control$kappa
            }else if(omega1<1e-2 & omega1>=1e-5){
              omega1 = omega1*(5e-2)
            }else if(omega1<1e-5){
              omega1 = omega1*(1e-5)
            }
            ii = ii+1
            if(ii>20){break}
            
            
          }
        }
        
      }
      
      
      
      #gamma
      gamma_old = gamma
      for(cov in 1:q){
        exp_zTg_quad = exp(z_t_quad %*% gamma)
        A_t_quad = c(h0_t_quad * exp_zTg_quad) * z_t_quad[,cov]
        df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* A_t_quad)
        A_t = c(quad.lambda) * as.matrix(aggregate( df[,2], list(df[,1]), FUN = sum )[,-1])
        #colnames(A_t) = NULL
        
        A_t_e = A_t[which(censor[,2]),]
        A_t_r = A_t[which(censor[,1]),]
        A_t_l = A_t[which(censor[,3]),]
        A_t_i2 = A_t[which(censor[,4]),]
        
        exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
        A_t_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * z_t_quad_i1[,cov]
        df = data.frame(id_quad = c(sapply(unique(id[which(censor_long[,4])]), rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* A_t_quad_i1)
        A_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2], list(df[,1]), FUN = sum )[,-1])
        #colnames(A_t_i1) = NULL
        
        
        
        if(sum(censor[,2]) > 0){
          gamma_score_e = apply(zt_e[,cov] - c(eXtB_e) * (A_t_e), 2,sum)
        }else{
          gamma_score_e = 0
        }
        gamma_score_r = t(A_t_r) %*% (-eXtB_r)
        gamma_score_l = t(A_t_l) %*% c(S_t_l * eXtB_l / (1 - S_t_l))
        gamma_score_i = t(A_t_i1) %*% (- c(eXtB_i * (S_t_i1) / (S_t_i1 - S_t_i2 ))) + 
          t(A_t_i2) %*% (c(eXtB_i * (S_t_i2) / (S_t_i1 - S_t_i2 )))
        
        gamma_score = gamma_score_e + gamma_score_r + gamma_score_l + gamma_score_i
        
        #second derivative of gamma
        A2_t_quad_long = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, control$gs.rules)) * A_t_quad
        #colnames(A2_t_quad_long) = NULL
        
        A_t_quad2_e = A2_t_quad_long[c(sapply(censor[,2], rep, control$gs.rules))]
        A_t_quad2_r = A2_t_quad_long[c(sapply(censor[,1], rep, control$gs.rules))]
        A_t_quad2_l = A2_t_quad_long[c(sapply(censor[,3], rep, control$gs.rules))]
        A_t_quad2_i2 = A2_t_quad_long[c(sapply(censor[,4], rep, control$gs.rules))]
        
        A_t_quad2_i1 = c(rep(quad.w, sum(censor[,4]))) * 
          c(sapply(quad.lambda_i1, rep, control$gs.rules)) * A_t_quad_i1
        #colnames(A_t_quad2_i1) = NULL
        #colnames(z_t_quad) = NULL
        #colnames(z_t_quad_i1) = NULL
        
        
        if(sum(censor[,2]) > 0){
          gamma_hess_e = t(-c(sapply(c(eXtB_e), rep, control$gs.rules)) * as.matrix(z_t_quad[c(sapply(censor[,2], rep, control$gs.rules)),cov])) %*% as.matrix(A_t_quad2_e)
        }else{
          gamma_hess_e = 0
        }
        
        gamma_hess_r = t(-c(sapply(c(eXtB_r), rep, control$gs.rules)) * as.matrix(z_t_quad[c(sapply(censor[,1], rep, control$gs.rules)),cov])) %*% as.matrix(A_t_quad2_r)
        #gamma_hess_l = t(S_t_l * eXtB_l/(1-S_t_l)) %*% A_t_quad2_l - 
        #  t(c( S_t_l * eXtB_l^2/(1-S_t_l)^2) * A_t_l) %*% A_t_l 
        #gamma_hess_i = t(-eXtB_i *S_t_i1 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i1 +
        #  t(eXtB_i *S_t_i2 / (S_t_i1 - S_t_i2)) %*% A_t_quad2_i2 -
        #  t(c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2 ) * (A_t_i1 - A_t_i2)) %*% (A_t_i1 - A_t_i2)
        
        gamma_hess_l = - t(c( S_t_l * (eXtB_l^2)/(1-S_t_l)^2) * A_t_l) %*% A_t_l 
        gamma_hess_i = - t(c(sapply(c(eXtB_i *S_t_i1 / (S_t_i1 - S_t_i2)), rep, control$gs.rules)) * as.matrix(z_t_quad_i1[,cov])) %*% as.matrix(A_t_quad2_i1) -
          t(c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2 ) * (A_t_i1 - A_t_i2)) %*% (A_t_i1- A_t_i2)
        
        
        gamma_hess = gamma_hess_e + gamma_hess_r + gamma_hess_l + gamma_hess_i
        gamma_hess_neg = -gamma_hess
        
        
        gamma[cov] = gamma_old[cov] + solve(gamma_hess_neg)%*%gamma_score
        #print(c("update gamma"))
        
        #step size for gamma
        #likelihood
        exp_zTg_quad = exp(z_t_quad %*% gamma)
        h0_t_star_quad = h0_t_quad * exp_zTg_quad
        h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
        H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
        
        H0_t_e = H0_t[which(censor[,2])]
        H0_t_r = H0_t[which(censor[,1])]
        H0_t_l = H0_t[which(censor[,3])]
        #H0_t_i2 = H0_t[dat.baseline$interval == 1]
        
        exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
        h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
        h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
        H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
        
        
        exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
        h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
        h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
        H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
        
        
        #event
        #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
        #zt_e = zt_e[which(censor[,2])]
        zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[x])))
        pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
        
        #right
        #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
        pl_r = - sum(eXtB_r * H0_t_r)
        
        #left
        S_t_l = exp(-eXtB_l * H0_t_l)
        diff_l = 1 - S_t_l
        diff_l[which(diff_l < 1e-5)] = 1e-5
        pl_l = sum(log(diff_l))
        
        #interval
        S_t_i1 = exp(-eXtB_i * H0_t_i1)
        S_t_i2 = exp(-eXtB_i * H0_t_i2)
        diff_i = S_t_i1 - S_t_i2
        diff_i[which(diff_i < 1e-5)] = 1e-5
        pl_i = sum(log(diff_i))
        
        log_lik_old = log_lik
        ll_old_surv = ll_surv
        #log likelihood
        ll_surv = pl_e + pl_r + pl_l + pl_i
        log_lik = ll_surv - pl_ls - 
          sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum)))) - 
          theta.lambda* t(theta)%*%theta.G%*%theta -
          alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
        
        if((((ll_surv - ll_old_surv) < (-1e-6)) & control$step_size == 1)){
          ii = 0
          omega1 = 1/control$kappa
          
          
          while((ll_surv - ll_old_surv) < (-1e-6)){
            gamma[cov] = gamma_old[cov] +  omega1 * solve(gamma_hess_neg)%*%gamma_score
            #print(c("gamma ss", ll_surv - ll_old_surv))
            #likelihood
            exp_zTg_quad = exp(z_t_quad %*% gamma)
            h0_t_star_quad = h0_t_quad * exp_zTg_quad
            h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
            H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
            
            H0_t_e = H0_t[which(censor[,2])]
            H0_t_r = H0_t[which(censor[,1])]
            H0_t_l = H0_t[which(censor[,3])]
            #H0_t_i2 = H0_t[dat.baseline$interval == 1]
            
            exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
            h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
            h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
            H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
            
            exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
            h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
            h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
            H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
            
            
            #event
            #zt_e = zt_e[which(censor[,2])]
            zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
            pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
            
            #right
            #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
            pl_r = - sum(eXtB_r * H0_t_r)
            
            #left
            S_t_l = exp(-eXtB_l * H0_t_l)
            diff_l = 1 - S_t_l
            diff_l[which(diff_l < 1e-5)] = 1e-5
            pl_l = sum(log(diff_l))
            
            #interval
            S_t_i1 = exp(-eXtB_i * H0_t_i1)
            S_t_i2 = exp(-eXtB_i * H0_t_i2)
            diff_i = S_t_i1 - S_t_i2
            diff_i[which(diff_i < 1e-5)] = 1e-5
            pl_i = sum(log(diff_i))
            
            #log_lik_old = log_lik
            #log likelihood
            ll_surv = pl_e + pl_r + pl_l + pl_i
            #alpha.lambda* t(alpha)%*%alpha.G%*%alpha
            
            #print(c(ii, ll_surv, ll_old_surv, as.numeric(ll_surv < ll_old_surv)))
            
            #update value of omega1
            if(omega1>=1e-2){
              omega1 = omega1/control$kappa
            }else if(omega1<1e-2 & omega1>=1e-5){
              omega1 = omega1*(5e-2)
            }else if(omega1<1e-5){
              omega1 = omega1*(1e-5)
            }
            ii = ii+1
            #print(c("gamma_step"))
            if(ii>50){break}
            
            
          }
        }
        
      }
      
      
      
      #rownames(gamma) = NULL
      #theta
      psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
      df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* psi_t_star_quad)
      Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
      
      Psi_t_e = Psi_t[which(censor[,2]),]
      Psi_t_r = Psi_t[which(censor[,1]),]
      Psi_t_l = Psi_t[which(censor[,3]),]
      Psi_t_i2 = Psi_t[which(censor[,4]),]
      
      psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
      df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* psi_t_star_quad_i1)
      Psi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
      
      #psi_t_star_quad_i2 = c(exp_zTg_quad_i2) * quad.psi.y_i2
      #df = data.frame(id_quad = c(sapply(dat.baseline$id[dat.baseline$interval == 1], rep, control$gs.rules)), c(rep(quad.w, sum(dat.baseline$interval == 1)))* psi_t_star_quad_i2)
      #Psi_t_i2 = Psi_t_i1 + c(quad.lambda_i2) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
      
      
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
      #print(c("update theta"))
      
      #step size for theta
      #likelihood
      h0_t_quad = quad.psi.event %*% theta
      mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
      nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
      z_t_quad = mu_t_quad + nu_t_quad
      exp_zTg_quad = exp(z_t_quad %*% gamma)
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[which(censor[,2])]
      H0_t_r = H0_t[which(censor[,1])]
      H0_t_l = H0_t[which(censor[,3])]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
      nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
      z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
      exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
      nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
      z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
      exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
      zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
      #zt_e = zt_e[which(censor[,2])]
      zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      diff_l = 1 - S_t_l
      diff_l[which(diff_l < 1e-5)] = 1e-5
      pl_l = sum(log(diff_l))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      diff_i = S_t_i1 - S_t_i2
      diff_i[which(diff_i < 1e-5)] = 1e-5
      pl_i = sum(log(diff_i))
      
      
      log_lik_old = log_lik
      ll_old_surv = ll_surv
      
      #log likelihood
      ll_surv = pl_e + pl_r + pl_l + pl_i
      log_lik = ll_surv - pl_ls - 
        sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum)))) - 
        theta.lambda* t(theta)%*%theta.G%*%theta -
        alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
      
      
      if(((ll_surv - ll_old_surv) < (-1e-6) & control$step_size == 1)){
        i = 0
        omega1 = 1/control$kappa
        
        
        while((ll_surv - ll_old_surv) < (-1e-6)){
          theta = theta_old +  omega1 * D_matrix*(theta_score)
          #print(c("theta ss"))
          #likelihood
          eXtB = exp(fixed %*% beta)
          
          h0_t_quad = quad.psi.event %*% theta
          mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
          nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
          z_t_quad = mu_t_quad + nu_t_quad
          exp_zTg_quad = exp(z_t_quad %*% gamma)
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[which(censor[,2])]
          H0_t_r = H0_t[which(censor[,1])]
          H0_t_l = H0_t[which(censor[,3])]
          #H0_t_i2 = H0_t[dat.baseline$interval == 1]
          
          mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
          nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
          z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
          exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
          h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
          h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
          H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
          
          h0_t_quad_i2 = quad.psi.y_i2 %*% theta
          mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
          nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
          z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
          exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
          h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
          h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
          H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
          zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
          #zt_e = zt_e[which(censor[,2])]
          zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #left
          eXtB_l = exp(fixed_l %*% beta)
          S_t_l = exp(-eXtB_l * H0_t_l)
          diff_l = 1 - S_t_l
          diff_l[which(diff_l < 1e-5)] = 1e-5
          pl_l = sum(log(diff_l))
          
          #interval
          eXtB_i = exp(fixed_i %*% beta)
          S_t_i1 = exp(-eXtB_i * H0_t_i1)
          S_t_i2 = exp(-eXtB_i * H0_t_i2)
          diff_i = S_t_i1 - S_t_i2
          diff_i[which(diff_i < 1e-5)] = 1e-5
          pl_i = sum(log(diff_i))
          
          
          #log_lik_old = log_lik
          #log likelihood
          ll_surv = pl_e + pl_r + pl_l + pl_i
          #alpha.lambda* t(alpha)%*%alpha.G%*%alpha
          
          #print(c(log_lik_old, log_lik, omega1))
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/control$kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>20){break}
          
          
        }
      }
      
      
      
      #beta
      
      if(sum(censor[,2]) > 0){
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
      #print(c("update beta"))
      #step size for beta
      #likelihood
      eXtB = exp(fixed %*% beta)
      
      h0_t_quad = quad.psi.event %*% theta
      mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
      nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
      z_t_quad = mu_t_quad + nu_t_quad
      exp_zTg_quad = exp(z_t_quad %*% gamma)
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      
      H0_t_e = H0_t[which(censor[,2])]
      H0_t_r = H0_t[which(censor[,1])]
      H0_t_l = H0_t[which(censor[,3])]
      #H0_t_i2 = H0_t[dat.baseline$interval == 1]
      
      h0_t_quad_i1 = quad.psi.y_i1 %*% theta
      mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
      nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
      z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
      exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
      h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
      h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
      H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
      
      h0_t_quad_i2 = quad.psi.y_i2 %*% theta
      mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
      nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
      z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
      exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
      h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
      h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
      H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
      
      
      #event
      h0_t_e = psi_e %*% theta
      eXtB_e = exp(fixed_e %*% beta)
      #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
      zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
      #zt_e = zt_e[which(censor[,2])]
      zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
      pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
      
      #right
      eXtB_r = exp(fixed_r %*% beta)
      #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #left
      eXtB_l = exp(fixed_l %*% beta)
      S_t_l = exp(-eXtB_l * H0_t_l)
      diff_l = 1 - S_t_l
      diff_l[which(diff_l < 1e-5)] = 1e-5
      pl_l = sum(log(diff_l))
      
      #interval
      eXtB_i = exp(fixed_i %*% beta)
      S_t_i1 = exp(-eXtB_i * H0_t_i1)
      S_t_i2 = exp(-eXtB_i * H0_t_i2)
      diff_i = S_t_i1 - S_t_i2
      diff_i[which(diff_i < 1e-5)] = 1e-5
      pl_i = sum(log(diff_i))
      
      log_lik_old = log_lik
      ll_old_surv = ll_surv
      #log likelihood
      ll_surv = pl_e + pl_r + pl_l + pl_i
      log_lik = ll_surv - pl_ls - 
        sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum))))- 
        theta.lambda* t(theta)%*%theta.G%*%theta -
        alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
      
      
      #print(c("beta", log_lik_old, log_lik))
 
      if(((ll_surv < ll_old_surv) < (-1e-6) & control$step_size == 1)){
        i = 0
        omega1 = 1/control$kappa
        
        
        while((ll_surv < ll_old_surv) < (-1e-6)){
          beta = beta_old +  omega1 * solve(beta_neg_hess)%*%beta_score
          #print(c("beta ss"))
          #likelihood
          #likelihood
          eXtB = exp(fixed %*% beta)
          
          h0_t_quad = quad.psi.event %*% theta
          mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
          nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
          z_t_quad = mu_t_quad + nu_t_quad
          exp_zTg_quad = exp(z_t_quad %*% gamma)
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
          H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
          
          H0_t_e = H0_t[which(censor[,2])]
          H0_t_r = H0_t[which(censor[,1])]
          H0_t_l = H0_t[which(censor[,3])]
          #H0_t_i2 = H0_t[dat.baseline$interval == 1]
          
          h0_t_quad_i1 = quad.psi.y_i1 %*% theta
          mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
          nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
          z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
          exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
          h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
          h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
          H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
          
          h0_t_quad_i2 = quad.psi.y_i2 %*% theta
          mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
          nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
          z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
          exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
          h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
          h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
          H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
          
          
          #event
          h0_t_e = psi_e %*% theta
          eXtB_e = exp(fixed_e %*% beta)
          #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
          zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
          #zt_e = zt_e[which(censor[,2])]
          zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
          pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
          
          #right
          eXtB_r = exp(fixed_r %*% beta)
          #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #left
          eXtB_l = exp(fixed_l %*% beta)
          S_t_l = exp(-eXtB_l * H0_t_l)
          diff_l = 1 - S_t_l
          diff_l[which(diff_l < 1e-5)] = 1e-5
          pl_l = sum(log(diff_l))
          
          #interval
          eXtB_i = exp(fixed_i %*% beta)
          S_t_i1 = exp(-eXtB_i * H0_t_i1)
          S_t_i2 = exp(-eXtB_i * H0_t_i2)
          diff_i = S_t_i1 - S_t_i2
          diff_i[which(diff_i < 1e-5)] = 1e-5
          pl_i = sum(log(diff_i))
          
          #least squares
        
          #log_lik_old = log_lik
          #log likelihood
          ll_surv = pl_e + pl_r + pl_l + pl_i
          #alpha.lambda* t(alpha)%*%alpha.G%*%alpha
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/control$kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>20){break}
          
          
        }
        #print(c(omega1))
        
      }
      
      log_lik = ll_surv - pl_ls - 
        sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum))))- 
        theta.lambda* t(theta)%*%theta.G%*%theta -
        alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
      
      a_re_cols_old = a_re_cols
      
      for(cov in 1:q){
        
        B_t_re_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event.re[[cov]]
        df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* B_t_re_quad)
        B_t_re = c(quad.lambda) * as.matrix(aggregate( df[,2:(w[[cov]]+1)], list(df[,1]), FUN = sum )[,-1])
        
        Bphi_t_re_e = B_t_re[which(censor[,2]),]
        Bphi_t_re_r = B_t_re[which(censor[,1]),]
        Bphi_t_re_l = B_t_re[which(censor[,3]),]
        Bphi_t_re_i2 = B_t_re[which(censor[,4]),]
        
        B_t_re_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * quad.phi.y_i1.re[[cov]]
        df = data.frame(id_quad = c(sapply(unique(id)[(censor[,4])], rep, control$gs.rules)), 
                        c(rep(quad.w, sum((censor[,4]))))* B_t_re_quad_i1)
        B_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(w[[cov]]+1)], list(df[,1]), FUN = sum )[,-1])
        
        a1i_score_e = as.numeric(gamma[cov]) * phi_e.re[[cov]] - as.numeric(gamma[cov]) * as.numeric(eXtB_e) * (Bphi_t_re_e)
        a1i_score_r = - as.numeric(gamma[cov]) * as.numeric(eXtB_r) * (Bphi_t_re_r)
        a1i_score_l = as.numeric(gamma[cov]) * as.numeric(eXtB_l * S_t_l/(1-S_t_l)) * (Bphi_t_re_l)
        a1i_score_i = - as.numeric(gamma[cov]) * as.numeric(eXtB_i * S_t_i1/(S_t_i1-S_t_i2)) * (B_t_re_i1) + 
          as.numeric(gamma[cov]) * as.numeric(eXtB_i * S_t_i2/(S_t_i1-S_t_i2)) * (Bphi_t_re_i2)
        
        id_censor_order = c(unique(id)[which(censor[,2])], 
                            unique(id)[which(censor[,1])], 
                            unique(id)[which(censor[,3])], 
                            unique(id)[which(censor[,4])])
        a1i_score_s = rbind(as.matrix(a1i_score_e), as.matrix(a1i_score_r), as.matrix(a1i_score_l), as.matrix(a1i_score_i))
        a1i_score_s = data.frame(id_censor_order, a1i_score_s)[order(id_censor_order),-1]
        
        a1i_score_ls = data.frame(id, (cont[[cov]]-z_t_ia[,cov])*phi.mu.re[[cov]])
        a1i_score_ls = as.matrix(aggregate( a1i_score_ls[,2:(w[[cov]]+1)], list(a1i_score_ls[,1]), FUN = sum )[,-1])
        
        a1i_score = a1i_score_s + (1/sigma2_Et[[cov]])*a1i_score_ls - a_re_cols[[cov]]/matrix(rep(sigma2_re[[cov]], n), ncol = w[[cov]], byrow = T)
        
        a_cols_update = matrix(0, nrow = n, ncol = w[[cov]])
        
        S_t = exp(-H0_t*eXtB)
        
        S_t_i1_long = rep(0, length(S_t))
        S_t_i1_long[which(censor[,4])] = S_t_i1
        
        B_t_re_quad_i1_long = matrix(0, nrow = nrow(as.matrix(B_t_re_quad)), ncol = ncol(as.matrix(B_t_re_quad)))
        B_t_re_quad_i1_long[i1_ind_quad,] = B_t_re_quad_i1
        
        B_t_re_i1_long = matrix(0, nrow = nrow(B_t_re), ncol = ncol(B_t_re))
        B_t_re_i1_long[which(censor[,4]),] = B_t_re_i1
        
        quad.phi.y_i1_long.re = matrix(0, nrow = nrow(quad.phi.event.re[[cov]]), ncol = ncol(quad.phi.event.re[[cov]]))
        quad.phi.y_i1_long.re[i1_ind_quad,] = quad.phi.y_i1.re[[cov]]
        
        phi.mu.re.cov = as.matrix(phi.mu.re[[cov]])
        #quad.phi.event = as.matrix(quad.phi.event)
        B_t_re_quad = as.matrix(B_t_re_quad)
        quad.lambda_i1_long = rep(0, n)
        quad.lambda_i1_long[which(censor[,4])] = quad.lambda_i1
        
      
        #second derivative
        B_t_quad_re = c(h0_t_quad * exp_zTg_quad) * data.frame(quad.phi.event.re[[cov]])
        B_t_quad2_re = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, control$gs.rules)) * B_t_quad_re 
    
        B_t_quad_re_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * data.frame(quad.phi.y_i1.re[[cov]])
        B_t_quad_re_i1_2 = c(rep(quad.w, sum(censor[,4]))) * 
          c(sapply(quad.lambda_i1, rep, control$gs.rules)) * B_t_quad_re_i1
        
        temp = (quad.phi.event.re[[cov]] * gamma[cov]^2)  * B_t_quad2_re
        df = data.frame(id_quad = rep(c(sapply(unique(id), rep, control$gs.rules)), w[[cov]]),
                     c(rep(quad.w, n)) * temp)
        D_t_re_re = c(quad.lambda) * as.matrix(aggregate( df[,2:(w[[cov]]+1)], list(df[,1]), FUN = sum )[,-c(1)])
    
        temp = (quad.phi.y_i1.re[[cov]] * gamma[cov]^2) * B_t_quad_re_i1_2
        df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)),
                    c(rep(quad.w, sum(censor[,4]))) * temp)
        D_t_re_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(w[[cov]]+1)], list(df[,1]), FUN = sum )[,-c(1)])
    
        D_t_rere_i1_long = matrix(0, nrow = nrow(D_t_re_re), ncol = ncol(D_t_re_re))
        D_t_rere_i1_long[rep(which(censor[,4]), w[[cov]]),] = D_t_re_re_i1
        a1i_hess_ls = data.frame(id, phi.mu.re[[cov]]*phi.mu.re[[cov]])
        a1i_hess_ls = as.matrix(aggregate( a1i_hess_ls[,2:(w[[cov]]+1)], list(a1i_hess_ls[,1]), FUN = sum )[,-1])
    
        #a1a1_hess = - rep((c(censor[,2] + censor[,1]) * eXtB), w[[cov]]) * D_t_re_re + 
        #  rep((c(censor[,3]) * eXtB * (S_t/(1-S_t))), w[[cov]]) * D_t_re_re - 
        #  (rep((c(censor[,3]) * (eXtB^2) * (S_t/(1-S_t)^2)), w[[cov]]) * B_t_re) * B_t_re -
        #  rep((c(censor[,4]) * eXtB * (S_t_i1_long/(S_t_i1_long-S_t)) ), w[[cov]]) * D_t_rere_i1_long + 
        #  rep((c(censor[,4]) * eXtB * (S_t/(S_t_i1_long-S_t))), w[[cov]]) * D_t_re_re - 
        #  (rep((c(censor[,4]) * (eXtB^2) * (S_t_i1_long*S_t/(S_t_i1_long-S_t)^2)), w[[cov]]) * (B_t_re_i1_long - B_t_re)) * (B_t_re_i1_long - B_t_re) - 
        #  c(1/sigma2_Et[[cov]]) * a1i_hess_ls - 
        #  matrix(rep(1/sigma2_re[[cov]], n), ncol = w[[cov]], byrow = TRUE)
    
        a1a1_hess = - rep((c(censor[,2] + censor[,1]) * eXtB), w[[cov]]) * D_t_re_re - 
          (rep((c(censor[,3]) * (eXtB^2) * gamma[cov]^2 * (S_t/(1-S_t)^2)), w[[cov]]) * B_t_re) * B_t_re -
          rep((c(censor[,4]) * eXtB * (S_t_i1_long/(S_t_i1_long-S_t)) ), w[[cov]]) * D_t_rere_i1_long - 
          (rep((c(censor[,4]) * (eXtB^2) * gamma[cov]^2 * (S_t_i1_long*S_t/(S_t_i1_long-S_t)^2)), w[[cov]]) * (B_t_re_i1_long - B_t_re)) * (B_t_re_i1_long - B_t_re) - 
          c(1/sigma2_Et[[cov]]) * a1i_hess_ls - 
          matrix(rep(1/sigma2_re[[cov]], n), ncol = w[[cov]], byrow = TRUE)
        
        
        a_re_cols[[cov]] = a_re_cols_old[[cov]] + as.matrix(a1i_score/(-a1a1_hess))
        #a_re_cols = a_re_cols_old + as.matrix(a1i_score)
        #print(c("update kappa"))
        
        a_re_long_old = a_re_long
        #a_re_cols_pad_old = a_re_cols_pad
        a_re_cols_quad_old = a_re_cols_quad
        
        a_re_long[[cov]] = unlist(c(a_re_cols[[cov]]))
        #a_re_cols_pad[,control$re_ind] = as.matrix(a_re_cols)
        a_re_cols_quad[[cov]] = matrix(sapply(a_re_cols[[cov]], rep, control$gs.rules), ncol = w[[cov]])
        
        
        #step size for random effects
        #likelihood
        h0_t_quad = quad.psi.event %*% theta
        mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
        nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
        z_t_quad = mu_t_quad + nu_t_quad
        exp_zTg_quad = exp(z_t_quad %*% gamma)
        h0_t_star_quad = h0_t_quad * exp_zTg_quad
        h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
        H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
        
        H0_t_e = H0_t[which(censor[,2])]
        H0_t_r = H0_t[which(censor[,1])]
        H0_t_l = H0_t[which(censor[,3])]
        #H0_t_i2 = H0_t[dat.baseline$interval == 1]
        
        h0_t_quad_i1 = quad.psi.y_i1 %*% theta
        mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
        nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
        z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
        exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
        h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
        h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
        H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
        
        h0_t_quad_i2 = quad.psi.y_i2 %*% theta
        mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
        nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
        z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
        exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
        h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
        h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
        H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
        
        #event
        h0_t_e = psi_e %*% theta
        eXtB_e = exp(fixed_e %*% beta)
        #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
        zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
        #zt_e = zt_e[which(censor[,2])]
        zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
        pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
        
        #right
        eXtB_r = exp(fixed_r %*% beta)
        #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
        pl_r = - sum(eXtB_r * H0_t_r)
        #print(c(pl_r))
        #left
        eXtB_l = exp(fixed_l %*% beta)
        S_t_l = exp(-eXtB_l * H0_t_l)
        diff_l = 1 - S_t_l
        diff_l[which(diff_l < 1e-5)] = 1e-5
        pl_l = sum(log(diff_l))
        
        #interval
        eXtB_i = exp(fixed_i %*% beta)
        S_t_i1 = exp(-eXtB_i * H0_t_i1)
        S_t_i2 = exp(-eXtB_i * H0_t_i2)
        diff_i = S_t_i1 - S_t_i2
        diff_i[which(diff_i < 1e-5)] = 1e-5
        pl_i = sum(log(diff_i))
        
        #least squares
        
        #mu_t_ia = sapply(seq_along(phi.mu), function(x) phi.mu[[x]] %*% alpha[[x]])
        
        nu_t_ia = sapply(seq_along(phi.mu.re), function(x) apply(phi.mu.re[[x]] * sapply(data.frame(a_re_cols[[x]]), rep, n.obs), 1, sum))
        
        z_t_ia = mu_t_ia + nu_t_ia
        
        pl_ls = sum((sapply(seq_along(cont), function(x) (1/(2*sigma2_Et[[x]]))*sum((cont[[x]] - z_t_ia[,x])^2))))
        
        log_lik_old = log_lik
        ll_old_surv = ll_surv
        #log likelihood
        ll_surv = pl_e + pl_r + pl_l + pl_i
        log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
          sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum)))) - 
          theta.lambda* t(theta)%*%theta.G%*%theta -
          alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
        
        #print(c("re", log_lik_old, log_lik))
        re_ss = 1
        
        if(((log_lik- log_lik_old) < (-1e-6) & re_ss == 1)){
          ii = 0
          omega1 = 1/control$kappa
          
          
          while((log_lik- log_lik_old) < (-1e-6)){
            
            a_re_cols[[cov]] = a_re_cols_old[[cov]] + omega1 * as.matrix(a1i_score/(-a1a1_hess))
            #print(c("kappa ss", omega1, c(log_lik - log_lik_old)))
            
            a_re_long[[cov]] = unlist(c(a_re_cols[[cov]]))
            a_re_cols_quad[[cov]] = matrix(sapply(a_re_cols[[cov]], rep, control$gs.rules), ncol = w[[cov]])
           
            #step size for random effects
            #likelihood
            #mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
            nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
            z_t_quad = mu_t_quad + nu_t_quad
            exp_zTg_quad = exp(z_t_quad %*% gamma)
            h0_t_star_quad = h0_t_quad * exp_zTg_quad
            h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
            H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
            
            
            H0_t_e = H0_t[which(censor[,2])]
            H0_t_r = H0_t[which(censor[,1])]
            H0_t_l = H0_t[which(censor[,3])]
            #H0_t_i2 = H0_t[dat.baseline$interval == 1]
            
            #mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
            nu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
            z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
            exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
            h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
            h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
            H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
            
            #mu_t_quad_i2 = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]][i1_ind_quad ==1,] %*% alpha[[x]])
            nu_t_quad_i2 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]][i1_ind_quad ==1,] * a_re_cols_quad[[x]][i1_ind_quad ==1,], 1, sum))
            z_t_quad_i2 = mu_t_quad_i2 + nu_t_quad_i2
            exp_zTg_quad_i2 = exp(z_t_quad_i2 %*% gamma)
            h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
            h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
            H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
            
            
            #event
            #egWtA_e = exp(c(gamma) * (W_short[dat.baseline$event == 1] %*% alpha_v))
            zt_e = data.frame(lapply(seq_along(phi), function(x) (phi[[x]] %*% alpha[[x]])[which(censor[,2]),] + apply((phi_e.re[[x]] * a_re_cols[[x]][which(censor[,2]),]), 1, sum)))
            #zt_e = zt_e[which(censor[,2])]
            zt_g_e = unlist(lapply(seq_along(zt_e), function(x) zt_e[[x]] * c(gamma[[x]])))
            pl_e = sum(log(h0_t_e) + fixed_e %*% beta + sum(zt_g_e) - eXtB_e*H0_t_e)
            
            #right
            #egWtA_r = exp(c(gamma) * (W_short[dat.baseline$right == 1] %*% alpha_v))
            pl_r = - sum(eXtB_r * H0_t_r)
            #print(c(pl_r)) ## SOMETHING IS WRONG WITH THE STEP SIZE HERE?? not getting back to previous value
            #left
            S_t_l = exp(-eXtB_l * H0_t_l)
            diff_l = 1 - S_t_l
            diff_l[which(diff_l < 1e-5)] = 1e-5
            pl_l = sum(log(diff_l))
            
            #interval
            S_t_i1 = exp(-eXtB_i * H0_t_i1)
            S_t_i2 = exp(-eXtB_i * H0_t_i2)
            diff_i = S_t_i1 - S_t_i2
            diff_i[which(diff_i < 1e-5)] = 1e-5
            pl_i = sum(log(diff_i))
            
            #least squares
            
            #mu_t_ia = sapply(seq_along(phi.mu), function(x) phi.mu[[x]] %*% alpha[[x]])
            
            nu_t_ia = sapply(seq_along(phi.mu.re), function(x) apply(phi.mu.re[[x]] * sapply(data.frame(a_re_cols[[x]]), rep, n.obs), 1, sum))
            
            z_t_ia = mu_t_ia + nu_t_ia
            
            pl_ls = sum((sapply(seq_along(cont), function(x) (1/(2*sigma2_Et[[x]]))*sum((cont[[x]] - z_t_ia[,x])^2))))
            
            #log_lik_old = log_lik
            #log likelihood
            ll_surv = pl_e + pl_r + pl_l + pl_i
            log_lik = pl_e + pl_r + pl_l + pl_i - pl_ls - 
              sum((sapply(seq_along(a_re_cols), function(x) (1/(2*sigma2_re[[x]])) * apply(a_re_cols[[x]] * a_re_cols[[x]], 2, sum)))) - 
              theta.lambda* t(theta)%*%theta.G%*%theta -
              alpha.lambda* t(alpha[[1]])%*%alpha.G%*%alpha[[1]]
            
            #update value of omega1
            if(omega1>=1e-2){
              omega1 = omega1/control$kappa
            }else if(omega1<1e-2 & omega1>=1e-5){
              omega1 = omega1*(5e-2)
            }else if(omega1<1e-5){
              omega1 = omega1*(1e-5)
            }
            ii = ii+1
            if(ii>50){break}
            
            
          }
          #print(c(omega1))
          
        }
        
        
        
      }
      
      
      #hist(a_re_cols[,1])
      #hist(a_re_cols[,2])
      
      #for(u in 1:m){
      #  if((theta[u]< (5e-4) & theta_score[u] > 1e-1)){
      #    #pos_temp[u] = 0
      #    theta[u] = 0.1
      #  }
      #}
      
      #pos = rep(1, p+q+m+r+w*n)
      #for(u in 1:m){
      #  if((theta[u]< (1e-1) & theta_score[u] < (-0.001))){
      #    pos[p+q+u] = 0
      #  }
      #}
      
      print(c(iter, beta, gamma, theta, unlist(alpha), log_lik))
      
      #print(c(iter, theta, theta_score_pos/theta_score_neg, log_lik))
      
      #all(abs(1 - (theta_score_pos/theta_score_neg))[pos[(p+q+1):(p+q+m)]==1] < 1e-2)
      
      ll_save = c(ll_save, log_lik)
      conv_record = 0
      if(all(abs(c(beta - beta_old, gamma - gamma_old, theta - theta_old, unlist(alpha) - unlist(alpha_old))) < 1e-6) & 
         all(abs((data.frame(a_re_cols) - data.frame(a_re_cols_old))) < 1e-5)  ){
        conv_record = 1
        break
      }
      
      if(it == 1 & iter == 100){
        break
      }
      
      
      
    } #end inner loop
    
    
    print(c(iter, beta, gamma, log_lik))
    
    r_all = sum(unlist(r))
    w_all = sum(unlist(w))
    
    Hessian = matrix(0, nrow = (p + q + m + r_all + w_all*n), ncol = (p + q + m + r_all + w_all*n))
    #Hessian = matrix(0, nrow = (p + q + m + r + n), ncol = (p + q + m + r + n))
    
    
    #least squares
    mu_t_ia = sapply(seq_along(phi.mu), function(x) phi.mu[[x]] %*% alpha[[x]])
    nu_t_ia = sapply(seq_along(phi.mu.re), function(x) apply(phi.mu.re[[x]] * sapply(data.frame(a_re_cols[[x]]), rep, n.obs), 1, sum))
    z_t_ia = mu_t_ia + nu_t_ia
    
    #phi.mu.summed  = data.frame(id, phi.mu.re)
    #phi.mu.summed = as.matrix(aggregate( phi.mu.summed[,2:(w+1)], list(phi.mu.summed[,1]), FUN = sum )[,-1])
    #phi.mu.summed_alpha  = data.frame(id, phi.mu)
    #phi.mu.summed_alpha = as.matrix(aggregate( phi.mu.summed_alpha[,2:(r+1)], list(phi.mu.summed_alpha[,1]), FUN = sum )[,-1])
    
    
    eXtB_e = exp(fixed_e %*% beta)
    eXtB_r = exp(fixed_r %*% beta)
    eXtB_l = exp(fixed_l %*% beta)
    eXtB_i = exp(fixed_i %*% beta)
    
    h0_t_e = psi_e %*% theta
    
    h0_t_quad = quad.psi.event %*% theta
    mu_t_quad = sapply(seq_along(quad.phi.event), function(x) quad.phi.event[[x]] %*% alpha[[x]])
    nu_t_quad = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.event.re[[x]] * a_re_cols_quad[[x]], 1, sum))
    z_t_quad = mu_t_quad + nu_t_quad
    exp_zTg_quad = exp(z_t_quad %*% gamma)
    h0_t_star_quad = h0_t_quad * exp_zTg_quad
    h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = control$gs.rules, byrow = FALSE)
    H0_t = quad.lambda * apply(matrix(rep(quad.w, n), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad, 2, sum)
    
    H0_t_e = H0_t[which(censor[,2])]
    H0_t_r = H0_t[which(censor[,1])]
    H0_t_l = H0_t[which(censor[,3])]
    H0_t_i2 = H0_t[which(censor[,4])]
    
    h0_t_quad_i1 = quad.psi.y_i1 %*% theta
    mu_t_quad_i1 = sapply(seq_along(quad.phi.y_i1), function(x) quad.phi.y_i1[[x]] %*% alpha[[x]])
    nu_t_quad_i1 = sapply(seq_along(quad.phi.event.re), function(x) apply(quad.phi.y_i1.re[[x]] * a_re_cols_quad[[x]][i1_ind_quad == 1,], 1, sum))
    z_t_quad_i1 = mu_t_quad_i1 + nu_t_quad_i1
    exp_zTg_quad_i1 = exp(z_t_quad_i1 %*% gamma)
    h0_t_star_quad_i1 = h0_t_quad_i1 * exp_zTg_quad_i1
    h0_t_star_quad_i1 = matrix(h0_t_star_quad_i1, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
    H0_t_i1 = quad.lambda_i1 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i1, 2, sum)
    
    #h0_t_quad_i2 = quad.psi.y_i2 %*% theta
    #z_t_quad_i2 = (quad.phi.event[i1_ind_quad ==1,] %*% alpha) + apply((quad.phi.event.re * a_re_cols_quad)[i1_ind_quad ==1,], 1, sum)
    #exp_zTg_quad_i2 = exp(c(gamma) * z_t_quad_i2)
    #h0_t_star_quad_i2 = h0_t_quad_i2 * exp_zTg_quad_i2
    #h0_t_star_quad_i2 = matrix(h0_t_star_quad_i2, ncol = sum(censor[,4]), nrow = control$gs.rules, byrow = FALSE)
    #H0_t_i2 = H0_t_i1 + quad.lambda_i2 * apply(matrix(rep(quad.w, sum(censor[,4])), nrow = control$gs.rules, byrow = FALSE) *  h0_t_star_quad_i2, 2, sum)
    
    H0_t_i1_long = rep(0, n)
    H0_t_i1_long[which(censor[,4])] = H0_t_i1
    
    S_t = exp(-eXtB * H0_t)
    S_t_i1_long = rep(0, n)
    
    S_t_l = exp(-eXtB_l * H0_t_l)
    S_t_i1 = exp(-eXtB_i * H0_t_i1)
    S_t_i2 = exp(-eXtB_i * H0_t_i2)
    
    S_t = exp(-eXtB * H0_t)
    S_t_i1_long = rep(0, n)
    S_t_i1_long[which(censor[,4])] = S_t_i1
    
    #first derivative of gamma
    A_t_quad = c(h0_t_quad * exp_zTg_quad) * z_t_quad
    df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* A_t_quad)
    A_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(q+1)], list(df[,1]), FUN = sum )[,-1])
    #colnames(A_t) = NULL
    
    A_t_e = A_t[which(censor[,2]),]
    A_t_r = A_t[which(censor[,1]),]
    A_t_l = A_t[which(censor[,3]),]
    A_t_i2 = A_t[which(censor[,4]),]
    
    A_t_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * z_t_quad_i1
    df = data.frame(id_quad = c(sapply(unique(id[which(censor_long[,4])]), rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* A_t_quad_i1)
    A_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(q+1)], list(df[,1]), FUN = sum )[,-1])
    #colnames(A_t_i1) = NULL
    
    A_t_il_long = matrix(0, nrow = n, ncol = q)
    A_t_il_long[which(censor[,4]),] = A_t_i1
    
    #second derivative of gamma
    A2_t_quad_long = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, control$gs.rules)) * A_t_quad
    #colnames(A2_t_quad_long) = NULL
    
    A_t_quad2_e = A2_t_quad_long[c(sapply(censor[,2], rep, control$gs.rules)),]
    A_t_quad2_r = A2_t_quad_long[c(sapply(censor[,1], rep, control$gs.rules)),]
    A_t_quad2_l = A2_t_quad_long[c(sapply(censor[,3], rep, control$gs.rules)),]
    A_t_quad2_i2 = A2_t_quad_long[c(sapply(censor[,4], rep, control$gs.rules)),]
    
    A_t_quad2_i1 = c(rep(quad.w, sum(censor[,4]))) * 
      c(sapply(quad.lambda_i1, rep, control$gs.rules)) * A_t_quad_i1
    #colnames(A_t_quad2_i1) = NULL
    #colnames(z_t_quad) = NULL
    #colnames(z_t_quad_i1) = NULL
    
    #A2_t_i1_long = matrix(0, nrow =n, ncol = q)
    #A2_t_i1_long[which(censor[,4]),] = A_t_quad2_i1
    
    #d RE d2 gamma - not done
    #A_dRE_t_quad2_long = c(c(gamma) * h0_t_quad * exp_zTg_quad * (z_t_quad)^2 + 2 * h0_t_quad * exp_zTg_quad * z_t_quad) * quad.phi.event.re
    #df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* A_dRE_t_quad2_long)
    #A2_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
    
    #A2_dRE_t_e = A2_dRE[which(censor[,2]),]
    #A2_dRE_t_r = A2_dRE[which(censor[,1]),]
    #A2_dRE_t_l = A2_dRE[which(censor[,3]),]
    #A2_dRE_t_i2 = A2_dRE[which(censor[,4]),]
    
    #A_dRE_t_quad_long = c(h0_t_quad * exp_zTg_quad + h0_t_quad * exp_zTg_quad * z_t_quad * c(gamma)) * quad.phi.event.re
    #df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* A_dRE_t_quad_long)
    #A_dRE = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
    
    #first derivative of theta
    psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
    df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* psi_t_star_quad)
    Psi_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
    
    Psi_t_e = Psi_t[which(censor[,2]),]
    Psi_t_r = Psi_t[which(censor[,1]),]
    Psi_t_l = Psi_t[which(censor[,3]),]
    Psi_t_i2 = Psi_t[which(censor[,4]),]
    
    psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
    df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* psi_t_star_quad_i1)
    Psi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-1])
    
    Psi_t_i1_long = matrix(0, nrow = nrow(Psi_t), ncol = ncol(Psi_t))
    Psi_t_i1_long[which(censor[,4]),] = Psi_t_i1
    
    #first derivative of alpha
    B_t_quad = c(h0_t_quad * exp_zTg_quad) * data.frame(lapply(seq_along(quad.phi.event), function(x) gamma[x] * quad.phi.event[[x]]))
    df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* B_t_quad)
    B_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r_all+1)], list(df[,1]), FUN = sum )[,-1])
    
    Bphi_t_e = B_t[which(censor[,2]),]
    Bphi_t_r = B_t[which(censor[,1]),]
    Bphi_t_l = B_t[which(censor[,3]),]
    Bphi_t_i2 = B_t[which(censor[,4]),]
    
    B_t_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * data.frame(lapply(seq_along(quad.phi.y_i1), function(x) gamma[x] * quad.phi.y_i1[[x]]))
    df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* B_t_quad_i1)
    Bphi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r_all+1)], list(df[,1]), FUN = sum )[,-1])
    
    #B_t_i1_long = matrix(0, nrow = nrow(B_t), ncol = ncol(B_t))
    #B_t_i1_long[which(censor[,4]),] = Bphi_t_i1
    
    B_t_i1_long = matrix(0, nrow = nrow(B_t), ncol = ncol(B_t))
    B_t_i1_long[which(censor[,4]),] = Bphi_t_i1
    
    
    #re version
    B_t_quad_re = c(h0_t_quad * exp_zTg_quad) * data.frame(lapply(seq_along(quad.phi.event.re), function(x) gamma[x] * quad.phi.event.re[[x]]))
    df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* B_t_quad_re)
    B_t_re = c(quad.lambda) * as.matrix(aggregate( df[,2:(w_all+1)], list(df[,1]), FUN = sum )[,-1])
    
    Bphi_t_re_e = B_t_re[which(censor[,2]),]
    Bphi_t_re_r = B_t_re[which(censor[,1]),]
    Bphi_t_re_l = B_t_re[which(censor[,3]),]
    Bphi_t_re_i2 = B_t_re[which(censor[,4]),]
    
    B_t_quad_re_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * data.frame(lapply(seq_along(quad.phi.y_i1.re), function(x) gamma[x] * quad.phi.y_i1.re[[x]]))
    df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* B_t_quad_i1)
    Bphi_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(w_all+1)], list(df[,1]), FUN = sum )[,-1])
    
    B_t_i1_re_long = matrix(0, nrow = nrow(B_t_re), ncol = ncol(B_t_re))
    B_t_i1_re_long[which(censor[,4]),] = Bphi_t_re_i1
    
    #second derivative of alpha
    B_t_quad = c(h0_t_quad * exp_zTg_quad) * data.frame(lapply(seq_along(quad.phi.event), function(x) gamma[x]*quad.phi.event[[x]]))
    B_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, control$gs.rules)) * B_t_quad 
    B_t_quad2_e = B_t_quad2[c(sapply(censor[,2], rep, control$gs.rules)),]
    B_t_quad2_r = B_t_quad2[c(sapply(censor[,1], rep, control$gs.rules)),]
    B_t_quad2_l = B_t_quad2[c(sapply(censor[,3], rep, control$gs.rules)),]
    B_t_quad2_i2 = B_t_quad2[c(sapply(censor[,4], rep, control$gs.rules)),]
    B_t_quad_i1_2 = c(rep(quad.w, sum(censor[,4]))) * 
      c(sapply(quad.lambda_i1, rep, control$gs.rules)) * B_t_quad_i1 
    
    # gamma alpha hessian
    
    
    #E_t_quad = c(gamma) *  c(h0_t_quad * exp_zTg_quad * z_t_quad ) * quad.phi.event + 
    #  c(h0_t_quad * exp_zTg_quad ) * quad.phi.event
    #df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* E_t_quad)
    #E_t = c(quad.lambda) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    #Ephi_t_e = E_t[which(censor[,2]),]
    #Ephi_t_r = E_t[which(censor[,1]),]
    #Ephi_t_l = E_t[which(censor[,3]),]
    #Ephi_t_i2 = E_t[which(censor[,4]),]
    
    #E_t_quad_i1 = c(gamma) *  c(h0_t_quad_i1 * exp_zTg_quad_i1 * z_t_quad_i1) * quad.phi.y_i1 + 
    #  c(h0_t_quad_i1 * exp_zTg_quad_i1 ) * quad.phi.y_i1
    #df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* E_t_quad_i1)
    #Ephi_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-1])
    
    #Ephi_t_i1_long = matrix(0, nrow = nrow(E_t), ncol = ncol(E_t))
    #Ephi_t_i1_long[which(censor[,4]),] = Ephi_t_i1
    
    
    # gamma alpha hessian: re version
    #E_t_quad_re = (h0_t_quad * exp_zTg_quad * z_t_quad ) * quad.phi.event.re + 
    #  c(h0_t_quad * exp_zTg_quad ) * quad.phi.event.re
    #df = data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* E_t_quad_re)
    #E_t_re = c(quad.lambda) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
    
    #Ephi_t_re_e = E_t_re[which(censor[,2]),]
    #Ephi_t_re_r = E_t_re[which(censor[,1]),]
    #Ephi_t_re_l = E_t_re[which(censor[,3]),]
    #Ephi_t_re_i2 = E_t_re[which(censor[,4]),]
    
    #E_t_quad_re_i1 = c(gamma) *  c(h0_t_quad_i1 * exp_zTg_quad_i1 * z_t_quad_i1) * quad.phi.y_i1.re + 
    #  c(h0_t_quad_i1 * exp_zTg_quad_i1 ) * quad.phi.y_i1.re
    #df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* E_t_quad_re_i1)
    #Ephi_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
    
    #Ephi_t_re_i1_long = matrix(0, nrow = nrow(E_t_re), ncol = ncol(E_t_re))
    #Ephi_t_re_i1_long[which(censor[,4]),] = Ephi_t_re_i1
    
    
    # theta gamma hessian (need to multiply this with z_quad)
    psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
    psi_t_quad2 = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, control$gs.rules)) * psi_t_star_quad 
    psi_t_quad2_e = psi_t_quad2[c(sapply(censor[,2], rep, control$gs.rules)),]
    psi_t_quad2_r = psi_t_quad2[c(sapply(censor[,1], rep, control$gs.rules)),]
    psi_t_quad2_l = psi_t_quad2[c(sapply(censor[,3], rep, control$gs.rules)),]
    psi_t_quad2_i2 = psi_t_quad2[c(sapply(censor[,4], rep, control$gs.rules)),]
    
    psi_t_star_quad_i1 = c(exp_zTg_quad_i1) * quad.psi.y_i1
    psi_t_quad2_i1 = c(rep(quad.w, sum(censor[,4]))) * 
      c(sapply(quad.lambda_i1, rep, control$gs.rules)) * psi_t_star_quad_i1 
    
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
    betagamma_hess_e = t(fixed_e) %*% (c(-eXtB_e) * A_t_e)
    betagamma_hess_r = t(fixed_r) %*% (c(-eXtB_r) * A_t_r)
    betagamma_hess_l = t(fixed_l) %*% (c(eXtB_l * (S_t_l)/ (1 - S_t_l)) * A_t_l) -
      t(fixed_l) %*% (c(eXtB_l^2 * S_t_l * H0_t_l/ (1 - S_t_l)^2) * A_t_l)
    betagamma_hess_i = t(fixed_i) %*% (c(-eXtB_i * S_t_i1 /(S_t_i1 - S_t_i2)) * A_t_i1) +
      t(fixed_i) %*% (c(eXtB_i * S_t_i2 /(S_t_i1 - S_t_i2)) * A_t_i2) -
      t(fixed_i) %*% (c(eXtB_i * S_t_i1 * S_t_i2 * (H0_t_i1 * eXtB_i - H0_t_i2 * eXtB_i)/(S_t_i1 - S_t_i2)^2) * (A_t_i1 - A_t_i2))
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
    betaalpha_hess_e = t(c(-eXtB_e) * fixed_e) %*% Bphi_t_e
    betaalpha_hess_r = t(c(-eXtB_r) * fixed_r) %*% Bphi_t_r
    betaalpha_hess_l = t(fixed_l * c(eXtB_l * (S_t_l - eXtB_l*H0_t_l*S_t_l - S_t_l^2)  / (1 - S_t_l)^2) ) %*% Bphi_t_l
    betaalpha_hess_i = t(fixed_i) %*% (c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)) * Bphi_t_i1) + 
      t(fixed_i) %*% (c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)) * Bphi_t_i2) -
      t(fixed_i * c(eXtB_i * S_t_i1 * S_t_i2 * (eXtB_i*H0_t_i1 - eXtB_i*H0_t_i2)/(S_t_i1 - S_t_i2)^2) ) %*% (Bphi_t_i1 - Bphi_t_i2)
    betaalpha_hess = betaalpha_hess_e + betaalpha_hess_r + betaalpha_hess_l + betaalpha_hess_i
    Hessian[1:p, (p+q+m+1):(p+q+m+r_all)] = betaalpha_hess
    Hessian[(p+q+m+1):(p+q+m+r_all), 1:p] = t(Hessian[1:p, (p+q+m+1):(p+q+m+r_all)])
    
    #beta kappa
    betaa1_hess = do.call(rbind, lapply(seq_len(w_all), function(X) fixed)) * 
      c(c(censor[,2] + censor[,1]) * c(eXtB) * B_t_re) +
      do.call(rbind, lapply(seq_len(w_all), function(X) fixed)) * 
      c(c(censor[,3]) * c(eXtB * S_t/(1-S_t)) * B_t_re) -
      do.call(rbind, lapply(seq_len(w_all), function(X) fixed)) * 
      c(c(censor[,3]) *  c(eXtB * S_t * (H0_t*eXtB)/(1-S_t)^2) * B_t_re) -
      do.call(rbind, lapply(seq_len(w_all), function(X) fixed)) * 
      c(c(censor[,4]) * c(eXtB * S_t_i1_long/(S_t_i1_long-S_t)) * B_t_i1_re_long) + 
      do.call(rbind, lapply(seq_len(w_all), function(X) fixed)) * 
      c(c(censor[,4]) * c(eXtB * S_t/(S_t_i1_long-S_t)) * B_t_re) - 
      do.call(rbind, lapply(seq_len(w_all), function(X) fixed)) * 
      c(c(censor[,4]) * c(eXtB * S_t_i1_long * S_t * (H0_t_i1_long*eXtB - H0_t*eXtB)/(S_t_i1_long-S_t)^2) * (B_t_i1_re_long - B_t_re))
    Hessian[1:p, (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)] = betaa1_hess
    Hessian[(p+q+m+r_all+1):(p+q+m+r_all+w_all*n), 1:p] = t(Hessian[1:p, (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)])
    
    #gamma
    #gamma gamma
    if(sum(censor[,2]) > 0){
      gamma_hess_e = t(-c(sapply(c(eXtB_e), rep, control$gs.rules)) * as.matrix(z_t_quad[c(sapply(censor[,2], rep, control$gs.rules)),])) %*% as.matrix(A_t_quad2_e)
    }else{
      gamma_hess_e = 0
    }
    
    gamma_hess_r = t(-c(sapply(c(eXtB_r), rep, control$gs.rules)) * as.matrix(z_t_quad[c(sapply(censor[,1], rep, control$gs.rules)),])) %*% as.matrix(A_t_quad2_r)
    gamma_hess_l = t(c(sapply(c(eXtB_l * (S_t_l)/(1-S_t_l)), rep, control$gs.rules)) * as.matrix(z_t_quad[c(sapply(censor[,3], rep, control$gs.rules)),])) %*% as.matrix(A_t_quad2_l) - 
      t(c( S_t_l * eXtB_l^2/(1-S_t_l)^2) * A_t_l) %*% A_t_l 
    gamma_hess_i = t(c(sapply(c(-eXtB_i *S_t_i1 / (S_t_i1 - S_t_i2)), rep, control$gs.rules)) * as.matrix(z_t_quad_i1)) %*% as.matrix(A_t_quad2_i1) +
      t(c(sapply(c(eXtB_i *S_t_i2 / (S_t_i1 - S_t_i2)), rep, control$gs.rules)) * as.matrix(z_t_quad[c(sapply(censor[,4], rep, control$gs.rules)),])) %*% as.matrix(A_t_quad2_i2) -
      t(c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2 ) * (A_t_i1 - A_t_i2)) %*% (A_t_i1 - A_t_i2)
    gamma_hess = gamma_hess_e + gamma_hess_r + gamma_hess_l + gamma_hess_i
    Hessian[(p+1):(p+q), (p+1):(p+q)] = gamma_hess
    
    #gamma theta
    if(sum(censor[,2]) > 0){
      gammatheta_hess_e = t(c(sapply(c(-eXtB_e), rep, control$gs.rules)) * psi_t_quad2_e)%*% as.matrix(z_t_quad)[c(sapply(censor[,2], rep, control$gs.rules)),]
    }else{
      gammatheta_hess_e = 0
    }
    gammatheta_hess_r = t(c(sapply(c(-eXtB_r), rep, control$gs.rules)) * psi_t_quad2_r) %*% as.matrix(z_t_quad)[c(sapply(censor[,1], rep, control$gs.rules)),]
    gammatheta_hess_l = t(c(sapply(c(eXtB_l * (S_t_l)/(1-S_t_l)), rep, control$gs.rules)) *psi_t_quad2_l) %*% as.matrix(z_t_quad)[c(sapply(censor[,3], rep, control$gs.rules)),] -
      t(Psi_t_l) %*% (c(eXtB_l^2 *(S_t_l)/(1-S_t_l)^2)* A_t_l) 
    gammatheta_hess_i = t(c(sapply(c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)), rep, control$gs.rules)) *psi_t_quad2_i1) %*% as.matrix(z_t_quad_i1) + 
      t(c(sapply(c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)), rep, control$gs.rules)) *psi_t_quad2_i2) %*% as.matrix(z_t_quad)[c(sapply(censor[,4], rep, control$gs.rules)),] -
      t(Psi_t_i1 - Psi_t_i2) %*% (c(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2) * (A_t_i1 - A_t_i2)) 
    gammatheta_hess = gammatheta_hess_e + gammatheta_hess_r + gammatheta_hess_l + gammatheta_hess_i
    Hessian[(p+1):(p+q), (p+q+1):(p+q+m)] = gammatheta_hess
    Hessian[(p+q+1):(p+q+m), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+1):(p+q+m)])
    
    #gamma alpha
    if(sum(censor[,2])>0){
      gammaalpha_hess_e =  t(c(-eXtB_e) *  as.matrix(A_t_quad2_e)) %*% as.matrix(B_t_quad[c(sapply(censor[,2], rep, control$gs.rules)),])
    }else{
      gammaalpha_hess_e = 0
    }
    
    
    gammaalpha_hess_r = t(c(sapply(c(-eXtB_r), rep, control$gs.rules)) *  as.matrix(A_t_quad2_r)) %*% as.matrix(B_t_quad[c(sapply(censor[,1], rep, control$gs.rules)),])
    gammaalpha_hess_l = t(c(sapply(c(eXtB_l * (S_t_l )/((1-S_t_l))), rep, control$gs.rules)) *  as.matrix(A_t_quad2_l)) %*% as.matrix(B_t_quad[c(sapply(censor[,3], rep, control$gs.rules)),]) -
      t(as.numeric(eXtB_l^2 * S_t_l/((1-S_t_l)^2)) * A_t_l) %*% Bphi_t_l
    gammaalpha_hess_i = t(c(sapply(-eXtB_i * S_t_i1/(S_t_i1 - S_t_i2), rep, control$gs.rules)) *  as.matrix(A_t_quad2_i1)) %*% as.matrix(B_t_quad_i1) + 
      t(c(sapply(-eXtB_i * S_t_i2/(S_t_i1 - S_t_i2), rep, control$gs.rules)) *  as.matrix(A_t_quad2_i2)) %*% as.matrix(B_t_quad[c(sapply(censor[,4], rep, control$gs.rules)),])-
      t(as.numeric(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2)* (A_t_i1 - A_t_i2)) %*% (Bphi_t_i1-Bphi_t_i2)
    gammaalpha_hess = gammaalpha_hess_e + gammaalpha_hess_r + gammaalpha_hess_l + gammaalpha_hess_i
    Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r_all)] = gammaalpha_hess
    Hessian[(p+q+m+1):(p+q+m+r_all), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+1):(p+q+m+r_all)])
    
    #gamma kappa = COME BACK TO THIS
    gamma_quad_phi_re = data.frame(lapply(seq_along(quad.phi.event.re), function(x) gamma[x]*quad.phi.event.re[[x]]))
    gamma_quad_phi_i1_re = data.frame(lapply(seq_along(quad.phi.y_i1.re), function(x) gamma[x]*quad.phi.y_i1.re[[x]]))
    
    ##gamma goes with the z not the B_t
    ## create a B_t without the gamma
    ## add z with gamma
    ## use this to create E_t, then use E_t or E_t + B_t
    
    temp2 = lapply(seq_along(quad.phi.event.re), function(x) c(h0_t_quad * exp_zTg_quad) * (quad.phi.event.re[[x]]))
    #df2 = lapply(seq_along(temp2), function(x) data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* temp2[[x]]))
    #B_t_re_g_l = lapply(seq_along(df2), function(x)  c(quad.lambda) * as.matrix(aggregate( df2[[x]][,2:(w[[x]]+1)], list(df2[[x]][,1]), FUN = sum )[,-1]))
    
    z_t_quad_gamma = data.frame(lapply(seq_along(gamma), function(x) gamma[x] * z_t_quad[,x]))
    temp = do.call(rbind, lapply(seq_len(w_all), function(x) z_t_quad_gamma)) * 
      unlist(temp2)
    df = data.frame(id_quad = rep(c(sapply(unique(id), rep, control$gs.rules)), w_all),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*n)), c(rep(quad.w, n)) * temp)
    E_t_re = c(quad.lambda) * as.matrix(aggregate( df[,3:(q+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    temp2_i1 = lapply(seq_along(quad.phi.y_i1.re), function(x) c(h0_t_quad_i1 * exp_zTg_quad_i1) * (quad.phi.y_i1.re[[x]]))
    z_t_quad_gamma_i1 = data.frame(lapply(seq_along(gamma), function(x) gamma[x] * z_t_quad_i1[,x]))
    temp = do.call(rbind, lapply(seq_len(w_all), function(X) z_t_quad_gamma_i1)) * 
      c(unlist(temp2_i1))
    df = data.frame(id_quad = rep(c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), w_all),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*sum(censor[,4]))), c(rep(quad.w, sum(censor[,4]))) * temp)
    E_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(q+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    E_t_re_i1_long = matrix(0, nrow = n*w_all, ncol = q)
    E_t_re_i1_long[c(sapply(1:w_all, function(x) x*which(censor[,4]))), ] = E_t_re_i1 
    
    B_t_re_g = B_t_re_g_i1 = matrix(0, nrow = nrow(E_t_re), ncol = ncol(E_t_re))
    temp2 = lapply(seq_along(quad.phi.event.re), function(x) c(h0_t_quad * exp_zTg_quad) * (quad.phi.event.re[[x]]))
    df2 = lapply(seq_along(temp2), function(x) data.frame(id_quad = c(sapply(unique(id), rep, control$gs.rules)), c(rep(quad.w, n))* temp2[[x]]))
    B_t_re_g_l = lapply(seq_along(df2), function(x)  c(quad.lambda) * as.matrix(aggregate( df2[[x]][,2:(w[[x]]+1)], list(df2[[x]][,1]), FUN = sum )[,-1]))
    
    temp2_i1 = lapply(seq_along(quad.phi.y_i1.re), function(x) c(h0_t_quad_i1 * exp_zTg_quad_i1) * (quad.phi.y_i1.re[[x]]))
    df2_i1 = lapply(seq_along(temp2_i1), function(x) data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* temp2_i1[[x]]))
    B_t_re_g_l_i1 = lapply(seq_along(df2_i1), function(x)  c(quad.lambda_i1) * as.matrix(aggregate( df2_i1[[x]][,2:(w[[x]]+1)], list(df2_i1[[x]][,1]), FUN = sum )[,-1]))
    B_t_re_g_l_i1_long = list() 
    
    phi.re_g = matrix(0, ncol=q, nrow = n*w_all)
    
    for(cov in 1:q){
      B_t_re_g_l_i1_long[[cov]] = matrix(0,nrow=n,ncol=w[[cov]])
      B_t_re_g_l_i1_long[[cov]][which(censor[,4]),] = B_t_re_g_l_i1[[cov]]
      
      if(cov==1){
        B_t_re_g[1:(n*w[[cov]]),cov] = c(B_t_re_g_l[[cov]])
        B_t_re_g_i1[1:(n*w[[cov]]),cov] = c(B_t_re_g_l_i1_long[[cov]])
        phi.re_g[1:(n*w[[cov]]),cov] = sapply(seq_along(phi.re), function(x) c(phi.re[[x]]))[[cov]]
        
      }else{
        B_t_re_g[(n*w[[cov-1]] + 1):(n*w[[cov-1]] + n*w[[cov]]),cov] = c(B_t_re_g_l[[cov]])
        B_t_re_g_i1[(n*w[[cov-1]] + 1):(n*w[[cov-1]] + n*w[[cov]]),cov] = c(B_t_re_g_l_i1_long[[cov]])
        phi.re_g[(n*w[[cov-1]] + 1):(n*w[[cov-1]] + n*w[[cov]]),cov] = sapply(seq_along(phi.re), function(x) c(phi.re[[x]]))[[cov]]
      }
    }
    
    E_t_re = E_t_re + B_t_re_g
    E_t_re_i1_long = E_t_re_i1_long + B_t_re_g_i1
    
    gammaa1_hess = (rep(censor[,2],w_all) * phi.re_g) - 
                       (censor[,2] + censor[,1]) * (c(eXtB) * E_t_re) +
      (c(censor[,3]) * (c(eXtB * S_t/(1-S_t))) * E_t_re) -
      (c(censor[,3]) * c(gamma) * (c((eXtB)^2 * S_t/(1-S_t)^2) * do.call(rbind, lapply(seq_len(w_all), function(X) A_t)) * c(B_t_re))) - 
      (c(censor[,4]) * (c(eXtB * S_t_i1_long/(S_t_i1_long-S_t)) * E_t_re_i1_long)) + 
      (c(censor[,4]) * (c(eXtB * S_t/(S_t_i1_long-S_t)) * E_t_re)) - 
      (c(censor[,4]) * c(gamma) * (c((eXtB)^2 * S_t_i1_long * S_t/(S_t_i1_long-S_t)^2) * do.call(rbind, lapply(seq_len(w_all), function(X) A_t_il_long - A_t)) * c(B_t_i1_re_long -B_t_re)))
    Hessian[(p+1):(p+q), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)] = gammaa1_hess
    Hessian[(p+q+m+r_all+1):(p+q+m+r_all+w_all*n), (p+1):(p+q)] = t(Hessian[(p+1):(p+q), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)])
    
    #theta
    theta_hess_e = - t(as.numeric(1/(h0_t_e^2)) * psi_e) %*% psi_e
    theta_hess_l = - t(as.numeric(eXtB_l^2 * S_t_l/(1 - S_t_l)^2) * Psi_t_l) %*% Psi_t_l
    theta_hess_i = - t(as.numeric(eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1 - S_t_i2)^2) * (Psi_t_i1 - Psi_t_i2)) %*% (Psi_t_i1 - Psi_t_i2)
    theta_hess = theta_hess_e + theta_hess_l + theta_hess_i
    Hessian[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = theta_hess
    
    #theta alpha
    gamma_quad_phi = data.frame(lapply(seq_along(quad.phi.event), function(x) gamma[x]*quad.phi.event[[x]]))
    gamma_quad_phi_i1 = data.frame(lapply(seq_along(quad.phi.y_i1), function(x) gamma[x]*quad.phi.y_i1[[x]]))
    
    if(sum(censor[,2])>0){
      thetaalpha_hess_e = t(c(sapply(c(-eXtB_e), rep, control$gs.rules)) *psi_t_quad2_e) %*% as.matrix(c(rep(quad.w, sum(censor[,1]))) * as.matrix(c(sapply(quad.lambda[which(censor[,2])], rep, control$gs.rules)) * gamma_quad_phi[c(sapply(censor[,2], rep, control$gs.rules)),]))
    }else{
      thetaalpha_hess_e = 0
    }
    
    thetaalpha_hess_r = t(c(sapply(c(-eXtB_r), rep, control$gs.rules)) * psi_t_quad2_r) %*% as.matrix( gamma_quad_phi[c(sapply(censor[,1], rep, control$gs.rules)),])
    thetaalpha_hess_l = t(c(sapply(c(eXtB_l * (S_t_l)/((1-S_t_l))), rep, control$gs.rules)) * psi_t_quad2_l) %*% as.matrix(gamma_quad_phi[c(sapply(censor[,3], rep, control$gs.rules)),]) -
      t(c(eXtB_l^2 * S_t_l/((1-S_t_l)^2)) * Psi_t_l) %*% Bphi_t_l
    thetaalpha_hess_i = t(c(sapply(c(-eXtB_i * (S_t_i1)/(S_t_i1 - S_t_i2)), rep, control$gs.rules)) * psi_t_quad2_i1) %*% as.matrix(gamma_quad_phi_i1) +
      t(c(sapply(c(eXtB_i * (S_t_i2)/(S_t_i1 - S_t_i2)), rep, control$gs.rules)) * psi_t_quad2_i2) %*% as.matrix(gamma_quad_phi[c(sapply(censor[,4], rep, control$gs.rules)),]) -
      t(c(eXtB_i^2 * S_t_i1 * S_t_i2/((S_t_i1-S_t_i2)^2)) * (Psi_t_i1 - Psi_t_i2)) %*% (Bphi_t_i1 - Bphi_t_i2)
    Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r_all)] = thetaalpha_hess_e + thetaalpha_hess_r + thetaalpha_hess_l + thetaalpha_hess_i
    Hessian[(p+q+m+1):(p+q+m+r_all), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+1):(p+q+m+r_all)])
    
    #theta a1
    gamma_quad_phi_re = data.frame(lapply(seq_along(quad.phi.event.re), function(x) gamma[x]*quad.phi.event.re[[x]]))
    gamma_quad_phi_i1_re = data.frame(lapply(seq_along(quad.phi.y_i1.re), function(x) gamma[x]*quad.phi.y_i1.re[[x]]))
    temp = (unlist(c(gamma_quad_phi_re))) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) psi_t_quad2))
    df = data.frame(id_quad = rep(c(sapply(unique(id), rep, control$gs.rules)), w_all),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*n)), c(rep(quad.w, n)) * temp)
    F_t = c(quad.lambda) * as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    temp = (unlist(c(gamma_quad_phi_i1_re))) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) psi_t_quad2_i1))
    df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*sum(censor[,4]))), c(rep(quad.w, sum(censor[,4]))) * temp)
    F_t_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(m+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    F_t_i1_long = matrix(0, nrow = nrow(F_t), ncol = ncol(F_t))
    F_t_i1_long[rep(which(censor[,4]), w_all),] = F_t_i1
    
    thetaa1_hess = rep((- c(censor[,2] + censor[,1]) * eXtB), w_all) * F_t + 
      rep((c(censor[,3]) * eXtB * S_t/(1-S_t)), w_all) * F_t -
      rep((c(censor[,3]) * (eXtB^2) * S_t/((1-S_t)^2)), w_all) * c(B_t_re) * do.call(rbind, lapply(seq_len(w_all), function(X) Psi_t)) -
      rep((c(censor[,4]) * eXtB * S_t_i1_long/(S_t_i1_long-S_t)), w_all) * F_t_i1_long + 
      rep((c(censor[,4]) * eXtB * S_t/(S_t_i1_long-S_t) ), w_all) * F_t -
      rep((c(censor[,4]) * (eXtB^2) * S_t * S_t_i1_long/((S_t_i1_long-S_t)^2)), w_all) * c(B_t_i1_re_long - B_t_re) * do.call(rbind, lapply(seq_len(w_all), function(X) (Psi_t_i1_long - Psi_t)))
    
    Hessian[(p+q+1):(p+q+m), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)] = t(thetaa1_hess)
    Hessian[(p+q+m+r_all+1):(p+q+m+r_all+w_all*n), (p+q+1):(p+q+m)] = t(Hessian[(p+q+1):(p+q+m), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)])
    
    
    #alpha
    if(sum(censor[,2])>0){
      alpha_hess_e = t(-c(sapply(c(eXtB_e), rep, control$gs.rules)) * B_t_quad2_e) %*% as.matrix((gamma_quad_phi)[c(sapply(censor[,2], rep, control$gs.rules)),])
    }else{
      alpha_hess_e = 0
    }
    ## gaussian quad stuff in both parts or only one???
    alpha_hess_r = t(-c(sapply(c(eXtB_r), rep, control$gs.rules)) * B_t_quad2_r) %*% as.matrix((gamma_quad_phi)[c(sapply(censor[,1], rep, control$gs.rules)),])
    alpha_hess_l = t(c(sapply(c(eXtB_l*(S_t_l)/(1-S_t_l)), rep, control$gs.rules)) * B_t_quad2_l) %*% as.matrix(gamma_quad_phi[c(sapply(censor[,3], rep, control$gs.rules)),]) -
      t(c((eXtB_l^2 * S_t_l/(1-S_t_l)^2)) * Bphi_t_l) %*% Bphi_t_l
    alpha_hess_i = t(-c(sapply(c(eXtB_i*S_t_i1 /(S_t_i1-S_t_i2)), rep, control$gs.rules)) * B_t_quad_i1_2) %*% as.matrix(gamma_quad_phi_i1)  + 
      t(c(sapply(c(eXtB_i*S_t_i2 /(S_t_i1-S_t_i2)), rep, control$gs.rules)) * B_t_quad2_i2) %*% as.matrix(gamma_quad_phi[c(sapply(censor[,4], rep, control$gs.rules)),]) -
      t(c((eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i1 + 
      t(c((eXtB_i^2 * S_t_i1 * S_t_i2/(S_t_i1-S_t_i2)^2)) * (Bphi_t_i1 - Bphi_t_i2)) %*% Bphi_t_i2
    Hessian[(p+q+m+1):(p+q+m+r_all),(p+q+m+1):(p+q+m+r_all)] = alpha_hess_e + alpha_hess_r + alpha_hess_l + alpha_hess_i
    
    #alpha a1
    gamma_quad_phi_re = data.frame(lapply(seq_along(quad.phi.event.re), function(x) gamma[x]*quad.phi.event.re[[x]]))
    gamma_quad_phi_i1_re = data.frame(lapply(seq_along(quad.phi.y_i1.re), function(x) gamma[x]*quad.phi.y_i1.re[[x]]))
    
    temp = (unlist(c(gamma_quad_phi_re))) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) B_t_quad2))
    df = data.frame(id_quad = rep(c(sapply(unique(id), rep, control$gs.rules)), w_all),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*n)), c(rep(quad.w, n)) * temp)
    D_t_re = c(quad.lambda) * as.matrix(aggregate( df[,3:(r_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    temp = (unlist(c(gamma_quad_phi_i1_re))) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) B_t_quad_i1_2))
    df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*sum(censor[,4]))), c(rep(quad.w, sum(censor[,4]))) * temp)
    D_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(r_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    D_t_re_i1_long = matrix(0, nrow = nrow(D_t_re), ncol = ncol(D_t_re))
    D_t_re_i1_long[rep(which(censor[,4]), w_all),] = D_t_re_i1
    
    alphaa1_hess = - rep((c(censor[,2] + censor[,1]) * eXtB), w_all) * D_t_re + 
      rep((c(censor[,3]) * eXtB * (S_t/(1-S_t)) ), w_all) * D_t_re - 
      c(rep((c(censor[,3]) * (eXtB^2) * (S_t/(1-S_t)^2)), w_all) * B_t_re) * do.call(rbind, lapply(seq_len(w_all), function(X) B_t)) -
      rep((c(censor[,4]) * eXtB * (S_t_i1_long/(S_t_i1_long-S_t))), w_all) * D_t_re_i1_long + 
      rep((c(censor[,4]) * eXtB * (S_t/(S_t_i1_long-S_t))), w_all) * D_t_re - 
      c(rep((c(censor[,4]) * (eXtB^2) * (S_t_i1_long*S_t/(S_t_i1_long-S_t)^2)), w_all) * (B_t_i1_re_long - B_t_re)) * do.call(rbind, lapply(seq_len(w_all), function(X) (B_t_i1_long - B_t)))
    
    Hessian[(p+q+m+1):(p+q+m+r_all), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)] = t(alphaa1_hess)
    Hessian[(p+q+m+r_all+1):(p+q+m+r_all+w_all*n), (p+q+m+1):(p+q+m+r_all)] = t(Hessian[(p+q+m+1):(p+q+m+r_all), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)])
    
    #a1 a1
    #temp = (c(quad.phi.event[,re_ind]) * 
    #          do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(gamma^2), rep, control$gs.rules)) * B_t_quad2[,re_ind])))
    #df = data.frame(id_quad = rep(c(sapply(dat.baseline$id, rep, control$gs.rules)), w),
    #                re_no = c(sapply(re_ind, rep, control$gs.rules*n)), temp)
    
    #second derivative of alpha
    B_t_quad_re = c(h0_t_quad * exp_zTg_quad) * data.frame(quad.phi.event.re)
    B_t_quad2_re = c(rep(quad.w, n)) * c(sapply(quad.lambda, rep, control$gs.rules)) * B_t_quad_re 
    
    B_t_quad_re_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * data.frame(quad.phi.y_i1.re)
    #df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)), c(rep(quad.w, sum(censor[,4])))* B_t_quad_i1)
    #Bphi_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-1])
    
    B_t_quad_re_i1_2 = c(rep(quad.w, sum(censor[,4]))) * 
      c(sapply(quad.lambda_i1, rep, control$gs.rules)) * B_t_quad_re_i1
    
    temp = (unlist(c(gamma_quad_phi_re))) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) B_t_quad2_re))
    df = data.frame(id_quad = rep(c(sapply(unique(id), rep, control$gs.rules)), w_all),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*n)), c(rep(quad.w, n)) * temp)
    D_t_re_re = c(quad.lambda) * as.matrix(aggregate( df[,3:(w_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    temp = (unlist(c(gamma_quad_phi_i1_re))) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) B_t_quad_re_i1_2))
    df = data.frame(id_quad = c(sapply(unique(id)[which(censor[,4])], rep, control$gs.rules)),
                    re_no = c(sapply(1:w_all, rep, control$gs.rules*sum(censor[,4]))), c(rep(quad.w, sum(censor[,4]))) * temp)
    D_t_re_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,3:(w_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    D_t_rere_i1_long = matrix(0, nrow = nrow(D_t_re_re), ncol = ncol(D_t_re_re))
    D_t_rere_i1_long[rep(which(censor[,4]), w_all),] = D_t_re_re_i1
    
    a1a1_hess = - rep((c(censor[,2] + censor[,1]) * eXtB), w_all) * D_t_re_re + 
      rep((c(censor[,3]) * eXtB * (S_t/(1-S_t))), w_all) * D_t_re_re - 
      c(rep((c(censor[,3]) * (eXtB^2) * (S_t/(1-S_t)^2)), w_all) * B_t_re) * do.call(rbind, lapply(seq_len(w_all), function(X) B_t_re)) -
      rep((c(censor[,4]) * eXtB * (S_t_i1_long/(S_t_i1_long-S_t)) ), w_all) * D_t_rere_i1_long + 
      rep((c(censor[,4]) * eXtB * (S_t/(S_t_i1_long-S_t))), w_all) * D_t_re_re - 
      c(rep((c(censor[,4]) * (eXtB^2) * (S_t_i1_long*S_t/(S_t_i1_long-S_t)^2)), w_all) * (B_t_i1_re_long - B_t_re)) * do.call(rbind, lapply(seq_len(w_all), function(X) (B_t_i1_re_long - B_t_re)))
    
    for(re in 1:w_all){
      ## check these dimensions??
      diag(Hessian[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = as.matrix(a1a1_hess)[((re-1)*n + 1):((re-1)*n + n),re]
      
      if(re > 1){
        diag(Hessian[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1a1_hess[((re-1)*n + 1):((re-1)*n + n), (re-1)]
        diag(Hessian[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n)]) = diag(Hessian[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)])
      }
      
    }
    
    #add penalty part 
    Q_theta = Q_alpha = Q_a1 = H_epsilon = matrix(0, nrow = (p+q+m+r_all+w_all*n), ncol = (p+q+m+r_all+w_all*n))
    
    #least squares part
    #H_epsilon[(p+q+m+1):(p+q+m+r_all),(p+q+m+1):(p+q+m+r_all)] = - c(1/sigma2_Et) * t(phi.mu) %*% phi.mu #alpha alpha
    
    H_epsilon[(p+q+m+1):(p+q+m+r_all),(p+q+m+1):(p+q+m+r_all)] = - t(as.matrix(data.frame(lapply(seq_along(phi.mu), function(x) (1/sigma2_Et[[x]]) * phi.mu[[x]])))) %*% as.matrix(data.frame(phi.mu))
    
    temp = unlist(data.frame(phi.mu.re)) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) data.frame(lapply(seq_along(phi.mu), function(x) c(1/sigma2_Et[[x]]) * phi.mu[[x]]))))
    df = data.frame(id, re_no = c(sapply(1:w_all, rep, n_long)), temp)
    alphaa1_hess_ls = -  as.matrix(aggregate( df[,3:(r_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    
    H_epsilon[(p+q+m+1):(p+q+m+r_all), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)] = t(alphaa1_hess_ls) #alpha re
    H_epsilon[ (p+q+m+r_all+1):(p+q+m+r_all+w_all*n), (p+q+m+1):(p+q+m+r_all)] = t(H_epsilon[(p+q+m+1):(p+q+m+r_all), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)])
    
    temp = unlist(data.frame(phi.mu.re)) * 
      do.call(rbind, lapply(seq_len(w_all), function(X) data.frame(lapply(seq_along(phi.mu.re), function(x) c(1/sigma2_Et[[x]]) * phi.mu.re[[x]]))))
    df = data.frame(id,
                    re_no = c(sapply(1:w_all, rep, n_long)), temp)
    a1a1_hess_ls = - as.matrix(aggregate( df[,3:(w_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
    ##dimensions issue??
    
    for(re in 1:w_all){
      diag(H_epsilon[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n),re]
      
      if(re > 1){
        diag(H_epsilon[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n), (re-1)]
        diag(H_epsilon[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n)]) = diag(H_epsilon[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)])
      }
      
    }
    
    #theta part
    Q_theta[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = - theta.lambda*theta.G
    
    #alpha part
    #Q_alpha[(p+q+m+1):(p+q+m+r_all),(p+q+m+1):(p+q+m+r_all)] = - alpha.lambda*alpha.G
    ## NEED TO FIX THIS
    
    #re part
    diag(Q_a1[(p+q+m+r_all+1):(p+q+m+r_all+w_all*n),(p+q+m+r_all+1):(p+q+m+r_all+w_all*n)]) = -c(1/c(sapply(unlist(sigma2_re), rep, n)))
    
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
    
    print("Updating variance components...")
    
    df.theta.old = df.theta
    #df.alpha.old = df.alpha
    df.epsilon.old = df.epsilon
    df.a_re.old = df.a_re
    
    df.theta = sum(diag(H_full_inv %*% (-Q_theta)))
    df.theta = m - df.theta
    theta.sigma2 = t(theta) %*% theta.G %*% theta / df.theta
    theta.lambda = c(1/(2*theta.sigma2))
    #theta.lambda = 0
    
    #df.alpha = sum(diag(H_full_inv %*% (-Q_alpha)))
    #df.alpha = r_all - df.alpha
    #alpha.sigma2 = t(alpha[[1]]) %*% alpha.G %*% alpha[[1]] / df.alpha
    #alpha.lambda = c(1/(2*alpha.sigma2))
    alpha.lambda = 0
    
    
    
    
    ## NEED TO WORK THIS OUT FOR q > 1
    #df.epsilon = sum(diag(H_full_inv %*% (-H_epsilon)))
    #df.epsilon = n_long - df.epsilon
    #sigma2_Et = as.numeric(sum((cont - z_t_ia)^2)/ df.epsilon)
    #sigma2_Et = 0.0025
    
    
    for(cov in 1:q){
      H_eps_temp = matrix(0, nrow = (p+q+m+r_all+w_all*n), ncol = (p+q+m+r_all+w_all*n))
      
      if(cov == 1){
        H_eps_temp[(p+q+m+1):(p+q+m+r[[cov]]),(p+q+m+1):(p+q+m+r[[cov]])] = -(1/sigma2_Et[[cov]]) * t(phi.mu[[cov]]) %*% phi.mu[[cov]]
        
        H_eps_temp[(p+q+m+1):(p+q+m+r[[cov]]), (p+q+m+r_all+1):(p+q+m+r_all+(n*w[[cov]]))] = alphaa1_hess_ls[1:(n*w[[cov]]),1:r[[cov]]]
        H_eps_temp[(p+q+m+r_all+1):(p+q+m+r_all+(n*w[[cov]])), (p+q+m+1):(p+q+m+r[[cov]])] = t(H_epsilon[(p+q+m+1):(p+q+m+r[[cov]]), (p+q+m+r_all+1):(p+q+m+r_all+(n*w[[cov]]))])
        
        for(re in 1:w[[cov]]){
          diag(H_eps_temp[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n),re]
          
          if(re > 1){
            diag(H_eps_temp[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n), (re-1)]
            diag(H_eps_temp[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n)]) = diag(H_epsilon[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)])
          }
          
        }
        
        
      }else{
        
        H_eps_temp[(p+q+m+cumsum(unlist(r))[cov-1]+1):(p+q+m+cumsum(unlist(r))[cov-1]+r[[cov]]),(p+q+m+cumsum(unlist(r))[cov-1]+1):(p+q+m+cumsum(unlist(r))[cov-1]+r[[cov]])] = -(1/sigma2_Et[[cov]]) * t(phi.mu[[cov]]) %*% phi.mu[[cov]]
        
        H_eps_temp[(p+q+m+cumsum(unlist(r))[cov-1]+1):(p+q+m+cumsum(unlist(r))[cov-1]+r[[cov]]), (p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(n*w[[cov]]))] = alphaa1_hess_ls[((n*cumsum(unlist(w))[cov-1])+1):((n*cumsum(unlist(w))[cov-1])+(n*w[[cov]])),(cumsum(unlist(r))[cov-1]+1):(cumsum(unlist(r))[cov-1] + r[[cov]])]
        H_eps_temp[(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(n*w[[cov]])),(p+q+m+cumsum(unlist(r))[cov-1]+1):(p+q+m+cumsum(unlist(r))[cov-1]+r[[cov]])] = t(H_epsilon[(p+q+m+cumsum(unlist(r))[cov-1]+1):(p+q+m+cumsum(unlist(r))[cov-1]+r[[cov]]), (p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(n*w[[cov]]))])
        
        for(re in 1:w[[cov]]){
          diag(H_eps_temp[(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + n), (p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + n)]) = a1a1_hess_ls[((n*cumsum(unlist(w))[cov-1])+(re-1)*n + 1):((n*cumsum(unlist(w))[cov-1])+(re-1)*n + n),(cumsum(unlist(w))[cov-1]+re)]
          
          if(re > 1){
            diag(H_eps_temp[(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-2)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-2)*n + n), (p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + n)]) = a1a1_hess_ls[((n*cumsum(unlist(w))[cov-1])+(re-1)*n + 1):((n*cumsum(unlist(w))[cov-1])+(re-1)*n + n), (cumsum(unlist(w))[cov-1]+re-1)]
            diag(H_eps_temp[(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + n), (p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-2)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-2)*n + n)]) = diag(H_epsilon[(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-2)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-2)*n + n), (p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + 1):(p+q+m+r_all+(n*cumsum(unlist(w))[cov-1])+(re-1)*n + n)])
          }
          
        }
        
      }
      
      df.epsilon[cov] = sum(diag(H_full_inv %*% (-H_eps_temp)))
      df.epsilon[cov] = n_long - df.epsilon[cov]
      sigma2_Et[[cov]] = as.numeric(sum((cont[[cov]] - z_t_ia[,cov])^2)/ df.epsilon[cov])
      
    }
    
    
    
    
    if(any(unlist(sigma2_Et) < 0)){
      for(sig in which(unlist(sigma2_Et) < 0)){
        sigma2_Et[[sig]] = as.numeric(sum((cont[[sig]] - z_t_ia[,sig])^2)/ n_long)
      }
    }
    
    sigma2_re_old = sigma2_re
    sigma2_save = NULL
    for(re in 1:w_all){
      sigma2_temp = sum((data.frame(a_re_cols)[,re])^2)/(n)
      sigma2_save = c(sigma2_save, sigma2_temp)
    }
    
    for(cov in 1:q){
      if(cov==1){
        sigma2_re[[cov]] = c(sigma2_save[1:w[[cov]]])
      }else{
        sigma2_re[[cov]] = c(sigma2_save[(cumsum(unlist(w))[cov-1]+1):(cumsum(unlist(w))[cov-1]+w[[cov]])])
      }
      
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
    
    #if(any(sigma2_re < 1e-5)){
    #  sigma2_re[which(sigma2_re < 1e-5)] = 1e-5
    
    #}
    
    #if(any(c(alpha.sigma2, theta.sigma2, sigma2_a0, sigma2_a1, sigma2_Et) < 0)){
    #  print(c(alpha.sigma2, theta.sigma2, sigma2_a0, sigma2_a1, sigma2_Et))
    #  print(c("Variance estimate is negative"))
    #  break
    #}
    
    print(c(it, unlist(sigma2_re), unlist(sigma2_Et)))
    
    if(all(abs(c(df.theta.old - df.theta )) < 1) & 
       all(abs(df.epsilon.old - df.epsilon) < 1) & all(abs(unlist(sigma2_re) - unlist(sigma2_re_old)) < 0.05) ){
      break
    }
    #if(abs(df.epsilon.old - df.epsilon) < 0.5 & all(abs(sigma2_re - sigma2_re_old) < 0.01) ){
    #  break
    #}
    
    
    
  } #end outer loop
  
  
  #add penalty part
  Q_theta = Q_alpha = Q_a1 = H_epsilon = matrix(0, nrow = (p+q+m+r_all+w_all*n), ncol = (p+q+m+r_all+w_all*n))
  
  #least squares part
  #H_epsilon[(p+q+m+1):(p+q+m+r_all),(p+q+m+1):(p+q+m+r_all)] = - c(1/sigma2_Et) * t(phi.mu) %*% phi.mu #alpha alpha
  
  H_epsilon[(p+q+m+1):(p+q+m+r_all),(p+q+m+1):(p+q+m+r_all)] = - t(as.matrix(data.frame(lapply(seq_along(phi.mu), function(x) (1/sigma2_Et[[x]]) * phi.mu[[x]])))) %*% as.matrix(data.frame(phi.mu))
  
  temp = unlist(data.frame(phi.mu.re)) * 
    do.call(rbind, lapply(seq_len(w_all), function(X) data.frame(lapply(seq_along(phi.mu), function(x) c(1/sigma2_Et[[x]]) * phi.mu[[x]]))))
  df = data.frame(id, re_no = c(sapply(1:w_all, rep, n_long)), temp)
  alphaa1_hess_ls = -  as.matrix(aggregate( df[,3:(r_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  
  H_epsilon[(p+q+m+1):(p+q+m+r_all), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)] = t(alphaa1_hess_ls) #alpha re
  H_epsilon[ (p+q+m+r_all+1):(p+q+m+r_all+w_all*n), (p+q+m+1):(p+q+m+r_all)] = t(H_epsilon[(p+q+m+1):(p+q+m+r_all), (p+q+m+r_all+1):(p+q+m+r_all+w_all*n)])
  
  temp = unlist(data.frame(phi.mu.re)) * 
    do.call(rbind, lapply(seq_len(w_all), function(X) data.frame(lapply(seq_along(phi.mu.re), function(x) c(1/sigma2_Et[[x]]) * phi.mu.re[[x]]))))
  df = data.frame(id,
                  re_no = c(sapply(1:w_all, rep, n_long)), temp)
  a1a1_hess_ls = - as.matrix(aggregate( df[,3:(w_all+2)], list(df[,1], df[,2]), FUN = sum )[,-c(1:2)])
  ##dimensions issue??
  
  for(re in 1:w_all){
    diag(H_epsilon[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n),re]
    
    if(re > 1){
      diag(H_epsilon[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1a1_hess_ls[((re-1)*n + 1):((re-1)*n + n), (re-1)]
      diag(H_epsilon[(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n), (p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n)]) = diag(H_epsilon[(p+q+m+r_all+(re-2)*n + 1):(p+q+m+r_all+(re-2)*n + n), (p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)])
    }
    
  }
  
  #theta part
  Q_theta[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = - theta.lambda*theta.G
  
  #alpha part
  #Q_alpha[(p+q+m+1):(p+q+m+r_all),(p+q+m+1):(p+q+m+r_all)] = - alpha.lambda*alpha.G
  
  ## NEED TO FIX THIS
  
  #re part
  diag(Q_a1[(p+q+m+r_all+1):(p+q+m+r_all+w_all*n),(p+q+m+r_all+1):(p+q+m+r_all+w_all*n)]) = -c(1/c(sapply(unlist(sigma2_re), rep, n)))
  
  H_full = Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1
  
  
  
   
  ##### Qu2 matrix for hessian
  Qu = matrix(0,nrow = n,ncol = (p+q+m+r_all+w_all*n))
  
  Qu[censor[,2],(1:p)] = as.matrix((fixed_e) * c(1- eXtB_e * H0_t_e))
  Qu[censor[,1],(1:p)] = as.matrix((fixed_r) * c(- eXtB_r * H0_t_r))
  Qu[censor[,3],(1:p)] = as.matrix((fixed_l) * c(S_t_l * eXtB_l * H0_t_l / (1 - S_t_l)))
  Qu[censor[,4],(1:p)] = as.matrix(- (fixed_i)  * c((S_t_i1 * eXtB_i * H0_t_i1 - S_t_i2 * eXtB_i * H0_t_i2) / (S_t_i1 - S_t_i2)))
  
  Qu[censor[,2],(p+1):(p+q)] = as.matrix(zt_e - c(eXtB_e) * A_t_e)
  Qu[censor[,1],(p+1):(p+q)] = as.matrix(- c(eXtB_r) * A_t_r)
  Qu[censor[,3],(p+1):(p+q)] = as.matrix(c(S_t_l * eXtB_l  / (1 - S_t_l)) * A_t_l)
  Qu[censor[,4],(p+1):(p+q)] = as.matrix(- c(eXtB_i * (S_t_i1) / (S_t_i1 - S_t_i2 ) ) * A_t_i1 + c(eXtB_i * (S_t_i2) / (S_t_i1 - S_t_i2 ) ) * A_t_i2)
  
  Qu[censor[,2],(p+q+1):(p+q+m)] = as.matrix(as.numeric(1/h0_t_e)* psi_e - as.numeric(eXtB_e)* Psi_t_e)
  Qu[censor[,1],(p+q+1):(p+q+m)] = as.matrix(as.numeric(- eXtB_r)* Psi_t_r)
  Qu[censor[,3],(p+q+1):(p+q+m)] = as.matrix(as.numeric(eXtB_l * S_t_l/(1-S_t_l))* Psi_t_l)
  Qu[censor[,4],(p+q+1):(p+q+m)] = as.matrix(- as.numeric(eXtB_i * S_t_i1/(S_t_i1-S_t_i2))* Psi_t_i1 + as.numeric(eXtB_i * S_t_i2/(S_t_i1-S_t_i2))* Psi_t_i2)
  
  
  Qu[censor[,2],(p+q+m+1):(p+q+m+r_all)] = as.matrix(data.frame(lapply(seq_along(phi_e), function(x) gamma[x] * phi_e[[x]])) - as.numeric(eXtB_e) * Bphi_t_e)
  Qu[censor[,1],(p+q+m+1):(p+q+m+r_all)] = as.matrix(- as.numeric(eXtB_r) * Bphi_t_r)
  Qu[censor[,3],(p+q+m+1):(p+q+m+r_all)] = as.matrix(as.numeric(eXtB_l * S_t_l/(1-S_t_l)) * Bphi_t_l)
  Qu[censor[,4],(p+q+m+1):(p+q+m+r_all)] = as.matrix(- as.numeric(eXtB_i * S_t_i1 /(S_t_i1-S_t_i2)) *  Bphi_t_i1 +  
                                                   as.numeric(eXtB_i * S_t_i2 /(S_t_i1-S_t_i2)) *  Bphi_t_i2)
  Qu[,(p+q+m+1):(p+q+m+r_all)] = as.matrix(data.frame(lapply(seq_along(cont), function(x) c(1/sigma2_Et[[x]]) * aggregate((c(cont[[x]] - z_t_ia[,x]) * phi.mu[[x]]), list(id), FUN = sum)[,-1])))
  
  B_t_re_quad_i1 = c(h0_t_quad_i1 * exp_zTg_quad_i1) * data.frame(quad.phi.y_i1.re)
  df = data.frame(id_quad = c(sapply(unique(id)[(censor[,4])], rep, control$gs.rules)), 
                  c(rep(quad.w, sum((censor[,4]))))* B_t_re_quad_i1)
  B_t_re_i1 = c(quad.lambda_i1) * as.matrix(aggregate( df[,2:(w_all+1)], list(df[,1]), FUN = sum )[,-1])
  
  a1i_score_e = data.frame(lapply(seq_along(phi_e.re), function(x) gamma[x] * phi_e.re[[x]])) - data.frame(as.numeric(eXtB_e) * (Bphi_t_re_e))
  a1i_score_r = - as.numeric(eXtB_r) * (Bphi_t_re_r)
  a1i_score_l = as.numeric(eXtB_l * S_t_l/(1-S_t_l)) * (Bphi_t_re_l)
  a1i_score_i = - as.numeric(eXtB_i * S_t_i1/(S_t_i1-S_t_i2)) * (B_t_re_i1) + 
    as.numeric(eXtB_i * S_t_i2/(S_t_i1-S_t_i2)) * (Bphi_t_re_i2)
  
  a1i_score_s = rbind(as.matrix(a1i_score_e), as.matrix(a1i_score_r), as.matrix(a1i_score_l), as.matrix(a1i_score_i))
  
  
  
  a1i_score_ls = lapply(seq_along(cont), function(x) (cont[[x]]-z_t_ia[,x])*phi.mu.re[[x]])
  a1i_score_ls = data.frame(lapply(seq_along(a1i_score_ls), function(x)  c(1/sigma2_Et[[x]]) * as.matrix(aggregate(a1i_score_ls[[x]][,1:(w[[x]])], list(id), FUN = sum )[,-1])))
  
  a1i_score = a1i_score_s + a1i_score_ls -data.frame(a_re_cols)/matrix(rep(unlist(sigma2_re), n), ncol = w_all, byrow = T)
  
  for(re in 1:w_all){
    diag(Qu[,(p+q+m+r_all+(re-1)*n + 1):(p+q+m+r_all+(re-1)*n + n)]) = a1i_score[,re]
  }
  
  #TwoLRtheta = 2*theta.lambda*( theta.G %*% theta)
  #TwoLRalpha = 2*alpha.lambda*( alpha.G %*% alpha)
  
  #check this part
  #Sp = Qu-matrix(rep(c(rep(0,p+q),TwoLRtheta, rep(0,r+w*n)),n),n,byrow=T)/n
  #Sp2 = Sp-matrix(rep(c(rep(0,p+q+m),TwoLRalpha, rep(0,w*n)),n),n,byrow=T)/n
  Qu = t(Qu)%*%Qu
  
  
  pos = rep(1, nrow(H_full))
  for(u in 1:m){
    if((theta[u]< (1e-1) & theta_score[u] <= (-0.001))){
      pos[p+q+u] = 0
    }
  }
  
  
  #calculate asymptotic variance
  #M2_trunc = (-(Hessian[1:(p+q+m+r),1:(p+q+m+r)] + H_epsilon[1:(p+q+m+r),1:(p+q+m+r)] + 
  #                Q_theta[1:(p+q+m+r),1:(p+q+m+r)] + Q_alpha[1:(p+q+m+r),1:(p+q+m+r)]))
  #
  
  print("Calculating asymptotic covariances...")
  
  corr = M2_inv = matrix(0, nrow(H_full), nrow(H_full))
  diag(corr) = pos
  #corr[!pos,]=0
  
  M2_inv[(pos==1), (pos==1)] = solve(-H_full[(pos==1), (pos==1)])
  A_omega = corr %*% M2_inv %*% t(corr)
  
  #cov_H = A_omega %*% (-((Hessian + H_epsilon))) %*% A_omega
  cov_H_RE = A_omega %*% (-((Hessian + H_epsilon + Q_a1))) %*% A_omega
  
  #check this part
  cov_Q = A_omega %*% (Qu) %*% A_omega
  
  #cov_H_RE_2 = solve((-((Hessian + H_epsilon + Q_theta + Q_alpha)))[1:(p+q+m+r_all), 1:(p+q+m+r_all)] - 
  #                     (-((Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1)))[1:(p+q+m+r_all), (p+q+m+r_all+1):(p+q+m+r_all+n*w_all)] %*%
  #                     solve((-((Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1)))[(p+q+m+r_all+1):(p+q+m+r_all+n*w_all), (p+q+m+r_all+1):(p+q+m+r_all+n*w_all)]) %*%
  #                     (-((Hessian + H_epsilon + Q_theta + Q_alpha + Q_a1)))[(p+q+m+r_all+1):(p+q+m+r_all+n*w_all), 1:(p+q+m+r_all)])
  
  
  #param = c(beta, gamma, theta, alpha)
  
  evalues = eigen(cov_H_RE, symmetric = TRUE)$values
  evectors = eigen(cov_H_RE, symmetric = TRUE)$vectors
  
  if(any(evalues < 0)){
    evalues[which(evalues < 0)] = 1e-10
    cov_H_RE_ee = evectors %*% diag(evalues) %*% t(evectors)
  }else{
    cov_H_RE_ee = NA
  }
  
  #output
  pars = list(beta = beta, gamma = gamma, theta = theta, alpha = c(alpha), 
              a_re_cols = a_re_cols)
  score = list(beta_score = beta_score, gamma_score = gamma_score, theta_score = theta_score,
               alpha_score = alpha_score, a1i_score = a1i_score)
  knots = list(theta.knots = event.kn)
  variance = list(sigma2_re = sigma2_re, sigma2_Et = sigma2_Et, 
                  theta.lambda = theta.lambda, alpha.lambda = alpha.lambda)
  penalty = list(R = theta.G)
  out = list(pars = pars, knots = knots, variance = variance, penalty = penalty, 
             score = score, M2_inv = M2_inv, cov_H_RE = cov_H_RE, 
             H_full = H_full, cov_Q = cov_Q,
             ll_save = ll_save, cov_H_RE_ee = cov_H_RE_ee,
             phi = phi, mc = mc, pll = log_lik, iterations = iter, conv_record = conv_record, data = data)
  return(out)
  
}


theta_penalty_f = function(order, IntKnt, bryKnt, norm = 2){
  if(norm == 2){
    
    ordSp = order
    dgrSp = ordSp - 1
    numIntKnt = length(IntKnt)
    numSp = numIntKnt+ordSp
    
    minTime = min(bryKnt)
    maxTime = max(bryKnt)
    
    R=matrix(0, nrow=numSp, ncol=numSp) 
    xknots = c(rep(minTime, ordSp), IntKnt, rep(maxTime, ordSp))
    for (ii in 1:numSp){
      for (jj in ii:numSp){
        if (jj - ii<ordSp){
          kntset = xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]]
          kntsum = 0
          for (kk in 1:(length(kntset)-1)){
            kntsum = kntsum + mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, Boundary.knots=bryKnt, 
                                      derivs=dgrSp)[ii]*mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, 
                                                                Boundary.knots=bryKnt,derivs=dgrSp)[jj]*(kntset[kk+1]-kntset[kk])
          }
          R[ii, jj] = kntsum
        }
      }
    }
    R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
    
  }
  
  
  return(R)
  
  
}




