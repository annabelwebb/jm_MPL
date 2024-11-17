
summary.jm_MPL=function(obj,se="penalised", full=FALSE, HR_display = TRUE){
  col.names = c("Estimate","Std. Error", "z-value", "Pr(>|z|)")
  col.names.HR = c("Hazard Ratio","Lower CI", "Upper CI", "Pr(>|z|)")
  p=length(obj$pars$beta)
  q=length(obj$pars$gamma)
  m=length(obj$pars$theta)
  r=length(obj$pars$alpha)
  
  if(se=="penalised"){
    seB=sqrt(diag(obj$cov_H_RE))[1:p]
    seG=sqrt(diag(obj$cov_H_RE))[(p+1):(p+q)]
    seT=sqrt(diag(obj$cov_H_RE))[(p+q+1):(p+q+m)]
    seA=sqrt(diag(obj$cov_H_RE))[(p+q+m+1):(p+q+m+r)]
    
    matxB = cbind(obj$pars$beta,seB,obj$pars$beta/seB,2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matxB)=col.names
    rownames(matxB)=row.names(try.MLE$pars$beta)
    matHRB = cbind(exp(obj$pars$beta), exp(obj$pars$beta - 1.96*seB),exp(obj$pars$beta + 1.96*seB), 
                   2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matHRB)=col.names.HR
    rownames(matHRB)=row.names(try.MLE$pars$beta)
    
    matxG = cbind(obj$pars$gamma,seG,obj$pars$gamma/seG,2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matxG)=col.names
    rownames(matxG)="Association"
    matHRG = cbind(exp(obj$pars$gamma), exp(obj$pars$gamma - 1.96*seG),exp(obj$pars$gamma + 1.96*seG), 
                   2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matHRG)=col.names.HR
    rownames(matHRG)="Association"
    
    matxReg = rbind(matxB, matxG)
    matHR = rbind(matHRB, matHRG)
    
    matxT = cbind(obj$pars$theta,seT,obj$pars$theta/seT,2*(1-pnorm(abs(obj$pars$theta/seT))))
    colnames(matxT)=col.names
    
    
    matxA = cbind(obj$pars$alpha,seA,obj$pars$alpha/seA,2*(1-pnorm(abs(obj$pars$alpha/seA))))
    colnames(matxA)=col.names
    rownames(matxA)=colnames(data.frame(try.MLE$phi))
    
  }
  if(se=="BLUP"){
    seB=sqrt(diag(obj$cov_H_RE_2))[1:p]
    seG=sqrt(diag(obj$cov_H_RE_2))[(p+1):(p+q)]
    seT=sqrt(diag(obj$cov_H_RE_2))[(p+q+1):(p+q+m)]
    seA=sqrt(diag(obj$cov_H_RE_2))[(p+q+m+1):(p+q+m+r)]
    
    
    matxB = cbind(obj$pars$beta,seB,obj$pars$beta/seB,2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matxB)=col.names
    rownames(matxB)=row.names(try.MLE$pars$beta)
    matHRB = cbind(exp(obj$pars$beta), exp(obj$pars$beta - 1.96*seB),exp(obj$pars$beta + 1.96*seB), 
                   2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matHRB)=col.names.HR
    rownames(matHRB)=row.names(try.MLE$pars$beta)
    
    matxG = cbind(obj$pars$gamma,seG,obj$pars$gamma/seG,2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matxG)=col.names
    rownames(matxG)="Association"
    matHRG = cbind(exp(obj$pars$gamma), exp(obj$pars$gamma - 1.96*seG),exp(obj$pars$gamma + 1.96*seG), 
                   2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matHRG)=col.names.HR
    rownames(matHRG)="Association"
    
    matxReg = rbind(matxB, matxG)
    matHR = rbind(matHRB, matHRG)
    
    matxT = cbind(obj$pars$theta,seT,obj$pars$theta/seT,2*(1-pnorm(abs(obj$pars$theta/seT))))
    colnames(matxT)=col.names
    
    
    matxA = cbind(obj$pars$alpha,seA,obj$pars$alpha/seA,2*(1-pnorm(abs(obj$pars$alpha/seA))))
    colnames(matxA)=col.names
    rownames(matxA)=colnames(data.frame(try.MLE$phi))
    
  }
  if(se=="Bayesian"){
    seB=sqrt(diag(obj$M2_inv))[1:p]
    seG=sqrt(diag(obj$M2_inv))[(p+1):(p+q)]
    seT=sqrt(diag(obj$M2_inv))[(p+q+1):(p+q+m)]
    seA=sqrt(diag(obj$M2_inv))[(p+q+m+1):(p+q+m+r)]
    
    
    matxB = cbind(obj$pars$beta,seB,obj$pars$beta/seB,2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matxB)=col.names
    rownames(matxB)=row.names(try.MLE$pars$beta)
    matHRB = cbind(exp(obj$pars$beta), exp(obj$pars$beta - 1.96*seB),exp(obj$pars$beta + 1.96*seB), 
                   2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matHRB)=col.names.HR
    rownames(matHRB)=row.names(try.MLE$pars$beta)
    
    matxG = cbind(obj$pars$gamma,seG,obj$pars$gamma/seG,2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matxG)=col.names
    rownames(matxG)="Association"
    matHRG = cbind(exp(obj$pars$gamma), exp(obj$pars$gamma - 1.96*seG),exp(obj$pars$gamma + 1.96*seG), 
                   2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matHRG)=col.names.HR
    rownames(matHRG)="Association"
    
    matxReg = rbind(matxB, matxG)
    matHR = rbind(matHRB, matHRG)
    
    matxT = cbind(obj$pars$theta,seT,obj$pars$theta/seT,2*(1-pnorm(abs(obj$pars$theta/seT))))
    colnames(matxT)=col.names
    
    
    matxA = cbind(obj$pars$alpha,seA,obj$pars$alpha/seA,2*(1-pnorm(abs(obj$pars$alpha/seA))))
    colnames(matxA)=col.names
    rownames(matxA)=colnames(data.frame(try.MLE$phi))
    
  }
  if(se=="Q"){
    seB=sqrt(diag(obj$cov_Q))[1:p]
    seG=sqrt(diag(obj$cov_Q))[(p+1):(p+q)]
    seT=sqrt(diag(obj$cov_Q))[(p+q+1):(p+q+m)]
    seA=sqrt(diag(obj$cov_Q))[(p+q+m+1):(p+q+m+r)]
    
    
    matxB = cbind(obj$pars$beta,seB,obj$pars$beta/seB,2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matxB)=col.names
    rownames(matxB)=row.names(try.MLE$pars$beta)
    matHRB = cbind(exp(obj$pars$beta), exp(obj$pars$beta - 1.96*seB),exp(obj$pars$beta + 1.96*seB), 
                   2*(1-pnorm(abs(obj$pars$beta/seB))))
    colnames(matHRB)=col.names.HR
    rownames(matHRB)=row.names(try.MLE$pars$beta)
    
    matxG = cbind(obj$pars$gamma,seG,obj$pars$gamma/seG,2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matxG)=col.names
    rownames(matxG)="Association"
    matHRG = cbind(exp(obj$pars$gamma), exp(obj$pars$gamma - 1.96*seG),exp(obj$pars$gamma + 1.96*seG), 
                   2*(1-pnorm(abs(obj$pars$gamma/seG))))
    colnames(matHRG)=col.names.HR
    rownames(matHRG)="Association"
    
    matxReg = rbind(matxB, matxG)
    matHR = rbind(matHRB, matHRG)
    
    matxT = cbind(obj$pars$theta,seT,obj$pars$theta/seT,2*(1-pnorm(abs(obj$pars$theta/seT))))
    colnames(matxT)=col.names
    
    
    matxA = cbind(obj$pars$alpha,seA,obj$pars$alpha/seA,2*(1-pnorm(abs(obj$pars$alpha/seA))))
    colnames(matxA)=col.names
    rownames(matxA)=colnames(data.frame(try.MLE$phi))
    
  }
  out=list(Pars=matxReg,Alpha = matxA, Theta=matxT, HR=matHR, 
           HR_display = HR_display, 
           inf=list(call=obj$mc,full=full,ploglik=obj$pll,iter=obj$iter,convergence=obj$conv_record,variance=obj$variance))
  class(out) = "summary.jm_MPL"
  out
  
}


print.jm_MPL=function(x){
  cat("\nLog-likelihood : ",x$ploglik,"\n",sep="")
  cat("\nHazard function smoothing value: ",x$variance$theta.lambda,"\n" ,sep = "")
  cat("\nLongitudinal function smoothing value: ",x$variance$alpha.lambda,"\n" ,sep = "")
  cat("\nConvergence: ",x$convergence,"\n",sep="")
  cat("\nProportional hazards regression parameters :\n")
  vectg=c(x$Pars)
  print(vectg)
  cat("\nLongitudinal model regression parameters :\n")
  print(x$Alpha)
  cat("\nBaseline hazard parameters : \n")
  print(x$Theta)
  cat("\n")
}

print.summary.jm_MPL=function(x,...){
  inf = x$inf
  print(inf$call)
  cat("\n-----\n")
  
  cat("MPL Joint Model for Time-to-Event and Longitudinal Data","\n")
  cat("Penalized log-likelihood  :  ",inf$ploglik,"\n",sep="")
  cat(ifelse(inf$convergence[1]==1,
             "Convergence : Yes, ",
             "Convergence : No, "),
      inf$iter," iterations\n",sep="")
  
  cat("\nHazard function smoothing value: ",inf$variance$theta.lambda,"\n" ,sep = "")
  cat("Longitudinal function smoothing value: ",inf$variance$alpha.lambda,"\n" ,sep = "")
  
  
  #cat("Data : ",inf$data$name,"\n",sep="")  
  #cat("No. of obs. : ",length(inf$data$time),"\n",sep="")  
  #cat("No. of events : ",sum(inf$data$censoring==1),"\n",sep="")  
  #cat("No. of right cens. : ",sum(inf$data$censoring==0),"\n",sep="") 
  #cat("No. of left cens. : ",sum(inf$data$censoring==2),"\n",sep="") 
  #cat("No. of interval cens. : ",sum(inf$data$censoring==3),"\n",sep="") 
  
  cat("\n-----\n")
  
  if(x$HR_display == TRUE){
    cat("\nHazard ratios : \n",sep="")
    printCoefmat(x$HR, P.values=TRUE, has.Pvalue=TRUE) 
    
  }else{
    cat("\nProportional hazards regression parameters : \n",sep="")
    printCoefmat(x$Pars, P.values=TRUE, has.Pvalue=TRUE) 
    
  }
  
  cat("\n-----\n")
  cat("\nLongitudinal model regression parameters : \n",sep="")
  printCoefmat(x$Alpha, P.values=TRUE, has.Pvalue=TRUE) 
  cat("\nEstimated measurement error variance: ",inf$variance$sigma2_Et,"\n" ,sep = "")
  cat("Estimated random effects variance(s): ",inf$variance$sigma2_re,"\n" ,sep = "")
  
  if(inf$full){
    cat("\n-----\n")
    cat("\nBaseline hazard parameter vector : \n",sep="")
    printCoefmat(x$Theta, P.values=TRUE, has.Pvalue=TRUE)
    
  }
}

plot.jm_MPL=function(x,se="penalised",ask=TRUE,which=1:3,range=NULL){
  which.plot=rep(TRUE,3)
  if(!is.null(which)){which.plot[-which]=FALSE}
  if(sum(which.plot)==1){ask=FALSE}
  if(ask){oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  }
  
  #control=x$control
  knots=x$knots$theta.knots
  pos=x$pars$theta
  n.x=1000
  min.t = ifelse(is.null(range), max(0, x$knots$theta.knots$bound[1]), range[1])
  max.t = ifelse(is.null(range), x$knots$theta.knots$bound[2], range[2])
  V_x_X = seq(min.t, max.t,length=n.x)
  pq=length(x$pars$beta) + length(x$pars$gamma)
  m=length(x$pars$theta)
  if(se == "penalised"){
    covar=x$cov_H_RE
  }else if(se == "BLUP"){
    covar=x$cov_H_RE_2
  }else if(se == "Bayesian"){
    covar=x$M2_inv
  }else{
    break("Unknown specification for covariance matrix.")
  }
  t.covar=covar[(pq+1):(pq+m),(pq+1):(pq+m)]
  
  anyplot=function(j,se,V_x_X,knots=x$knots$theta.knots){
    if(j==1){
      Ppsi=mSpline(V_x_X, knots = knots$int, Boundary.knots = knots$bound)
    }
    if(j>1){
      Ppsi=mSpline(V_x_X, knots = knots$int, Boundary.knots = knots$bound, integral = TRUE)
    }
    #var.Hh0=(Ppsi%*%(t.covar)%*%t(Ppsi))
    #pos.var=var.Hh0>0
    Hh0=c(Ppsi%*%matrix(x$pars$theta,ncol=1))
    V_x_X=V_x_X
    if(j < 3){
      sd.Hh0=sqrt(diag((c(1/Hh0) * Ppsi) %*% t.covar %*% t(c(1/Hh0) * Ppsi)))
      upper=log(Hh0)+1.96*sd.Hh0
      lower=log(Hh0)-1.96*sd.Hh0
      upper=exp(upper)
      lower=exp(lower)
      
    }else if(j==3){
      S0 = exp(-Hh0)
      logitS0 =log(S0/(1-S0))
      sd.S0=sqrt(diag((c(-1/(1-S0)) * Ppsi) %*% t.covar %*% t(c(-1/(1-S0)) * Ppsi)))
      upper=logitS0+1.96*sd.S0
      lower=logitS0-1.96*sd.S0
      upper=exp(upper)/(1 + exp(upper))
      lower=exp(lower)/(1 + exp(lower))
      Hh0 = S0
    }
    xlim=range(V_x_X)
    ylim=c(0,ifelse(j<3,max(upper[V_x_X<xlim[2]]),1))
    plot(1,1,xlim=xlim,ylim=ylim,main=paste("Estimate of the",c(" baseline hazard"," cumulative baseline hazard"," baseline survival")[j]," function",sep=""),
         xlab="Survival time",ylab=c(expression(h[0]*(t)),expression(H[0]*(t)),expression(S[0]*(t)))[j],type="n")
    polygon(x = c(V_x_X, rev(V_x_X)), y = c(lower, rev(upper)), 
            border = NA, col = adjustcolor("grey", alpha.f=0.5) )
    lines(V_x_X,Hh0,lwd=1.1)
    
  }
  
  if(which.plot[1]){anyplot(1,se=se,V_x_X)}
  if(which.plot[2]){anyplot(2,se=se,V_x_X)}
  if(which.plot[3]){anyplot(3,se=se,V_x_X)}
  
}


#JM predict: 
# hazard or survival function
# population mean (with CIs) or individual (no CIs)
# both - need baseline covariate values (for cox model and for longitudinal model if relevant)
# individual - need to give an i to take from the dataset to get random effects, will use baseline covariates from that


predict.jm_MPL = function(object,se="penalised",type="hazard",population=TRUE,base.cov=NULL,time_func,RE=NULL,
                           i=NULL,time=NULL,range=NULL,prob=0.95){
  
  beta=object$pars$beta
  gamma=object$pars$gamma
  theta=object$pars$theta
  alpha=object$pars$alpha[[1]]
  p=length(beta)
  q=length(gamma)
  m=length(theta)
  r=length(alpha[[1]])
  
  mc = object$mc
  mm=match(c("time.formula","long.formula","data", "id"), names(mc),0)
  mc = mc[c(1,mm)]
  form=lapply(list(mc$time.formula,mc$long.formula),as.formula)
  time_mf=model.frame(form[[1]][-2], object$data)
  time_mf = cbind(object$data$id, time_mf)
  long_mf=model.frame(form[[2]][[1]], object$data)
  long_mf = cbind(object$data$id, long_mf)
  
  if(se == "penalised"){
    covar=object$cov_H_RE
  }else if(se == "BLUP"){
    covar=object$cov_H_RE_2
  }else if(se == "Bayesian"){
    covar=object$M2_inv
  }else if(se == "Q"){
    covar=object$cov_Q
  }else{
    break("Unknown specification for covariance matrix.")
  }
  
  if(length(i)>1){warning("Only the first observation will be considered.\n",call.=FALSE)}
  if(is.null(time)){
    n.x=1000
    min.t = ifelse(is.null(range), max(0, object$knots$theta.knots$bound[1]), range[1])
    max.t = ifelse(is.null(range), object$knots$theta.knots$bound[2], range[2])
    V_x_X = seq(min.t, max.t,length=n.x)
  }else if(length(time) == 2){
    n.x=1000
    V_x_X = seq(time[1], time[2],length=n.x)
  }else{
    n.x=length(time)
    V_x_X = time
  }
  
  if(population==FALSE){
    if(!is.null(i) & is.null(RE) & is.null(base.cov)){
      long.i.ind = which(object$data$id == i)
      xT = time_mf[long.i.ind[1],-c(1)]
      W = long_mf[long.i.ind[1],-c(1:3)]
      re = object$pars$a_re_cols[[1]][i,]
    }else if(is.null(i) & !is.null(RE) & !is.null(base.cov)){
      xT = base.cov[[1]]
      W = base.cov[[2]]
      re = RE
    }else if(!is.null(i) & !is.null(base.cov)){
      long.i.ind = which(object$data$id == i)
      xT = base.cov[[1]]
      W = base.cov[[2]]
      re = object$pars$a_re_cols[[1]][i,]
    }else if(is.null(i) & is.null(RE) & !is.null(base.cov)){
      xT = base.cov[[1]]
      W = base.cov[[2]]
      re = as.matrix(0, ncol = ncol(object$pars$a_re_cols[[1]]))
    }
    
    phi = time_func(cbind(V_x_X, sapply(W, rep, length(V_x_X))))[[1]]
    phi.re = time_func(cbind(V_x_X, sapply(W, rep, length(V_x_X))))[[2]]
    z_t_est = phi %*% alpha + phi.re %*% re
    Mu=c(exp(as.matrix(xT)%*%beta))
    
    if(type=="hazard"){
      psi=mSpline(V_x_X, knots=object$knots$theta.knots$int, Boundary.knots = object$knots$theta.knots$bound_knots)
      h0=psi%*%theta
      h_t = h0 * c(Mu) * exp(c(gamma)*z_t_est)
      #no confidence intervals for an individual
      out=data.frame(time=V_x_X, z_t = z_t_est, mid=h_t, se=NA, low=NA, high=NA)
    }else{
      quad.lambda = ((V_x_X) - 0)/2
      quad.mu = ((V_x_X) + 0)/2
      quad.y = t(as.matrix(quad.lambda) %*% rules[[15]]$x + quad.mu) #one row per time, one column per i
      quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = object$knots$theta.knots$int, Boundary.knots = object$knots$theta.knots$bound))
      quad.phi.event = time_func(cbind(c(quad.y), c(sapply(W, rep, 15))))[[1]]
      quad.phi.event.re = time_func(cbind(c(quad.y), c(sapply(W, rep, 15))))[[2]]
      quad.w = rules[[15]]$w
      
      h0_t_quad = quad.psi.event %*% theta
      z_t_quad = (quad.phi.event %*% alpha) + quad.phi.event.re %*% re
      exp_zTg_quad = exp(c(gamma) * z_t_quad)
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      h0_t_star_quad = matrix(h0_t_star_quad, ncol = length(V_x_X), nrow = 15, byrow = FALSE)
      H0_t = quad.lambda * apply(matrix(rep(quad.w, length(V_x_X)), nrow = 15, byrow = FALSE) *  h0_t_star_quad, 2, sum)
      H_t = Mu * H0_t
      S_t = exp(-H_t)
      #no confidence intervals for an individual
      out=data.frame(time=V_x_X, z_t = z_t_est, mid=S_t, se=NA, low=NA, high=NA)
    }

  }else if(population == TRUE){
    if(!is.null(i) & is.null(base.cov)){
      long.i.ind = which(object$data$id == i)
      xT = time_mf[long.i.ind[1],-c(1:2)]
      W = long_mf[long.i.ind[1],-c(1:3)]
    }else if(is.null(i) & !is.null(base.cov)){
      xT = base.cov[[1]]
      W = base.cov[[2]]
    }
    
    phi = time_func(cbind(V_x_X, sapply(W, rep, length(V_x_X))))[[1]]
    z_t_est = phi %*% alpha
    #ci
    se_z_t = sqrt(diag((phi) %*% covar[(p+q+m+1):(p+q+m+r), (p+q+m+1):(p+q+m+r)] %*% t(phi)))
    z_t_est_ll = z_t_est - 1.96*se_z_t
    z_t_est_ul = z_t_est + 1.96*se_z_t
    
    Mu=c(exp(t(as.matrix(xT))%*%beta))
    
    if(type=="hazard"){
      psi=mSpline(V_x_X, knots=object$knots$theta.knots$int, Boundary.knots = object$knots$theta.knots$bound_knots)
      h0=psi%*%theta
      h_t = h0 * c(Mu) * exp(c(gamma)*z_t_est)
      
      log_se = sqrt(diag(as.matrix(cbind(xT, z_t_est, c(1/h0) * psi, c(gamma) * phi)) %*% covar[1:(p+q+m+r), 1:(p+q+m+r)] %*% t(as.matrix(cbind(xT, z_t_est, c(1/h0) * psi, c(gamma) * phi)))))
      log_ll = log(h_t) - 1.96*log_se
      log_ul = log(h_t) + 1.96*log_se
      h_t_ll = exp(log_ll)
      h_t_ul = exp(log_ul)
      
      out=data.frame(time=V_x_X, z_t = z_t_est, zt_ll = z_t_est_ll, zt_ul = z_t_est_ul, 
                     mid=h_t, low=h_t_ll, high=h_t_ul)
      
    }else{
      
      S_t_save = S_t_ll_save = S_t_ul_save = NULL
      
      for(t in 1:n.x){
        quad.lambda = ((V_x_X[t]) - 0)/2
        quad.mu = ((V_x_X[t]) + 0)/2
        quad.y = t(as.matrix(quad.lambda) %*% rules[[15]]$x + quad.mu) #one row per time, one column per i
        quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = object$knots$theta.knots$int, Boundary.knots = object$knots$theta.knots$bound))
        quad.phi.event = time_func(cbind(c(quad.y), c(sapply(W, rep, 15))))[[1]]
        quad.w = rules[[15]]$w
        
        h0_t_quad = quad.psi.event %*% theta
        z_t_quad = (quad.phi.event %*% alpha)
        exp_zTg_quad = exp(c(gamma) * z_t_quad)
        h0_t_star_quad = h0_t_quad * exp_zTg_quad
        H0_t = quad.lambda * sum(quad.w * h0_t_star_quad)
        H_t = Mu * H0_t
        S_t = exp(-H_t)
        
        A_t = quad.lambda * sum(c(quad.w) *  h0_t_quad * exp_zTg_quad * z_t_quad)
        
        Psi_t_star_quad = quad.lambda * apply(c(quad.w) * c(exp_zTg_quad) * quad.psi.event, 2, sum)
        
        D_t_star_quad = quad.lambda * apply(c(gamma) * c(quad.w) * c(exp_zTg_quad) * c(h0_t_quad) * quad.phi.event, 2, sum)
        
        logit_se = sqrt(c(-H_t/(1-S_t)*unlist(xT), -A_t*Mu/(1-S_t), -Psi_t_star_quad*Mu/(1-S_t), -D_t_star_quad*Mu/(1-S_t)) %*% covar[1:(p+q+m+r), 1:(p+q+m+r)] %*% 
               (as.matrix(c(-H_t/(1-S_t)*unlist(xT), -A_t*Mu/(1-S_t), -Psi_t_star_quad*Mu/(1-S_t), -D_t_star_quad*Mu/(1-S_t)))))
        
        logit_ll = log(S_t/(1-S_t)) - 1.96*logit_se
        logit_ul = log(S_t/(1-S_t)) + 1.96*logit_se
        
        S_t_ll = exp(logit_ll)/(1+exp(logit_ll))
        S_t_ul = exp(logit_ul)/(1+exp(logit_ul))
        
        S_t_save = c(S_t_save, S_t)
        S_t_ll_save = c(S_t_ll_save, S_t_ll)
        S_t_ul_save = c(S_t_ul_save, S_t_ul)
        
        
      }
      
      out=data.frame(time=V_x_X, z_t = z_t_est, zt_ll = z_t_est_ll, zt_ul = z_t_est_ul, 
                     mid=S_t_save, low=S_t_ll_save, high=S_t_ul_save)
      
      
    }
    
  }
  
  
  #times=c(object$data$time[,1])
  #attributes(out)$inf=list(i=i[1],user.time=!is.null(time),prob=prob,upper.value=quantile(times,prob),max=max(times),m=m,risk=(type=="hazard"))
  #colnames(out)[2]=type
  class(out)=c("predict.tvc_mpl","data.frame")
  out
}




