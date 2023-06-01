library(nlme)
library(lme4)
library(expm)
library(Matrix)
library(CompQuadForm)
library(lmerTest)

Qbeta0_den = function(Q, ksi){
  d = dchisq(Q/ksi,df=1)/ksi
  d[which(d==Inf)]=0
  return(d)
}

binary_search_quantile = function(eigen_rho,per,start,end){
  low.val = start
  high.val = end
  low.per = 1-davies(start, eigen_rho)$Qq
  high.per = 1-davies(end, eigen_rho)$Qq
  iter = 0
  
  while(high.val-low.val>1e-10 & iter <= 500 ){
    iter = iter + 1
    
    mid.val = 0.5 *(high.val+low.val)
    mid.per = 1-davies(mid.val,eigen_rho)$Qq
    
    if (mid.per<per){
      low.val = mid.val
      high.val = high.val
    }
    if (mid.per>per){
      low.val = low.val
      high.val = mid.val
    }
  }
  ans = 0.5 *(high.val+low.val)
  return(ans)
}
numeric_pval = function(Qbeta0, ksi, pval_obs, rho, M, Q_beta0, Q_tau){
  ## find Qbeta0 < Qbeta0, (1-pval_obs) quantile
  L = length(Qbeta0)
  a = as.numeric(Qbeta0/ksi < qchisq(1-pval_obs,1))
  a = a[which(a==1)]
  
  
  b = matrix(0,length(Qbeta0),length(rho)-1) # rho neq 1
  for(t in 1:(length(rho)-1)){
    Q_rho = rho[t]*Q_beta0+ (1-rho[t])*Q_tau
    eigen_rho = c(Re(eigen(M)$values)*(1-rho[t]), ksi*rho[t])
    mixchi_quantile = binary_search_quantile(eigen_rho,1-pval_obs,0,1e6)
    b[,t]=(mixchi_quantile-rho[t]*Qbeta0)/(1-rho[t])
  }
  
  c = rep(0,L)
  d = rep(0,L)
  if (length(a)>0){
    for(j in 1:length(a)){
      c[j] = 1- davies(min(b[j,]),Re(eigen(M)$values))$Qq
      d[j]= Qbeta0_den(Qbeta0[j], ksi)
    }
  }
  return(c*d)
}

overall_effect_test = function(X, Z, study, y, G, K = NULL){
  ## test whether X is a matrix or a vector
  if (!is.null(dim(X))){
    p = ncol(X)
    Jp = rep(1,p); X1p = as.matrix(X %*% Jp)
  } else {
    X1p = as.matrix(X)
    study.dm = model.matrix(~study-1)
    X = study.dm * as.numeric(X1p)
  }
  
  ## test whether intercept is included in Z
  Z = as.matrix(Z); N = nrow(Z)
  ind = sapply(1:ncol(Z), function(i){all(rep(1,N) == Z[,i])})
  if (!any(ind)) Z = cbind(1,Z)
  
  ## study
  p = length(table(study))
  Jp = rep(1,p); 
  ni = table(study)
  
  ## burden
  dat = data.frame(X1p, Z, study, y)
  fit = lmer(y ~ Z + (1|study))
  result = summary(fit)
  gamma_hat = result$coefficients[,1]
  mu_hat = as.vector(Z %*% gamma_hat)
  #mu_hat = as.vector(fitted(fit))
  var_est = data.frame(VarCorr(fit))
  sigma_hhat = var_est$sdcor[1]
  sigma_ehat = var_est$sdcor[2]
  
  JJt = lapply(1:p, function(i){J = rep(1,ni[i]); J %*% t(J)})
  H = bdiag(JJt)
  IN = bdiag(replicate(N,1,simplify = F))
  sigma0_hat= sigma_hhat^2 * H + sigma_ehat^2 * IN
  sigma0_hat_inv = solve(sigma0_hat)
  
  Q_beta0 = as.numeric(t(y-mu_hat) %*% sigma0_hat_inv %*% X %*% Jp %*% t(Jp) %*% t(X) %*% 
                         sigma0_hat_inv %*% (y-mu_hat))
  
  P_hat = sigma0_hat_inv - sigma0_hat_inv %*% Z %*% solve(t(Z)%*%sigma0_hat_inv %*% Z) %*% 
    t(Z) %*% sigma0_hat_inv
  
  ksi = as.numeric(t(Jp) %*% t(X) %*% P_hat %*% X %*% Jp)
  
  scaled_Qbeta0= as.numeric(Q_beta0/ksi)
  
  
  ## UCSKAT
  dat = data.frame(X1p, Z, study, y)
  fit = lmer(y ~ -1 + X1p+ Z + (1|study))
  result = summary(fit)
  #print(result)
  gamma_hat = result$coefficients[,1]
  mu_hat = as.vector(cbind(X1p, Z) %*% gamma_hat)
  #mu_hat = as.vector(fitted(fit))
  var_est = data.frame(VarCorr(fit))
  sigma_hhat = var_est$sdcor[1]
  sigma_ehat = var_est$sdcor[2]
  
  JJt = lapply(1:p, function(i){J = rep(1,ni[i]); J %*% t(J)})
  H = bdiag(JJt)
  IN = bdiag(replicate(N,1,simplify = F))
  sigma_hat= sigma_hhat^2 * H + sigma_ehat^2 * IN
  sigma_hat_inv = solve(sigma_hat)
  
  XZ  = cbind(X1p,Z)
  P_hat = sigma_hat_inv - sigma_hat_inv %*% XZ %*% solve(t(XZ)%*%sigma_hat_inv %*% XZ) %*% 
    t(XZ) %*% sigma_hat_inv
  
  if (!is.numeric(G)) stop("Study specific design matrix is not numeric!")
  if (is.null(K)){K = G %*% t(G)} # default: linear kernel
  
  
  Q_tau =  as.numeric(t(y-mu_hat) %*% sigma_hat_inv %*% X %*% K %*% t(X) %*% 
                        sigma_hat_inv %*% (y-mu_hat))
  
  K_sqrt = Re(sqrtm(K))
  M = as.matrix(K_sqrt %*% t(X) %*% P_hat %*% X %*% K_sqrt)
  
  ## optimal joint test
  rho = c(0,0.1^2,0.2^2,0.3^2,0.4^2,0.5^2,0.5,1)
  pval_rho = rep(0,length(rho))
  
  ## observed p value
  for(t in 1:length(rho)){
    Q_rho = rho[t]*Q_beta0+ (1-rho[t])*Q_tau
    eigen_rho = c(Re(eigen(M)$values)*(1-rho[t]), ksi*rho[t])
    pval_rho[t] = davies(Q_rho,eigen_rho)$Qq
  }
  pval_obs = min(pval_rho)
  if (pval_obs < 0){
    pval = 0
  } else{
    pval = 1-integrate(numeric_pval,ksi = ksi, pval_obs = pval_obs, rho = rho, M = M, Q_beta0 = Q_beta0, Q_tau = Q_tau,
                       lower = 0, upper = Inf)$value
  }
  return(pval)
}
