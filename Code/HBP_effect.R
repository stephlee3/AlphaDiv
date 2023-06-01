library(nlme)
library(lme4)
library(expm)
library(Matrix)
library(CompQuadForm)
library(lmerTest)

HBP_test = function(X, Z, study, y, G, K = NULL){
  
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
  
  dat = data.frame(X1p, Z, study, y)
  fit = lmer(y ~ -1 + X1p+ Z + (1|study))
  result = summary(fit)
  gamma_hat = result$coefficients[,1]
  
  mu_hat = as.vector(cbind(X1p, Z) %*% gamma_hat)
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
  pval = davies(Q_tau,Re(eigen(M)$values))$Qq
  return(pval)
}
