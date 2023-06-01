library(nlme)
library(lme4)
library(expm)
library(Matrix)
library(CompQuadForm)
library(lmerTest)

common_effect_test = function(X, Z, study, y, G, K = NULL){
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
  
  if (is.null(K)){K = G %*% t(G)} # default: linear kernel
  eK = eigen(K,symmetric = T)
  eK.value = eK$values
  eK.vector = eK$vectors[,which(eK.value>1e-10)]
  Lambda_sqrt = as.matrix(diag(eK.value[which(eK.value>1e-10)]^(1/2)))
  XQlambda = as.matrix( X %*% eK.vector %*% Lambda_sqrt)
  
  group = rep(1,N)
  fit1.lme = try(lme(
    y~ -1 + Z,
    random = list(
      group= pdIdent(~-1+XQlambda),
      study= ~1
    )
  ))
  
  if (class(fit1.lme) != "try-error"){
    var.lme = VarCorr(fit1.lme)
    sigma_elme = as.numeric(var.lme["Residual",2])
    sigma_hlme = as.numeric(var.lme["(Intercept)",2])
    tau_lme = as.numeric(var.lme["XQlambda1",2])
    
    gamma_lme = as.matrix(fit1.lme$coefficients$fixed)
    mu_lme = as.vector(Z %*% gamma_lme)
    
    JJt = lapply(1:p, function(i){J = rep(1,ni[i]); J %*% t(J)})
    H = bdiag(JJt)
    IN = bdiag(replicate(N,1,simplify = F))
    sigma0_lme= sigma_hlme^2 * H + sigma_elme^2 * IN+ tau_lme^2 * X%*% K %*% t(X)
    sigma0_lme_inv = solve(sigma0_lme)

    
    Q_beta0_lme = as.numeric(t(y-mu_lme) %*% sigma0_lme_inv %*% X %*% Jp %*% t(Jp) %*% t(X) %*% 
                               sigma0_lme_inv %*% (y-mu_lme))
    
    P_hat_lme = sigma0_lme_inv - sigma0_lme_inv %*% Z %*% solve(t(Z)%*%sigma0_lme_inv %*% Z) %*% 
      t(Z) %*% sigma0_lme_inv
    
    ksi_lme = as.numeric(t(Jp) %*% t(X) %*% P_hat_lme %*% X %*% Jp)
    
    scaled_Qbeta0_lme= as.numeric(Q_beta0_lme/ksi_lme)
    pval = 1-pchisq(scaled_Qbeta0_lme,df=1)
  }
  else{
    pval = NA
  }
  return(pval)
}