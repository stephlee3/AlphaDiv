## library
library(tidyverse)
library(nlme)
library(lme4)
library(expm)
library(Matrix)
library(CompQuadForm)

source("./code/common_effect.R")
source("./code/HBP_effect.R")
source("./code/overall_effect.R")

nsim = 10000
level= 0.05 

# params
p_list = c(10, 30, 50)
ni_list = c(10, 30, 50)
r_list = c(3, 6, 12)
beta_list = c(0, 0.1, 0.2)
tau_list = c(0, 0.1, 0.2)
eg = expand.grid(p = p_list, r = r_list, ni = ni_list,
                 beta = beta_list, tau = tau_list,
                 stringsAsFactors = FALSE)

## You can specify job.id to change param settings
job.id = 1

# params
p = eg$p[job.id]
ni = eg$ni[job.id]
N = p *ni
r = eg$r[job.id]
q = 3
beta0 = eg$beta[job.id]
tau = eg$tau[job.id] 

# design matrix 
X = matrix(0,N,p)
Z = matrix(0,N,q)
h = rep(0,N)
eps = rep(0,N)
G = matrix(0,p,r)
alpha = rep(0,r)
gamma = rep(0,q)
sigma_h = 1
sigma_e = 1


simu_func = function(i){
  print(i)
  ## stage one 
  X = bdiag(replicate(p,rbinom(ni,1,0.7),simplify = F)) # hiv status with prob 0.7
  Z = cbind(1,rbinom(N,1,0.7),rbinom(N,1,0.4)) # intercept 1, MSM wp 0.7, gender wp 0.4
  colnames(Z)=c("Z0","MSM","gender")
  hi = rnorm(p,0,sigma_h)
  h = rep(hi, each = ni)
  eps = rnorm(N,0,sigma_e)
  study = factor(rep(c(1:p),each=ni))
  gamma = rep(1,q)
  
  if (r== 3){
    G1 = dummy(factor(sample(0:2,p,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G3 = dummy(factor(sample(0:1,p,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G1,G3)
    K = G %*% t(G)
  }
  
  if (r== 6){
    G1 = dummy(factor(sample(0:2,p,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G2 = dummy(factor(sample(0:3,p,prob = c(2,3,3,2), replace = T), levels= c(0,1,2,3)))
    G3 = dummy(factor(sample(0:1,p,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G1,G2,G3)
    K = G %*% t(G)
  }
  
  if (r== 12){
    G1 = dummy(factor(sample(0:2,p,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G2 = dummy(factor(sample(0:3,p,prob = c(2,3,3,2), replace = T), levels= c(0,1,2,3)))
    G3 = dummy(factor(sample(0:1,p,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G1,G2,G3)
    
    G1 = dummy(factor(sample(0:2,p,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G2 = dummy(factor(sample(0:3,p,prob = c(2,3,3,2), replace = T), levels= c(0,1,2,3)))
    G3 = dummy(factor(sample(0:1,p,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G,G1,G2,G3)
  }
  
  
  alpha = rnorm(r,0,1)
  beta = beta0+ G %*% alpha * tau
  y =rep(0,N)
  y = as.vector(X %*% beta + Z %*% gamma + h + eps)
  
  Jp = rep(1,p)
  X1p= as.matrix(X %*% Jp)
  colnames(X1p) = "HIV"
  
  pval = rep(0,3)
  pval[1] = try(common_effect_test(X,Z,study, y, G = G))
  pval[2] = try(HBP_test(X,Z,study, y, G = G))
  pval[3] = try(overall_effect_test(X, Z, study, y, G))
  names(pval) = c('common','HBP','overall')
  if (class(pval)!='numeric') pval = as.numeric(pval)
  print(pval)
  print("----------------")
  return(pval)
}

pval = simu_func(1)


