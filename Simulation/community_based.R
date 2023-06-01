## library
library(tidyverse)
library(MiSPU)
library(dirmult)
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


## params
m_list = c(10,30,50)
ni_list = c(10,30,50)
r_list = c(3,6,12)
alpha_div_list = c("observed species","shannon")
theta_list = c(-0.1, -0.3, -0.5) 
phi_list = c(0.1, 0.3, 0.5) 

eg = expand.grid( m = m_list,r = r_list, ni = ni_list,
                  alpha_div_method = alpha_div_list, phi = phi_list,theta = theta_list,
                  stringsAsFactors = FALSE) 

## You can specify job.id to change param settings
job.id = 300

# dimensions
m = eg$m[job.id]
ni = eg$ni[job.id]
N = m *ni
r = eg$r[job.id]
alpha_div_method = eg$alpha_div_method[job.id]
phi = eg$phi[job.id]
theta = eg$theta[job.id]

## simulate OTU count based on throat data
simuOTU <- function(nSam=100, s=12,ncluster = 20,mu = 1000, size = 25) {
  
  data(throat.tree,envir = environment())
  data(dd,envir = environment())
  
  tree <- get("throat.tree", envir  = environment())
  dd1 = get("dd",envir  = environment())
  tree.dist <- cophenetic(tree)
  obj <- pam(as.dist(tree.dist), ncluster,diss =TRUE)
  clustering <- obj$clustering
  otu.ids <- tree$tip.label
  
  p.est = dd1$pi
  names(p.est) <- names(dd1$pi)
  theta <- dd1$theta
  gplus <- (1 - theta) / theta
  p.est <- p.est[otu.ids]
  g.est <- p.est * gplus
  p.clus <- sort(tapply(p.est, clustering, sum), decreasing=T)
  scale2 = function(x)as.numeric(scale(x))
  
  
  comm <- matrix(0, nSam, length(g.est))
  rownames(comm) <- 1:nrow(comm)
  colnames(comm) <- names(g.est)
  # comm.p hold the underlying proportions
  comm.p <- comm
  nSeq <- rnbinom(nSam, mu = mu, size = size)
  for (i in 1:nSam) {
    comm.p[i, ] <- MiSPU::rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
  }
  
  otu.ids <- names(which(clustering == s))
  
  # No additional covariates in this case.
  OTU = comm[, otu.ids]
  return(list(informative.OTU = OTU,whole.OTU = comm, prob = comm.p, informative.prob = comm.p[,otu.ids]))
}


get_alpha_div = function(OTU.tab, method = c("observed species","shannon","simpson")){
  OTU.prop = t(apply(OTU.tab,1, function(x){x/sum(x)}))
  if (method[1] == "observed species")
    ans = rowSums(OTU.tab > 0)
  if (method[1] == "shannon")
    ans = apply(OTU.prop,1, function(x){sum(-log(x)*x, na.rm = T)})
  if (method[1] == "simpson")
    ans = apply(OTU.prop,1, function(x){sum(-x^2)})
  return(ans)
}

simu_alpha_div = function(pseudocount = 0){
  
  OTU.simu = simuOTU(nSam = N, mu = 1e5, size = 2500)
  OTU.tab = OTU.simu$whole.OTU + pseudocount
  p = ncol(OTU.tab)
  rc = rowSums(OTU.tab)
  true.OTU.prop= t(apply(OTU.tab,1, function(x){x/sum(x)}))
  sampling.OTU.prob = OTU.simu$prob
  
  ## study
  study = factor(rep(c(1:m),each=ni))
  
  ## alpha diversity
  true_alpha_div = get_alpha_div(OTU.tab,method = alpha_div_method)
  true_y = scale(true_alpha_div)
  
  ## simulate HIV status and adjusting covariates
  Z = cbind(1,rbinom(N,1,0.7),rbinom(N,1,0.4))
  
  beta = rnorm(m,theta,phi)
  x = sapply(1:N, function(i){
    linear_term = true_y[i] * beta[study[i]]
    pr = exp(linear_term)/(1+exp(linear_term))
    rbinom(1,1,prob = pr)
  }) 
  
  if (r== 3){
    G1 = dummy(factor(sample(0:2,m,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G3 = dummy(factor(sample(0:1,m,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G1,G3)
  }
  
  if (r== 6){
    G1 = dummy(factor(sample(0:2,m,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G2 = dummy(factor(sample(0:3,m,prob = c(2,3,3,2), replace = T), levels= c(0,1,2,3)))
    G3 = dummy(factor(sample(0:1,m,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G1,G2,G3)
  }
  
  if (r== 12){
    G1 = dummy(factor(sample(0:2,m,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G2 = dummy(factor(sample(0:3,m,prob = c(2,3,3,2), replace = T), levels= c(0,1,2,3)))
    G3 = dummy(factor(sample(0:1,m,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G1,G2,G3)
    
    G1 = dummy(factor(sample(0:2,m,prob = c(3,3,4), replace = T), levels= c(0,1,2)))
    G2 = dummy(factor(sample(0:3,m,prob = c(2,3,3,2), replace = T), levels= c(0,1,2,3)))
    G3 = dummy(factor(sample(0:1,m,prob = c(4,6), replace = T), levels= c(0,1)))
    G = cbind(G,G1,G2,G3)
  }
  
  ## quad kernel
  bk = rbinom(p, 1, 0.5)
  OTU.prop.bias = t(sapply(1:N, function(i){
    Gi = G[study[i],]
    logP = log(sampling.OTU.prob[i,]) + (1+sum(Gi))^2 * bk * eta
    Pstar = exp(logP); Pstar = Pstar/sum(Pstar)
  }))
  
  SeqDepth = runif(m,1000,5000)
  nSeq <- sapply(1:N, function(i){
    nSeq <- rnbinom(1, mu = SeqDepth[study[i]], size = 25)
  })
  
  obs.OTU.tab = t(sapply(1:N,function(i){
    rmultinom(1, nSeq[i], prob=OTU.prop.bias[i,])[, 1]
  }))
  
  obs_alpha_div = get_alpha_div(obs.OTU.tab, method = alpha_div_method)
  y = scale(obs_alpha_div)
  
  dat = data.frame(x,Z,study,y)
  ans = list(dat = dat, G= G)
  return(ans)
}


simu_alpha_div_test = function(i){
  print(i)
  dat0 = simu_alpha_div()
  dat = dat0$dat
  X = dat$x
  Z = dat[,c("X1","X2","X3")]
  study = dat$study
  y = dat$y
  
  G = dat0$G
  K = (1 + G %*% t(G))^2

  pval = rep(0,3)
  pval[1] = try(common_effect_test(X, Z, study, y, G, K = K))
  pval[2] = try(HBP_test(X, Z, study, y, G, K = K))
  pval[3] = try(overall_effect_test(X, Z, study, y, G, K = K))
  names(pval) = c('common','HBP','overall')
  
  if (class(pval)!='numeric') pval = as.numeric(pval)
  print(pval)
  print("-------------")
  return(pval)
}

pval = simu_alpha_div_test(1)



