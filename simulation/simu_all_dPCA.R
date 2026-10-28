rm(list=ls())


library(ggplot2)
library(dplyr)
library(tidyverse)
library(MASS)
library(Matrix)


arrID <- 1 ## local running
# arrID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) ## HPC running
# source('YOUR_DIRECTORY/utility/utility_function_dPCA.R')
source('./R/github_all/utility/utility_function_dPCA.R')


N <- 1000; p <- 5; K <- 3; r <- 3
params <- generate_params(p=p, K=K, r=r, betaELE=1, alpha=0.1)
## use generate_params_beta_diff() to generate differnt beta

time <- c(1:K)
A = rep(1,K) %x% diag(1,r)
B = time %x% diag(1,r)
C = cbind(A,B)

data_all <- generate_data(N=N, p=p, K=K, r=r, params=params)
data_all %>% names

data_all_long_wide <- pivot_longer_wider(Tt=data_all$Tt, X=data_all$X,
                                         K=K, p=p, r=r)
data_all_long_wide %>% names
Tt_long <- data_all_long_wide$Tt_long
X_long <- data_all_long_wide$X_long
Tt_wide <- data_all_long_wide$Tt_wide
X_wide <- data_all_long_wide$X_wide

OneNInr <- kronecker(rep(1, N), Diagonal(K * r))
OneNInr <- as.matrix(OneNInr)


## start point
W <- params$W ## fixed due to set.seed()
beta <- params$beta
Sigb <- params$Sigb
SigG <- params$SigG
sig2E <- params$sig2E

sigma.e <- sqrt(sig2E)
sigma.g <- diag(SigG)

uid <- unique(data_all_long_wide$X_long$ID)

x <- X_long$x

niter <- 1e4; tol <- 1e-3
likelihood <- NULL; beta.est <- NULL
Lambda.sum.record <- NULL; diff.LL.record <- NULL
Sigb.est <- NULL
sigma.e.est <- NULL
sigma.g.est <- NULL

InW <- diag(1,K) %x% W
mu.t <- C %*% t(beta)
mu.x <- InW %*% mu.t
Sig.t <- A %*% Sigb %*% t(A) + diag(1, K) %x% SigG
Sig.x <- InW %*% Sig.t %*% t(InW) + diag(1, K) %x% (sig2E * diag(p))

Sig.x.list <- list(Sig.x) ## give the first list

LL <- 0*0:niter ## length=niter+1, to make LL[1]=0
## otherwise NaN for LL[1], cannot compute LL[2]-LL
beta.tol <- 0*0:niter

Lambda.sum <- 0*0:niter
diff.LL.all <- 0*0:niter


for (v in 1:niter) {
  
  
  if (v > 1 && abs(diff.LL.all[v]) < tol && abs(diff.LL.all[v])>0) {
    break
  }
  
  InW <- diag(1,K) %x% W
  
  t.exp <- NULL; t.var <- NULL; tt.exp <- NULL
  b.exp <- NULL; b.var <- NULL; bb.exp <- NULL; bt.exp <- NULL; tb.var <- NULL
  
  sum.Sigb.norm <- 0; sum.SigG.norm <- 0
  
  sum.t.prod.x <- 0; sum.sigma.e.norm <- 0
  t_minus_Ab <- NULL
  sum.tt.diag <- matrix(0, r, r)
  t_minus_Ab_old <- NULL; t.exp_old <- NULL
  
  mu.t <- C %*% t(beta)
  mu.x <- InW %*% mu.t
  
  Sig.t <- A %*% Sigb %*% t(A) + diag(1, K) %x% SigG
  Sig.x <- InW %*% Sig.t %*% t(InW) + diag(1, K) %x% (sig2E * diag(p))
  
  Sig.t.x <- Sig.t %*% t(InW) ## Cov(ti, xi)
  Sig.b.x <- Sigb %*% t(A) %*% t(InW) ## Cov(bi, xi)
  
  Sig.x.inv <- chol2inv(chol(Sig.x))
  Sig.x.det <- det(Sig.x) ## determinant
  
  Sig.g.x <- (diag(K) %x% SigG) %*% t(InW) ## Cov(gi, xi)
  for (i in uid) {
    row.i <- ((i-1)*K*p+1):(i*K*p)
    
    ####================= E-step =================####
    t.exp.i <- mu.t + Sig.t.x %*% Sig.x.inv %*% (x[row.i] - mu.x) ## E(t_i | x_i)
    t.var.i <- Sig.t - Sig.t.x %*% Sig.x.inv %*% t(Sig.t.x) ## Var(t_i | x_i)
    tt.exp.i <- t.var.i + tcrossprod(t.exp.i) ## E(t_i t_i' | x_i)
    
    b.exp.i<- Sig.b.x %*% Sig.x.inv %*% (x[row.i] - mu.x) ## E(b_i | x_i)
    b.var.i <- Sigb - Sig.b.x %*% Sig.x.inv %*% t(Sig.b.x) ## Var(b_i | x_i)
    bb.exp.i <- b.var.i + tcrossprod(b.exp.i) ## E(b_i b_i' | x_i)
    
    ## bt
    bt.exp.i <- bb.exp.i %*% t(A)
    tb.var.i <- A %*% b.var.i
    
    ## g
    g.exp.i <- Sig.g.x %*% Sig.x.inv %*% (x[row.i] - mu.x)
    g.var.i <- (diag(K) %x% SigG) - Sig.g.x %*% Sig.x.inv %*% t(Sig.g.x)
    gg.exp.i <- g.var.i + tcrossprod(g.exp.i) ## E(g_i g_i' | x_i)
    
    
    ####================= M-step =================####
    sum.t.prod.x.i <- 0 ## for each i, different
    for (k in 1:K) {
      row.k <- ((k-1)*r+1):(k*r)
      xik <- X_long[X_long$ID==i & X_long$time==k, 'x'] %>% as.matrix()
      sum.t.prod.x.i <- sum.t.prod.x.i + xik %*% t(t.exp.i[row.k])  ## p x r
      ## extract (k,k) diagonal r*r block of tt.exp.i and accumulate
      tt.block.kk <- tt.exp.i[row.k, row.k]   ## r*r block E[t_ik t_ik' | x_i]
      sum.tt.diag <- sum.tt.diag + tt.block.kk
    }
    sum.t.prod.x <- sum.t.prod.x + sum.t.prod.x.i
    
    sum.sigma.e.norm.i <- crossprod(x[row.i]) - 
      2*t(t.exp.i) %*% t(InW) %*% x[row.i] +
      t(t.exp.i) %*% crossprod(InW) %*% t.exp.i +
      sum(diag(crossprod(InW) %*% t.var.i))
    
    sum.sigma.e.norm <- sum.sigma.e.norm + sum.sigma.e.norm.i
    
    ## to compute beta
    t_minus_Ab.i <- as.data.frame(t.exp.i - A %*% b.exp.i)
    t_minus_Ab <- rbind(t_minus_Ab, t_minus_Ab.i) %>% as.matrix
    t.exp <- rbind(t.exp, t.exp.i)
    
    sum.Sigb.norm.i <- Diag(bb.exp.i)
    sum.Sigb.norm <- sum.Sigb.norm + sum.Sigb.norm.i
    
    sum.SigG.norm.i <- Diag(gg.exp.i)
    sum.SigG.norm <- sum.SigG.norm + sum.SigG.norm.i
    
  } ## end i loop
  
  #### ================= update parameters ================= ####
  W_raw <- sum.t.prod.x %*% solve(sum.tt.diag)   ## p x r
  W <- orth(W_raw, type = 'SVD')
  sigma.e <- sqrt(sum.sigma.e.norm / (N*K*p)) %>% as.numeric()
  sig2E <- sigma.e^2
  
  Sigb <- sum.Sigb.norm/N
  Sigb <- Diag(Sigb)
  
  if (r==1) {
    SigG <- diag(block_average(sum.SigG.norm, K, r)/N %>% as.matrix) %>% as.matrix()
  } else {
    SigG <- diag(block_average(sum.SigG.norm, K, r)/N %>% as.vector) %>% as.matrix()
  }
  
  sigma.g <- diag(SigG)
  
  beta_1 <- t(ginv(C) %*% average_by_stride(t_minus_Ab, N, (K*r)))[-(1:r)]
  beta_0_all <- t.exp - matrix(1, nrow=N, ncol=1) %x% B %*% beta_1
  beta_0 <- average_by_stride(
    average_by_stride(beta_0_all, N, (K*r)), K,r)
  beta <- matrix(c(beta_0, beta_1), nrow=1)
  
  
  #### likelihood, AIC, BIC
  LL.loop <- -0.5*N*log(2*pi)
  Lambda.sum.loop <- 0
  
  for (i in uid) {
    row.i <- ((i-1)*K*p+1):(i*K*p)
    LL.loop <- LL.loop - 0.5*determinant(Sig.x)$modulus - 0.5*(
      (t(x[row.i]) - t(mu.x)) %*% Sig.x.inv %*% (x[row.i] - mu.x)
    )
    
    Lambda.sum.loop <- Lambda.sum.loop +
      (t(x[row.i]) - t(mu.x)) %*% Sig.x.inv %*% (x[row.i] - mu.x)
  }
  LL[v+1] <- LL.loop
  Lambda.sum[v+1] <- Lambda.sum.loop
  
  #### for p large, use log(L[v] / L[v-1]) to check convergence
  
  #### RECORD
  Sig.x.list <- append(Sig.x.list, list(Sig.x))
  
  Lambda.sum.iter <- tibble(iter = v, Lambda.sum = Lambda.sum[v])
  Lambda.sum.record <- bind_rows(Lambda.sum.record, Lambda.sum.iter)
  
  likelihood.iter <- tibble(iter = v, LL = LL[v])
  likelihood <- bind_rows(likelihood, likelihood.iter)
  
  beta.est.iter <- tibble(iter = v, beta=beta)
  beta.est <- bind_rows(beta.est, beta.est.iter)
  
  Sigb.est.iter <- tibble(iter = v, 
                          Sigb = matrix(diag(Sigb), nrow=1, ncol=r))
  Sigb.est <- bind_rows(Sigb.est, Sigb.est.iter)
  
  sigma.e.est.iter <- tibble(iter = v, sigma.e = sigma.e)
  sigma.e.est <- bind_rows(sigma.e.est, sigma.e.est.iter)
  
  sigma.g.rec <- matrix(sigma.g, nrow=1, ncol=r)
  sigma.g.est.iter <- tibble(iter = v, sigma.g = sigma.g.rec)
  sigma.g.est <- bind_rows(sigma.g.est, sigma.g.est.iter)
  
  ## use beta.tol to record the difference
  beta.tol[v+1] <- beta.est$beta[v, (r+1)]
  
  if (v > 1) {
    diff.LL <- - 0.5*N*(
      determinant(Sig.x.list[[v+1]], logarithm = T)$modulus 
      - determinant(Sig.x.list[[v]], logarithm = T)$modulus
    ) - 0.5*(
      Lambda.sum[v+1] - Lambda.sum[v]
    )
    diff.LL.all[v+1] <- diff.LL
  }
  
  diff.LL.iter <- tibble(iter = v, diff.LL.all = diff.LL.all[v+1])
  diff.LL.record <- bind_rows(diff.LL.record, diff.LL.iter)
  
  
  cat("================================================", "\n")
  cat("iteration =", v, "\n")
  cat("sigma.e =", sigma.e, "\n")
  cat("sigma.g =", sigma.g, "\n")
  cat("beta =", beta, "\n")
  cat("Sigb =", Sigb, "\n")
  cat("LL_diff=", abs(LL[v+1] - LL[v]), "\n")
  cat("diff.LL=", diff.LL.all[v+1], "\n")
  cat("beta_diff=", beta.tol[v+1] - beta.tol[v], "\n")
  cat("likelihood is ", LL[v+1], "\n")
  cat("W[1,1] is ", W[1,1], "\n")
}

##  to show variance explained, need to store t.big in wide format (nrow=N*K, ncol=r)
t.exp.long <- data.frame(
  PublicID = rep(uid, each=(K*r)),
  visit = rep(rep(1:K, each=r), times=N),
  comp = rep(1:r, times=N*K),
  t = t.exp
)
t.exp.wide <- t.exp.long %>% 
  pivot_wider(
    id_cols = c(PublicID, visit),
    names_from = comp,
    values_from = t,
    names_prefix = "PC"
  ) %>%
  arrange(PublicID, visit)

save(t.exp.wide, file=paste0("t.exp.wide_", "ID_", arrID, ".RData"))

save(likelihood, file=paste0("likelihood_", "ID_", arrID, ".RData"))
save(beta.est, file=paste0("beta.est_", "ID_", arrID, ".RData"))
save(Sigb.est, file=paste0("Sigb.est_", "ID_", arrID, ".RData"))
save(sigma.e.est, file=paste0("sigma.e.est_", "ID_", arrID, ".RData"))
save(sigma.g.est, file=paste0("sigma.g.est_", "ID_", arrID, ".RData"))
save(W, file=paste0("W_", "ID_", arrID, ".RData"))

save(Sig.x.list, file=paste0("Sig.x.list_", "ID_", arrID, ".RData"))
save(Lambda.sum.record, file=paste0("Lambda.sum.record_", "ID_", arrID, ".RData"))
save(diff.LL.record, file=paste0("diff.LL.record_", "ID_", arrID, ".RData"))




