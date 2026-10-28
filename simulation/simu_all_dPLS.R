rm(list=ls())


library(dplyr)
library(tidyverse)
library(MASS)
library(Matrix)


arrID <- 1 ## local running
# arrID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) ## HPC running
# source('YOUR_DIRECTORY/utility/utility_function_dPLS.R')
source('./R/github_all/utility/utility_function_dPLS.R')


## -----------------------------------------------------------------------
## Data generation
## -----------------------------------------------------------------------
N <- 1000; p <- 5; q <- 4; K <- 3; r <- 3

params <- generate_params(p=p, q=q, K=K, r=r, betaELE_t=1, betaELE_u=1)
params %>% names

time <- c(1:K)
A  <- rep(1,K) %x% diag(1,r)
B  <- time %x% diag(1,r)
C  <- cbind(A, B)

data_all <- generate_data_2ds(N, p, q, K, r, params)

data_all_long_wide <- pivot_longer_wider_2ds(
  Tt=data_all$Tt, U=data_all$U,
  X=data_all$X,   Y=data_all$Y,
  K=K, p=p, q=q, r=r)

Tt_long <- data_all_long_wide$Tt_long
U_long  <- data_all_long_wide$U_long
X_long  <- data_all_long_wide$X_long
Y_long  <- data_all_long_wide$Y_long

uid <- unique(X_long$ID)
x   <- X_long$x
y   <- Y_long$y

## -----------------------------------------------------------------------
## Starting values (true params — for simulation study)
## -----------------------------------------------------------------------
W      <- params$W;      W.true <- params$W
V      <- params$V;      V.true <- params$V
beta_t <- params$beta_t; beta_u <- params$beta_u
Sigb   <- params$Sigb
Siga   <- params$Siga
Sigd   <- params$Sigd
SigG   <- params$SigG
SigH   <- params$SigH
sig2E  <- params$sig2E
sig2F  <- params$sig2F

sigma.e <- sqrt(sig2E); sigma.f <- sqrt(sig2F)
sigma.g <- diag(SigG);  sigma.h <- diag(SigH)

## -----------------------------------------------------------------------
## EM algorithm setup
## -----------------------------------------------------------------------
niter <- 1e4; tol <- 1e-3

likelihood    <- NULL
beta_t.est    <- NULL;  beta_u.est  <- NULL
Sigb.est      <- NULL;  Siga.est    <- NULL;  Sigd.est   <- NULL
sigma.e.est   <- NULL;  sigma.f.est <- NULL
sigma.g.est   <- NULL;  sigma.h.est <- NULL
Lambda.sum.record <- NULL
diff.LL.record    <- NULL

## Initial joint covariance (for list tracking)
InW <- diag(1,K) %x% W
InV <- diag(1,K) %x% V
Sig.t   <- A %*% (Sigb+Siga) %*% t(A) + diag(1,K) %x% SigG
Sig.u   <- A %*% (Sigb+Sigd) %*% t(A) + diag(1,K) %x% SigH
Sig.t.u <- A %*% Sigb %*% t(A)
Sig.x   <- InW %*% Sig.t %*% t(InW) + diag(1,K) %x% (sig2E*diag(p))
Sig.y   <- InV %*% Sig.u %*% t(InV) + diag(1,K) %x% (sig2F*diag(q))
Sig.x.y <- InW %*% Sig.t.u %*% t(InV)
Sig.xy.all <- blockm(Sig.x, Sig.x.y, Sig.y)
Sig.xy.all.list <- list(Sig.xy.all)

LL          <- 0*0:niter
Lambda.sum  <- 0*0:niter
diff.LL.all <- 0*0:niter
beta.tol    <- 0*0:niter

## -----------------------------------------------------------------------
## EM Loop
## -----------------------------------------------------------------------
for (v in 1:niter) {
  
  if (v > 1 && abs(diff.LL.all[v]) < tol && abs(diff.LL.all[v]) > 0) break
  
  
  sum.t.prod.x  <- 0;  sum.tt.diag.W  <- matrix(0, r, r)
  sum.u.prod.y  <- 0;  sum.uu.diag.V  <- matrix(0, r, r)
  sum.sigma.e.norm <- 0;  sum.sigma.f.norm <- 0
  sum.Sigb.norm    <- 0;  sum.Siga.norm    <- 0;  sum.Sigd.norm <- 0
  sum.SigG.norm    <- 0;  sum.SigH.norm    <- 0
  t_minus_Ab_Aa    <- NULL
  u_minus_Ab_Ad    <- NULL
  t.exp.joint.all  <- NULL
  u.exp.joint.all  <- NULL
  
  
  InW <- diag(1,K) %x% W
  InV <- diag(1,K) %x% V
  
  mu.t <- C %*% t(beta_t);  mu.u <- C %*% t(beta_u)
  mu.x <- InW %*% mu.t;     mu.y <- InV %*% mu.u
  
  Sig.t   <- A %*% (Sigb+Siga) %*% t(A) + diag(1,K) %x% SigG
  Sig.u   <- A %*% (Sigb+Sigd) %*% t(A) + diag(1,K) %x% SigH
  Sig.t.u <- A %*% Sigb %*% t(A)           ## Cov(t_i, u_i) = A Sigb A'
  
  Sig.x   <- InW %*% Sig.t   %*% t(InW) + diag(1,K) %x% (sig2E*diag(p))
  Sig.y   <- InV %*% Sig.u   %*% t(InV) + diag(1,K) %x% (sig2F*diag(q))
  Sig.x.y <- InW %*% Sig.t.u %*% t(InV)   ## Cov(x_i, y_i)
  
  Sig.xy.all     <- blockm(Sig.x, Sig.x.y, Sig.y)
  Sig.xy.all.inv <- chol2inv(chol(Sig.xy.all))
  Sig.x.inv      <- chol2inv(chol(Sig.x))
  Sig.y.inv      <- chol2inv(chol(Sig.y))
  
  Sig.t.x  <- Sig.t   %*% t(InW)
  Sig.t.y  <- Sig.t.u %*% t(InV)
  Sig.t.xy <- cbind(Sig.t.x, Sig.t.y)
  
  Sig.u.x  <- t(Sig.t.u) %*% t(InW)
  Sig.u.y  <- Sig.u %*% t(InV)
  Sig.u.xy <- cbind(Sig.u.x, Sig.u.y)
  
  Sig.a.x  <- Siga %*% t(A) %*% t(InW)
  Sig.d.y  <- Sigd %*% t(A) %*% t(InV)
  
  Sig.g.x  <- (diag(K) %x% SigG) %*% t(InW)
  Sig.h.y  <- (diag(K) %x% SigH) %*% t(InV)
  
  Sig.b.x  <- Sigb %*% t(A) %*% t(InW)
  Sig.b.y  <- Sigb %*% t(A) %*% t(InV)
  
  t.var.joint <- Sig.t - Sig.t.xy %*% Sig.xy.all.inv %*% t(Sig.t.xy)
  u.var.joint <- Sig.u - Sig.u.xy %*% Sig.xy.all.inv %*% t(Sig.u.xy)
  
  ## -----------------------------------------------------------------------
  ## E-step + M-step accumulators (loop over individuals)
  ## -----------------------------------------------------------------------
  for (i in uid) {
    
    row.i.x <- ((i-1)*K*p+1):(i*K*p)
    row.i.y <- ((i-1)*K*q+1):(i*K*q)
    
    dx <- x[row.i.x] - mu.x    # Kp * 1
    dy <- y[row.i.y] - mu.y    # Kq * 1
    dz <- rbind(dx, dy)        # K(p+q) * 1
    
    ## ====== E-step ======================================================
    
    t.exp.i  <- mu.t + Sig.t.xy %*% Sig.xy.all.inv %*% dz
    tt.exp.i <- t.var.joint + tcrossprod(t.exp.i)
    
    u.exp.i  <- mu.u + Sig.u.xy %*% Sig.xy.all.inv %*% dz
    uu.exp.i <- u.var.joint + tcrossprod(u.exp.i)
    
    b.exp.i <- cbind(Sig.b.x, Sig.b.y) %*% Sig.xy.all.inv %*% dz
    b.var.i <- Sigb - cbind(Sig.b.x, Sig.b.y) %*%
      Sig.xy.all.inv %*% t(cbind(Sig.b.x, Sig.b.y))
    bb.exp.i <- b.var.i + tcrossprod(b.exp.i)
    
    a.exp.i <- Sig.a.x %*% Sig.x.inv %*% dx
    a.var.i <- Siga - Sig.a.x %*% Sig.x.inv %*% t(Sig.a.x)
    aa.exp.i <- a.var.i + tcrossprod(a.exp.i)
    
    d.exp.i <- Sig.d.y %*% Sig.y.inv %*% dy
    d.var.i <- Sigd - Sig.d.y %*% Sig.y.inv %*% t(Sig.d.y)
    dd.exp.i <- d.var.i + tcrossprod(d.exp.i)
    
    g.exp.i <- Sig.g.x %*% Sig.x.inv %*% dx
    g.var.i <- (diag(K) %x% SigG) - Sig.g.x %*% Sig.x.inv %*% t(Sig.g.x)
    gg.exp.i <- g.var.i + tcrossprod(g.exp.i)
    
    h.exp.i <- Sig.h.y %*% Sig.y.inv %*% dy
    h.var.i <- (diag(K) %x% SigH) - Sig.h.y %*% Sig.y.inv %*% t(Sig.h.y)
    hh.exp.i <- h.var.i + tcrossprod(h.exp.i)
    
    ## ====== M-step accumulators =========================================
    
    sum.t.prod.x.i <- 0
    sum.u.prod.y.i <- 0
    
    for (k in 1:K) {
      row.k <- ((k-1)*r+1):(k*r)
      
      xik <- X_long[X_long$ID==i & X_long$time==k, 'x'] %>% as.matrix()
      yik <- Y_long[Y_long$ID==i & Y_long$time==k, 'y'] %>% as.matrix()
      
      ## SUM_TX: x_ik E[t_ik|x_i,y_i]'   (p * r)
      sum.t.prod.x.i <- sum.t.prod.x.i + xik %*% t(t.exp.i[row.k])
      
      ## SUM_UY: y_ik E[u_ik|x_i,y_i]'   (q * r)
      sum.u.prod.y.i <- sum.u.prod.y.i + yik %*% t(u.exp.i[row.k])
      
      ## SUM_TT_W: (k,k) diagonal block of tt.exp.i.joint  (r * r)
      sum.tt.diag.W <- sum.tt.diag.W + tt.exp.i[row.k, row.k]
      
      ## SUM_UU_V: (k,k) diagonal block of uu.exp.i.joint  (r * r)
      sum.uu.diag.V <- sum.uu.diag.V + uu.exp.i[row.k, row.k]
    }
    
    sum.t.prod.x <- sum.t.prod.x + sum.t.prod.x.i
    sum.u.prod.y <- sum.u.prod.y + sum.u.prod.y.i
    
    sum.sigma.e.norm.i <- crossprod(dx + mu.x) -
      2*t(t.exp.i) %*% t(InW) %*% x[row.i.x] +
      t(t.exp.i) %*% crossprod(InW) %*% t.exp.i +
      sum(diag(crossprod(InW) %*% t.var.joint))
    sum.sigma.e.norm <- sum.sigma.e.norm + sum.sigma.e.norm.i
    
    sum.sigma.f.norm.i <- crossprod(dy + mu.y) -
      2*t(u.exp.i) %*% t(InV) %*% y[row.i.y] +
      t(u.exp.i) %*% crossprod(InV) %*% u.exp.i +
      sum(diag(crossprod(InV) %*% u.var.joint))
    sum.sigma.f.norm <- sum.sigma.f.norm + sum.sigma.f.norm.i
    
    t_minus_Ab_Aa <- rbind(t_minus_Ab_Aa,
                           as.data.frame(t.exp.i - A %*% b.exp.i - A %*% a.exp.i))
    u_minus_Ab_Ad <- rbind(u_minus_Ab_Ad,
                           as.data.frame(u.exp.i - A %*% b.exp.i - A %*% d.exp.i))
    
    t.exp.joint.all <- rbind(t.exp.joint.all, t.exp.i)
    u.exp.joint.all <- rbind(u.exp.joint.all, u.exp.i)
    
    sum.Sigb.norm <- sum.Sigb.norm + Diag(bb.exp.i)
    sum.Siga.norm <- sum.Siga.norm + Diag(aa.exp.i)
    sum.Sigd.norm <- sum.Sigd.norm + Diag(dd.exp.i)
    
    sum.SigG.norm <- sum.SigG.norm + Diag(gg.exp.i)
    sum.SigH.norm <- sum.SigH.norm + Diag(hh.exp.i)
    
  } ## end i loop
  
  ## ======================================================================
  ## M-step: update parameters
  ## ======================================================================
  
  W_raw <- sum.t.prod.x %*% solve(sum.tt.diag.W)
  W     <- orth(W_raw)
  
  V_raw <- sum.u.prod.y %*% solve(sum.uu.diag.V)
  V     <- orth(V_raw)
  
  sigma.e <- sqrt(sum.sigma.e.norm / (N*K*p)) %>% as.numeric()
  sigma.f <- sqrt(sum.sigma.f.norm / (N*K*q)) %>% as.numeric()
  sig2E   <- sigma.e^2;  sig2F <- sigma.f^2
  
  beta_t_1 <- t(ginv(C) %*% average_by_stride(
    as.matrix(t_minus_Ab_Aa), N, (K*r)))[-(1:r)]
  beta_t_0_all <- t.exp.joint.all -
    matrix(1, N, 1) %x% B %*% beta_t_1
  beta_t_0 <- average_by_stride(
    average_by_stride(beta_t_0_all, N, (K*r)), K, r)
  beta_t <- matrix(c(beta_t_0, beta_t_1), nrow=1)
  
  beta_u_1 <- t(ginv(C) %*% average_by_stride(
    as.matrix(u_minus_Ab_Ad), N, (K*r)))[-(1:r)]
  beta_u_0_all <- u.exp.joint.all -
    matrix(1, N, 1) %x% B %*% beta_u_1
  beta_u_0 <- average_by_stride(
    average_by_stride(beta_u_0_all, N, (K*r)), K, r)
  beta_u <- matrix(c(beta_u_0, beta_u_1), nrow=1)
  
  Sigb <- Diag(sum.Sigb.norm / N)
  Siga <- Diag(sum.Siga.norm / N)
  Sigd <- Diag(sum.Sigd.norm / N)
  
  SigG <- diag(block_average(sum.SigG.norm, K, r) / N) %>% as.matrix()
  SigH <- diag(block_average(sum.SigH.norm, K, r) / N) %>% as.matrix()
  sigma.g <- diag(SigG);  sigma.h <- diag(SigH)
  
  LL.loop     <- -0.5*N*(p+q)*log(2*pi)
  Lambda.loop <- 0
  
  for (i in uid) {
    row.i.x <- ((i-1)*K*p+1):(i*K*p)
    row.i.y <- ((i-1)*K*q+1):(i*K*q)
    dz.i <- rbind(x[row.i.x]-mu.x, y[row.i.y]-mu.y)
    
    LL.loop     <- LL.loop - 0.5*determinant(Sig.xy.all)$modulus -
      0.5*(t(dz.i) %*% Sig.xy.all.inv %*% dz.i)
    Lambda.loop <- Lambda.loop + t(dz.i) %*% Sig.xy.all.inv %*% dz.i
  }
  
  LL[v+1]         <- LL.loop
  Lambda.sum[v+1] <- Lambda.loop
  
  ## ======================================================================
  ## Record and convergence
  ## ======================================================================
  Sig.xy.all.list <- append(Sig.xy.all.list, list(Sig.xy.all))
  
  if (v > 1) {
    diff.LL <- -0.5*N*(
      determinant(Sig.xy.all.list[[v+1]], logarithm=T)$modulus -
        determinant(Sig.xy.all.list[[v]],   logarithm=T)$modulus
    ) - 0.5*(Lambda.sum[v+1] - Lambda.sum[v])
    diff.LL.all[v+1] <- diff.LL
  }
  
  likelihood    <- bind_rows(likelihood,
                             tibble(iter=v, LL=LL[v]))
  beta_t.est    <- bind_rows(beta_t.est,
                             tibble(iter=v, beta_t=beta_t))
  beta_u.est    <- bind_rows(beta_u.est,
                             tibble(iter=v, beta_u=beta_u))
  Sigb.est      <- bind_rows(Sigb.est,
                             tibble(iter=v, Sigb=matrix(diag(Sigb), 1, r)))
  Siga.est      <- bind_rows(Siga.est,
                             tibble(iter=v, Siga=matrix(diag(Siga), 1, r)))
  Sigd.est      <- bind_rows(Sigd.est,
                             tibble(iter=v, Sigd=matrix(diag(Sigd), 1, r)))
  sigma.e.est   <- bind_rows(sigma.e.est, tibble(iter=v, sigma.e=sigma.e))
  sigma.f.est   <- bind_rows(sigma.f.est, tibble(iter=v, sigma.f=sigma.f))
  sigma.g.est   <- bind_rows(sigma.g.est,
                             tibble(iter=v, sigma.g=matrix(sigma.g, 1, r)))
  sigma.h.est   <- bind_rows(sigma.h.est,
                             tibble(iter=v, sigma.h=matrix(sigma.h, 1, r)))
  Lambda.sum.record <- bind_rows(Lambda.sum.record,
                                 tibble(iter=v, Lambda.sum=Lambda.sum[v]))
  diff.LL.record    <- bind_rows(diff.LL.record,
                                 tibble(iter=v, diff.LL=diff.LL.all[v+1]))
  
  beta.tol[v+1] <- beta_t.est$beta_t[v, r+1]
  
  cat("================================================\n")
  cat("iter =", v, "\n")
  cat("sigma.e =", sigma.e, "  sigma.f =", sigma.f, "\n")
  cat("beta_t  =", beta_t,  "\n")
  cat("beta_u  =", beta_u,  "\n")
  cat("Sigb    =", diag(Sigb), "\n")
  cat("LL_diff =", LL[v+1]-LL[v], "  diff.LL =", diff.LL.all[v+1], "\n")
  cat("W[1,1]  =", W[1,1], "  W.true[1,1] =", W.true[1,1], "\n")
  cat("V[1,1]  =", V[1,1], "  V.true[1,1] =", V.true[1,1], "\n")
  
} ## end EM loop

## -----------------------------------------------------------------------
## Save results
## -----------------------------------------------------------------------
save(likelihood,    file=paste0("likelihood_2ds_ID_",    arrID, ".RData"))
save(beta_t.est,    file=paste0("beta_t.est_2ds_ID_",   arrID, ".RData"))
save(beta_u.est,    file=paste0("beta_u.est_2ds_ID_",   arrID, ".RData"))
save(Sigb.est,      file=paste0("Sigb.est_2ds_ID_",     arrID, ".RData"))
save(Siga.est,      file=paste0("Siga.est_2ds_ID_",     arrID, ".RData"))
save(Sigd.est,      file=paste0("Sigd.est_2ds_ID_",     arrID, ".RData"))
save(sigma.e.est,   file=paste0("sigma.e.est_2ds_ID_",  arrID, ".RData"))
save(sigma.f.est,   file=paste0("sigma.f.est_2ds_ID_",  arrID, ".RData"))
save(sigma.g.est,   file=paste0("sigma.g.est_2ds_ID_",  arrID, ".RData"))
save(sigma.h.est,   file=paste0("sigma.h.est_2ds_ID_",  arrID, ".RData"))
save(W,             file=paste0("W_2ds_ID_",             arrID, ".RData"))
save(V,             file=paste0("V_2ds_ID_",             arrID, ".RData"))
save(diff.LL.record,file=paste0("diff.LL.record_2ds_ID_",arrID, ".RData"))
