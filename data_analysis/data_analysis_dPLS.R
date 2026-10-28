rm(list=ls())


library(dplyr)
library(tidyverse)
library(MASS)
library(Matrix)
library(lme4)
library(OmicsPLS)

##### PLEASE REPLACE THE DIRECTORY YOURSELF
# source('YOUR_DIRECTORY/utility/utility_function_dPLS.R')
source('./R/github_all/utility/utility_function_dPLS.R')
load('./R/github_all/data_analysis/data/X_example.RData')
load('./R/github_all/data_analysis/data/Y_example.RData')
out_dir <- './R/github_all/data_analysis/results_dPLS'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

niter <- 1e4; tol <- 1e-3
## also test the time slopes of each component separately
## otherwise please set as FALSE
test_per_component <- TRUE



## -----------------------------------------------------------------------
## Data
## -----------------------------------------------------------------------
X_wide <- X_example %>% arrange(ID, time)
Y_wide <- Y_example %>% arrange(ID, time)
stopifnot(identical(X_wide$ID, Y_wide$ID), identical(X_wide$time, Y_wide$time))
features_x <- setdiff(colnames(X_wide), c("ID", "time"))
features_y <- setdiff(colnames(Y_wide), c("ID", "time"))

N <- length(unique(X_wide$ID))
p <- length(features_x)
q <- length(features_y)
K <- length(unique(X_wide$time))
stopifnot(nrow(X_wide) == N*K, all(sort(unique(X_wide$time)) == 1:K))

## x: stacked as (x_1', ..., x_N')', with x_i = (x_i1', ..., x_iK')' and x_ik of length p
## y: the same, with y_ik of length q
x <- as.vector(t(as.matrix(X_wide[, features_x])))
y <- as.vector(t(as.matrix(Y_wide[, features_y])))
uid <- unique(X_wide$ID)


EM_dPLS <- function(x, y, N, p, q, K, r, init,
                    fix_beta_t1 = FALSE, fix_beta_u1 = FALSE,
                    niter = 1e4, tol = 1e-3, label = "H1") {

  fix_beta_t1 <- rep_len(fix_beta_t1, r)
  fix_beta_u1 <- rep_len(fix_beta_u1, r)

  time <- c(1:K)
  A <- rep(1,K) %x% diag(1,r)
  B <- time %x% diag(1,r)
  C <- cbind(A,B)
  C.ginv <- ginv(C)

  W <- init$W; V <- init$V
  beta_t <- init$beta_t; beta_u <- init$beta_u
  Sigb <- init$Sigb; Siga <- init$Siga; Sigd <- init$Sigd
  SigG <- init$SigG; SigH <- init$SigH
  sig2E <- init$sig2E; sig2F <- init$sig2F
  beta_t[r + which(fix_beta_t1)] <- 0
  beta_u[r + which(fix_beta_u1)] <- 0

  X.mat <- matrix(x, nrow=K*p)
  Y.mat <- matrix(y, nrow=K*q)
  sum.xx <- sum(X.mat^2); sum.yy <- sum(Y.mat^2)

  ## joint covariance of (x_i, y_i)
  Sig_xy <- function(W, V, Sigb, Siga, Sigd, SigG, SigH, sig2E, sig2F) {
    InW <- diag(1,K) %x% W
    InV <- diag(1,K) %x% V
    Sig.t   <- A %*% (Sigb+Siga) %*% t(A) + diag(1,K) %x% SigG
    Sig.u   <- A %*% (Sigb+Sigd) %*% t(A) + diag(1,K) %x% SigH
    Sig.t.u <- A %*% Sigb %*% t(A)
    blockm(InW %*% Sig.t %*% t(InW) + diag(1,K) %x% (sig2E*diag(p)),
           InW %*% Sig.t.u %*% t(InV),
           InV %*% Sig.u %*% t(InV) + diag(1,K) %x% (sig2F*diag(q)))
  }

  
  loglik <- function(W, V, beta_t, beta_u, Sigb, Siga, Sigd, SigG, SigH, sig2E, sig2F) {
    mu.x <- (diag(1,K) %x% W) %*% C %*% t(beta_t)
    mu.y <- (diag(1,K) %x% V) %*% C %*% t(beta_u)
    Sig.xy.all <- Sig_xy(W, V, Sigb, Siga, Sigd, SigG, SigH, sig2E, sig2F)
    Dz <- rbind(X.mat - as.vector(mu.x), Y.mat - as.vector(mu.y))
    -0.5*N*K*(p+q)*log(2*pi) - 0.5*N*determinant(Sig.xy.all)$modulus[1] -
      0.5*sum(Dz * (chol2inv(chol(Sig.xy.all)) %*% Dz))
  }

  LL <- loglik(W, V, beta_t, beta_u, Sigb, Siga, Sigd, SigG, SigH, sig2E, sig2F)
  trace <- NULL
  converged <- FALSE

  for (v in 1:niter) {

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

    Sig.t.xy <- cbind(Sig.t %*% t(InW), Sig.t.u %*% t(InV))         ## Cov(t_i, [x_i; y_i])
    Sig.u.xy <- cbind(t(Sig.t.u) %*% t(InW), Sig.u %*% t(InV))      ## Cov(u_i, [x_i; y_i])
    Sig.b.xy <- cbind(Sigb %*% t(A) %*% t(InW), Sigb %*% t(A) %*% t(InV)) ## Cov(b_i, [x_i; y_i])
    Sig.a.x  <- Siga %*% t(A) %*% t(InW)
    Sig.d.y  <- Sigd %*% t(A) %*% t(InV)
    Sig.g.x  <- (diag(K) %x% SigG) %*% t(InW)
    Sig.h.y  <- (diag(K) %x% SigH) %*% t(InV)

    ## ====== E-step (all subjects at once) ===============================
    Dx <- X.mat - as.vector(mu.x)
    Dy <- Y.mat - as.vector(mu.y)
    Dz <- rbind(Dx, Dy)

    t.exp <- as.vector(mu.t) + Sig.t.xy %*% Sig.xy.all.inv %*% Dz
    u.exp <- as.vector(mu.u) + Sig.u.xy %*% Sig.xy.all.inv %*% Dz
    b.exp <- Sig.b.xy %*% Sig.xy.all.inv %*% Dz
    a.exp <- Sig.a.x %*% Sig.x.inv %*% Dx
    d.exp <- Sig.d.y %*% Sig.y.inv %*% Dy
    g.exp <- Sig.g.x %*% Sig.x.inv %*% Dx
    h.exp <- Sig.h.y %*% Sig.y.inv %*% Dy

    t.var.joint <- Sig.t - Sig.t.xy %*% Sig.xy.all.inv %*% t(Sig.t.xy)
    u.var.joint <- Sig.u - Sig.u.xy %*% Sig.xy.all.inv %*% t(Sig.u.xy)
    b.var <- Sigb - Sig.b.xy %*% Sig.xy.all.inv %*% t(Sig.b.xy)
    a.var <- Siga - Sig.a.x %*% Sig.x.inv %*% t(Sig.a.x)
    d.var <- Sigd - Sig.d.y %*% Sig.y.inv %*% t(Sig.d.y)
    g.var <- (diag(K) %x% SigG) - Sig.g.x %*% Sig.x.inv %*% t(Sig.g.x)
    h.var <- (diag(K) %x% SigH) - Sig.h.y %*% Sig.y.inv %*% t(Sig.h.y)

    ## ====== M-step ======================================================
    sum.t.prod.x <- matrix(0, p, r);  sum.tt.diag.W <- matrix(0, r, r)
    sum.u.prod.y <- matrix(0, q, r);  sum.uu.diag.V <- matrix(0, r, r)
    for (k in 1:K) {
      row.k <- ((k-1)*r+1):(k*r)
      row.x.k <- ((k-1)*p+1):(k*p)
      row.y.k <- ((k-1)*q+1):(k*q)
      sum.t.prod.x  <- sum.t.prod.x + X.mat[row.x.k, , drop=FALSE] %*% t(t.exp[row.k, , drop=FALSE])
      sum.u.prod.y  <- sum.u.prod.y + Y.mat[row.y.k, , drop=FALSE] %*% t(u.exp[row.k, , drop=FALSE])
      sum.tt.diag.W <- sum.tt.diag.W + N*t.var.joint[row.k, row.k] + tcrossprod(t.exp[row.k, , drop=FALSE])
      sum.uu.diag.V <- sum.uu.diag.V + N*u.var.joint[row.k, row.k] + tcrossprod(u.exp[row.k, , drop=FALSE])
    }

    WtW <- crossprod(InW); VtV <- crossprod(InV)
    sum.sigma.e.norm <- sum.xx - 2*sum(t.exp * (t(InW) %*% X.mat)) +
      sum(t.exp * (WtW %*% t.exp)) + N*sum(diag(WtW %*% t.var.joint))
    sum.sigma.f.norm <- sum.yy - 2*sum(u.exp * (t(InV) %*% Y.mat)) +
      sum(u.exp * (VtV %*% u.exp)) + N*sum(diag(VtV %*% u.var.joint))

    W <- orth(sum.t.prod.x %*% solve(sum.tt.diag.W))
    V <- orth(sum.u.prod.y %*% solve(sum.uu.diag.V))
    sig2E <- as.numeric(sum.sigma.e.norm / (N*K*p))
    sig2F <- as.numeric(sum.sigma.f.norm / (N*K*q))

    beta_t_1 <- (C.ginv %*% rowMeans(t.exp - A %*% b.exp - A %*% a.exp))[-(1:r)]
    beta_t_1[fix_beta_t1] <- 0
    beta_t_0 <- average_by_stride(rowMeans(t.exp) - B %*% beta_t_1, K, r)
    beta_t <- matrix(c(beta_t_0, beta_t_1), nrow=1)

    beta_u_1 <- (C.ginv %*% rowMeans(u.exp - A %*% b.exp - A %*% d.exp))[-(1:r)]
    beta_u_1[fix_beta_u1] <- 0
    beta_u_0 <- average_by_stride(rowMeans(u.exp) - B %*% beta_u_1, K, r)
    beta_u <- matrix(c(beta_u_0, beta_u_1), nrow=1)

    Sigb <- diag(diag(b.var) + rowSums(b.exp^2)/N, nrow=r)
    Siga <- diag(diag(a.var) + rowSums(a.exp^2)/N, nrow=r)
    Sigd <- diag(diag(d.var) + rowSums(d.exp^2)/N, nrow=r)
    SigG <- diag(average_by_stride(diag(g.var) + rowSums(g.exp^2)/N, K, r), nrow=r)
    SigH <- diag(average_by_stride(diag(h.var) + rowSums(h.exp^2)/N, K, r), nrow=r)

    ## ====== log-likelihood and convergence ==============================
    LL.new <- loglik(W, V, beta_t, beta_u, Sigb, Siga, Sigd, SigG, SigH, sig2E, sig2F)
    diff.LL <- LL.new - LL
    LL <- LL.new

    trace <- bind_rows(trace, tibble(iter = v, LL = LL, diff.LL = diff.LL,
                                     beta_t = beta_t, beta_u = beta_u,
                                     Sigb = matrix(diag(Sigb), nrow=1),
                                     Siga = matrix(diag(Siga), nrow=1),
                                     Sigd = matrix(diag(Sigd), nrow=1),
                                     sigma.g2 = matrix(diag(SigG), nrow=1),
                                     sigma.h2 = matrix(diag(SigH), nrow=1),
                                     sigma.e = sqrt(sig2E), sigma.f = sqrt(sig2F)))

    cat(label, "| iter =", v, "| LL =", LL, "| diff.LL =", diff.LL, "\n")

    if (abs(diff.LL) < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) warning(label, ": EM did not converge in ", niter, " iterations")

  InW <- diag(1,K) %x% W
  InV <- diag(1,K) %x% V
  mu.t <- C %*% t(beta_t);  mu.u <- C %*% t(beta_u)
  Sig.t   <- A %*% (Sigb+Siga) %*% t(A) + diag(1,K) %x% SigG
  Sig.u   <- A %*% (Sigb+Sigd) %*% t(A) + diag(1,K) %x% SigH
  Sig.t.u <- A %*% Sigb %*% t(A)
  Sig.xy.all.inv <- chol2inv(chol(Sig_xy(W, V, Sigb, Siga, Sigd, SigG, SigH, sig2E, sig2F)))
  Dz <- rbind(X.mat - as.vector(InW %*% mu.t), Y.mat - as.vector(InV %*% mu.u))
  t.exp <- as.vector(mu.t) +
    cbind(Sig.t %*% t(InW), Sig.t.u %*% t(InV)) %*% Sig.xy.all.inv %*% Dz
  u.exp <- as.vector(mu.u) +
    cbind(t(Sig.t.u) %*% t(InW), Sig.u %*% t(InV)) %*% Sig.xy.all.inv %*% Dz
  t.exp.wide <- matrix(t.exp, ncol=r, byrow=TRUE)
  u.exp.wide <- matrix(u.exp, ncol=r, byrow=TRUE)
  colnames(t.exp.wide) <- colnames(u.exp.wide) <- paste0("PC", 1:r)

  list(W = W, V = V, beta_t = beta_t, beta_u = beta_u,
       Sigb = Sigb, Siga = Siga, Sigd = Sigd, SigG = SigG, SigH = SigH,
       sig2E = sig2E, sig2F = sig2F,
       LL = LL, iter = v, converged = converged, trace = trace,
       t.exp.wide = t.exp.wide, u.exp.wide = u.exp.wide,
       fix_beta_t1 = fix_beta_t1, fix_beta_u1 = fix_beta_u1)
}


############################# stage 0: set r ###################################

pca_x <- prcomp(X_wide[, features_x], center = TRUE, scale. = TRUE)
pca_y <- prcomp(Y_wide[, features_y], center = TRUE, scale. = TRUE)
cum_var_explained_x <- cumsum((pca_x$sdev)^2 / sum((pca_x$sdev)^2))
cum_var_explained_y <- cumsum((pca_y$sdev)^2 / sum((pca_y$sdev)^2))

## define r = max{r1, r2} (or set r manually here)
r1 <- which(cum_var_explained_x >= 0.8)[1]
r2 <- which(cum_var_explained_y >= 0.8)[1]
r <- max(r1, r2)
print(paste0("selected r = ", r))
n_show <- min(r+2, p, q)
cum_var_explained_x_y <- data.frame(
  PC = 1:n_show,
  CumulativeVarianceX = cum_var_explained_x[1:n_show],
  CumulativeVarianceY = cum_var_explained_y[1:n_show]
)
print(cum_var_explained_x_y)


####################### stage 1: initial theta by O2PLS + LMM #######################

## O2PLS on the two stacked matrices, joint components only
fit.o2m <- o2m(X=as.matrix(X_wide[, features_x]), Y=as.matrix(Y_wide[, features_y]),
               n = r, nx = 0, ny = 0)
summary(fit.o2m)
t. <- as.matrix(fit.o2m$Tt)
u. <- as.matrix(fit.o2m$U)
W. <- as.matrix(fit.o2m$W.)
V. <- as.matrix(fit.o2m$C.)

W <- align_sign(W.); V <- align_sign(V.)
t <- t. %*% diag(sign(colSums(W * W.)), nrow=r)
u <- u. %*% diag(sign(colSums(V * V.)), nrow=r)
colnames(t) <- colnames(u) <- paste0("PC", 1:r)

t.wide <- bind_cols(X_wide[, c('ID', 'time')], as_tibble(t))
u.wide <- bind_cols(Y_wide[, c('ID', 'time')], as_tibble(u))

## LMM per component: PC_j ~ time + (1|ID), on t and on u
pc_cols <- paste0("PC", 1:r)
lmm_results_t <- vector("list", r); names(lmm_results_t) <- pc_cols
lmm_results_u <- vector("list", r); names(lmm_results_u) <- pc_cols
for (pc in pc_cols) {
  fmla <- as.formula(paste(pc, "~ time + (1|ID)"))
  lmm_results_t[[pc]] <- lmer(fmla, data = t.wide)
  lmm_results_u[[pc]] <- lmer(fmla, data = u.wide)
}

## initial of sigma.e and sigma.f
e.matrix <- scale(as.matrix(X_wide[, features_x]), center = TRUE, scale = FALSE) - t %*% t(W)
f.matrix <- scale(as.matrix(Y_wide[, features_y]), center = TRUE, scale = FALSE) - u %*% t(V)

Sigb <- diag(rep(10, r), nrow=r)
init_H1 <- list(
  W = W, V = V,
  beta_t = t(sapply(lmm_results_t, fixef)) %>% as.vector() %>% matrix(nrow=1), ## (beta_t0', beta_t1')
  beta_u = t(sapply(lmm_results_u, fixef)) %>% as.vector() %>% matrix(nrow=1), ## (beta_u0', beta_u1')
  Sigb = Sigb, Siga = Sigb, Sigd = 5*Sigb,
  SigG = diag(sapply(lmm_results_t, function(m) attr(VarCorr(m), "sc")^2), nrow=r), ## residual variance
  SigH = diag(sapply(lmm_results_u, function(m) attr(VarCorr(m), "sc")^2), nrow=r),
  sig2E = var(as.vector(e.matrix)),
  sig2F = var(as.vector(f.matrix))
)


####################### stage 2: EM for the full model (H1) #######################

fit_H1 <- EM_dPLS(x, y, N, p, q, K, r, init = init_H1, niter = niter, tol = tol, label = "H1")

cat("\n=== H1: LL =", fit_H1$LL, "| iterations =", fit_H1$iter, "\n")
cat("beta_t0 =", fit_H1$beta_t[1:r], "  beta_t1 =", fit_H1$beta_t[(r+1):(2*r)], "\n")
cat("beta_u0 =", fit_H1$beta_u[1:r], "  beta_u1 =", fit_H1$beta_u[(r+1):(2*r)], "\n")


####################### stage 3: LRT on beta_t1 and beta_u1 #######################

## null models start from the H1 estimates
init_H0 <- fit_H1[c("W", "V", "beta_t", "beta_u", "Sigb", "Siga", "Sigd",
                    "SigG", "SigH", "sig2E", "sig2F")]

## each null model: which slopes are fixed at 0 (fix_beta_t1, fix_beta_u1)
null_models <- list(
  "H0_t: beta_t1 = 0"  = list(t = rep(TRUE, r),  u = rep(FALSE, r)),
  "H0_u: beta_u1 = 0"  = list(t = rep(FALSE, r), u = rep(TRUE, r)),
  "H0_tu: both = 0"    = list(t = rep(TRUE, r),  u = rep(TRUE, r))
)
if (test_per_component && r > 1) {
  for (j in 1:r) {
    null_models[[paste0("H0_t: beta_t1", j, " = 0 (PC", j, ")")]] <-
      list(t = (1:r) == j, u = rep(FALSE, r))
    null_models[[paste0("H0_u: beta_u1", j, " = 0 (PC", j, ")")]] <-
      list(t = rep(FALSE, r), u = (1:r) == j)
  }
}

fit_H0 <- lapply(names(null_models), function(h) {
  EM_dPLS(x, y, N, p, q, K, r, init = init_H0,
          fix_beta_t1 = null_models[[h]]$t, fix_beta_u1 = null_models[[h]]$u,
          niter = niter, tol = tol, label = h)
})
names(fit_H0) <- names(null_models)

## LRT statistic: -2 (LL_H0 - LL_H1) ~ chi^2(df), df = number of slopes fixed at 0
lrt_table_dPLS <- data.frame(
  Hypothesis = names(null_models),
  LL_null = sapply(fit_H0, function(f) f$LL),
  LL_full = fit_H1$LL,
  df = sapply(null_models, function(h) sum(h$t) + sum(h$u)),
  row.names = NULL
)
lrt_table_dPLS$LRT_stat <- pmax(-2*(lrt_table_dPLS$LL_null - lrt_table_dPLS$LL_full), 0)
lrt_table_dPLS$p_value <- pchisq(lrt_table_dPLS$LRT_stat, df = lrt_table_dPLS$df,
                                 lower.tail = FALSE)

cat("\n============================================\n")
cat("Likelihood Ratio Tests (dPLS)\n")
cat("--------------------------------------------\n")
print(lrt_table_dPLS)
cat("============================================\n")


####################### save #######################

t.exp.wide <- bind_cols(X_wide[, c('ID', 'time')], as_tibble(fit_H1$t.exp.wide))
u.exp.wide <- bind_cols(Y_wide[, c('ID', 'time')], as_tibble(fit_H1$u.exp.wide))
W <- fit_H1$W; rownames(W) <- features_x
V <- fit_H1$V; rownames(V) <- features_y

save(r, cum_var_explained_x_y, file=file.path(out_dir, "r.RData"))
save(fit_H1, file=file.path(out_dir, "fit_H1_dPLS.RData"))
save(fit_H0, file=file.path(out_dir, "fit_H0_dPLS.RData"))
save(lrt_table_dPLS, file=file.path(out_dir, "lrt_table_dPLS.RData"))
save(t.exp.wide, u.exp.wide, file=file.path(out_dir, "t_u.exp.wide.RData"))
save(W, V, file=file.path(out_dir, "W_V.RData"))
