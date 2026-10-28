rm(list=ls())


library(dplyr)
library(tidyverse)
library(MASS)
library(Matrix)
library(lme4)


## Real data analysis with dPCA, illustrated on the example data X_example.
##
## Steps:
##   0. choose the number of components r (cumulative variance explained >= 80%)
##   1. initial values by PCA + linear mixed model (LMM) on the PC scores
##   2. fit the full model (H1) by the EM algorithm
##   3. likelihood ratio test (LRT) on the time slopes beta_1:
##        H0:   beta_1 = 0                          (df = r)
##        H0_j: beta_1j = 0 for component j = 1..r  (df = 1)
##      each null model is re-fitted by EM, started from the H1 estimates
##
## Input data format (same as X_example, see generate_example_data.R):
##   one row per subject and time point; column 'ID', column 'time' coded 1,...,K,
##   all other columns are features. Every subject is observed at all K time points.

# source('YOUR_DIRECTORY/utility/utility_function_dPCA.R')
source('./R/github_all/utility/utility_function_dPCA.R')
load('./R/github_all/data_analysis/data/X_example.RData')
out_dir <- './R/github_all/data_analysis/results_dPCA'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

niter <- 1e4; tol <- 1e-3        ## EM: stop when the increase of log-likelihood < tol
test_per_component <- TRUE       ## also test beta_1j = 0 for each component separately


## -----------------------------------------------------------------------
## Data
## -----------------------------------------------------------------------
X_wide <- X_example %>% arrange(ID, time)
features <- setdiff(colnames(X_wide), c("ID", "time"))

N <- length(unique(X_wide$ID))
p <- length(features)
K <- length(unique(X_wide$time))
stopifnot(nrow(X_wide) == N*K, all(sort(unique(X_wide$time)) == 1:K))

## x: stacked as (x_1', ..., x_N')', with x_i = (x_i1', ..., x_iK')' and x_ik of length p
x <- as.vector(t(as.matrix(X_wide[, features])))
uid <- unique(X_wide$ID)


## -----------------------------------------------------------------------
## EM algorithm for dPCA
##
## x_ik = W t_ik + e_ik,  t_ik = beta_0 + beta_1 * k + b_i + g_ik
## init      : list(W, beta, Sigb, SigG, sig2E), beta = (beta_0', beta_1')
## fix_beta1 : logical (length 1 or r), TRUE fixes the time slope of that component at 0
## The log-likelihood is evaluated at the updated parameters of each iteration.
## -----------------------------------------------------------------------
EM_dPCA <- function(x, N, p, K, r, init, fix_beta1 = FALSE,
                    niter = 1e4, tol = 1e-3, label = "H1") {

  fix_beta1 <- rep_len(fix_beta1, r)

  time <- c(1:K)
  A <- rep(1,K) %x% diag(1,r)
  B <- time %x% diag(1,r)
  C <- cbind(A,B)
  C.ginv <- ginv(C)

  W <- init$W; beta <- init$beta; Sigb <- init$Sigb; SigG <- init$SigG
  sig2E <- init$sig2E
  beta[r + which(fix_beta1)] <- 0

  X.mat <- matrix(x, nrow=K*p)   ## column i is x_i (Kp x N)
  sum.xx <- sum(X.mat^2)

  ## observed-data log-likelihood
  loglik <- function(W, beta, Sigb, SigG, sig2E) {
    InW <- diag(1,K) %x% W
    mu.x <- InW %*% C %*% t(beta)
    Sig.t <- A %*% Sigb %*% t(A) + diag(1,K) %x% SigG
    Sig.x <- InW %*% Sig.t %*% t(InW) + diag(1,K) %x% (sig2E * diag(p))
    D <- X.mat - as.vector(mu.x)
    -0.5*N*K*p*log(2*pi) - 0.5*N*determinant(Sig.x)$modulus[1] -
      0.5*sum(D * (chol2inv(chol(Sig.x)) %*% D))
  }

  LL <- loglik(W, beta, Sigb, SigG, sig2E)
  trace <- NULL
  converged <- FALSE

  for (v in 1:niter) {

    InW <- diag(1,K) %x% W
    mu.t <- C %*% t(beta)
    mu.x <- InW %*% mu.t

    Sig.t <- A %*% Sigb %*% t(A) + diag(1,K) %x% SigG
    Sig.x <- InW %*% Sig.t %*% t(InW) + diag(1,K) %x% (sig2E * diag(p))
    Sig.x.inv <- chol2inv(chol(Sig.x))

    Sig.t.x <- Sig.t %*% t(InW)                   ## Cov(t_i, x_i)
    Sig.b.x <- Sigb %*% t(A) %*% t(InW)           ## Cov(b_i, x_i)
    Sig.g.x <- (diag(K) %x% SigG) %*% t(InW)      ## Cov(g_i, x_i)

    ####================= E-step (all subjects at once) =================####
    D <- X.mat - as.vector(mu.x)                  ## x_i - mu.x, Kp x N
    t.exp <- as.vector(mu.t) + Sig.t.x %*% Sig.x.inv %*% D   ## E(t_i | x_i), Kr x N
    b.exp <- Sig.b.x %*% Sig.x.inv %*% D                     ## E(b_i | x_i), r x N
    g.exp <- Sig.g.x %*% Sig.x.inv %*% D                     ## E(g_i | x_i), Kr x N

    ## conditional variances are the same for all subjects
    t.var <- Sig.t - Sig.t.x %*% Sig.x.inv %*% t(Sig.t.x)
    b.var <- Sigb - Sig.b.x %*% Sig.x.inv %*% t(Sig.b.x)
    g.var <- (diag(K) %x% SigG) - Sig.g.x %*% Sig.x.inv %*% t(Sig.g.x)

    ####================= M-step =================####
    ## W: W_raw = SUM_TX %*% solve(SUM_TT), then orth()
    ## SUM_TX = sum_i sum_k x_ik E[t_ik|x_i]',  SUM_TT = sum_i sum_k E[t_ik t_ik'|x_i]
    sum.t.prod.x <- matrix(0, p, r)
    sum.tt.diag <- matrix(0, r, r)
    for (k in 1:K) {
      row.k <- ((k-1)*r+1):(k*r)
      row.x.k <- ((k-1)*p+1):(k*p)
      sum.t.prod.x <- sum.t.prod.x + X.mat[row.x.k, , drop=FALSE] %*% t(t.exp[row.k, , drop=FALSE])
      sum.tt.diag <- sum.tt.diag + N*t.var[row.k, row.k] + tcrossprod(t.exp[row.k, , drop=FALSE])
    }

    ## sigma.e: sum_i E[||x_i - (I_K x W) t_i||^2 | x_i], using the current W
    WtW <- crossprod(InW)
    sum.sigma.e.norm <- sum.xx - 2*sum(t.exp * (t(InW) %*% X.mat)) +
      sum(t.exp * (WtW %*% t.exp)) + N*sum(diag(WtW %*% t.var))

    W <- orth(sum.t.prod.x %*% solve(sum.tt.diag), type = 'SVD')
    sig2E <- as.numeric(sum.sigma.e.norm / (N*K*p))

    Sigb <- diag(diag(b.var) + rowSums(b.exp^2)/N, nrow=r)
    SigG <- diag(average_by_stride(diag(g.var) + rowSums(g.exp^2)/N, K, r), nrow=r)

    ## beta_1 from E[t_i - A b_i], then beta_0 given beta_1;
    ## under H0 the constrained slopes stay at 0
    beta_1 <- (C.ginv %*% rowMeans(t.exp - A %*% b.exp))[-(1:r)]
    beta_1[fix_beta1] <- 0
    beta_0 <- average_by_stride(rowMeans(t.exp) - B %*% beta_1, K, r)
    beta <- matrix(c(beta_0, beta_1), nrow=1)

    #### log-likelihood and convergence
    LL.new <- loglik(W, beta, Sigb, SigG, sig2E)
    diff.LL <- LL.new - LL
    LL <- LL.new

    trace <- bind_rows(trace, tibble(iter = v, LL = LL, diff.LL = diff.LL,
                                     beta = beta,
                                     Sigb = matrix(diag(Sigb), nrow=1),
                                     sigma.g2 = matrix(diag(SigG), nrow=1),
                                     sigma.e = sqrt(sig2E)))

    cat(label, "| iter =", v, "| LL =", LL, "| diff.LL =", diff.LL, "\n")

    if (abs(diff.LL) < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) warning(label, ": EM did not converge in ", niter, " iterations")

  ## E(t_i | x_i) in wide format (nrow = N*K, ncol = r), with the final parameters
  InW <- diag(1,K) %x% W
  mu.t <- C %*% t(beta)
  Sig.t <- A %*% Sigb %*% t(A) + diag(1,K) %x% SigG
  Sig.x <- InW %*% Sig.t %*% t(InW) + diag(1,K) %x% (sig2E * diag(p))
  t.exp <- as.vector(mu.t) + Sig.t %*% t(InW) %*% chol2inv(chol(Sig.x)) %*%
    (X.mat - as.vector(InW %*% mu.t))
  t.exp.wide <- matrix(t.exp, ncol=r, byrow=TRUE)
  colnames(t.exp.wide) <- paste0("PC", 1:r)

  list(W = W, beta = beta, Sigb = Sigb, SigG = SigG, sig2E = sig2E,
       LL = LL, iter = v, converged = converged, trace = trace,
       t.exp.wide = t.exp.wide, fix_beta1 = fix_beta1)
}


############################# stage 0: set r ###################################

pca.all <- prcomp(X_wide[, features], center = TRUE, scale. = FALSE)
var_explained <- (pca.all$sdev)^2 / sum((pca.all$sdev)^2)
cum_var_explained <- cumsum(var_explained)

## set r data-driven (or set r manually here)
r <- which(cum_var_explained >= 0.8)[1]
print(paste0("selected r = ", r))
cum_var_explained <- data.frame(
  PC = 1:min(r+2, length(cum_var_explained)),
  CumulativeVariance = cum_var_explained[1:min(r+2, length(cum_var_explained))]
)
print(cum_var_explained)


####################### stage 1: initial theta by PCA + LMM #######################

t. <- pca.all$x[, 1:r, drop=FALSE]   ## dim = (N*K) x r
W. <- pca.all$rotation[, 1:r, drop=FALSE]

## make the largest (in absolute value) element of each loading positive,
## flip the corresponding scores accordingly
W <- align_sign(W.)
t <- t. %*% diag(sign(colSums(W * W.)), nrow=r)
colnames(t) <- paste0("PC", 1:r)

t_wide <- bind_cols(X_wide[, c('ID', 'time')], as_tibble(t))

## LMM per component: PC_j ~ time + (1|ID)
pc_cols <- paste0("PC", 1:r)
lmm_results <- vector("list", r)
names(lmm_results) <- pc_cols
for (pc in pc_cols) {
  lmm_results[[pc]] <- lmer(as.formula(paste(pc, "~ time + (1|ID)")), data = t_wide)
}

## initial of sigma.e
e.matrix <- scale(as.matrix(X_wide[, features]), center = TRUE, scale = FALSE) - t %*% t(W)
sig2E <- var(as.vector(e.matrix))

init_H1 <- list(
  W = W,
  beta = t(sapply(lmm_results, fixef)) %>% as.vector() %>% matrix(nrow=1), ## (beta_0', beta_1')
  Sigb = diag(sapply(lmm_results, function(m) VarCorr(m)$ID[1]), nrow=r),    ## random intercept variance
  SigG = diag(sapply(lmm_results, function(m) attr(VarCorr(m), "sc")^2), nrow=r), ## residual variance
  sig2E = sig2E
)


####################### stage 2: EM for the full model (H1) #######################

fit_H1 <- EM_dPCA(x, N, p, K, r, init = init_H1, niter = niter, tol = tol, label = "H1")

cat("\n=== H1: LL =", fit_H1$LL, "| iterations =", fit_H1$iter, "\n")
cat("beta_0 =", fit_H1$beta[1:r], "\n")
cat("beta_1 =", fit_H1$beta[(r+1):(2*r)], "\n")


####################### stage 3: LRT on beta_1 #######################

## null models start from the H1 estimates
init_H0 <- fit_H1[c("W", "beta", "Sigb", "SigG", "sig2E")]

null_models <- list("H0: beta_1 = 0" = rep(TRUE, r))
if (test_per_component && r > 1) {
  for (j in 1:r) {
    null_models[[paste0("H0: beta_1", j, " = 0 (PC", j, ")")]] <- (1:r) == j
  }
}

fit_H0 <- lapply(names(null_models), function(h) {
  EM_dPCA(x, N, p, K, r, init = init_H0, fix_beta1 = null_models[[h]],
          niter = niter, tol = tol, label = h)
})
names(fit_H0) <- names(null_models)

## LRT statistic: -2 (LL_H0 - LL_H1) ~ chi^2(df), df = number of slopes fixed at 0
lrt_table_dPCA <- data.frame(
  Hypothesis = names(null_models),
  LL_null = sapply(fit_H0, function(f) f$LL),
  LL_full = fit_H1$LL,
  df = sapply(null_models, sum),
  row.names = NULL
)
lrt_table_dPCA$LRT_stat <- pmax(-2*(lrt_table_dPCA$LL_null - lrt_table_dPCA$LL_full), 0)
lrt_table_dPCA$p_value <- pchisq(lrt_table_dPCA$LRT_stat, df = lrt_table_dPCA$df,
                                 lower.tail = FALSE)

cat("\n============================================\n")
cat("Likelihood Ratio Tests (dPCA)\n")
cat("--------------------------------------------\n")
print(lrt_table_dPCA)
cat("============================================\n")


####################### save #######################

t.exp.wide <- bind_cols(X_wide[, c('ID', 'time')], as_tibble(fit_H1$t.exp.wide))
W <- fit_H1$W
rownames(W) <- features

save(r, cum_var_explained, file=file.path(out_dir, "r.RData"))
save(fit_H1, file=file.path(out_dir, "fit_H1_dPCA.RData"))
save(fit_H0, file=file.path(out_dir, "fit_H0_dPCA.RData"))
save(lrt_table_dPCA, file=file.path(out_dir, "lrt_table_dPCA.RData"))
save(t.exp.wide, file=file.path(out_dir, "t.exp.wide.RData"))
save(W, file=file.path(out_dir, "W.RData"))
