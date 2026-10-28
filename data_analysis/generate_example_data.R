rm(list=ls())


library(dplyr)
library(tidyverse)
library(MASS)
library(Matrix)


## Generate the example datasets X_example and Y_example used by
## data_analysis_dPCA.R (X only) and data_analysis_dPLS.R (X and Y).
## Both datasets are simulated from the dPLS model, so that X and Y share the
## joint random effect b_i and are correlated. X alone follows a dPCA model.
##
## Format of the output (same as the real data):
##   one row per subject and time point, arranged by ID and time;
##   column 'ID'   : subject identifier
##   column 'time' : time point, coded 1,...,K
##   other columns : features (X_1,...,X_p for X; Y_1,...,Y_q for Y)

# source('YOUR_DIRECTORY/utility/utility_function_dPLS.R')
source('./R/github_all/utility/utility_function_dPLS.R')
out_dir <- './R/github_all/data_analysis/data'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## -----------------------------------------------------------------------
## Data configuration
## -----------------------------------------------------------------------
N <- 100; p <- 5; q <- 4; K <- 3; r <- 2

## beta = (beta_0 for comp 1..r, beta_1 for comp 1..r)
## component 1 changes over time (beta_1 = 0.5), component 2 does not (beta_1 = 0),
## so that the LRT on the time slopes has something to detect
params <- generate_params(p=p, q=q, K=K, r=r,
                          betaELE_t=c(1, 1, 0.5, 0),
                          betaELE_u=c(1, 1, 0.5, 0))
## generate_params() resets the seed, set it after generating the parameters
set.seed(2026)

data_all <- generate_data_2ds(N, p, q, K, r, params)


## -----------------------------------------------------------------------
## Wide format: one row per (ID, time)
## -----------------------------------------------------------------------
## data_all$X is N x (K*p) with columns (x_i1, ..., x_iK), each x_ik of length p
to_wide <- function(M, n_feature, prefix) {
  wide <- matrix(as.vector(t(M)), ncol = n_feature, byrow = TRUE)
  colnames(wide) <- paste0(prefix, 1:n_feature)
  bind_cols(tibble(ID = rep(1:N, each = K), time = rep(1:K, times = N)),
            as_tibble(wide))
}

X_example <- to_wide(data_all$X, p, "X_")
Y_example <- to_wide(data_all$Y, q, "Y_")

head(X_example)
head(Y_example)

save(X_example, file=file.path(out_dir, "X_example.RData"))
save(Y_example, file=file.path(out_dir, "Y_example.RData"))
