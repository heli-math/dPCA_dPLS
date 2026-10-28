rm(list=ls())


library(dplyr)
library(tidyverse)
library(MASS)
library(Matrix)

###### PLEASE REPLACE THE DIRECTORY YOURSELF
# source('YOUR_DIRECTORY/utility/utility_function_dPLS.R')
source('./R/github_all/utility/utility_function_dPLS.R')
out_dir <- './R/github_all/data_analysis/data'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## -----------------------------------------------------------------------
## Data configuration
## -----------------------------------------------------------------------
N <- 100; p <- 5; q <- 4; K <- 3; r <- 2

params <- generate_params(p=p, q=q, K=K, r=r,
                          betaELE_t=c(1, 1, 0.5, 0),
                          betaELE_u=c(1, 1, 0.5, 0))
set.seed(2026)

data_all <- generate_data_2ds(N, p, q, K, r, params)


## -----------------------------------------------------------------------
## long format: one row per (ID, time)
## -----------------------------------------------------------------------

to_long <- function(M, n_feature, prefix) {
  wide <- matrix(as.vector(t(M)), ncol = n_feature, byrow = TRUE)
  colnames(wide) <- paste0(prefix, 1:n_feature)
  bind_cols(tibble(ID = rep(1:N, each = K), time = rep(1:K, times = N)),
            as_tibble(wide))
}

X_example <- to_long(data_all$X, p, "X_")
Y_example <- to_long(data_all$Y, q, "Y_")

head(X_example)
head(Y_example)

save(X_example, file=file.path(out_dir, "X_example.RData"))
save(Y_example, file=file.path(out_dir, "Y_example.RData"))
