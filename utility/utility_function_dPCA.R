
# Fun: prepare ------------------------------------------------------------

## to align the sign of loading matrix W
align_sign <- function(W) {
  for (i in 1:ncol(W)) {
    max_idx <- which.max(abs(W[, i]))
    if (W[max_idx, i] < 0) {
      W[, i] <- -W[, i]
    }
  }
  return(W)
}



block_diag <- function(W_list) {
  if (!is.list(W_list)) stop("Input must be a list of matrices")
  if (!all(sapply(W_list, is.matrix))) stop("All elements must be matrices")
  
  rows <- sum(sapply(W_list, nrow))
  cols <- sum(sapply(W_list, ncol))
  
  out <- matrix(0, nrow = rows, ncol = cols)
  
  r_start <- 1
  c_start <- 1
  for (W in W_list) {
    r_idx <- r_start:(r_start + nrow(W) - 1)
    c_idx <- c_start:(c_start + ncol(W) - 1)
    out[r_idx, c_idx] <- W
    r_start <- r_start + nrow(W)
    c_start <- c_start + ncol(W)
  }
  
  return(out)
}




Diag <- function(mat) {
  if (is.matrix(mat) && nrow(mat) == 1 && ncol(mat) == 1) {
    return(mat)
  } else {
    return(mat %>% diag %>% diag)
  }
}



## transfer vector into n*r matrix
## then calculate mean by column, obtain length=r vector
average_by_stride <- function(vector, n, r) {
  matrix_form <- matrix(vector, nrow = n, ncol = r, byrow = T)
  colMeans(matrix_form)
}


group_variances <- function(data, r) {
  groups <- split(data, (seq_along(data) - 1) %% r + 1)
  variances <- sapply(groups, var, na.rm = TRUE)
  return(variances)
}


## calculate average by \otimes
block_average <- function(matrix, K, r) {
  vec <- diag(matrix)
  average <- average_by_stride(vec, K, r)
  return(average)
}



#' Construct a block-diagonal matrix
#'
#' @param A Numerical matrix. Upper diagonal block.
#' @param B Numerical matrix. Off diagonal block.
#' @param C Numerical matrix. Lower diagonal block
#'
#' @return A block diagonal matrix of the form \eqn{\begin{bmatrix} A & B \\ B' & C \end{bmatrix}}.
#'
#' @details This function is typically called for constructing the covariance matrix of $(x,y)$.
#'
#' @export
blockm<-function(A,B,C) {
  #input: Matrices A,B,C
  #output: the block matrix
  # A    B
  #t(B)  C
  M = rbind(cbind(A,B),cbind(t(B),C))
  return(M)
}


#' Check if matrices satisfy input conditions
#'
#' @param X Should be numeric matrix.
#' @param Y Should be numeric matrix.
#' @return NULL
#' @details This function throws an error if any of the elements is \code{NA}, \code{Inf}, \code{NaN} or \code{nrow(X)} doesn't match \code{nrow(Y)}.
#' 
#' @keywords internal
#' @export
input_checker <- function(X, Y = NULL) {
  if(!is.numeric(X)) stop("Input is not numeric, but of mode ",mode(X),"\n")
  if(!is.matrix(X)) stop("Input is not a matrix, but of class ",class(X),"\n")
  if(any(is.na(X))) stop("Input contains NA's or NaN's, consider imputing with impute_matrix","\n")
  if(any(is.infinite(X))) stop("Input contains non-finite elements, consider imputing with impute_matrix","\n")
  
  if (!is.null(Y)) {
    if(!is.numeric(Y)) stop("Input is not numeric, but of mode ",mode(Y),"\n")
    if(!is.matrix(Y)) stop("Input is not a matrix, but of class ",class(Y),"\n")
    if(any(is.na(Y))) stop("Input contains NA's or NaN's, consider imputing with impute_matrix","\n")
    if(any(is.infinite(Y))) stop("Input contains non-finite elements, consider imputing with impute_matrix","\n")
    if(nrow(X) != nrow(Y)) stop("# rows don't match: ",nrow(X)," versus ",nrow(Y),"\n")
    if(!identical(rownames(X), rownames(Y))) warning("Caution: Rownames don't match!","\n")
  }
  NULL
}

#' Orthogonalize a matrix
#'
#' @param X Numeric vector or matrix.
#' @param X_true (optional) A 'true' matrix/vector. Used to correct the sign of the orthonormalized X if QR is used. Only the first column is corrected.
#' @param type A character or numeric. Should be one of "QR" or "SVD".
#' @return An orthogonalized representation of \eqn{X}
#' @details Choosing type='QR' uses a QR decomposition of X to produce orthonormal columns. For type=='SVD' it uses an SVD decomposition.
#' The columns are corrected for sign.
#' @examples
#' orth(c(3,4))
#' round(crossprod(orth(matrix(rnorm(500),100,5))),4)
#' orth(matrix(1:9,3,3),type='QR')[,1] - orth(1:3); orth(matrix(1:9,3,3),type='SVD')[,1] - orth(1:3);
#' @export
orth <- function(X, X_true = NULL, type = c("QR", "SVD")) {
  X <- as.matrix(X)
  if (!is.null(X_true)) {
    X_true <- as.matrix(X_true)
    input_checker(X,X_true)
    if(ncol(X) != ncol(X_true)) stop("# columns don't match:",ncol(X),"versus",ncol(X_true),"\n")
  }else {
    input_checker(X)
  }
  
  type <- match.arg(type)
  if (type == "SVD") {
    e <- svd(X)
    e <- tcrossprod(e$u, e$v)
  } else {
    e <- qr.Q(qr(X))
  }
  if(is.null(X_true)) {
    sign_e <- sign(crossprod(e,X)) * diag(1,ncol(e))
  } else {
    sign_e <- sign(crossprod(e,X_true)) * diag(1,ncol(e))
  }
  if(any(diag(sign_e)==0)){warning("Orthogonalization made some columns orthogonal to original columns","\n")}
  
  return(e %*% sign_e)
}


## Calculate group mean
group_mean <- function(vec, p) {
  K <- length(vec) / p
  mat <- matrix(vec, ncol = p, byrow = FALSE)
  colMeans(mat)
}

group_var <- function(vec, p) {
  K <- length(vec) / p
  mat <- matrix(vec, ncol = p, byrow = FALSE)
  apply(mat, 2, var)
}



# Fun: generate parameters -------------------------------------------

## Generate parameter values of a PO2PLS model
generate_params <- function(X, p, K, r, betaELE, alpha = 0.1) {
  
  if(length(alpha) == 1) alpha <- rep(alpha, 2)
  if(!(length(alpha) %in% c(1,2))) stop("length alpha should be 1 or 2")
  
  time = c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  # pve = c(0.9, 0.05, 0.03, 0.015, 0.005)  # must sum to 1
  # pve.true = pve[1:r]/sum(pve[1:r])
  outp <- list(
    set.seed(123), ## for W
    W = orth(matrix(rnorm(p*r), p, r)+1),
    set.seed(NULL), ## reset seed, finish W
    beta = matrix(betaELE, ncol = 2*r),
    Sigb = diag(r:1,r)
    ### make the 1st pve most
    # pve = c(0.9, 0.05, 0.03, 0.015, 0.005),  # must sum to 1
    # pve.true = pve[1:r]/sum(pve[1:r]),
    # Sigb = diag(r*pve.true,r)
  )
  
  outp$SigG = alpha[2]/(1-alpha[2]) * (
    diag(group_var(C %*% t(outp$beta), r)) + outp$Sigb)
  
  return(with(outp, {
    c(outp,
      sig2E = alpha[1]/(1-alpha[1]) * (
        mean(1/alpha[2]*SigG))
    )
  }))
}



## generate data using different beta
generate_params_beta_diff <- function(X, p, K, r, betaELE, alpha = 0.1) {
  
  if(length(alpha) == 1) alpha <- rep(alpha, 2)
  if(!(length(alpha) %in% c(1,2))) stop("length alpha should be 1 or 2")
  if(length(betaELE) == 1) betaELE <- rep(betaELE, 2*r)
  if(!(length(betaELE) %in% c(1,2*r))) stop("length alpha should be 1 or 2")
  
  time = c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  outp <- list(
    set.seed(123), ## for W
    W = orth(matrix(rnorm(p*r), p, r)+1),
    set.seed(NULL), ## reset seed, finish W
    beta = matrix(betaELE, ncol = 2*r),
    Sigb = diag(r:1,r)
  )
  
  outp$SigG = alpha[2]/(1-alpha[2]) * (
    diag(group_var(C %*% t(outp$beta), r)) + outp$Sigb)
  
  return(with(outp, {
    c(outp,
      sig2E = alpha[1]/(1-alpha[1]) * (
        mean(1/alpha[2]*SigG))
    )
  }))
}


## Generate parameters with different W's at different time points
generate_params_KW <- function(X, p, K, r, betaELE, alpha = 0.1) {
  
  if (length(alpha) == 1) alpha <- rep(alpha, 2)
  if (!(length(alpha) %in% c(1,2))) stop("length alpha should be 1 or 2")
  
  time <- 1:K
  A <- rep(1, K) %x% diag(1, r)
  B <- time %x% diag(1, r)
  C <- cbind(A, B)
  
  set.seed(123)
  W_list <- lapply(1:K, function(k) {
    orth(matrix(rnorm(p * r), p, r) + 1)
  })
  set.seed(NULL)
  
  beta <- matrix(betaELE, ncol = 2 * r)
  Sigb = diag(r:1,r)
  
  SigG <- alpha[2]/(1 - alpha[2]) * (
    diag(group_var(C %*% t(beta), r)) + Sigb
  )
  
  sig2E <- alpha[1]/(1 - alpha[1]) * (
    mean(1/alpha[2] * SigG)
  )
  
  return(list(
    W_list = W_list,
    beta = beta,
    Sigb = Sigb,
    SigG = SigG,
    sig2E = sig2E
  ))
}





# Fun: generate data and latent variable, and record ---------------------------


#' Generate data for a Dynamic Components Model
#' 
#' @param N Positive integer. Sample size to be simulated
#' @param p Positive integer. Number of features to be simulated
#' @param K positive integer. Number of time points to be simulated
#' @param r Positive integer. Number of principal components to be simulated
#' @param distr Functions. One of the random number generator functions (e.g. rnorm) for generating latent variables
#' 
#' @return A list containing the X, t, and b dataframe, in long-format
generate_data <- function(N, p, K, r, params, distr = rnorm) {
  W = params$W
  beta = params$beta
  sig2E = params$sig2E
  SigG = params$SigG
  Sigb = params$Sigb
  
  p = nrow(W)
  r = ncol(W)
  
  Gamma = diag(K) %x% W %*% cbind(diag(r*K), matrix(0, nrow=r*K, ncol=r))
  
  time <- c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  VarZ = blockm(
    A %*% Sigb %*% t(A) + diag(K) %x% SigG,
    A %*% Sigb,
    Sigb
  )
  
  Z <- MASS::mvrnorm(N, mu = rbind(C%*%t(beta), matrix(0, nrow=r, ncol=1)), Sigma = VarZ)
  
  E <- scale(matrix(distr(N*p*K), ncol = 1)) * sqrt(sig2E)
  
  dat <- Z %*% t(Gamma) + matrix(E, nrow=N)
  return(data_all = list(X = dat,
                         Tt = Z[, 1:(K*r)],
                         b = Z[, (K*r+1):(K*r+r)]))
}


## Generate data using K loading matrices W1,...,WK
## Only need to change Gamma
generate_data_KW <- function(N, p, K, r, params, distr = rnorm) {
  W_list <- params$W_list
  beta   <- params$beta
  sig2E  <- params$sig2E
  SigG   <- params$SigG
  Sigb   <- params$Sigb
  
  stopifnot(length(W_list) == K)
  
  Gamma <- block_diag(W_list) %*% cbind(diag(r*K), matrix(0, nrow=r*K, ncol=r))
  time <- c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  VarZ = blockm(
    A %*% Sigb %*% t(A) + diag(K) %x% SigG,
    A %*% Sigb,
    Sigb
  )
  
  Z <- MASS::mvrnorm(N, mu = rbind(C%*%t(beta), matrix(0, nrow=r, ncol=1)), Sigma = VarZ)
  
  E <- scale(matrix(distr(N*p*K), ncol = 1)) * sqrt(sig2E)
  
  dat <- Z %*% t(Gamma) + matrix(E, nrow=N)
  return(data_all = list(X = dat,
                         Tt = Z[, 1:(K*r)],
                         b = Z[, (K*r+1):(K*r+r)]))
}



# Fun: long and wide-format data ------------------------------------------

#' Pivot long and wide-format data
#' 
#' @param Tt A numeric matrix. The latent variable matrix of size \eqn{N \times (K \times r)}.
#' @param X A numeric matrix. The feature matrix of size \eqn{N \times (K \times p)}.
#' @param K Positive integer. Number of time points.
#' @param p Positive integer. Number of features.
#' @param r Positive integer. Number of components.
#' 
#' @return A list containing the long-format and wide-format data of \eqn{Tt} and \eqn{X}. The long-format data has four columns: ID, time, comp (or feature), and value. The wide-format data has one row for each ID and time, with columns for each component (or feature) in the format "comp_1", "comp_2", ..., "feature_1", "feature_2", ..., etc.
pivot_longer_wider <- function(Tt, X, K, p, r) {
  
  ## sub-functions
  pivot_longer_comp <- function(Tt = data$Tt, K, r) {
    ID <- rep(1:N, each = K*r)
    time <- rep(rep(1:K, each = r), times = N)
    comp <- rep(rep(1:r, times = K), times = N)
    value <- as.vector(t(Tt))
    
    Tt_long <- data.frame(ID = ID, time = time, comp = comp, t = value)
    return(Tt_long)
  }
  
  pivot_longer_data <- function(X = data$X, K, p) {
    ID <- rep(1:N, each = K*p)
    time <- rep(rep(1:K, each = p), times = N)
    feature <- rep(rep(1:p, times = K), times = N)
    value <- as.vector(t(X))
    
    X_long <- data.frame(ID = ID, time = time, feature = feature, x = value)
    return(X_long)
  }
  
  pivot_wider_comp <- function(Tt_long, r) {
    Tt_long <- Tt_long %>%
      mutate(comp = paste0("t_", comp))
    
    Tt_wide <- Tt_long %>%
      pivot_wider(names_from = comp, values_from = t) %>%
      arrange(ID, time)
    
    return(Tt_wide)
  }
  
  pivot_wider_data <- function(X_long, p) {
    X_long <- X_long %>%
      mutate(feature = paste0("x_", feature))
    
    X_wide <- X_long %>%
      pivot_wider(names_from = feature, values_from = x) %>%
      arrange(ID, time)
    
    return(X_wide)
  }
  
  ## main function
  Tt_long <- pivot_longer_comp(Tt, K = K, r = r)
  X_long <- pivot_longer_data(X, K = K, p = p)
  
  Tt_wide <- pivot_wider_comp(Tt_long, r)
  X_wide <- pivot_wider_data(X_long, p)
  
  return(list(Tt_long = Tt_long, X_long = X_long,
              Tt_wide = Tt_wide, X_wide = X_wide))
}

