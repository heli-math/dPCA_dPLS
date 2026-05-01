
# Fun: prepare ------------------------------------------------------------

align_sign <- function(w) {
  for (i in 1:ncol(w)) {
    max_idx <- which.max(abs(w[, i]))
    if (w[max_idx, i] < 0) {
      w[, i] <- -w[, i]
    }
  }
  return(w)
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


average_by_stride <- function(vector, n, r) {
  matrix_form <- matrix(vector, nrow = n, ncol = r, byrow = T)
  colMeans(matrix_form)
}


group_variances <- function(data, r) {
  groups <- split(data, (seq_along(data) - 1) %% r + 1)
  variances <- sapply(groups, var, na.rm = TRUE)
  return(variances)
}



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


## For 2-datasets: Construct a block-diagonal matrix with a list of matrices
blockmm <- function(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15) {
  #input: Matrices A1~A15
  #output: the block matrix
  # A1      A2      A3      A4      A5
  #t(A2)    A6      A7      A8      A9
  #t(A3)    t(A7)   A10     A11     A12
  #t(A4)    t(A8)   t(A11)  A13     A14
  #t(A5)    t(A9)   t(A12)  t(A14)  A15
  # Construct rows of the block matrix
  row1 <- cbind( A1,   A2,   A3,   A4,   A5)
  row2 <- cbind(t(A2), A6,   A7,   A8,   A9)
  row3 <- cbind(t(A3), t(A7), A10,  A11,  A12)
  row4 <- cbind(t(A4), t(A8), t(A11), A13, A14)
  row5 <- cbind(t(A5), t(A9), t(A12), t(A14), A15)
  
  M <- rbind(row1, row2, row3, row4, row5)
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

generate_params <- function(X, Y, p, q, K, r, betaELE_t, betaELE_u, alpha = 0.1) {
  if(length(alpha) == 1) alpha <- rep(alpha, 4)
  if(!(length(alpha) %in% c(1,4))) stop("length alpha should be 1 or 4")
  
  time = c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  outp <- list(
    set.seed(123),
    W = orth(matrix(rnorm(p*r), p, r)+1),
    V = orth(matrix(rnorm(q*r), q, r)+1),
    set.seed(NULL),
    beta_t = matrix(betaELE_t, ncol = 2*r),
    beta_u = matrix(betaELE_u, ncol = 2*r),
    Sigb = diag(r:1,r),
    Siga = diag(r:1,r)/2,
    Sigd = diag(r:1, r)/3
  )
  
  outp$SigG = alpha[3]/(1-alpha[3]) * (
    diag(group_var(C %*% t(outp$beta_t), r)) + outp$Sigb + outp$Siga)
  outp$SigH = alpha[4]/(1-alpha[4]) * (
    diag(group_var(C %*% t(outp$beta_u), r)) + outp$Sigb + outp$Sigd)
  
  return(with(outp, {
    c(outp,
      sig2E = alpha[1]/(1-alpha[1]) * (
        mean(1/alpha[3]*SigG))
      ,
      sig2F = alpha[2]/(1-alpha[2]) * (
        mean(1/alpha[4]*SigH))
    )
  }))
}




generate_params_KW <- function(X, Y, p, q, K, r, betaELE_t, betaELE_u, alpha = 0.1) {
  if(length(alpha) == 1) alpha <- rep(alpha, 4)
  if(!(length(alpha) %in% c(1,4))) stop("length alpha should be 1 or 4")
  
  time = c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  set.seed(123)
  W_list <- lapply(1:K, function(k) {
    orth(matrix(rnorm(p * r), p, r) + 1)
  })
  V_list <- lapply(1:K, function(k) {
    orth(matrix(rnorm(q * r), q, r) + 1)
  })
  set.seed(NULL)
  
  beta_t = matrix(betaELE_t, ncol = 2*r)
  beta_u = matrix(betaELE_u, ncol = 2*r)
  Sigb = diag(r:1,r)
  Siga = diag(r:1,r)/2
  Sigd = diag(r:1, r)/3
  
  SigG = alpha[3]/(1-alpha[3]) * (
    diag(group_var(C %*% t(beta_t), r)) + Sigb + Siga)
  SigH = alpha[4]/(1-alpha[4]) * (
    diag(group_var(C %*% t(beta_u), r)) + Sigb + Sigd)
  
  sig2E = alpha[1]/(1-alpha[1]) * (
    mean(1/alpha[3]*SigG))
  sig2F = alpha[2]/(1-alpha[2]) * (
    mean(1/alpha[4]*SigH))
  
  return(list(
    W_list = W_list,
    V_list = V_list,
    beta_t = beta_t,
    beta_u = beta_u,
    Sigb = Sigb,
    Siga = Siga,
    Sigd = Sigd,
    SigG = SigG,
    SigH = SigH,
    sig2E = sig2E,
    sig2F = sig2F
  ))
}



generate_params_ad <- function(X, Y, p, q, K, r, betaELE_t, betaELE_u, alpha = 0.1) {
  if(length(alpha) == 1) alpha <- rep(alpha, 4)
  if(!(length(alpha) %in% c(1,4))) stop("length alpha should be 1 or 4")
  
  time = c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  outp <- list(
    set.seed(123),
    W = orth(matrix(rnorm(p*r), p, r)+1),
    V = orth(matrix(rnorm(q*r), q, r)+1),
    set.seed(NULL),
    beta_t = matrix(betaELE_t, ncol = 2*r),
    beta_u = matrix(betaELE_u, ncol = 2*r),
    
    Sigb = diag(r:1,r)/2,
    Siga = diag(r:1,r)*2,
    Sigd = diag(r:1, r)*10
  )
  
  outp$SigG = alpha[3]/(1-alpha[3]) * (
    diag(group_var(C %*% t(outp$beta_t), r)) + outp$Sigb + outp$Siga)
  outp$SigH = alpha[4]/(1-alpha[4]) * (
    diag(group_var(C %*% t(outp$beta_u), r)) + outp$Sigb + outp$Sigd)
  
  return(with(outp, {
    c(outp,
      sig2E = alpha[1]/(1-alpha[1]) * (
        mean(1/alpha[3]*SigG))
      ,
      sig2F = alpha[2]/(1-alpha[2]) * (
        mean(1/alpha[4]*SigH))
    )
  }))
}





# Fun: generate data and record -------------------------------------------

#' Generate data for a Dynamic Components Model
#' 
#' @param N Positive integer. Sample size to be simulated
#' @param p Positive integer. Number of features to be simulated
#' @param K positive integer. Number of time points to be simulated
#' @param r Positive integer. Number of principal components to be simulated
#' @param distr Functions. One of the random number generator functions (e.g. rnorm) for generating latent variables
#' 
#' @return A list containing observed X, Y; latent T, U, b, in long-format
generate_data_2ds <- function(N, p, q, K, r, params, distr = rnorm) {
  W = params$W
  V = params$V
  beta_t = params$beta_t
  beta_u = params$beta_u
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigG = params$SigG
  SigH = params$SigH
  Sigb = params$Sigb
  Siga = params$Siga
  Sigd = params$Sigd
  
  p = nrow(W)
  q = nrow(V)
  r = ncol(W)
  
  Gamma = Gamma <- rbind( ## also can use 'sparse'
    (diag(1,K)%x%W) %*% cbind(diag(1,(r*K)), matrix(0, nrow=r*K, ncol=(3*r+K*r))),
    (diag(1,K)%x%V) %*% cbind(matrix(0, nrow=(r*K), ncol=(2*r+K*r)), diag(1,(r*K)), 
                              matrix(0, nrow=(r*K), ncol=r))
  )
  
  time <- c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  A1 <- A %*% (Sigb+Siga) %*% t(A) + diag(1,K) %x% SigG
  A2 <- A %*% Sigb; A3 <- A %*% Siga; A4 <- A %*% Sigb %*% t(A)
  A5 <- matrix(0, nrow = r*K, ncol = r)
  A6 <- Sigb; A7 <- matrix(0, nrow = r, ncol = r); A8 <- Sigb %*% t(A)
  A9 <- matrix(0, nrow = r, ncol = r); A10 <- Siga
  A11 <- matrix(0, nrow = r, ncol = r*K); A12 <- matrix(0, nrow = r, ncol = r)
  A13 <- A %*% (Sigb+Sigd) %*% t(A) + diag(1,K) %x% SigH
  A14 <- A %*% Sigd; A15 <- Sigd
  
  VarZ = blockmm(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15)
  
  
  # Z <- matrix(distr(N*p*K, mean = C%*%t(beta)), ncol = 1) %*% chol(VarZ)
  Z <- MASS::mvrnorm(N, mu = rbind(C%*%t(beta_t), matrix(0, nrow=r, ncol=1),
                                   matrix(0, nrow=r, ncol=1), C%*%t(beta_u), 
                                   matrix(0, nrow=r, ncol=1)),
                     Sigma = VarZ)
  
  E <- scale(matrix(distr(N*p*K), ncol = 1)) * sqrt(sig2E)
  Ff <- scale(matrix(distr(N*q*K), ncol = 1)) * sqrt(sig2F)
  
  dat <- Z %*% t(Gamma) + cbind(matrix(E, nrow=N),matrix(Ff, nrow=N))
  return(list(X = dat[, 1:(K*p)],
              Y = dat[, (K*p+1):(K*p+K*q)],
              Tt = Z[, 1:(K*r)],
              b = Z[, (K*r+1):(K*r+r)],
              a = Z[, (K*r+r+1):(K*r+2*r)],
              U = Z[, (K*r+2*r+1):(K*r+2*r+K*r)],
              d = Z[, (K*r+2*r+K*r+1):(K*r+2*r+K*r+r)]
  ))
}




generate_data_2ds_KW <- function(N, p, q, K, r, params, distr = rnorm) {
  W_list = params$W_list
  V_list = params$V_list
  beta_t = params$beta_t
  beta_u = params$beta_u
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigG = params$SigG
  SigH = params$SigH
  Sigb = params$Sigb
  Siga = params$Siga
  Sigd = params$Sigd
  
  stopifnot(length(W_list) == K)
  
  Gamma = Gamma <- rbind(
    block_diag(W_list) %*% cbind(diag(1,(r*K)), matrix(0, nrow=r*K, ncol=(3*r+K*r))),
    block_diag(V_list) %*% cbind(matrix(0, nrow=(r*K), ncol=(2*r+K*r)), diag(1,(r*K)), 
                                 matrix(0, nrow=(r*K), ncol=r))
  )
  
  time <- c(1:K)
  A = rep(1,K) %x% diag(1,r)
  B = time %x% diag(1,r)
  C = cbind(A,B)
  
  A1 <- A %*% (Sigb+Siga) %*% t(A) + diag(1,K) %x% SigG
  A2 <- A %*% Sigb; A3 <- A %*% Siga; A4 <- A %*% Sigb %*% t(A)
  A5 <- matrix(0, nrow = r*K, ncol = r)
  A6 <- Sigb; A7 <- matrix(0, nrow = r, ncol = r); A8 <- Sigb %*% t(A)
  A9 <- matrix(0, nrow = r, ncol = r); A10 <- Siga
  A11 <- matrix(0, nrow = r, ncol = r*K); A12 <- matrix(0, nrow = r, ncol = r)
  A13 <- A %*% (Sigb+Sigd) %*% t(A) + diag(1,K) %x% SigH
  A14 <- A %*% Sigd; A15 <- Sigd
  
  VarZ = blockmm(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15)
  
  Z <- MASS::mvrnorm(N, mu = rbind(C%*%t(beta_t), matrix(0, nrow=r, ncol=1),
                                   matrix(0, nrow=r, ncol=1), C%*%t(beta_u), 
                                   matrix(0, nrow=r, ncol=1)),
                     Sigma = VarZ)
  
  E <- scale(matrix(distr(N*p*K), ncol = 1)) * sqrt(sig2E)
  Ff <- scale(matrix(distr(N*q*K), ncol = 1)) * sqrt(sig2F)
  
  dat <- Z %*% t(Gamma) + cbind(matrix(E, nrow=N),matrix(Ff, nrow=N))
  return(list(X = dat[, 1:(K*p)],
              Y = dat[, (K*p+1):(K*p+K*q)],
              Tt = Z[, 1:(K*r)],
              b = Z[, (K*r+1):(K*r+r)],
              a = Z[, (K*r+r+1):(K*r+2*r)],
              U = Z[, (K*r+2*r+1):(K*r+2*r+K*r)],
              d = Z[, (K*r+2*r+K*r+1):(K*r+2*r+K*r+r)]
  ))
}




# Fun: long and wide-format data ------------------------------------------

pivot_longer_wider_2ds <- function(Tt, U, X, Y, K, p, q, r) {
  pivot_longer_comp_T <- function(Tt = data$Tt, K, r) {
    ID <- rep(1:N, each = K*r)
    time <- rep(rep(1:K, each = r), times = N)
    comp <- rep(rep(1:r, times = K), times = N)
    value <- as.vector(t(Tt))
    
    Tt_long <- data.frame(ID = ID, time = time, comp = comp, t = value)
    return(Tt_long)
  }
  
  pivot_longer_comp_U <- function(U = data$U, K, r) {
    ID <- rep(1:N, each = K*r)
    time <- rep(rep(1:K, each = r), times = N)
    comp <- rep(rep(1:r, times = K), times = N)
    value <- as.vector(t(U))
    
    U_long <- data.frame(ID = ID, time = time, comp = comp, u = value)
    return(U_long)
  }
  
  pivot_longer_data_X <- function(X = data$X, K, p) {
    ID <- rep(1:N, each = K*p)
    time <- rep(rep(1:K, each = p), times = N)
    feature <- rep(rep(1:p, times = K), times = N)
    value <- as.vector(t(X))
    
    X_long <- data.frame(ID = ID, time = time, feature = feature, x = value)
    return(X_long)
  }
  
  pivot_longer_data_Y <- function(Y = data$Y, K, q) {
    ID <- rep(1:N, each = K*q)
    time <- rep(rep(1:K, each = q), times = N)
    feature <- rep(rep(1:q, times = K), times = N)
    value <- as.vector(t(Y))
    
    Y_long <- data.frame(ID = ID, time = time, feature = feature, y = value)
    return(Y_long)
  }
  
  pivot_wider_comp_T <- function(Tt_long, r) {
    Tt_long <- Tt_long %>%
      mutate(comp = paste0("t_", comp))
    
    Tt_wide <- Tt_long %>%
      pivot_wider(names_from = comp, values_from = t) %>%
      arrange(ID, time)
    
    return(Tt_wide)
  }
  
  pivot_wider_comp_U <- function(U_long, r) {
    U_long <- U_long %>%
      mutate(comp = paste0("u_", comp))
    
    U_wide <- U_long %>%
      pivot_wider(names_from = comp, values_from = u) %>%
      arrange(ID, time)
    
    return(U_wide)
  }
  
  pivot_wider_data_X <- function(X_long, p) {
    X_long <- X_long %>%
      mutate(feature = paste0("x_", feature))
    
    X_wide <- X_long %>%
      pivot_wider(names_from = feature, values_from = x) %>%
      arrange(ID, time)
    
    return(X_wide)
  }
  
  pivot_wider_data_Y <- function(Y_long, q) {
    Y_long <- Y_long %>%
      mutate(feature = paste0("y_", feature))
    
    Y_wide <- Y_long %>%
      pivot_wider(names_from = feature, values_from = y) %>%
      arrange(ID, time)
    
    return(Y_wide)
  }
  
  Tt_long <- pivot_longer_comp_T(Tt, K = K, r = r)
  X_long <- pivot_longer_data_X(X, K = K, p = p)
  
  Tt_wide <- pivot_wider_comp_T(Tt_long, r)
  X_wide <- pivot_wider_data_X(X_long, p)
  
  U_long <- pivot_longer_comp_U(U, K = K, r = r)
  Y_long <- pivot_longer_data_Y(Y, K = K, q = q)
  
  U_wide <- pivot_wider_comp_U(U_long, r)
  Y_wide <- pivot_wider_data_Y(Y_long, q)
  
  return(list(Tt_long = Tt_long, X_long = X_long,
              Tt_wide = Tt_wide, X_wide = X_wide,
              U_long = U_long, Y_long = Y_long,
              U_wide = U_wide, Y_wide = Y_wide))
}

