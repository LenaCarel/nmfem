#' Extract log-likelihood from a mixture of multinomials
#'
#' @param X a matrix of dimension \code{N} (number of observation) \code{x M} (number of variables) containing multinomials observations.
#' @param Theta matrix of dimension \code{M x H}.
#' @param Lambda matrix of dimension \code{H x K}. Can be \code{NULL}.
#' @param p vector containing the proportions of each cluster. Must be of dimension \code{K} (or \code{H} if \code{Lambda} is \code{NULL}).
#' @return The function returns the log-likelihood of the data to the model
#' @examples
#' travelers <- travelers[ ,-1]
#' M <- ncol(travelers)
#' K <- 5
#'
#' Theta0    <- t(dplyr::sample_n(travelers, K))
#' Theta0    <- Theta0 / matrix(rep(apply(Theta0, 2, sum), M), nrow = M, ncol = K, byrow = TRUE)
#' travelers <- as.matrix(travelers)
#' p0        <- rep(1 / K, K)
#'
#' llh <- loglik_mult(travelers, Theta0, p = p0)
#' llh
#'
#' @export

loglik_mult <- function(X, Theta, Lambda = NULL, p){
  
  H <- ncol(Theta)
  if(is.null(Lambda)) Lambda <- diag(nrow = H, ncol = H)
  
  lTL <- log(Theta %*% Lambda)
  lTL[lTL == -Inf] <- 0
  
  Mat = X %*% lTL
  rowmax_Mat = apply(Mat, 1, max)
  
  Mat = exp(Mat - rowmax_Mat)
  out_mat <- log( Mat %*% p)
  out_mat[out_mat == -Inf] <- 0
  output <- sum(out_mat) + sum(rowmax_Mat)
  
  return(output)
}

#' Logarithm of a factorial
#'
#' @param x an integer
#' @return log(factorial(x))
#' @keywords internal
lnfact <- function(x){
  if(x == 0){
    output <- 0
  }
  else{
    output <- sum(log(c(1:x)))
  }
  return(output)
}

#' Log-likelihood function from M-step, when t is estimated
#'
#' @param Xtt matrix. Correspond to the matrix multiplication of t(X) and t.
#' @param Theta matrix of dimension \code{M x H}.
#' @param Lambda matrix of dimension \code{H x K}.
#' @return Returns the log-likelihood function from M-step, when t is estimated.
#' @keywords internal

Q <- function(Xtt, Theta, Lambda){
  tl <- Theta %*% Lambda
  tl[tl == 0] <- 10^(-12)
  tl[is.na(tl)] <- 10^(-12)

  output <- sum(Xtt * log(tl))

  return(output)
}

#' Update of lambda
#'
#' @param Xtt matrix. Correspond to the matrix multiplication of t(X) and t.
#' @param Theta matrix of dimension \code{M x H}.
#' @param Lambda matrix of dimension \code{H x K}.
#' @return Returns the updated value of Lambda
#' @keywords internal

lambda_update <- function(Xtt, Theta, Lambda){

  H <- ncol(Theta)
  K <- ncol(Lambda)
  M <- nrow(Theta)

  mpow <- floor(log10(Theta))
  mpow[mpow == -Inf] <- 0
  pow <- max(mpow)
  if(pow > 200) Theta <- Theta*10^(-pow)

  tmp <- Xtt / (Theta %*% Lambda)
  tmp <- ifelse(is.na(tmp), 0, tmp)
  tmp <- ifelse(tmp == Inf, 0, tmp)

  mpow <- floor(log10(tmp))
  mpow[mpow == -Inf] <- 0
  pow <- max(mpow)
  if(pow > 200) tmp <- tmp * 10 ^ (- pow)
  sum_Theta <- matrix(rep(apply(Theta, 2, sum),K), nrow = H, ncol = K)

  Lambda <- Lambda * (t(Theta) %*% tmp) / sum_Theta
  Lambda <- Lambda / matrix(rep(apply(Lambda, 2, sum), H), nrow = H, ncol = K, byrow = TRUE)
  return(Lambda)
}

#' Update of theta
#'
#' @param Xtt matrix. Correspond to the matrix multiplication of t(X) and t.
#' @param Theta matrix of dimension \code{M x H}.
#' @param Lambda matrix of dimension \code{H x K}.
#' @return Returns the updated value of Theta
#' @keywords internal

theta_update <- function(Xtt, Theta, Lambda){

  H <- ncol(Theta)
  K <- ncol(Lambda)
  M <- nrow(Theta)

  tmp <- Xtt / (Theta %*% Lambda)
  tmp <- ifelse(is.na(tmp), 0, tmp)
  tmp <- ifelse(tmp == Inf, 0, tmp)
  sum_Lambda <- matrix(rep(apply(Lambda, 1, sum), M), nrow = M, ncol = H, byrow = TRUE)

  Theta <- Theta * (tmp %*% t(Lambda)) / sum_Lambda
  return(Theta)
}

#' Update of t
#'
#' @param X a matrix of dimension \code{N} (number of observation) \code{x M} (number of variables) containing multinomials observations.
#' @param Theta matrix of dimension \code{M x H}.
#' @param Lambda matrix of dimension \code{H x K}.
#' @param p vector containing the proportions of each cluster. Must be of dimension \code{K}.
#' @return Returns the updated value of Theta
#' @keywords internal

t_update <- function(X, Theta, Lambda, p){
  n <- nrow(X)
  K <- ncol(Lambda)
  t <- matrix(data = NA, nrow = n, ncol = K)

  logTL <- log(Theta %*% Lambda)
  logTL[logTL == -Inf] <- log(1e-60)
  TL <- Theta %*% Lambda
  for(i in 1:n){
    logvec <- log(p) + X[i, ] %*% logTL
    Expdifflogvec <- exp(matrix(data = logvec, nrow = K, ncol = K, byrow=TRUE) - matrix(data = logvec, nrow = K, ncol = K, byrow = FALSE))
    t[i, ] <- 1 / rowSums(Expdifflogvec)
  }
  return(t)
}



