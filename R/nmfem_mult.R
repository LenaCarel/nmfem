#' NMF-EM algorithm for mixture of multinomials
#'
#' Proceed to an NMF-EM algorithm on mixture of multinomials dataset. In comparison to the classical EM algorithm, the number of parameters to estimate is lower. For more explanation, see pre-print of Carel and Alquier (2017) <arXiv:1709.03346>.
#'
#' @param X a matrix containing multinomials observations of dimension \code{N} (number of observation) \code{x M} (number of variables).
#' @param H number of words.
#' @param K number of clusters.
#' @param path path to the directory to save the initialization or to load it. NULL by default, won't save or load it.
#' @param eps_init convergence criterion on the initialization. Default value is 1e-3.
#' @param eps_M convergence criterion on the Maximization step. Default value is 1e-8.
#' @param eps_llh convergence criterion on the log-likelihood. Default value is 1e-5.
#' @return A list with the elements:
#' \item{Theta}{matrix of dimension \code{M x H}. Contains a dictionnary of redundant components.}
#' \item{Lambda}{matrix of dimension \code{H x K}. Contains the expression of the \code{K} clusters in the dictionnary.}
#' \item{llh}{log-likelihood of the model.}
#' \item{p}{vector containing the proportions of each cluster.}
#' \item{posterior}{matrix containing for each observation the posterior probability to belong to each cluster.}
#' @examples
#'
#' # Example on a data sample
#' x <- dplyr::sample_n(travelers[,-1],900)
#' out <- nmfem_mult(x, H = 4, K = 7)
#' # Display first cluster profile
#' display_profile(t((out$Theta %*% out$Lambda)[ ,1]))
#' # Display first word profile
#' display_profile(t(out$Theta[ ,1]), color = "Greens")
#'
#' # Example on the complete data - it needs a few minutes to run
#' \dontrun{
#' nmfem_travelers <- nmfem_mult(travelers[ ,-1], H = 5, K = 10)
#' Theta <- nmfem_travelers$Theta
#' Lambda <- nmfem_travelers$Lambda
#'
#' # Display first cluster profile
#' display_profile(t((Theta %*% Lambda)[ ,1]))
#'
#' # Display first word profile
#' display_profile(t(Theta[ ,1]), color = "Greens")}
#'
#' @importFrom stats dist runif
#' @export

nmfem_mult <- function(X, H, K, path = NULL, eps_init = 1e-3, eps_M = 1e-8, eps_llh = 1e-5){

  #### Entry parameters verification ----
  if(!is.null(path)){
    path <- ifelse(is.character(path), path, NULL)
  }

  #### Parameters ----
  M <- ncol(X)
  n <- nrow(X)

  Theta0     <- t(dplyr::sample_n(X, K))
  Theta0     <- Theta0 / matrix(rep(apply(Theta0, 2, sum), M), nrow = M, ncol = K, byrow = TRUE)
  X          <- as.matrix(X)
  Lambda0    <- diag(nrow = K, ncol = K)
  p0         <- rep(1 / K, K)
  posterior0 <- matrix(data = NA, nrow = n, ncol = K)

  crit     <- Inf

  restart_init <- TRUE

  fichier <- paste0("HK", K, ".RData")

  #### Initialization ----
  if(!is.null(path) && H != K){
    load(paste(path, fichier, sep = "/"))
  }else{

    llh <- loglik_mult(X, Theta0, Lambda0, p0)

    while(restart_init){
      initEM       <- mixtools::multmixEM(y = X, k = K, epsilon = eps_init)
      Theta0       <- t(initEM$theta)
      restart_init <- any(dist(t(Theta0)) < 1e-4)
    }

    p0               <- initEM$lambda
    posterior0       <- initEM$posterior
    rownames(Theta0) <- colnames(X)

    if(!is.null(path)) save(Theta0, Lambda0, p0, posterior0, file = paste(path, fichier, sep = "/"))
  }

  llh <- loglik_mult(X, Theta0, Lambda0, p0)

  if(H != K){
    #### First Maximization step ----
    Xtt       <- t(X) %*% posterior0
    Theta     <- matrix(data <- runif(M*H), nrow <- M, ncol <- H)
    Lambda    <- matrix(data <- runif(H*K), nrow <- H, ncol <- K)
    p         <- p0
    posterior <- posterior0
    restart   <- TRUE

    q      <- Inf
    Expect <- Q(Xtt, Theta, Lambda)
    crit2  <- Inf

    while(crit2 > eps_M){
      q <- Expect

      Lambda <- lambda_update(Xtt, Theta, Lambda)
      Theta <- theta_update(Xtt, Theta, Lambda)

      Expect <- Q(Xtt, Theta, Lambda)
      crit2 <- abs((Expect - q) / Expect)
    }

    Theta <- Theta / matrix(rep(apply(Theta, 2, sum), M), nrow = M, ncol = H, byrow = TRUE)
    Theta <- ifelse(is.na(Theta), 0, Theta)
    llh <- loglik_mult(X, Theta, Lambda, p)

    #### Main Expectation-Maximization loop ----
    crit <- Inf
    while(crit > eps_llh){
      prev_llh <- llh

      posterior <- t_update(X, Theta, Lambda, p)
      p         <- apply(posterior, 2, sum) / sum(posterior)
      Xtt       <- t(X) %*% posterior
      theta     <- Theta
      lambda    <- Lambda

      # Expectation step
      Expect <- Q(Xtt, Theta, Lambda)
      q      <- Inf
      crit2  <- Inf

      # Maximization step
      while(crit2 > eps_M){
        q      <- Expect

        lambda <- lambda_update(Xtt, theta, lambda)
        theta  <- theta_update(Xtt, theta, lambda)

        Expect <- Q(Xtt, Theta, Lambda)
        crit2  <- abs((Expect-q)/Expect)
      }
      Theta  <- theta
      Theta <- Theta / matrix(rep(apply(Theta, 2, sum), M), nrow = M, ncol = H, byrow = TRUE)
      Lambda <- lambda
      llh    <- loglik_mult(X, Theta, Lambda, p)
      crit   <- abs((llh - prev_llh) / llh)
    }
  } else {
    Theta     <- Theta0
    Lambda    <- Lambda0
    p         <- p0
    posterior <- posterior0
  }
  out <- list()
  out$Theta <- Theta
  out$Lambda <- Lambda
  out$llh <- llh
  out$p <- p
  out$posterior <- posterior

  return(out)
}
