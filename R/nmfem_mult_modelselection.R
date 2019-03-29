#' Model selection in NMF-EM algorithm for mixture of multinomials
#'
#' The function proceed to a model selection with NMF-EM algorithm on mixture of multinomials dataset.
#' First, the function plots the log-likelihood in function of K.
#' Second, log-likelihood is plotted in function of H.
#' We recommend the user to choose K and H by slope heuristic method.
#'
#' @param X numeric matrix containing multinomials observations of dimension \code{N} (number of observation) \code{x M} (number of variables).
#' @param maxK integer. Maximum number of clusters to be tested. By default, function tests from 2 to 30 clusters.
#' @param save logical. Whether the result of each parameter couple \code{(H,K)} tested got to be saved.
#' @param path path to save the results if \code{save = TRUE}. By default, it is the working directory. Three directories are created to save the results. Directory "Initializations" contains the initialization of the algorithm for each value of K. Matrices are saved in directory "Matrices", and plots in directory "Results".
#' @examples
#'
#' # Example on the complete data - needs around an hour to run
#' \dontrun{
#' nmfem_mult_modelselection(travelers[ ,-1])
#' }
#' @importFrom grDevices dev.copy dev.off png
#' @importFrom graphics plot
#' @export

nmfem_mult_modelselection <- function(X, maxK = 30, save = FALSE, path = "."){

  if(save){
    init <- paste0(path,"/Initializations")
    mat  <- paste0(path,"/Matrices")
    res  <- paste0(path,"/Results")
    dir.create(init, showWarnings = FALSE)
    dir.create(mat, showWarnings = FALSE)
    dir.create(res, showWarnings = FALSE)
  }else{
    init <- NULL
    mat  <- NULL
    res  <- NULL
  }

  mllh <- as.data.frame(matrix(nrow = maxK, ncol = maxK))

  # Selection of K
  for(K in 2:maxK){
    print(paste("Test for K =", K))
    H <- K
    if(is.na(mllh[H,K])){
      tmp <- nmfem_mult(X, H, K, init)
      mllh[H,K] <- tmp$llh
      if(save){
        Theta <- tmp$Theta
        Lambda <- tmp$Lambda
        p <- tmp$p
        t <- tmp$t
        save(Theta,Lambda,p,t, file=paste0(mat,"/Matrices_H",H,"_K",K,".RData"))
      }
    }
  }

  plot(x = c(2:maxK), y = diag(as.matrix(mllh[-1,-1])), type = "b",
       main = "EM likelihood", xlab = "K number of clusters", ylab = "Log-likelihood")
  if(save){
    dev.copy(png, paste0(res,"/llh_EM.png"))
    dev.off()
  }

  K <- NULL
  while(is.null(K)){
    K <- suppressWarnings(as.integer(readline("\n By looking at the plot, what value of K do you choose ? ")))
    if(is.na(K)) K <- NULL
  }

  for(H in (K-1):2){
    print(paste("H =", H, ", K =", K))
    if(is.na(mllh[H,K])){
      tmp       <- nmfem_mult(X, H, K, init)
      mllh[H,K] <- tmp$llh
      if(save){
        Theta     <- tmp$Theta
        Lambda    <- tmp$Lambda
        p <- tmp$p
        t <- tmp$t
        save(Theta, Lambda, p, t, file = paste0(mat,"/Matrices_H", H, "_K", K, ".RData"))
      }
    }
  }

  plot(x = c(2:K), y = mllh[2:K,K], type = "b", main = paste0("NMF-EM likelihood for K = ", K),
       ylab = "Log-likelihood", xlab = "H number of words")
  if(save){
    dev.copy(png, paste0(res,"/llh_K", K, ".png"))
    dev.off()
  }


  H <- NA
  while(is.na(H)){
    H <- suppressWarnings(as.integer(readline("\n By looking at the plot, what value of H do you choose ? ")))
    if(H>=K){
      answer <- readline("\n You chose H >= K! \n Code will run, but you won't reduce the number of parameters. \n Are you sure you want to keep it that way ? (y/n)")
      if(toupper(answer)=="N") H <- NA
    }
  }

  cat(paste("\n For your data, you chose to select the model with K =",K,"and H =",H,"."))
  cat( "\n" )
  if(save) cat(paste0('You can access corresponding matrices by using load("',mat,'/Matrices_H', H, '_K', K, '.RData")'))
}
