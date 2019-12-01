#' Model selection in NMF-EM algorithm for mixture of multinomials
#'
#' The function proceed to a model selection with NMF-EM algorithm on mixture of multinomials dataset.
#' The function use the automated dimension jump approach implemented in the capushe package.
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
#' @importFrom capushe Djump
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
  comp <- as.data.frame(matrix(nrow = maxK, ncol = maxK))

  # Selection of K
  for(K in 2:maxK){
    print(paste("Test for K =", K))
    H <- K
    if(is.na(mllh[H,K])){
      tmp <- nmfem_mult(X, H, K, init)
      mllh[H,K] <- tmp$llh
      comp[H,K] <- (length(X) - 1) * H + (H - 1) * K
      if(save){
        Theta <- tmp$Theta
        Lambda <- tmp$Lambda
        p <- tmp$p
        t <- tmp$t
        save(Theta,Lambda,p,t, file=paste0(mat,"/Matrices_H",H,"_K",K,".RData"))
      }
    }
  }

  dividor <- max(diag(as.matrix(comp[-1,-1])), na.rm = TRUE) * 2
  cap_data_K <- as.data.frame(cbind(paste0("K = ", 2:maxK), # Model name
                                    diag(as.matrix(comp[-1,-1]))/dividor, # Penalization
                                    diag(as.matrix(comp[-1,-1])), # Complexity
                                    -diag(as.matrix(mllh[-1,-1])) # Contrast
  ))
  colnames(cap_data_K) <- c("ModName", "Pen", "Comp", "Contrast")
  cap_data_K <- cap_data_K %>%
    mutate(
      Pen = as.numeric(as.character(Pen)),
      Comp = as.numeric(as.character(Comp)),
      Contrast = as.numeric(as.character(Contrast))
    )

  K <- suppressWarnings(as.numeric(gsub("K = ", "", capushe::Djump(cap_data_K)@model)))

  if(K < 11){
    Kp = 11
  }else{Kp = K}

  for(H in Kp:2){
    print(paste("H =", H, ", K =", K))
    if(is.na(mllh[H,K])){
      tmp       <- nmfem_mult(X, H, K, init)
      mllh[H,K] <- tmp$llh
      comp[H,K] <- (length(X) - 1) * H + (H - 1) * K
      if(save){
        Theta     <- tmp$Theta
        Lambda    <- tmp$Lambda
        p <- tmp$p
        t <- tmp$t
        save(Theta, Lambda, p, t, file = paste0(mat,"/Matrices_H", H, "_K", K, ".RData"))
      }
    }
  }

  dividor <- max(comp[2:Kp,K], na.rm = TRUE) * 2
  cap_data_H <- as.data.frame(cbind(paste0("H = ", 2:Kp), # Model name
                                    comp[2:Kp,K]/dividor, # Penalization
                                    comp[2:Kp,K], # Complexity
                                    -mllh[2:Kp,K] # Contrast
  ))
  colnames(cap_data_H) <- c("ModName", "Pen", "Comp", "Contrast")
  cap_data_H <- cap_data_H %>%
    mutate(
      Pen = as.numeric(as.character(Pen)),
      Comp = as.numeric(as.character(Comp)),
      Contrast = as.numeric(as.character(Contrast))
    )

  H <- suppressWarnings(as.numeric(gsub("H = ", "", capushe::Djump(cap_data_H)@model)))

  cat(paste("\n For your data, we advise you to select the model with K =",K,"and H =",H,"."))
  cat( "\n" )
  if(save) cat(paste0('You can access corresponding matrices by using load("',mat,'/Matrices_H', H, '_K', K, '.RData")'))
}
