# ------------------------------------------------------------------------------
# This script containes ad-hoc functions to implement the analysis of
# spatial correlation patterns of water quality variables in fluvial networks. 
#
# Paper: Analyzing Spatial Correlation Patterns of Water Quality Variables in Fluvial Networks
# Authors: N. Pronello, S. Castiglia, V. Frontuto, N. Golini, R. Ignaccolo, L. Ippoliti
# Journal: Journal of Agricultural, Biological, and Environmental Statistics 
# Year: 2026
# DOI: 10.1007/s13253-026-00726-9
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# cv_fun_cov_holes(): Cross-validation function for selecting the smoothing
# parameter gamma
#
# nn: number of water-bodies in the network domain
# gamma_values: grid of candidate values for gamma
# L: laplacian matrix
# cov_star: list of rempirical covariance matreices
# id_sample: id of sampled vertices in the network
# ------------------------------------------------------------------------------
cv_fun_cov_holes <- function(nn,gamma_values,L,cov_star,id_sample){
  
  LOOCV <- numeric(length(gamma_values))
  nvo <- length(id_sample)
  id_r <- 1:nvo
  
  for (i in seq_along(gamma_values)) {
    gamma <- gamma_values[i]
    K <- solve(diag(nn)+gamma*L)
    K_r_temp <- matrix(0,nvo,nvo)
    for (j in 1:nvo) {
      id_temp <- id_sample[j]
      weights <- K[id_temp,id_sample] 
      K_r_temp[j,] <- weights/sum(weights )
    }
    Sigma_k <- vector("list", nvo)
    
    for (j in 1:nvo) {
      weights <- diag(K_r_temp[j, ]) 
      weights <- weights[-j, -j] 
      weights <- weights / sum(weights)  
      S <- matrix(0,dim(cov_star[[1]])[1],dim(cov_star[[1]])[2])
      for(k in 1:(nvo-1)){
        S <- S+weights[k]*cov_star[[id_r[-j][k]]]
      }
      Sigma_k[[j]] <- S
    }
    
    loocv_error1 <- numeric(nvo)
    loocv_error2 <- numeric(nvo)
    for (j in 1:nvo) {
      Theta <- pseudoinverse(Sigma_k[[j]])  
      loocv_error1[j] <- norm(Sigma_k[[j]]-cov_star[[j]], type = "F")^2
      loocv_error2[j] <- sum(diag(Theta%*%cov_star[[j]]))
    }
    LOOCV[i] <- sqrt(sum(loocv_error1)*sum(loocv_error2))
  }
  return(cbind(gamma_values,LOOCV))
}


# ------------------------------------------------------------------------------
# betti_0_filtration_mst(): function to obtain Betti_0 curves from correlation
# matrices. This function follows from Lee H, Kang H, Chung MK, 
# Kim BN, Lee DS. Persistent brain network homology from the perspective of dendrogram. 
# IEEE Trans Med Imaging. 2012 Dec;31(12):2267-77. 
# doi: 10.1109/TMI.2012.2219590. Epub 2012 Sep 19. PMID: 23008247.
# adj_matrix: is a matrix that obtained as 1-cor().
# The output of this function are the filtration values and Betti_0 values.
# ------------------------------------------------------------------------------
betti_0_filtration_mst <- function(adj_matrix) {
  diag(adj_matrix) <- 0
  p <- dim(adj_matrix)[1]
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = T)
  g_mst <- mst(g)
  return(data_frame(eta=1-c(0,sort((E(g_mst)$weight))),c=1:p))
}

# ------------------------------------------------------------------------------
# ks_distance(): function to calculate the KS-distance between two Betti curves
# curve1, curve2: first vectors resulting from betti_0_filtration_mst()
# q: number of discretization points to consider
# ------------------------------------------------------------------------------
ks_distance <- function(curve1,curve2,q=100){
  x1 <- curve1
  x2 <- curve2
  y <- length(curve1):1
  x_common <- seq(0, 1, length.out = q)
  y1_interp <- approx(x1, y, xout = x_common, method = "constant", f = 0, rule = 1)$y
  y2_interp <- approx(x2, y, xout = x_common, method = "constant", f = 0, rule = 1)$y
  df_temp <- cbind(y1_interp,y2_interp)
  df_temp <- df_temp[!rowSums(!is.finite(df_temp)),]
  max(abs(df_temp[,1] - df_temp[,2]))
}
