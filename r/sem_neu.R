library(MASS)

rm(list = ls())

source("r/em_neu.R")

source("r/utils.R")

source("r/data_simulation.R")

source("r/dm_matrix_neu.R")


# compute_r_ij <- function(X, param_df, i, j, tol) {
#  
#   r_vec <- NULL
#   
#   # Get last parameter element
#   theta_final <- unlist(param_df[nrow(param_df), ])# get theta star
#   
#   # Compute rate until convergence
#   t <- 1
#   
#   repeat { # TODO: In while schleife umbauen mit abbruchs variable (siehe EM function)
#       
#     # Get current state of theta from em computation
#     theta_t <- unlist(param_df[t, ])
#     
#     # Define theta_i as the final theta and replace i-th value with current state at iteration t
#     # to replace the correct parameter, we convert the vector from len 6 to 5 and back again, since em function requires len 6 param vector
#     theta_t_i <- param_6_to_5(theta_final)
#     theta_t_i[i] <- param_6_to_5(theta_t)[i]
#     
#     theta_t_i = param_5_to_6(theta_t_i)
#     
#     theta_t1_i_tilde <- unlist(norm_em(X, max_iters = 1, initial_param_vec = theta_t_i))
#    
#     
#     theta_t1_i_tilde5 = unlist(param_6_to_5(theta_t1_i_tilde))
#     theta_t1_i_tilde5_trans = stabilizing_transformation(theta_t1_i_tilde5)
#     
#     theta_final5 = unlist(param_6_to_5(theta_final))
#     theta_final5_trans = stabilizing_transformation(theta_final5)
# 
#     theta_t5 = unlist(param_6_to_5(theta_t))
#     theta_t5_trans = stabilizing_transformation(theta_t5)
#     
#     
#     # Calculate ratio
#     r_vec <- append(r_vec, (theta_t1_i_tilde5_trans[j] - theta_final5_trans[j])/(theta_t5_trans[i] - theta_final5_trans[i]))
#     
#     # Increase iteration
#     t <- t + 1
#     
#     # Check if convergence criterion is hit or we're running out of original estimations
#     len <- length(r_vec)
#     if((len >= 20 && abs(r_vec[len] - r_vec[len-1]) < tol)) {
#       
#       break
#       
#     }
#   }
#   print(r_vec)
#   return(r_vec[len])
#   
# }
# 
# compute_DM <- function(X, param_df6, tol) {
#     
#     # hard coding: number of params that are affected by contaminated data
#     d = 3
#     
#     # Compute r_ij for all i and j
#     # set up DM_star
#     DM_star <- matrix(nrow = d, ncol = d)
#     for(i in 1:d) {
#         for(j in 1:d) {
#             DM_star[i,j] <- compute_r_ij(X, param_df6, i+2, j+2, tol) # +2 to ensure that only contaminated params are used but DM_star indices i, j = 1,2,3 remain
#         }
#         
#     }
#     return(DM_star)
# }
# 


sem <- function(X, param_df6, tol) {
    
    #param_df5 = apply(param_df6, 1, param_6_to_5)
  
    
    param_vec5 = param_6_to_5(unlist(param_df6[nrow(param_df6), ]))

  n <- length(X[,1])
  sig1 = sqrt(param_vec5[2])
  sig2 = sqrt(param_vec5[4])
  rho = param_vec5[5]
  
  # compute the DM matrix
  DM_star <- compute_DM2(X, param_df6, tol = 0.0001)
 
  # obtain final cov matrix from em algorithm
  #cov_final <- param_vec_to_list(unlist(param_df6[nrow(param_df6), ]))[[2]]

  # Compute G11 from the I_oc
  #print(c(n * sig1^(-2)/(1 - rho^2), (n/4) * (2 - rho^2)*(1/(1 - rho^2))))
  G11 <- diag(c( n * sig1^(-2) / (1 - rho^2), (n / 4) * (2 - rho^2)* (1 / (1 - rho^2)))) 
  
  # Compute G22 from the I_oc
  G21 = matrix(c( -n * sig1^(-1) *  sig2^(-1) * rho * (1 / (1 - rho^2)), 0,
            0, -(n / 4) * rho^2 * (1 / (1 - rho^2)), 
            0, - 0.5 * n * rho), nrow = 3, ncol = 2, byrow = TRUE)
  
  # Compute G12 from the I_oc
  G12 = t(G21) 
  
  # Compute G22 from the I_oc
  G22 = matrix(c( n * sig2^(-2) *  (1 - rho^2)^(-1), 0, 0, 
                  0, (n / 4) * (2 - rho^2) * (1 - rho^2)^(-1), -0.5 * n * rho, 
                  0, - 0.5 * n * rho, n *(1 + rho^2)), nrow = 3, ncol = 3, byrow = TRUE)


  
  # Compute Delta V*
  A = (G22 - G21 %*% solve(G11) %*% G12)
  DV_22 <- solve(diag(3)- t(DM_star)) %*% t(DM_star) %*% A
  #DV_22 <- (G22 - t(G12)%*%solve(G11)%*%G12)%*%DM_star%*%solve(diag(3)-DM_star)
  
  #setup 5x5 matrix
  I_oc = matrix(1:25, ncol = 5)
  delta_V = matrix(rep(0, times = 25), ncol = 5)
  delta_V[3:5, 3:5] = DV_22
  
  I_oc[1:2, 1:2] = G11
  
  I_oc[3:5, 1:2] = G21
  
  I_oc[1:2, 3:5] = G12
  
  I_oc[3:5, 3:5] = G22
  
  I_oc_inv = solve(I_oc)
  
  I_inv_final = I_oc_inv + delta_V
  
  return(I_inv_final)
  
}


# Test

data = simulate_data(18, missings = 0.4,  mu = c(1, 2), sigma= matrix(c(1,.5,.5,1),2,2))

x <- c(8,6,11,22,14,17,18,24,19,23,26,40,4,4,5,6,8,10)
y <- c(59,58,56,53,50,45,43,42,39,38,30,27,NA,NA,NA,NA,NA,NA)
data2 = as.matrix(data.frame(x = x, y = y))

epsilon_em = 0.000000001
param_df6 = norm_em(data2, max_iters = 1000, epsilon = epsilon_em, initial_param_vec = NULL)


V <- sem(data2, param_df6, tol = sqrt(epsilon_em))
