library(MASS)

library(mixtools)

library(plyr)

source("r/em_neu.R")

source("r/utils.R")




compute_r_ij <- function(X, param_df, i, j, tol) {
 
  r_vec <- NULL
  
  # Get last parameter element
  theta_final <- unlist(param_df[nrow(param_df), ])# get theta star
  
  # Compute rate until convergence
  t <- 1
  
  repeat { # TODO: In while schleife umbauen mit abbruchs variable (siehe EM function)
      
    # Get current state of theta from em computation
    theta_t <- unlist(param[t, ])
    
    # Define theta_i as the final theta and replace i-th value with current state at iteration t
    theta_t_i <- theta_final
    theta_t_i[i] <- theta_t[i]
    theta_t1_i_tilde <- norm_em(X, max_iters = 1, initial_param_vec = theta_t_i)
    
    # Calculate ratio
    r_vec <- append(r_vec, (theta_t1_i[j] - theta_final[j])/(theta_t[i] - theta_final[i]))
    
    # Increase iteration
    t <- t + 1
    
    # Check if convergence criterion is hit or we're running out of original estimations
    len <- length(r_vec)
    if((len >= 2 && abs(r_vec[last] - r_vec[len-1]) < tol)) {
      
      break
      
    }
  }
  return(r_vec[len])
  
}

compute_DM <- function(X, param_df6, tol) {
    
    # hard coding: number of params that are affected by contaminated data
    d = 3
    
    # convert param df to version of length 5 (see above for explanation)
    param_df5 = apply(param_df6, 1, param_6_to_5)
    
    param_df5_trans = apply(param_df5, 1, stabilizing_transformation)
    
    # select only the params that are affected by contaminated data
    param_df3 = param_df5_trans[, 3:5]
  
  
  # Compute r_ij for all i and j
  # set up DM_star
  DM_star <- matrix(nrow = d, ncol = d)
  for(i in 1:d) {
    for(j in 1:d) {
      DM_star[i,j] <- compute_r_ij(X, param_df3, i, j, tol)
    }
    
  }
  return(DM_star)
}

# Ab hier verstehn wir es jetzt

sem <- function(X, param_df6, tol) {
    
    param_df5 = apply(param_df6, 1, param_6_to_5)
  
    
    param_vec5 = param_df5[nrow(param_df5), ]

  n <- length(X[,1])
  sig1 = sqrt(param_vec5[2])
  sig2 = sqrt(param_vec5[4])
  rho = param_vec5[5]
  
  # compute the DM matrix
  DM_star <- compute_DM(X, param_df6, tol)
 
  # obtain final cov matrix from em algorithm
  cov_final <- param_vec_to_list(unlist(param_df[nrow(param_df), ]))[[2]]

  # Compute G11 from the I_oc
  G11 <- diag(c( n * sig1^(-2) / (1 - rho^2), (n / 4) * (2 - rho^2) (1 / (1 - rho^2)))) 
  
  # Compute G22 from the I_oc
  G21 = martix(c( -n * sig1^(-1) *  sig2^(-1) * rho * (1 / (1 - rho^2)), 0,
            0, -(n / 4) * rho^2 * (1 / (1 - rho^2)), 
            0, - 0.5 * n * rho), nrow = 3, ncol = 2, byrow = TRUE)
  
  # Compute G12 from the I_oc
  G12 = t(G21) 
  
  # Compute G22 from the I_oc
  G22 = martix(c( n * sig2^(-2) *  (1 - rho^2)^(-1), 0, 0, 
                  0, (n / 4) * (2 - rho^2) * (1 - rho^2)^(-1), -0.5 * n * rho, 
                  0, - 0.5 * n * rho, n *(1 + rho^2)), nrow = 3, ncol = 2, byrow = TRUE)

  
  # Compute Delta V*
  A = (G22 - G21 %*% solve(G11) %*% G12)
  DV_22 <- solve(diag(3)- t(DM_star)) %*% t(DV) %*% A
  
  #setup 5x5 matrix
  I_oc_inv = matrix(1:25, ncol = 5)
  
  I_oc_inv[1:2, 1:2] = G11
  
  I_oc_inv[3:5, 1:2] = G21
  
  I_oc_inv[1:2, 3:5] = G12
  
  I_oc_inv[3:5, 3:5] = G22 + DV_22
  
  
  
}

simulate_data = function(n, missings = 0.0, mu = c(0, 0), Sigma= matrix(c(1,2,2,1),2,2), seed = 12345)
{
  set.seed(seed)
  
  # Simulate bivariate normaldistibuted data
  data = mvrnorm(n = n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  # simulating missing data
  n_miss = round(n * missings)
  data[1 : n_miss, 2] = NA
  
  return(data)
}

# Test

data = simulate_data(30, 0.3,  mu = c(0, 0), Sigma= matrix(c(1,0,0,1),2,2))
data_complete = data[complete.cases(data),]
EM = normalmixEM(data_complete) 
params = list(mu = EM$mu, cov = matrix(c(EM$sigma[1], EM$lambda[1], EM$lambda[2], EM$sigma[2]),2,2))
EM$lambda
EM$sigma

V <- sem_algorithm(data, params, 10^-4)
