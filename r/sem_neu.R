library(MASS)

library(mixtools)

library(plyr)


compute_d <- function(params_vec){
    return(- 0.5 + sqrt(1/ 4 + length(params_vec))) # p/q-formula to compute dimension
    
}

param_vec_to_list <- function(params_vec){
    d = compute_d(params_vec)
    
    mu = params_vec[1:d]
    sigma = matrix(params_vec[(d+1):(d^2)], nrow = d , ncol = d )
    
    return(list(mu, sigma))
}

param_list_to_vec <- function(params_list){
    return(unlist(params_list))
}

p_list = param_vec_to_list(params_vec)

p_vec = param_list_to_vec(p_list)


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

calculateDM <- function(X, param_df, tol) {
  
  # compute dimension d
  d <- compute_d(unlist(param_df[1, ]))
  
  # Compute r_ij for all i and j
  DM <- matrix(nrow = d, ncol = d)
  for(i in 1:d) {
    for(j in 1:d) {
      DM[i,j] <- compute_r_ij(X, param_df, i, j, tol)
    }
    
  }
  return(DM)
}

sem_algorithm <- function(data, params, tolerance) {
  
  n <- length(data[,1])
  
  # params <- result$params
  
  
  # Get DM* matrix
  
  DM <- calculateDM(data, params, tolerance)
 
  
  # Get covariance matrix of MLE estimate (last step from em algorithm)
  
  cov <- tail(params, 1)[[1]]$cov
  
  # Preparation of c
  
  c <- -cov[1,2]^6+3*cov[1,1]*cov[2,2]*cov[1,2]^4-3*(cov[1,1]*cov[2,2]*cov[1,2])^2+(cov[1,1]*cov[2,2])^3
  
  # Calculate G1 if I_oc
  G1_11 <- cov[1,1]
  G1_12 <- 0
  G1_22 <- 2*cov[1,1]^2
  G1 <- (1/n)*matrix(c(G1_11, G1_12, G1_12, G1_22), nrow = 2, ncol = 2)
  
  # Calculate G2 if I_oc
  G2_11 <- cov[1,2]
  G2_12 <- 0
  G2_13 <- 0
  G2_21 <- 0
  G2_22 <- 2*cov[1,1]*cov[2,2]
  G2_23 <- 2*det(cov)^2*(cov[1,1]*cov[2,2]*cov[1,2]^2 - cov[1,2]^4)/c
  G2 <- (1/n)*matrix(c(G2_11, G2_12, G2_13, G2_21, G2_22, G2_23), nrow = 2, ncol = 3, byrow = TRUE)
  
  
  # Calculate G3 of I_oc
  G3_11 <- cov[2,2]
  G3_12 <- 0
  G3_13 <- 0
  G3_22 <- det(cov)^2*((cov[1,1]*cov[2,2])^2 - cov[1,2]^4)/c
  G3_23 <- 2*cov[2,2]*cov[1,2]
  G3_33 <- 2*cov[2,2]^2
  G3 <- (1/n)*matrix(c(G3_11, G3_12, G3_13, G3_12, G3_22, G3_23, G3_13, G3_23, G3_33), nrow = 3, ncol = 3)
  
  # Compute Delta V*
  DV <- (G3 - t(G2)%*%solve(G1)%*%G2)%*%DM%*%solve(diag(3)-DM)
  
  return(G3 + DV)
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
