library(MASS)
library(mixtools)
library(plyr)

#####################
### SEM Algorithm ###
#####################

calculateRatio <- function(data, params, i, j, tolerance) {
  
  indicator <- is.na(data[,2])
  vecR <- NULL
  
  # Get last parameter element (MLE)
  mle <- unlist(tail(params, 1)[[1]])
  
  # Compute rate until convergence
  t <- 1
  repeat {
    # Get current state of theta from em computation
    theta_t <- unlist(params[[t]])
    
    # Define theta_i as the final mle and replace i-th value with current state
    theta_t_i <- mle
    theta_t_i[[i]] <- theta_t[[i]]
    
    # Hier ist falsch
    lambda = c(theta_t_i[1,2], theta_t_i[2,1] )
    sigma = c(theta_t_i[1,1], theta_t_i[2,2] )
    print(lambda)
    theta_t_1_i <- unlist(normalmixEM(data, lambda = lambda, sigma = sigma )) # EM einbauen
    
    # Calculate ratio
    vecR <- append(vecR, (theta_t_1_i[[j]] - mle[[j]])/(theta_t[[i]] - mle[[i]]))
    
    # Increase iteration
    t <- t+1
    
    # Check if convergence criterion is hit or we're running out of original estimations
    last <- length(vecR)
    if((last >= 2 && abs(vecR[last] - vecR[last-1]) < tolerance)) {
      break
    }
    
    # Check if there is still a parameter from EM to calculate next iteration
    if(t >= length(params)) {
      warning("SEM did not converge for one component.")
      break
    }
    
  }
  # Just return last rate after convergence
  return(vecR[length(vecR)])
  
}



calculateDM <- function(data, params, tolerance) {
  
  # Parameters to estimate in DM*
  # 1: mu1, 2: mu2, 3: s_xx, 4/5: s_xy, 6: s_yy
  estimates <- c(2,4,6)
  
  # Number of parameters to calculate variance for
  d <- length(estimates)
  
  # Define empty DM matrix
  DM <- matrix(nrow = d, ncol = d)
  
  # Calculate any r_ij and store in DM
  for(i in 1:d) {
    for(j in 1:d) {
      DM[i,j] <- calculateRatio(data, params, estimates[i], estimates[j], tolerance)
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
