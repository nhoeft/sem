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
    
    param_df5 = t(apply(param_df6, 1, param_6_to_5))
    
    
    param_vec5 = unlist(param_df5[nrow(param_df5), ])
    
    n <- length(X[,1])
    sig1 = sqrt(param_vec5[2])
    sig2 = sqrt(param_vec5[4])
    rho = param_vec5[5]
    
    # compute the DM matrix
    DM_star <- compute_DM2(X, param_df6, tol = 0.0001)
    
    # obtain final cov matrix from em algorithm
    cov <- param_vec_to_list(unlist(param_df6[nrow(param_df6), ]))[[2]]
    
    # Preparation of c
    c <- -cov[1,2]^6+3*cov[1,1]*cov[2,2]*cov[1,2]^4-3*(cov[1,1]*cov[2,2]*cov[1,2])^2+(cov[1,1]*cov[2,2])^3
    
    # Calculate G1 of I_oc^-1
    G1_11 <- cov[1,1]
    G1_12 <- 0
    G1_22 <- 2*cov[1,1]^2
    
    G1 <- (1/n)*matrix(c(G1_11, G1_12, G1_12, G1_22), nrow = 2, ncol = 2)
    
    # Calculate G2 of I_oc^-1
    G2_11 <- cov[1,2]
    G2_12 <- 0
    G2_13 <- 0
    G2_21 <- 0
    G2_22 <- 2*cov[1,1]*cov[2,2]
    G2_23 <- 2*det(cov)^2*(cov[1,1]*cov[2,2]*cov[1,2]^2 - cov[1,2]^4)/c
    
    G2 <- (1/n)*matrix(c(G2_11, G2_12, G2_13, G2_21, G2_22, G2_23),
                       nrow = 2, ncol = 3, byrow = TRUE)
    
    # Calculate G3 of I_oc^-1
    G3_11 <- cov[2,2]
    G3_12 <- 0
    G3_13 <- 0
    G3_22 <- det(cov)^2*((cov[1,1]*cov[2,2])^2 - cov[1,2]^4)/c
    G3_23 <- 2*cov[2,2]*cov[1,2]
    G3_33 <- 2*cov[2,2]^2
    
    G3 <- (1/n)*matrix(c(G3_11, G3_12, G3_13, G3_12, G3_22, G3_23,
                         G3_13, G3_23, G3_33), nrow = 3, ncol = 3)
           
    
    # Compute Delta V*
    DV22 <- (G3 - t(G2)%*%solve(G1)%*%G2)%*%DM_star%*%solve(diag(3)-DM_star)
                
    # # # 
    
    
    
    # Compute G11 from the I_oc
   # print(c(n * sig1^(-2)/(1 - rho^2), (n/4) * (2 - rho^2)*(1/(1 - rho^2))))
   # G11 <- diag(c( n * sig1^(-2) / (1 - rho^2), (n / 4) * (2 - rho^2)* (1 / (1 - rho^2)))) 
    
   # # Compute G22 from the I_oc
   # G21 = matrix(c( -n * sig1^(-1) *  sig2^(-1) * rho * (1 / (1 - rho^2)), 0,
    #                0, -(n / 4) * rho^2 * (1 / (1 - rho^2)), 
   #                 0, - 0.5 * n * rho), nrow = 3, ncol = 2, byrow = TRUE)
    
   # # Compute G12 from the I_oc
   # G12 = t(G21) 
    
    # Compute G22 from the I_oc
   # G22 = matrix(c( n * sig2^(-2) *  (1 - rho^2)^(-1), 0, 0, 
   #                 0, (n / 4) * (2 - rho^2) * (1 - rho^2)^(-1), -0.5 * n * rho, 
   #                 0, - 0.5 * n * rho, n *(1 + rho^2)), nrow = 3, ncol = 3, byrow = TRUE)
    
   # print(dim(G11))
   # print(dim(G12))
  #  print(dim(G21))
   # print(dim(G22))
    
    # Compute Delta V*
  #  A = (G22 - G21 %*% solve(G11) %*% G12)
    #DV_22 <- A %*% t(DM_star) %*% solve(diag(3) - t(DM_star)) # corrected?
  #  DV_22 <- solve(diag(3)- t(DM_star)) %*% t(DM_star) %*% A
    
    #setup 5x5 matrix
    I_oc_inv = matrix(1:25, ncol = 5)
    
    I_oc_inv[1:2, 1:2] = G1
    
    I_oc_inv[3:5, 1:2] = G2
    
    I_oc_inv[1:2, 3:5] = t(G2)
    
    I_oc_inv[3:5, 3:5] = G3 + DV22
    
    #Compute difference as a measure of accuracy
    #I_oc = matrix(c(G1, G2, t(G2), G3), ncol = 5)
    
    #diff = I_oc - I_oc_inv
    
    #return(diff)
    return(I_oc_inv)
}


# Test
#with data from paper (doesnt work)
#x <- c(8,6,11,22,14,17,18,24,19,23,26,40,4,4,5,6,8,10)
#y <- c(59,58,56,53,50,45,43,42,39,38,30,27,NA,NA,NA,NA,NA,NA)
#data = matrix(c(x,y), ncol = 2)

#epsilon_em = 0.000000001
#param_df6 = norm_em(data, max_iters = 1000, epsilon = epsilon_em, initial_param_vec = NULL)

#V <- sem(data, param_df6, tol = sqrt(epsilon_em))


#Simulate data for different number of observations and different fraction of NAs
#Then calculate complete data VCOV and take difference
n_obs = c(20, 50, 100)
p_frac = c(0.1, 0.25, 0.45)

#generate empty lists to store the results in
V_c <- list()
diff <- list()


for (p in p_frac){
    for(n in n_obs){
        #empty lists to store results for each p in
        V_c_helper <- list()
        diff_helper <- list()
        #simulate data
        data = simulate_data(n, missings = p,  mu = c(1, 1), sigma= matrix(c(1,.5,.5,1),2,2))
        #run EM
        epsilon_em = 0.000000001
        param_df6 = norm_em(data, max_iters = 1000, epsilon = epsilon_em, initial_param_vec = NULL)
        #run SEM
        V <- sem(data, param_df6, tol = sqrt(epsilon_em))

        ###Calculating the complete-data variance covariance matrix 
        data_c <- data[!is.na(data[,2]), ]
        cov_c <- cov(data_c)

        n_data <- length(data_c[,1])

        # Preparation of c
        c <- -cov_c[1,2]^6+3*cov_c[1,1]*cov_c[2,2]*cov_c[1,2]^4-3*(cov_c[1,1]*cov_c[2,2]*cov_c[1,2])^2+(cov_c[1,1]*cov_c[2,2])^3

        # Calculate G1 of I_o^-1
        G1_11 <- cov_c[1,1]
        G1_12 <- 0
        G1_22 <- 2*cov_c[1,1]^2
        G1 <- (1/n_data)*matrix(c(G1_11, G1_12, G1_12, G1_22), nrow = 2, ncol = 2)

        # Calculate G2 of I_o^-1
        G2_11 <- cov_c[1,2]
        G2_12 <- 0
        G2_13 <- 0
        G2_21 <- 0
        G2_22 <- 2*cov_c[1,1]*cov_c[2,2]
        G2_23 <- 2*det(cov_c)^2*(cov_c[1,1]*cov_c[2,2]*cov_c[1,2]^2 - cov_c[1,2]^4)/c
        G2 <- (1/n_data)*matrix(c(G2_11, G2_12, G2_13, G2_21, G2_22, G2_23),
                   nrow = 2, ncol = 3, byrow = TRUE)

        # Calculate G3 of I_o^-1
        G3_11 <- cov_c[2,2]
        G3_12 <- 0
        G3_13 <- 0
        G3_22 <- det(cov_c)^2*((cov_c[1,1]*cov_c[2,2])^2 - cov_c[1,2]^4)/c
        G3_23 <- 2*cov_c[2,2]*cov_c[1,2]
        G3_33 <- 2*cov_c[2,2]^2
        G3 <- (1/n_data)*matrix(c(G3_11, G3_12, G3_13, G3_12, G3_22, G3_23,
                     G3_13, G3_23, G3_33), nrow = 3, ncol = 3)

        # Complete-data VCOV matrix = Inverse of observed Information matrix I_o
        V_c_helper[[n]] = matrix(c(G1, G2, t(G2), G3), ncol = 5)


        # Calculate difference between V and V_c 
        diff_helper[[n]] <- V - V_c_helper[[n]]; diff_helper[[n]]
    }
    #V_c[[p]] <- list(V_c_helper[[n]])
    #diff[[p]] <- list(diff_helper[[n]])
}

