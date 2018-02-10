##### Implementation of the SEM-Algorithm #####

# Source funktions needed
source("r/em_function.R")
source("r/utils.R")
source("r/data_simulation.R")


# Function for comuputing the DM matrix
compute_DM = function(data, theta_final, tol = 0.0001){
    
    theta_t = theta_final + c(0.2, 0.3, -0.3, 0.2, 0.2) # Start theta near the ML estiamtor
    use_params = c(3,4,5) # indexes of the contaminated parameter

    DM = matrix(rep(0, times = 9), nrow = 3, ncol = 3) # initialise matrix for dm values
    t = 0 # initialise t variable
    repeat{
        
        DM_old = DM
        
        # compute the r_ij
        for(i in 1:3){
            for(j in 1:3){
                
                theta_t_i = theta_final 
                theta_t_i[use_params[i]] = theta_t[use_params[i]] # replace the parameter with the estimated parameter from the one stem em
                theta_t_i_tilde = unlist(estimate_em(data, max_iters = 1, epsilon = 0.00001, initial_param_vec = theta_t_i)) # estimate parameters with initial parameter theta_t_i
                DM[i, j] = (theta_t_i_tilde[[use_params[j]]] - theta_final[use_params[j]]) / (theta_t[use_params[i]] -  theta_final[use_params[i]]) # compute r_ij and insert into the DM matrix
            }
        }
        
        theta_t = unlist(estimate_em(data, max_iters = 1, epsilon = 0.00001, initial_param_vec = theta_t)) # get next stem of the em algorithm
        t = t+1
        
        if(all(abs(DM - DM_old) < tol)){ # stopping criteria for convergence
            break()
        }
    }
    return(DM)
}

# Function for computing the parameter Covariances with the SEM-algorithm
sem <- function(X, param_vec, tol) {
    
    # Get n as length of the data
    n <- length(X[,1])
    
    # Extract parameters from parameter vector
    sig1 = sqrt(param_vec[2])
    sig2 = sqrt(param_vec[4])
    rho = param_vec[5]
    
    # compute the DM matrix
    DM_star <- compute_DM(data = X, theta_final = param_vec, tol = 0.0001)

    # Compute the components of G
    
    # Compute G11 from the I_oc
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
    
    #setup 5x5 matrix
    I_oc = matrix(1:25, ncol = 5)
    delta_V = matrix(rep(0, times = 25), ncol = 5)
    delta_V[3:5, 3:5] = DV_22
    
    I_oc[1:2, 1:2] = G11
    
    I_oc[3:5, 1:2] = G21
    
    I_oc[1:2, 3:5] = G12
    
    I_oc[3:5, 3:5] = G22
    
    I_oc_inv = solve(I_oc)
    
    I_inv_corrected = I_oc_inv + delta_V # correct Informationmatrix with delta_v
    
    return(I_inv_corrected)
}
