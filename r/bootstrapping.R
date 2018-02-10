##### Functions estimating the asymptotic VCOV matrix of the parameters via bootstrapping #####


# VCOV matrix of incomplete data problem
param_variance_boot = function(data, n_boots = 1000, epsilon = 0.0001){

    # create empty matrix for bootstrap parameters
    param_df = as.data.frame(matrix(ncol = 5, nrow = n_boots))
    
    # compute estimates for bootstrap samples
    for(i in 1:n_boots){
        data_boot = data[sample(1:nrow(data), nrow(data), replace = TRUE),]
        param_df[i, ] = unlist(tail(estimate_em(data_boot, max_iters = 1000, epsilon = epsilon, initial_param_vec = NULL), 1))
    }

    
    return(var(param_df))
}

# VCOV matrix of complete data problem
param_variance_real = function(n_runs = 1000, epsilon = 0.0001, n = 100, missings = 0.3,  mu = c(12, 52), sigma= matrix(c(36,35,35,81),2,2)){
    
    # create empty matrix for bootstrap parameters
    param_df = as.data.frame(matrix(ncol = 5, nrow = n_runs))
    seed = 12345
    
    # compute estimates for bootstrap samples
    for(i in 1:n_runs){
        data_real = simulate_data(n, missings = missings,  mu = mu, sigma= sigma, seed = seed + i)[[2]]
        param_df[i, ] = unlist(tail(estimate_em(data_real, max_iters = 1000, epsilon = epsilon, initial_param_vec = NULL), 1))
    }
    
    
    return(var(param_df))
}

