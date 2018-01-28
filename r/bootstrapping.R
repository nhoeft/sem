# bootstram

param_variance_boot = function(data, n_boots = 1000, epsilon = 0.0001){

    # create empty matrix for bootstrap parameters
    param_df = as.data.frame(matrix(ncol = 5, nrow = n_boots))
    
    # compute estimates for bootstrap samples
    for(i in 1:n_boots){
        data_boot = data[sample(1:nrow(data), nrow(data), replace = TRUE),]
        param_df[i, ] = unlist(tail(estimate_em(data_boot, max_iters = 1000, epsilon = epsilon, initial_param_vec = NULL), 1))
    }
    print(param_df)
    
    return(var(param_df))
}


data = simulate_data(1000, missings = 0.2,  mu = c(1, 2), sigma= matrix(c(1,.5,.5,1),2,2))
param_variance_boot(data)
#estimate_em(data, max_iters = 1000, epsilon = epsilon, initial_param_vec = NULL)