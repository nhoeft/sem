### EM Algorithm for bivariate normal data with missing values ###


source("r/utils.R")

# Note: EM function can be used for bivariate normal data with missing in one column
# initial_param_vec mustbe of the following form: (mu1, sig1_squared, mu2, sigma2_squared, rho)

# EM function for multivariate normal data 
estimate_em <- function (X, max_iters = 1000, epsilon = 0.0001, initial_param_vec = NULL){
    
    
    # prepare df to capture parameter estimates at each step
    param_df <- as.data.frame(matrix(1:5, nrow = 1))
    
    # initialize stopping variable and iteration counter
    continue_iterating = TRUE
    iter = 0
    
    n <- nrow(X)
    p <- ncol(X)
    
    # If not given explicitly to the function, we use sample means / covariances 
    # excluding the cases with missing for the initial estimates 
    
    if (is.null(initial_param_vec)){
        mu <- colMeans(X,na.rm=TRUE)
        variance <- diag(apply(X,2,var,na.rm=TRUE)) # compute variances with all non-na values
        sigma <- var(X, na.rm = TRUE) # compute covariance matrix with complete cases
        diag(sigma) = diag(variance) # use the variances that are computed with all non-na values
        
        
    } else{
        
        # extract values from starting parameter vector
        mu1 = initial_param_vec[1]
        mu2 = initial_param_vec[3]
        sig1_sq = initial_param_vec[2]
        sig2_sq = initial_param_vec[4]
        cov = sqrt(sig1_sq)*sqrt(sig2_sq)*initial_param_vec[5]
        
        # here we convert the parameter vector s.t. it has 6 elements, i.e. we replace
        # the correlation coefficient rho by the covariance (twice)
        initial_param_vec = c(mu1, mu2, sig1_sq, cov, cov, sig2_sq)
        
        # convert parameter vector into a mean vector and a covariance matrix
        mu <- param_vec_to_list(initial_param_vec)[[1]]
        sigma <- param_vec_to_list(initial_param_vec)[[2]]
        
    }
    
    
    
    # Update the parameter estimates with iterations of the EM algorithm.
    
    while (continue_iterating == TRUE){
        
        # E Step:
        #
        # Impute missing values in X by their conditional mean given the values 
        # in the complete column and current parameter estimates
        
        iter = iter + 1
        
        # set up data matrix for EM imputation
        X_imputed <- X
        
        for (i in 1:n){
            
            # select i-th row of the data
            x <- X[i,]
            
            # na is a boolean variable specifying for each of the values in x
            # whether it is missing
            na <- is.na(x)
            
            # only update the rows that contain a missing value
            if (any(na)){
                
                # take column that does not contain missings and invert variance
                sigma_inv <- solve(sigma[!na,!na]) 
                
                cov_vec <- sigma[na,!na]
                
                # replace missing value by conditional mean
                X_imputed[i,na] <- mu[na] + cov_vec %*% sigma_inv %*% (x[!na] - mu[!na])
            }
        }
        
        # M Step:
        #
        # Update the parameters by estimating them from the data with the imputed 
        # missing values at the current iteration
        
        mu_new <- colMeans(X_imputed)
        sigma_new <- var(X_imputed)
        
        # stopping criterion: break if EM has converged or if maximum number of 
        # iterations is reached
        if (iter < max_iters){
            continue_iterating <- any(abs(sigma-sigma_new) > epsilon) | any(abs(mu-mu_new) > epsilon)
        } else{
            continue_iterating <- FALSE
        }
        
        # update parameters (assign to to the variables of the previous iterations)
        mu <- mu_new
        sigma <- sigma_new
        
        # compute components for the parameter vector in the above mentioned form:
        # (mu1, sig1_squared, mu2, sigma2_squared, rho)

        mu1 = mu[1]
        mu2 = mu[2]
        sig1_sq = sigma[1, 1]
        sig2_sq = sigma[2, 2]
        rho = sigma[1, 2] / (sqrt(sig1_sq) * sqrt(sig2_sq))
        
        # store the parameter estimates at the current iteration
        param_df[iter, ] <- c(mu1, sig1_sq, mu2, sig2_sq, rho)
    }
    
    return(param_df)
}








