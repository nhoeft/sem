# data simulation

library(MASS)

# variables for simulation
# Todo simulations function von JM einfügen
sigma_true = matrix(c(1, 0.3, 0.3, 2), ncol = 2)
mu_true = c(1,2)
                
X_complete = mvrnorm(1000, mu = mu_true, Sigma = sigma_true)

X = X_complete
X[1:200, 2] = NA



# EM function for multivariate normal data 
norm_em <- function (X, max_iters = 1000, epsilon = 0.0001, initial_param_vec = NULL)
{
    #if (any(rowSums(is.na(X))>1))
     #   stop("This function can handle at most one missing value per case")
    
    # prepare df to capture parameter estimates at each step
    param_df <- as.data.frame(matrix(initial_param_vec, nrow = 1))
    
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
        
        mu <- param_vec_to_list(initial_param_vec)[[1]]
        sigma <- param_vec_to_list(initial_param_vec)[[2]]
    }


    # Update the parameter estimates with iterations of the EM algorithm.
    
    while (continue_iterating == TRUE) { # TODO: implement stopping criterium
        
        # E Step:
        #
        # Compute a filled in version of X, with missing values set to their
        # conditional mean given other values and current parameter estimates,
        # along with the total extra variance due to uncertainty in the
        # missing values.
        
        iter = iter + 1

        X_imputed <- X
        var_adjustment <- numeric(p) # empty vector of dimension 1xp
        for (i in 1:n) {
            x <- X[i,]
            na <- is.na(x)
            if (any(na)) {
                sigma_inv <- solve(sigma[!na,!na]) # take column that does not contain missings and invert variance
                cov_vec <- sigma[na,!na] # get covariance
                X_imputed[i,na] <- mu[na] + cov_vec %*% sigma_inv %*% (x[!na] - mu[!na])
                var_adjustment[na] <- var_adjustment[na] + sigma[na,na] - cov_vec %*% sigma_inv %*% cov_vec
            }
        }
        
        # M Step:
        #
        # Find new parameter estimates from the filled in version of X, by
        # taking simple averages over the filled in data.
        
        mu_new <- colMeans(X_imputed)
        X_imputed_centered <- t(t(X_imputed) - mu)
        sigma_new <- (t(X_imputed_centered) %*% X_imputed_centered + diag(var_adjustment)) / (n - 1)
        #sigma_new <- var(X_imputed)
        
        if (iter < max_iters){
            continue_iterating <- any(abs(sigma-sigma_new) > epsilon) | any(abs(mu-mu_new) > epsilon)
        } else{
            continue_iterating <- FALSE
        }
        
        mu <- mu_new
        sigma <- sigma_new
        
        # change iter to iter+1 if the initial values (i.e. t = 0) should be stored in the df as well
        param_df[iter, ] <- param_list_to_vec(mu, sigma)
    }
    
    return(param_df)
}



# Testing

cat("\nTrue covariance matrix:\n")
print(round(sigma_true,3))


cat("\nSample covariance of complete data:\n")
print(round(cov(X_complete),3))


cat("\nCovariance estimate from pairwise complete observations:\n")
print(round(cov(X,use="pairwise.complete.obs"),3))

cat("\nResult of maximum likelihood estimation with EM:\n")
print(norm_em(X, max_iters = 1000, epsilon = 0.0001)$sigma)



