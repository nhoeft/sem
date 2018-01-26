
# EM function for multivariate normal data 

norm_em <- function (X, iters, debug=0){

    
    n <- nrow(X)
    p <- ncol(X)
    
    # For the initial estimates we use sample means / covariances excluding the cases with missing
    
    mu <- colMeans(X,na.rm=TRUE)
    variance <- diag(apply(X,2,var,na.rm=TRUE)) # compute variances with all non-na values
    sigma <- var(X, na.rm = TRUE) # compute covariance matrix with complete cases
    diag(sigma) = diag(variance) # use the variances that are computed with all non-na values

    # Update the parameter estimates with iterations of the EM algorithm.
    
    for (t in 1:iters) { # TODO: implement stopping criterium
        
        # E Step:
        #
        # Compute a filled in version of X, with missing values set to their
        # conditional mean given other values and current parameter estimates,
        # along with the total extra variance due to uncertainty in the
        # missing values.
        
        X_imputed <- X
        
        for (i in 1:n) {
        
        #var_adjustment <- numeric(p) # empty vector of dimension 1xp
            x <- X[i,]
            na <- is.na(x)
            if (any(na)) {
                sigma_inv <- solve(sigma[!na,!na]) # take column that does not contain missings and invert variance
                cov_vec <- sigma[na,!na] # get covariance
                X_imputed[i,na] <- mu[na] + cov_vec %*% sigma_inv %*% (x[!na] - mu[!na])
                #var_adjustment[na] <- var_adjustment[na] + sigma[na,na] - cov_vec %*% sigma_inv %*% cov_vec
            }
        }
        
        # M Step:
        #
        # Find new parameter estimates from the filled in version of X, by
        # taking simple averages over the filled in data.
        
        mu <- colMeans(X_imputed)
        X_imputed_centered <- t(X_imputed) - mu
        #sigma <- (X_imputed_centered %*% t(X_imputed_centered) + diag(var_adjustment)) / n
        sigma <- var(X_imputed)
    }
    
    list(mu=mu,sigma=sigma)
}



cat("\nResult of maximum likelihood estimation with EM:\n")
print(norm_em(data, 30 ,debug=1))

