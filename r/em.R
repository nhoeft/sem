# data simulation

library(MASS)

V_true = matrix(c(1, 0.3, 0.3, 2), ncol = 2)
mu_true = c(1,2)
                
full_X = mvrnorm(100, mu = mu_true, Sigma = V_true)

X = full_X
X[1:20, 2] = NA



# EM function for multivariate normal data 

norm_em <- function (X, iters, debug=0)
{
    if (any(rowSums(is.na(X))>1))
        stop("This function can handle at most one missing value per case")
    
    n <- nrow(X)
    p <- ncol(X)
    
    # Find initial estimates, using simple averages with missing observations
    # omitted.  The initial covariance estimate is diagonal.
    
    mu <- colMeans(X,na.rm=TRUE)
    sigma <- diag(apply(X,2,var,na.rm=TRUE))
    
    if (debug>0) {
        cat("\nINITIAL VALUES\n\n")
        if (debug>1) { cat("X:\n"); print(X) }
        cat("mu:   ",round(mu,3),"\n")
        cat("sigma:\n")
        print(round(sigma,3))
    }
    
    # Update the parameter estimates with iterations of the EM algorithm.
    
    for (t in 1:iters) {
        
        # E Step:
        #
        # Compute a filled in version of X, with missing values set to their
        # conditional mean given other values and current parameter estimates,
        # along with the total extra variance due to uncertainty in the
        # missing values.
        
        filled_in_X <- X
        extra_var <- numeric(p)
        for (i in 1:n) {
            x <- X[i,]
            na <- is.na(x)
            if (any(na)) {
                sigma_inv <- solve(sigma[!na,!na])
                cov_vec <- sigma[na,!na]
                filled_in_X[i,na] <- 
                    mu[na] + cov_vec %*% sigma_inv %*% (x[!na] - mu[!na])
                extra_var[na] <- extra_var[na] + 
                    sigma[na,na] - cov_vec %*% sigma_inv %*% cov_vec
            }
        }
        
        # M Step:
        #
        # Find new parameter estimates from the filled in version of X, by
        # taking simple averages over the filled in data.
        
        mu <- colMeans(filled_in_X)
        Y <- t(filled_in_X) - mu
        sigma <- (Y %*% t(Y) + diag(extra_var)) / n
        
        if (debug>0) {
            cat("\nITERATION",t,"\n\n")
            if (debug>1) { cat("filled in X:\n"); print(round(filled_in_X,3)) }
            cat("mu:   ",round(mu,3),"\n")
            cat("sigma:\n")
            print(round(sigma,3))
            cat("\n")
        }
    }
    
    list(mu=mu,sigma=sigma)
}



# Testing

cat("\nTrue covariance matrix:\n")
print(round(V_true,3))


cat("\nSample covariance of complete data:\n")
print(round(cov(full_X),3))


cat("\nCovariance estimate from pairwise complete observations:\n")
print(round(cov(X,use="pairwise.complete.obs"),3))

cat("\nResult of maximum likelihood estimation with EM:\n")
print(norm_em(X,30,debug=1))

