# Function to simulate data.
# Generates missing at random pattern (MNAR) in  variable Y1 and Y2
# n = number of observations, missings = proportion of missing values in variable Y1 and in Y2 each,
# mu = vector of means, sigma = VCOV matrix.



simulate_data = function(n, missings = 0.0, mu = c(0, 0), sigma= matrix(c(1,2,2,1),2,2), seed = 12345)
{
    # Load function and set seed
    require(MASS) 
    set.seed(seed)
    
    # Simulate bivariate normaldistibuted data
    data = mvrnorm(n = n, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    colnames(data) = c("Y1", "Y2")
    data_without_miss = data
    
    # Simulating missing data in Y2 (MAR)
    n_miss_y2 = ceiling(n * missings)
    probabilities = ifelse(data[, "Y2"] < mu[2], 0, data[, "Y2"] / max(data[, "Y2"]))
    dropouts = sample(1:length(data[,"Y2"]), size = n_miss_y2, replace = FALSE, prob = probabilities)
    data[dropouts , "Y2"] = NA
    
    return(list(data, data_without_miss))
}
