# Function to simulate data.
# Generates missing at random pattern (MNAR) in  variable Y1 and Y2
# n = number of observations, missings = proportion of missing values in variable Y1 and in Y2 each,
# mu = vector of means, sigma = VCOV matrix.



simulate_data = function(n, missings = 0.0, mu = c(0, 0), sigma= matrix(c(1,2,2,1),2,2), seed = 12345)
{
    set.seed(seed)
    require(MASS)
    # Simulate bivariate normaldistibuted data
    data = mvrnorm(n = n, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    colnames(data) = c("Y1", "Y2")
    
    # Simulating missing data in Y2 (MAR)
    n_miss_y2 = ceiling(n * missings)
    probabilities = ifelse(data[, "Y2"] < mu[2], 0, data[, "Y2"] / max(data[, "Y2"]))
    dropouts = sample(1:length(data[,"Y2"]), size = n_miss_y2, replace = FALSE, prob = probabilities)
    data[dropouts , "Y2"] = NA
    
    # Simulating missing data in Y1 (MAR)
    #n_miss_y1 = ceiling(n * missings)
    #probabilities = ifelse(data[, "Y1"] < mu[1], 0, data[, "Y1"] / max(data[, "Y1"]))
    #probabilities = ifelse(is.na(data[, "Y2"]) == TRUE, 0, probabilities)
    #dropouts = sample(1:length(data[,"Y1"]), size = n_miss_y1, replace = FALSE, prob = probabilities)
    #data[dropouts , "Y1"] = NA
    
    return(data)
}

# Test

data = simulate_data(1000, missings = 0.2,  mu = c(1, 2), sigma= matrix(c(1,.5,.5,1),2,2))
data_complete = data[complete.cases(data),]
data[order(data[, "Y1"]),] # The larger the positive values in Y1, the more likely they were dropped
data[order(data[, "Y2"]),] # The larger the positive values in Y2, the more likely they were dropped

