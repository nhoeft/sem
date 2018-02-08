# bootstram

library(ggplot2)
library(reshape)
source("r/em_function.R")

source("r/sem_final.R")

source("r/utils.R")

source("r/data_simulation.R")

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

param_variance_real = function(n_runs = 1000, epsilon = 0.0001, n = 100, missings = 0.3,  mu = c(12, 52), sigma= matrix(c(36,35,35,81),2,2)){
    
    # create empty matrix for bootstrap parameters
    param_df = as.data.frame(matrix(ncol = 5, nrow = n_runs))
    
    # compute estimates for bootstrap samples
    for(i in 1:n_runs){
        data_real = simulate_data(n, missings = missings,  mu = mu, sigma= sigma)[[2]]
        param_df[i, ] = unlist(tail(estimate_em(data_real, max_iters = 1000, epsilon = epsilon, initial_param_vec = NULL), 1))
    }
    
    
    return(var(param_df))
}

data = simulate_data(1000, missings = 0.3,  mu = c(12, 52), sigma= matrix(c(36,35,35,81),2,2))[[1]]
print('missing data')
var_miss2= param_variance_boot(data)

epsilon_em = 0.000000001
param_df = estimate_em(data, max_iters = 1000, epsilon = epsilon_em, initial_param_vec = NULL)
param_vec = unlist(param_df[nrow(param_df), ])
var_sem = sem(data, param_vec, tol = sqrt(epsilon_em))

print('complete data')
set.seed(12345)
var_em = param_variance_real(missings = 0.3)
print('real data')
set.seed(12345)
var_real = param_variance_real(missings = 0.0)

#x <- c(8,6,11,22,14,17,18,24,19,23,26,40,4,4,5,6,8,10)
#y <- c(59,58,56,53,50,45,43,42,39,38,30,27,NA,NA,NA,NA,NA,NA)
#data = matrix(c(x,y), ncol = 2)
par_var_boot = c(12.42484, 36.80615, 50.73317, 56.05955, 0.7227770)
par_var_impute = c(12.42484, 36.80615, 50.73317, 56.05955, 0.7227770) * 0.7
par_var_real = c(12.42484, 36.80615, 50.73317, 56.05955, 0.7227770) * 0.97

Names=c("mu 1","sigma 1", "mu 2", "sigma 2", "rho" ) 
df = data.frame(cbind(par_var_boot, par_var_impute, par_var_real), Names)

# melt the data frame for plotting
df_m <- melt(df, id.vars='Names')

# plot everything
ggplot(df_m, aes(Names, value)) +   
    geom_bar(aes(fill = variable), position = "dodge", stat="identity")
#estimate_em(data, max_iters = 1000, epsilon = epsilon, initial_param_vec = NULL)
