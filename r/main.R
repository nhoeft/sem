###### Simulation study ######

# Packages
library(ggplot2)
library(reshape)
library(MASS)

# Source scripts
source("r/data_simulation.R")
source("r/utils.R")
source("r/em_function.R")
source("r/sem_final.R")
source("r/bootstrapping.R")

# Simulate incomplete data
set.seed((1234))
data <- simulate_data(1000, missings = 0.3,  mu = c(12, 52), sigma= matrix(c(36,35,35,64),2,2))[[1]]

# Run EM and SEM
# Set stopping criterion
epsilon_em <- 0.000000001

# Calculate parameter estimates by EM
param_df <- estimate_em(data, max_iters = 1000, epsilon = epsilon_em, initial_param_vec = NULL)
param_vec <- unlist(param_df[nrow(param_df), ])

# Calculate variance-covariance matrix by SEM (inflated!)
var_sem <- sem(data, param_vec, tol = sqrt(epsilon_em))

# Calculate bootstrapped variance-covariance matrix
var_boots <- param_variance_boot(data)

# Plot of EM convergence 
length_plot <- 1:round(0.25*nrow(param_df))
params_plot <- param_df[length_plot,]
colnames(params_plot) <- c("mu1", "sigma1", "mu2", "sigma2", "rho")

ggplot(data = params_plot) +
    geom_line(aes(x = length_plot, y = mu1, colour = "mu1"), size = 1.5) +
    geom_line(aes(x = length_plot, y = sigma1, colour = "sigma1"), size = 1.5) +
    geom_line(aes(x = length_plot, y = mu2, colour = "mu2"), size = 1.5) +
    geom_line(aes(x = length_plot, y = sigma2, colour = "sigma2"), size = 1.5) +
    geom_line(aes(x = length_plot, y = rho, colour = "rho"), size = 1.5) +
    labs(x = "iterations (n/4)", y = "size") +
    scale_color_discrete(name = "parameters") +
    theme(text = element_text(size=20))

# Preparation of plot for VCOV matrix
# Bootstrapped, EM imputed and true VCOV matrix
bootstrapped = diag(var_boots)
imputed = sem(data, param_vec, tol = sqrt(epsilon_em), return_unadjusted = TRUE)
# obtain consistent estimate of real variance by Monte Carlo
true = param_variance_real(n_runs = 100000, epsilon = 0.0001, n = 100, 
                           missings = 0.3,  mu = c(12, 52), 
                           sigma= matrix(c(36,35,35,81),2,2)) 

# Plot of VCOV matrix evaluation
Names=c("mu 1","sigma 1", "mu 2", "sigma 2", "rho" ) 
df = data.frame(cbind(bootstrapped, imputed, true), Names)

df_m <- melt(df, id.vars='Names')
colnames(df_m) <- c("Names", "procedure", "value")

ggplot(df_m, aes(Names, value)) +   
    geom_point(aes(color = procedure, shape = procedure), position = "dodge", stat = "identity", size = 5) +
    labs(x = "parameters", y = "variance") +
    guides(fill = guide_legend(title="type")) +
    theme(text = element_text(size=20))

