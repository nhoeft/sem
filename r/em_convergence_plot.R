#####Convergence of parameters plot
library("ggplot2")

#prepare em-parameter-estimates for plot
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
    scale_color_discrete(name = "parameters") 


