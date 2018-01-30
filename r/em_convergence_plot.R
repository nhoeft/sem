#####Convergence of parameters plot
library("ggplot2")

#prepare em-parameter-estimates for plot
length_plot <- 1:round(0.25*nrow(param_df6))
params_plot <- param_df6[length_plot,]

ggplot(data = params_plot) +
    geom_line(aes(x = length_plot, y = V1, colour = "V1")) +
    geom_line(aes(x = length_plot, y = V2, colour = "V2")) +
    geom_line(aes(x = length_plot, y = V3, colour = "V3")) +
    geom_line(aes(x = length_plot, y = V4, colour = "V4")) +
    geom_line(aes(x = length_plot, y = V5, colour = "V5")) +
    geom_line(aes(x = length_plot, y = V6, colour = "V6")) +
    labs(x = "iterations (n/4)", y = "size") +
    scale_color_discrete(name = "parameters") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          panel.background = element_blank())


