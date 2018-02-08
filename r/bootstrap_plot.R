bootstrapped = c(12.42484, 36.80615, 50.73317, 56.05955, 0.7227770)
imputed = c(12.42484, 36.80615, 50.73317, 56.05955, 0.7227770) * 0.7
true = c(12.42484, 36.80615, 50.73317, 56.05955, 0.7227770) * 0.97

Names=c("mu 1","sigma 1", "mu 2", "sigma 2", "rho" ) 
df = data.frame(cbind(bootstrapped, imputed, true), Names)

# melt the data frame for plotting
df_m <- melt(df, id.vars='Names')
colnames(df_m) <- c("Names", "procedure", "value")
# plot everything
ggplot(df_m, aes(Names, value)) +   
    geom_point(aes(color = procedure, shape = procedure), position = "dodge", stat = "identity", size = 5) +
    labs(x = "parameters", y = "variance") +
   guides(fill = guide_legend(title="type")) 
