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