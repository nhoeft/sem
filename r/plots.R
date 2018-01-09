library(ggplot2)
library(dplyr)
library(reshape2)


m = as.matrix(matrix(rnorm(9), ncol = 3))




heatmap_custom <- function(C){
    
    C = as.matrix(C)
    rownames(C) = colnames(C) = paste("X", 1:nrow(C))
    C_melted <- melt(C)
    
    
    g1 <- ggplot(data = C_melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile() + 
        scale_fill_continuous(low = "white", high = "black", name = "") + 
        #guides(fill=guide_legend(title="Correlation")) +
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank()) + 
        ggtitle("Some Title") + 
        theme_bw()  + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = "black"), plot.title = element_text(face="bold"),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
        theme(axis.text.y = element_text(size = 12, hjust = 0))
    
    
    return(g1)
}

heatmap_custom(m)
