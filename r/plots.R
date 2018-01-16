library(ggplot2)
library(dplyr)
library(reshape2)


heatmap_custom = function(C, param_names = NULL, title = NULL, legend_title = NULL){
    
    C = as.matrix(C)
    
    if(is.null(param_names)){
        rownames(C) = colnames(C) = paste("X", 1:nrow(C))
    } else {
        rownames(C) = colnames(C) = param_names
    }
    
    C_melted <- melt(C)
    
    
    g1  = ggplot(data = C_melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile() + 
        scale_fill_continuous(low = "#fef0d9", high = "#a50f15", name = "") + 
        #guides(fill=guide_legend(title=legend_title)) +
        ggtitle(title) + 
        theme_bw()  + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = "black"), 
              plot.title = element_text(face="bold", hjust = .5, size = 15),
              panel.border = element_rect(colour = "black", fill=NA, size=.3), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x = element_text(size = 15, face = "bold"),
              axis.text.y = element_text(size = 15, face = "bold")) 
        
    return(g1)
}

m = matrix(rnorm(9), ncol = 3)

heatmap_custom(m, c("Mu", "Sigma11", "Sigma12"), title = "Plot Title")
