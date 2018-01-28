
compute_DM2 = function(data, param_df6, tol = 0.0001){
    
    
    
    #data = simulate_data(1000, missings = 0.2,  mu = c(1, 2), sigma= matrix(c(1,.5,.5,1),2,2))
    
    
    theta_final = unlist(param_df6[nrow(param_df6), ])
    
    theta_t = theta_final #+ c(0.2, 0.3, -0.3, 0.2, 0.2, -0.5) ##???
    
    t = 0
    vecR <- NULL
    use_params = c(3,4,5)
    
    DM = matrix(rep(0, times = 9), nrow = 3, ncol = 3)
    
    repeat{
        
        DM_old = DM
        
        
        for(i in 1:3){
            
            for(j in 1:3){
                # Get current state of theta from em computation
                theta_t <- param_df6[t,]
                
                # Define theta_i as the final mle and replace i-th value with current state
                theta_t_i <- theta_final # evtl. theta_t nehmen
                theta_t_i[use_params[i]] <- theta_t[,use_params[i]]
                #5 to 6 translation
                theta_t_i = param_5_to_6(theta_t_i)
                
                # Do EM step using theta_i: convert parameters to list, then do step
               # paramList <- list(mu = theta_t_i[1:2],
               #                   cov = matrix(theta_t_i[3:6], nrow = 2, ncol = 2)) 
                theta_t_1_i = unlist(norm_em(data, max_iters = 1, epsilon = 0.00001, initial_param_vec = theta_t_i))
                
                #6 to 5 translation
                theta_t_1_i = param_6_to_5(theta_t_1_i)
                
                # Calculate ratio rij
                vecR <- append(vecR, (theta_t_1_i[[j]] - theta_final[[j]])/(theta_t[[i]] - theta_final[[i]]))
                
                DM[i, j] <- (theta_t_1_i[[use_params[j]]] - theta_final[use_params[j]])/(theta_t[use_params[i+1]] - theta_final[use_params[i+1]])
                # # # 
                
                
                #take last theta from last EM-iteration
                #theta_t_i = theta_final
                #theta_t_i[use_params[i]] = theta_t[use_params[i]]
                
                #theta_t_i = param_5_to_6(theta_t_i)
                
                #theta_t_i_tilde = unlist(norm_em(data, max_iters = 1, epsilon = 0.00001, initial_param_vec = theta_t_i))
                #theta_t_i_tilde = param_6_to_5(theta_t_i_tilde)
                
                #DM[i, j] = (theta_t_i_tilde[[use_params[j]]] - theta_final[use_params[j]]) / (theta_t[use_params[i]] -  theta_final[use_params[i]])
                
            }
            
        }
        
        print(DM-DM_old)
        
        theta_t = param_5_to_6(theta_t)
        theta_t = unlist(norm_em(data, max_iters = 1, epsilon = 0.00001, initial_param_vec = theta_t))
        theta_t = param_6_to_5(theta_t)
        
        t = t+1
        
        if(all(abs(DM - DM_old) < tol)){
            break()
        }
        
        
    }
    
    
    
    
    
    return(DM)
}
