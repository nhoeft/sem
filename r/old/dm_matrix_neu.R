
compute_DM2 = function(data, param_df6, tol = 0.0001){
    
    
    
    #data = simulate_data(1000, missings = 0.2,  mu = c(1, 2), sigma= matrix(c(1,.5,.5,1),2,2))
    
    
    theta_final = unlist(param_df6[nrow(param_df6), ])
    
    theta_t = theta_final + c(0.2, 0.3, -0.3, 0.2, 0.2, -0.5)
    
    
    t = 0
    use_params = c(3,4,5)
    
    DM = matrix(rep(0, times = 9), nrow = 3, ncol = 3)
    
    repeat{
        
        DM_old = DM
        
        
        for(i in 1:3){
            
            for(j in 1:3){
                
                theta_t_i = param_6_to_5(theta_final)
                theta_t_i[use_params[i]] = param_6_to_5(theta_t)[use_params[i]]
                
                theta_t_i = param_5_to_6(theta_t_i)
                
                theta_t_i_tilde = unlist(norm_em(data, max_iters = 1, epsilon = 0.00001, initial_param_vec = theta_t_i))
                theta_t_i_tilde = param_6_to_5(theta_t_i_tilde)
                
                DM[i, j] = (theta_t_i_tilde[[use_params[j]]] - param_6_to_5(theta_final)[use_params[j]]) / (param_6_to_5(theta_t)[use_params[i]] -  param_6_to_5(theta_final)[use_params[i]])
                
            }
            
        }
        
        #print(DM-DM_old)
        
        #theta_t = param_5_to_6(theta_t)
        theta_t = unlist(norm_em(data, max_iters = 1, epsilon = 0.00001, initial_param_vec = theta_t))
        #theta_t = param_6_to_5(theta_t)
        
        t = t+1
        
        if(all(abs(DM - DM_old) < tol)){
            break()
        }
        
        
    }
    
    
    
    
    
    return(DM)
}
