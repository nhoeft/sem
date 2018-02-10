### utils for sem project ###


# convert param vec of length 6 (with covariances) to param vec of length 5 (with correlation instead)
param_6_to_5 = function(param_vec6){
    
    mu1 = param_vec6[1]
    mu2 = param_vec6[2]
    sig1_sq = param_vec6[3]
    sig2_sq = param_vec6[6]
    rho = param_vec6[4] / (sqrt(sig1_sq)*sqrt(sig2_sq))

    
    # adjust order according to contanimated data (see EM Algorithm and extensions book p. 128)
    param_vec5 = c(mu1, sig1_sq, mu2, sig2_sq, rho)
    
    return(param_vec5)
}



# convert param vec of length 5 to param vec of length 6 
param_5_to_6 = function(param_vec5){
    
    mu1 = param_vec5[1]
    mu2 = param_vec5[3]
    sig1_sq = param_vec5[2]
    sig2_sq = param_vec5[4]
    cov = sqrt(sig1_sq)*sqrt(sig2_sq)*param_vec5[5]
    
    
    # adjust order according to contanimated data (see EM Algorithm and extensions book p. 128)
    param_vec6 = c(mu1, mu2, sig1_sq, cov, cov, sig2_sq)
    
    return(param_vec6)
}





# log / fisher transformation of variances and correlation for higher numerical stability
stabilizing_transformation = function(param_vec5){
    
    rho = param_vec5[5]
    
    # for increased numerical stability use log variances and fisher transformation of rho (see EM Algorith and extensions book p.129)
    z = 0.5*log((1+rho) / (1-rho))
    
    # adjust order according to contanimated data (see EM Algorithm and extensions book p. 128)
    param_vec5_trans = c(param_vec5[1], log(param_vec5[2]), param_vec5[3], log(param_vec5[4]), z)
    
    return(param_vec5_trans)
    
    
}



compute_d = function(param_vec){
    return(- 0.5 + sqrt(1/ 4 + length(param_vec))) # p/q-formula to compute dimension
    
}

param_vec_to_list = function(param_vec){
    #d = compute_d(param_vec)
    
    mu = param_vec[1:2]
    sigma = matrix(param_vec[3:6], nrow = 2 , ncol = 2)
    
    return(list(mu, sigma))
}

param_list_to_vec = function(param_list){
    return(unlist(param_list))
}



