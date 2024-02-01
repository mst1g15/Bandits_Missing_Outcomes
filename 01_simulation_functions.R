
allocate_missing <- function(x, k, y, alpha){
  
  des_vec <- c(1, x, k, y) #one row of design matrix 
  f <- alpha %*% des_vec  #linear predictor
  indic <- rbinom(1, 1, 1/(1+exp(-f)))
  
  return(indic)
}

generate_response <- function(x, k, beta){
  
  des_vec <- c(1, x, k, x*k) #one row of design matrix 
  f <- beta %*% des_vec  #linear predictor
  indic <- rbinom(1, 1, 1/(1+exp(-f)))
  
  return(indic)
}


compute_prob <- function(x, k, beta){
  
  des_vec <- c(1, x, k, x*k) #one row of design matrix 
  f <- beta %*% des_vec  #linear predictor
  prob= 1/(1+exp(-f))
  
  return(prob)
}

comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY=FALSE)
}


run_sim <- function(alpha, beta, N, Nt, method_name, missing_approach, setting_i){
  
  #if(missing_approach =="mean_impute_restrict"){
  #  restrict <- TRUE
  #  missing_approach <- "mean_impute"
  #}else if (missing_approach=="MI_restrict"){
  #  restrict <- TRUE
  #  missing_approach <- "MI"
  #}else if (missing_approach =="MI_covariate_restrict"){
  #  restrict <- TRUE
  #  missing_approach <- "MI_covariate"
  #}else{
  #  restrict <- FALSE
  #}
  
  
  
  #Priors for Bayesian methods 
  prior0 <- c(1,1)
  prior1 <- c(1,1)
  #N=2
  #Nt=2
  #method_name="FR"
  
  #record responses 
  y <- matrix(0, N, Nt)
  pi <- matrix(0, N, Nt)
  #record responses and retain missingness (for mean impute, imputed values go to y)
  y_obs <- matrix(0, N, Nt)   
  #record missing observations in each arm 
  N_miss0 = N_miss1 = matrix(0, N, Nt)
  
  #number of observatinos in each arm 
  N_obs0 = N_obs1 = matrix(0, N, Nt)
  
  #pstar (unconditional) and pstar in the observed
  pstar = pstar_obs = matrix(0, N, Nt)
  
  #record estimated probability of success and bias
  p_0 = p_1 = bias_0 = bias_1 = matrix(NA, N, Nt)
  
  #record observed number of successes
  ONS = matrix(0, N, Nt)
  
  #record assigned arm 
  arm_choice = matrix(0, N, Nt)
  
  obs_vec <- rep(0, Nt)
  
  for (i in 1:N){
    
    s_obs0=s_obs1 = 0
    f_obs0=f_obs1 = 0
    
    for (t in 1:Nt){
     # print(t)
      #this is for GI in Python code 
      seed=sample(1:100000, 1)
      
      #ratio = N_obs0[i, (t-1)]/N_obs1[i, (t-1)] 
      #print(ratio)
      
      #if (length(ratio)==0){
      #  ratio = 0.5
      #}else if(is.infinite(ratio) | is.nan(ratio)){
      #  ratio=0.5
      #}
      
      #if(restrict==TRUE & (t>1) & ratio < 0.1){
      #    arm=1
      #  }else if(restrict==TRUE & (t>1) & ratio > 0.9){
      #    arm=0
      #  }else{
        
      
        #compute (in Python) the index for arm 0
        #k=as.integer(0)
        sk0=as.integer(prior0[1])
        skt= as.integer(s_obs0)
        fk0=as.integer(prior0[2])
        fkt=as.integer(f_obs0)
        #reticulate::source_python("indices_output.py")
        #cat("t=", t, sk0, skt, fk0, fkt, "\n")
        index0 <-match.fun(method_name)(sk0, skt, fk0, fkt, Nt, t, K)
        #compute the index for arm1
        #k=as.integer(1)
        sk0=as.integer(prior1[1])
        skt=as.integer(s_obs1)
        fk0=as.integer(prior1[2])
        fkt=as.integer(f_obs1)
        #cat(sk0, skt, fk0, fkt, "\n")
        
        index1 <-match.fun(method_name)(sk0, skt, fk0, fkt, Nt, t, K)
        
        pi[i, t]=index1/(index0+index1)
        
        #select best arm, randomize in case of a tie 
        
        arm = which(c(index0, index1) == max(c(index0, index1)))-1
        if(length(arm)==2){
          arm=sample(c(0,1), 1)
        }
        
       # }
      #add arm choice
      arm_choice[i, t]=arm
      
      # generate response
      res = generate_response(x=X[t], k=arm, beta)
      
      #add response 
      y[i,t] = res
      
      # generate missingness
      if (length(alpha)==1){
        miss=0
      }else{
        miss=allocate_missing(X[t], arm, res, alpha)
      }
      
      #add response 
      y[i,t] = ifelse(miss==0, res, NA)
      y_obs[i,t] = ifelse(miss==0, res, NA)
      
      
      obs_vec[t] <- ifelse(miss==0, 1, 0)
      

      #add total observed number of successes before any imputation happens
      ONS[i, t] <- sum(y[i,1:t], na.rm=T)
      
      
        if(t==1){
          p0_t=NA
          p1_t=NA
        }else{
          p0_t=p_0[i,(t-1)]
          p1_t=p_1[i,(t-1)]
        }
      
      
      if(missing_approach=="MI_covariate"){
      
        est <- handle_missing(y_t=y[i,1:t], arm_t=arm_choice[i,1:t], 
                              p0_t, p1_t, approach=missing_approach, impute_after_first=FALSE, X[1:t])
      }else{
        est <- handle_missing(y_t=y[i,1:t], arm_t=arm_choice[i,1:t], 
                              p0_t, p1_t, approach=missing_approach, impute_after_first=FALSE)
      }
      
      
      
      

      #response after missingness approach applied 
      y[i,1:t]  <- est$y
      p_0[i,t] = est$p0 
      bias_0[i,t] = est$p0 - compute_prob(X[t], 0, beta)
      p_1[i,t] =  est$p1
      bias_1[i,t]=est$p1 - compute_prob(X[t], 1, beta)
      pstar[i,t] <- est$pstar
      
      s_obs0 = as.integer(est$s_obs0) 
      f_obs0= as.integer(est$f_obs0) 
      s_obs1 = as.integer(est$s_obs1) 
      f_obs1 = as.integer(est$f_obs1) 
      
      
      #unconditional pstar 

      N_obs0[i,t] <- sum(arm_choice[i,1:t]==0 & !is.na(y[i,1:t]))
      N_obs1[i,t] <- sum(arm_choice[i,1:t]==1 & !is.na(y[i,1:t]))
      
      pstar_obs[i,t] <- sum(arm_choice[i,1:t]==1 & obs_vec[1:t]==1) / sum(obs_vec[1:t]==1)
      
      
    }
    
  }
  
  
  
  
  bias0_mean <- colMeans(bias_0, na.rm=T)
  bias0_sd <- colSds(bias_0, na.rm=T)
  
  bias1_mean <- colMeans(bias_1, na.rm=T)
  bias1_sd <- colSds(bias_1, na.rm=T)
  
  
  N0_300 <- N_obs0[,Nt]
  p0_300 <- p_0[, Nt]
  
  N1_300 <- N_obs1[,Nt]
  p1_300 <- p_1[, Nt]
  
  
  #calculate expected bias according to Bowden and Trippa's formula 
  expbias_p0 <- rep(0, Nt)
  expbias_p1 <- rep(0, Nt)
  
  for (t in 1:Nt){
    expbias_p0[t] <- -cov(p_0[,t], N_obs0[,t])/mean(N_obs0[,t], na.rm=T)
    expbias_p1[t] <- -cov(p_1[,t], N_obs1[,t])/mean(N_obs1[,t], na.rm=T)
    
  }
  
  N1_mean<- colMeans(t(apply(arm_choice, 1, cumsum)), na.rm=T)
  N1_sd <- colSds(t(apply(arm_choice, 1, cumsum)), na.rm=T)
   
  pstar_mean <- colMeans(pstar, na.rm=T)
  pstar_sd <- colSds(pstar, na.rm=T)
  pstar_obs_mean <-  colMeans(pstar_obs, na.rm=T)
  pstar_obs_sd <- colSds(pstar_obs, na.rm=T)
  
  ONS_mean<- colMeans(ONS, na.rm=T)
  ONS_sd <- colSds(ONS, na.rm=T)
  
  results_all <- 
    data.frame(method=method_name, beta_setting=as.character(setting_i$beta), 
               alpha_setting=as.character(setting_i$alpha),
               missing_approach=missing_approach, 
               sample_size=1:(Nt),  
               N1_mean=N1_mean, N1_sd=N1_sd,
               pstar_mean=pstar_mean, pstar_sd=pstar_sd,
               pstar_obs_mean = pstar_obs_mean, pstar_obs_sd = pstar_obs_sd,
               ONS_mean=ONS_mean, ONS_sd=ONS_sd, 
               bias0_mean = bias0_mean, bias0_sd = bias0_sd,
               bias1_mean = bias1_mean, bias1_sd = bias1_sd,  
               expbias_p0 = expbias_p0, expbias_p1=expbias_p1)
  
  
  
  
  
  
  #p_300_all <- data.frame(method=method_name, beta_setting=as.character(setting_i$beta), 
  #                        alpha_setting=as.character(setting_i$alpha), missing_approach=missing_approach, 
  #                        p0_300=p0_300, N0_300=N0_300, p1_300=p1_300, N1_300=N1_300)
  
  
  list(results_all=results_all)
  
}




calc_offline_est <- function(y_obs, pi, arm_choice){
  
  
  #cc, weighted
  rowMeans((y_obs==1)*(arm_choice==1)*(1/pi), na.rm=T)
  rowMeans((y_obs==1)*(arm_choice==0)*(1/(1-pi)), na.rm=T)
  
  
  
}

