
create_plot <- function(method_name, alpha, beta, missing_approach, Nt){
  
#Priors for Bayesian methods 
prior0 <- c(1,1)
prior1 <- c(1,1)
N=1
K=2
#method_name="FR"

#record responses 
y <- matrix(0, N, Nt)  #may include imputed
y_obs <- matrix(0, N, Nt)   #observed only, missing values indicated with NA

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

index1_all = index0_all =matrix(0, N, Nt)



X <- c(rep(0, Nt/2), rep(1, Nt/2))
#randomize
X <- sample(X)

i=1

s_obs0=s_obs1 = 0
f_obs0=f_obs1 = 0

for (t in 1:Nt){
  #this is for GI in Python code 
  seed=sample(1:100000, 1)
  
 
    #compute (in Python) the index for arm 0
    #k=as.integer(0)
    sk0=as.integer(prior0[1])
    skt= as.integer(s_obs0)
    fk0=as.integer(prior0[2])
    fkt=as.integer(f_obs0)
    #reticulate::source_python("indices_output.py")
    
    index0 <-match.fun(method_name)(sk0, skt, fk0, fkt, Nt, t, K)
    index0_all[i, t] <- index0
    #compute the index for arm1
    #k=as.integer(1)
    sk0=as.integer(prior1[1])
    skt=as.integer(s_obs1)
    fk0=as.integer(prior1[2])
    fkt=as.integer(f_obs1)
    index1 <-match.fun(method_name)(sk0, skt, fk0, fkt, Nt, t, K)
    index1_all[i, t] <- index1
    
    #select best arm, randomize in case of a tie 
    
    arm = which(c(index0, index1) == max(c(index0, index1)))-1
    if(length(arm)==2){
      arm=sample(c(0,1), 1)
    }
    
  
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
  
  y_obs[i,t] = y[i,t]
  
  #add total observed number of successes before any imputation happens
  ONS[i, t] <- sum(y[i,1:t], na.rm=T)
  
  
  if(t==1){
    p0_t=NA
    p1_t=NA
  }else{
    p0_t=p_0[i,(t-1)]
    p1_t=p_1[i,(t-1)]
  }
  
  
  if(missing_approach=="MI_covariates" | missing_approach=="MI_covariate_current"){
    
    est <- handle_missing(y_t=y[i,1:t], arm_t=arm_choice[i,1:t], 
                          p0_t, p1_t, approach=missing_approach, X[1:t])
  }else{
    est <- handle_missing(y_t=y[i,1:t], arm_t=arm_choice[i,1:t], 
                          p0_t, p1_t, approach=missing_approach)
  }
  
  
  
  
  
  
  #response after missingness approach applied 
  y[i,1:t]  <- est$y
  p_0[i,t] = est$p0 
  bias_0[i,t] = est$p0 - compute_prob(0, 0, beta)
  p_1[i,t] =  est$p1
  bias_1[i,t]=est$p1 - compute_prob(0, 1, beta)
  
  #unconditional pstar 
  sum(arm_choice[1:t])/length(arm_choice[1:t])
  
  N_obs0[i,t] <- sum(arm_choice[i,1:t]==0 * !is.na(y[i,1:t]))
  N_obs1[i,t] <- sum(arm_choice[i,1:t]==1 * !is.na(y[i,1:t]))
  
  pstar_obs[i,t] <- sum(arm_choice[i,1:t]==1 * !is.na(y[i,1:t])) / sum(!is.na(y[i,1:t]))
  pstar[i,t] <-sum(arm_choice[i,1:t]==1) / t
  
  
  if(!is.na(y[i,t])){
    
    if(arm==0){
      s_obs0 = s_obs0 + y[i,t]
      f_obs0= f_obs0 + 1 - y[i,t]
    }else{
      s_obs1 = s_obs1 + y[i,t]
      f_obs1 = f_obs1 + 1 - y[i,t]
    }
    
  }
  
  
}




plot_table <- tibble(index0=as.numeric(index0_all), 
       index1=as.numeric(index1_all),
       arm = as.numeric(arm_choice), 
       y_obs = as.numeric(y_obs),
       y=as.numeric(y), 
       pstar = as.numeric(pstar), 
       pstar_observed=as.numeric(pstar_obs))[1:(Nt-1),]  #remove last row as it may not have a response


plot_table$arm <- as.factor(plot_table$arm)
plot_table$y[is.na(plot_table$y)] <-"Null"
plot_table$y <- as.factor(plot_table$y)
plot_table$y_obs[is.na(plot_table$y_obs)] <-"Null"
plot_table$y_obs <- as.factor(plot_table$y_obs)


index_plot <- ggplot(plot_table) + geom_point(aes(x=1:(Nt-1), y=index0, color="0"), size=2) + 
  geom_point(aes(x=1:(Nt-1), y=index1, color="1"), size=2) +
  scale_color_manual(values=c("blue", "red")) +ylab("Index") + xlab("") + 
  ylim(min(plot_table$index0, plot_table$index1)-0.5, max(plot_table$index0, plot_table$index1)+0.5)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

arm_plot <- ggplot(plot_table) + geom_point(aes(x=1:(Nt-1), y="------", color=arm), size=10, shape=15) +
  scale_shape_manual(values=c(15, 15)) + 
  scale_color_manual(values=c("blue", "red"))+ylab("Arm") + xlab("") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


obs_response_plot <-  ggplot(plot_table) +
  geom_point(aes(x=1:(Nt-1), y="------", color=y_obs, shape=y_obs), size=10)+ 
  scale_shape_manual(values=c(15, 15, 4))+
  scale_color_manual(values=c("lightblue", "pink", "black"))+
  ylab("Observed  Response") + xlab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


response_plot <-  ggplot(plot_table) +
  geom_point(aes(x=1:(Nt-1), y="------", color=y, shape=y), size=10)+ 
  scale_shape_manual(values=c(15, 15, 4))+
  scale_color_manual(values=c("lightblue", "pink", "black"))+
  ylab("Response") + xlab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


pstar_plot <-  ggplot(plot_table) +
  geom_line(aes(x=1:(Nt-1), y=pstar), linewidth=1)+
  geom_line(aes(x=1:(Nt-1), y=pstar_observed), linetype="longdash", linewidth=1)+ 
  ylab("Pstar") + xlab("Sample Size") + ylim(-0.05,1.2)

if(missing_approach=="cc"){
  
  p <- ggarrange(index_plot, NULL, arm_plot, NULL, response_plot, NULL, pstar_plot, 
                 ncol=1,
                 heights=c(1, -0.3, 1, -0.3, 1, -0.3, 1),
                 legend="none")
  
}else{
  p <- ggarrange(index_plot, NULL, arm_plot, NULL, obs_response_plot, NULL, 
                 response_plot, NULL, pstar_plot, 
                 ncol=1,
            heights=c(1, -0.3, 1, -0.3, 1, -0.3, 1, -0.3, 1),
            legend="none")
}

  p

}
