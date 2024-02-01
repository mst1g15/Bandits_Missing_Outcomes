##19 Jan 2023: 



handle_missing <- function(y_t, arm_t, p0_t, p1_t, approach, impute_after_first=FALSE, ...){
  

    #if no missing values, just do complete cases. 
    if(sum(is.na(y_t))==0 ){
      approach="cc"
    }
  
    
    #meanimpute 
    if(approach=="mean_impute_current" | approach=="mean_impute"){
      
      if((impute_after_first==TRUE & is.na(p0_t)) | (impute_after_first==TRUE & is.na(p1_t))){  #if there is not yet an estimate, just do complete case. 
        approach="cc"
      }else{ 
  
      if(is.na(p0_t)){  #if there is not yet an estimate, just do complete case. 
       p0_t=0.5
      }
        
      if(is.na(p1_t)){  #if there is not yet an estimate, just do complete case. 
          p1_t=0.5
      }
        
        new_y <- y_t
        
        
        num_na_arm1 <- length(new_y[is.na(new_y)& arm_t==1])   #number of missing values arm 1
        new_y[is.na(new_y)& arm_t==1] <- rbinom(num_na_arm1, 1, p1_t) #impute with mean 
        
        num_na_arm0 <- length(new_y[is.na(new_y)& arm_t==0])  #number of missing values arm0
        new_y[is.na(new_y)& arm_t==0] <- rbinom(num_na_arm0, 1, p0_t)  #impute with mean 
        
        model <- glm(new_y~arm_t, family = "binomial")
        new_p0 <- 1/(1+exp(-model$coefficients[1]))
        new_p1 <- 1/(1+exp(-sum(model$coefficients)))
        
        pstar <- sum(arm_t==1) / length(arm_t)
        
        s_obs0 = sum(new_y[arm_t==0]==1)
        f_obs0 = sum(new_y[arm_t==0]==0)
        s_obs1 = sum(new_y[arm_t==1]==1)
        f_obs1 = sum(new_y[arm_t==1]==0)
      
        if(approach=="mean_impute"){
          new_y <- y_t
        }
      }
    }
  
   
    
    if(approach=="MI"){
      
      if (sum(is.na(y_t))==length(y_t) ){
        new_y <- y_t
        new_p0=p0_t
        new_p1=p1_t
        pstar=NA
        
        s_obs0 = f_obs0 = s_obs1 = f_obs1 = 0
        
      }else if (sum(arm_t==1)==length(arm_t) | sum(arm_t==0)==length(arm_t)){
        approach="cc"
       
        }else{ error_handle <- try(mice(cbind(y_t, arm_t), m=1, print=F))
        
          if(length(error_handle)==1){
            approach="cc"
          }else{
          
          impute_step <- mice(cbind(y_t, arm_t), m=10, print=F)
          #cat("y_t", y_t, "\n", "arm_t=", arm_t)
          fit_model <- with(impute_step, glm(y_t~arm_t, family = "binomial"))
          
          pooled_model <- summary(pool(fit_model))
        
          new_p0 <- 1/(1+exp(-pooled_model$estimate[1]))
          new_p1 <- 1/(1+exp(-sum(pooled_model$estimate)))
          
          cat("MI estimates", new_p0, new_p1, "\n")

                  
          pstar <- sum(arm_t==1) / length(arm_t)
          
          
          
          s_obs0 = new_p0 * sum(arm_t==0)
          f_obs0 = (1-new_p0) * sum(arm_t==0)
          s_obs1 = new_p1 * sum(arm_t==1)
          f_obs1 = (1-new_p1) * sum(arm_t==1)
          
          
          new_y <-  y_t
          
          
          }
      }
    }
  
    if(approach=="MI_covariate"){
      
      if(sum(is.na(y_t))==length(y_t)){
        new_y <- y_t
        new_p0=p0_t
        new_p1=p1_t
        pstar=NA
        
        s_obs0 = f_obs0 = s_obs1 = f_obs1 = 0
        
      }else if (sum(arm_t==1)==length(arm_t) | sum(arm_t==0)==length(arm_t)){
        approach="cc"
        
      }else{
        #check you don't have collinearity in either of these 
        error_handle1 <- try(mice(cbind(y_t, arm_t, ... ), m=1, print=F))
        error_handle2 <- try(mice(cbind(y_t, arm_t), m=1, print=F))
        
        #sometimes, you can still get collinearity that slips through 
        #check through logged events
        if(length(error_handle1)!=1){
          MI <- mice(cbind(y_t, arm_t, ... ), m=1, print=F)
          check <- "collinear" %in% MI$loggedEvents$meth
        }else{
          check=FALSE
        }
        
        
          if(length(error_handle1)==1 | length(error_handle2)==1 | check ){
            approach="cc"
          }else{
            
          #cat("y_t=", y_t, "\n", "arm_t=", arm_t, "\n", "x=", ... , "\n")
         
          impute_step <- mice(cbind(y_t, arm_t, ...), m=10, print=F)
            
          fit_model <- with(impute_step, glm(y_t~arm_t, family = "binomial"))
          
          pooled_model <- summary(pool(fit_model))

          new_p0 <- 1/(1+exp(-pooled_model$estimate[1]))
          new_p1 <-  1/(1+exp(-sum(pooled_model$estimate)))
          

          pstar <- sum(arm_t==1) / length(arm_t)
          
          
          
          
          s_obs0 = new_p0 * sum(arm_t==0)
          f_obs0 = (1-new_p0) * sum(arm_t==0)
          s_obs1 = new_p1 * sum(arm_t==1)
          f_obs1 = (1-new_p1) * sum(arm_t==1)
          
          
          new_y <-  y_t
        }
      }
    }
  
      
      if(approach=="impute_zero"){
        
        new_y <- y_t
        new_y[is.na(new_y)] <- 0
        
        model <- glm(new_y~arm_t, family = "binomial")
        new_p0 <- 1/(1+exp(-model$coefficients[1]))
        new_p1 <- 1/(1+exp(-sum(model$coefficients)))
        
        pstar <- sum(arm_t==1) / length(arm_t)
        
        s_obs0 = sum(new_y[arm_t==0]==1, na.rm=T)
        f_obs0 = sum(new_y[arm_t==0]==0, na.rm=T)
        s_obs1 = sum(new_y[arm_t==1]==1, na.rm=T)
        f_obs1 = sum(new_y[arm_t==1]==0, na.rm=T)
        
        
      }
  

      if(approach=="cc"){
        
        if(sum(is.na(y_t))==length(y_t)){
          new_y <- y_t
          new_p0=p0_t
          new_p1=p1_t
          pstar = NA
          
          
          s_obs0 = f_obs0 = s_obs1 = f_obs1 = 0
          
          
        }else{
          #complete case analysis - no handling of missing data.  
          model <- glm(y_t~arm_t, family = "binomial")
          new_y <- y_t
          new_p0 <- 1/(1+exp(-model$coefficients[1]))
          new_p1 <- 1/(1+exp(-sum(model$coefficients)))
          
          
          s_obs0 = sum(new_y[arm_t==0]==1, na.rm=T)
          f_obs0 = sum(new_y[arm_t==0]==0, na.rm=T)
          s_obs1 = sum(new_y[arm_t==1]==1, na.rm=T)
          f_obs1 = sum(new_y[arm_t==1]==0, na.rm=T)
          
          pstar = sum(arm_t==1) / length(arm_t)
        }
      }
      
  
return(list(y=new_y, p0=as.numeric(new_p0), 
                     p1=as.numeric(new_p1), 
            pstar=pstar, 
            s_obs0=s_obs0, s_obs1=s_obs1, 
            f_obs0=f_obs0, f_obs1=f_obs1))

}




