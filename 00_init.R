library(reticulate)
library(tidyverse)
library(matrixStats)
library(gridExtra)
library(grid)
library(ggpubr)
library(mice)
library(foreach)
library(doParallel)
library(shiny)


#set up to work with Python----------------------------------------------------------------------------
#setwd("C:/Users/lshmt4/Desktop/plos-one-main/simulation")
#
use_python("C:/Users/Mia/miniconda3")

#load useful information for simulation ------------------------------------------------------------------

Method= c("FR", "RTS", "CB", "UCB", "RGI", "GI", "RBI")

beta_names <- c("Null", "Null_low", "Null_high",
                "Alternative", "Alternative_low","Alternative_high",
                "Coveffect", "Synergistic","Antagonistic")

#beta parameters for model for the response
beta_list <- list("Null" =c(0, 0, 0, 0), #null case - pk=0.5
                  "Null_low" =c(-1.5, 0, 0, 0), #null case - pk=0.182
                  "Null_high" =c(1.5, 0, 0, 0), #null case - pk=0.818
                  "Alternative"= c(0, 0, 2, 0), # alternative, no effect of covariate
                  "Alternative_low" = c(0, 0, 0.7, 0),  
                  "Alternative_low2" = c(0, 0, 0.3, 0),  
                  "Alternative_high" = c(-1, 0, 2, 0),  #more contrast between arms 
                  "Coveffect" = c(0, -2, 2, 0),      #alternative, effect of covariate 
                  "Synergistic" = c(0, 1, 1, 0),  # covariate makes choice for arm 1 harder 
                  "Antagonistic"= c(0, -1, 1, -2)) # choice of best arm is different 



alpha_names <- c("Complete", "MCAR",
                 "MAR_X1", "MAR_X0",
                 "MAR_T1_high", "MAR_T1_low",
                 "MNAR_Y0", "MNAR_Y1", "MNAR_T1Y0", 
                 "MNAR_X1Y1","MNAR_X0Y0")

#alpha parameters for the model for missingness mechanism------------------------
alpha_list <- list("Complete"= "FALSE",                     
                   "MCAR" = c(-1, 0, 0, 0), 
                   "MAR_X1" =  c(-5, 5, 0, 0), 
                   #"MAR_X0" =  c(0, -5, 0, 0),
                   "MAR_T1_high"=c(-5, 0, 5, 0),
                   "MAR_T1_low"=c(-5, 0, 3, 0),
                   "MNAR_Y0" = c(0, 0, 0, -5), 
                   #"MNAR_Y1" = c(-5, 0, 0, 5), 
                   "MNAR_T1Y0" =  c(0, 0, 0.5, -4.5), 
                   #"MNAR_X1Y1" =  c(-5, 4.5, 0, 0.5), 
                   "MNAR_X0Y0" =  c(0, -4.5, 0, -0.5))

missing_names <- c("cc", "mean_impute_current", 
                   "mean_impute", "MI", "MI_covariate", "impute_zero")
#load common functions ---------------------------------------------------------------


prob_success <- function(beta, x, t){ 
  #given parameters of the data-generating mechanism, 
  #and a value for the covariate (x) and treatment arm (k), output probabilty 
  #of success. Allows for treatment-covariate interactions 
  prob <- 1/(1 + exp(-(beta[1]+ beta[2]*x+beta[3]*t+beta[4]*x*t)))   
  
  return(prob)
}


prob_missing<- function(alpha, beta, x, t){
  
  y <- ifelse(prob_success(beta, x, t) >= 0.5, 1, 0)
  
  
  prob <- 1/(1+exp(-(alpha[1]+alpha[2]*x +alpha[3]*t+ alpha[4]*y)))
  
  return(prob)
  
}
#load common functions for the simullation study 



prob_missing_linear <- function(alpha, beta, x, t){
  #given the parameters of the model for the response (beta)
  # and the parameters of the model for missingness (alpha), 
  # return probability of successful outcome given x (covariate) and t (treatment)
  
  prob <- 1/(1+exp(-((alpha[1]+alpha[4]*beta[1])+(alpha[2]+alpha[4]*beta[2])*x +
                       (alpha[3]+alpha[4]*beta[3])*t  + alpha[4]*beta[4]*t*x)))
  
  return(prob)
  
}



setting_summary <- function(alpha=c(0, 0, 0, 0), beta){ 
  #create table showing success probabilities for all combinations of t and x 
  
  arm0 <- c(prob_success(beta, 0, 0), prob_success(beta, 1, 0)) 
  arm1 <- c(prob_success(beta, 0, 1), prob_success(beta, 1, 1))
  
  if(is.character(alpha)==TRUE){
    p_arm0=p_arm1=c(0, 0)
  }else{
    
    p_arm0 <- c(prob_missing(alpha,beta, 0, 0), prob_missing(alpha,beta, 1, 0)) 
    p_arm1 <- c(prob_missing(alpha,beta, 0, 1), prob_missing(alpha,beta, 1, 1)) 
    
  }
  

    
  result <- rbind(c(arm0, mean(arm0)), 
                  c(arm1, mean(arm1)), 
                  c(p_arm0, mean(p_arm0)), 
                  c(p_arm1, mean(p_arm1)))
  colnames(result) <- c("conditional on X=x0", "conditional on X=x1", "marginal")
  rownames(result) <- c("p(success for arm0)",
                        "p(success for arm1)",
                        "p(missing value for arm0)",
                        "p(missing value for arm1)")
  
  
  
  return(result)
  
}

