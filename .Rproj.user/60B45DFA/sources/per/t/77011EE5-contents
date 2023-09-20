##### Functions
library(secr)
library(tidyverse)
library(rlang)

##### HHN encounter rate function

lambda_hhn <- function(d,lambda0,sigma){
  return(lambda0*exp(-d^2/(2*(sigma^2))))
}

##### HEX encounter rate function

lambda_hex <- function(d,lambda0,sigma){
  return(lambda0*exp(-d/sigma))
}


##### Introduce ghosts as 
##### 1. As some probability of detected individuals
##### 2. As some proportion of detected individuals
##### 3. As some fixed number n

# alternatively, introduce as probability or proportion of detections
# by substituting nrow(capthist) with sum(capthist)

introduce_ghost <- function(capthist,n,type = 'proportion'){
  ###convert secr capture history file to dataframe
  capt <- data.frame(capthist)                
  
  ###number of ghosts to generate
  
  if(type == 'probability')newg <- rbinom(1,nrow(capthist),n)
  if(type == 'proportion')newg <- round(n*nrow(capthist),0)
  if(type == 'fixed') newg <- n
  
  ###Remove single captures for which can't be made to ghosts to append back later
  capt_sing <- capt %>% group_by(ID) %>% filter(n() == 1)
  
  ###subset data to all individuals with more that 2 captures
  capt_plus <- capt %>% group_by(ID) %>% filter(n() > 1)
  
  ###Change individual ID randomly for the number of ghosts to be created
  capt_plus$ID[sample(1:nrow(capt_plus),newg)] <- paste0('ghost',1:newg)
  
  ###Add back the single captures
  capt <- bind_rows(capt_sing,capt_plus) %>% data.frame()
  
  ###convert back to a capthist object
  capt <- make.capthist(captures = capt,traps = traps(capthist),fmt = 'trapID')
  return(capt)
}


###### remove single captures

remove_single_captures <- function(capthist){
  return(capthist %>% data.frame() %>% 
           group_by(ID) %>% 
           filter(n()>1) %>% 
           data.frame() %>% 
           make.capthist(captures = .,traps = traps(capthist),fmt = 'trapID')) 
}



modifymodel <- function(model, user_model){
  modifyList(model,user_model[intersect(names(model),names(user_model))]) 
}


scr2_lik <- function(caphist, mesh, detectfn = 'HHN', simulations = FALSE, model = NULL, startparams = NULL,method = 'Nelder-Mead',hessian = F){
  
  ## Define the standard formula
  correct_model = list(D~1,lambda0~1,sigma~1)
  ##Name the list with the formulae
  names(correct_model) <- lapply(correct_model, f_lhs)
  
  ## Use initial parameters if provided
  if(!is.null(model)) {
    names(model) <- lapply(model, f_lhs)
    model = modifymodel(correct_model,model)
  }else{
    model = correct_model
  }
  
  ## Remove single encounters
  caphist <- remove_single_captures(caphist)    
  
  if(detectfn == 'HHN') lambda_x <- lambda_hhn  ### Hazard half normal encounter rate
  if(detectfn == 'HEX') lambda_x <- lambda_hex  ### Hazard exponential encounter rate
  
  ## Create design matrix
  desmat <- lapply(model,model.matrix,data = data.frame(cbind(D = 1,lambda0 = 1,sigma = 1, covariates(mesh))))
  
  ## Calculate number of parameters
  npars = lapply(desmat,ncol)
  
  ### extract number of individuals detected
  n <- nrow(caphist)  
  
  ### Area of each pixel
  a <- attr(mesh,'spacing')^2/100^2  
  
  
  ## Extract/Define usage
  if (simulations == FALSE) {
    usage = usage(traps(caphist))
  }else{
    usage = rep(1, times = nrow(traps(caphist)))
  }
  
  ## Extract traps object
  traps = traps(caphist)
  
  ## distance from traps to each mesh point
  distmat = edist(traps,mesh) 
  
  ## Set up the parameter vector
  if(!is.null(startparams)){  
    params = startparams
  }else{
    autopars = autoini(caphist,mesh)
    D = numeric(length = npars$D)
    D[1] = log(autopars$D)
    lambda0 = numeric(length = npars$lambda0)
    lambda0[1] = log(1 - exp(-autopars$g0))
    sigma = numeric(length = npars$sigma)
    sigma[1] = log(autopars$sigma)
    params = c(D,lambda0,sigma)
  }
  
  ## Vector of the index of each parameter
  args.index = split(1:length(params),rep(c('D','lambda0','sigma'),times = c(npars$D,npars$lambda0,npars$sigma)))
  
  ## Remove occasions 
  caphist = t(caphist[,1,])   ###Remove a dimension as not using occasions
  
  ## Define the lambda matrix
  lambda_x_mat = matrix(nrow = nrow(traps), ncol = nrow(mesh))
  
  ## Define the p. vector
  prob_capt = numeric(length = nrow(mesh))
  
  ## Local function to optimise 
  optimiser <- local ({
    lik_opt <- function(params){
      D = exp(desmat$D %*% params[args.index$D])
      lambda0 = exp(desmat$lambda0 %*% params[args.index$lambda0])
      sigma = exp(desmat$sigma %*% params[args.index$sigma])  
      
      lambda_x_mat <<- t(apply(t(distmat), 2,lambda_x ,lambda0 = lambda0,sigma = sigma)%*%diag(as.numeric(usage)))
      
      lik_all_ind <- sapply(1:ncol(caphist),function(x) log(mean(exp(colSums(dpois(caphist[,x],lambda_x_mat,log = T)))*D*a)))
      
      lam <- colSums(lambda_x_mat) 
      prob_no_capt <- exp(-lam) ### Probability of 0 encounters
      prob_sing_capt <- lam*exp(-lam) ### Probability of 1 encounter 
      prob_capt <<- 1 - prob_no_capt - prob_sing_capt  ### probability of atleast 2 encounters
      
      ### the full likelihood
      loglik = - sum(D * a * prob_capt) - lfactorial(n) + sum(lik_all_ind)  
      
      return(-loglik)
    }
  })
  
  args = list(f = optimiser,p = params,hessian = hessian, print.level = 1)
  
  ## Run optimiser 
  mod = do.call(nlm,args)
  
  ## Output optimised parameters, lambda matrix and p.
  
  return(mod)
  
}
