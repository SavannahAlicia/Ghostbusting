##### Functions
library(secr)
library(tidyverse)
library(rlang)

##### Half normal encounter rate function

lambda_hhn <- function(d,lambda0,sigma){
  return(lambda0*exp(-d^2/(2*(sigma^2))))
}

##### Exponential encounter rate function

lambda_hex <- function(d,lambda0,sigma){
  return(lambda0*exp(-d/sigma))
}


##### half normal detection function
g_hn <- function(d,g0,sigma){
  return(g0*exp(-d^2/(2*(sigma^2))))
}


##### Introduce ghosts as 
##### 1. As some probability of detected individuals
##### 2. As some proportion of detected individuals
##### 3. As some fixed number n

# alternatively, introduce as probability or proportion of detections
# by substituting nrow(capthist) with sum(capthist)

introduce_ghost <- function(capthist,n,type = 'proportion',noccasions = 1){
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
  capt <- make.capthist(captures = capt,traps = traps(capthist),fmt = 'trapID',noccasions = noccasions)
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



scr2_lik <- function(capthist, mesh, detectfn = 'HHN', simulations = FALSE, model = NULL, startparams = NULL,hessian = F){
  
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
  capthist <- remove_single_captures(capthist)    
  
  if(detectfn == 'HHN') lambda_x <- lambda_hhn  ### Hazard half normal encounter rate
  if(detectfn == 'HEX') lambda_x <- lambda_hex  ### Hazard exponential encounter rate
  
  ## Create design matrix
  desmat <- lapply(model,model.matrix,data = data.frame(cbind(D = 1,lambda0 = 1,sigma = 1, covariates(mesh))))
  
  ## Calculate number of parameters
  npars = lapply(desmat,ncol)
  
  ### extract number of individuals detected
  n <- nrow(capthist)  
  
  ### Area of each pixel
  a <- attr(mesh,'spacing')^2/100^2  
  
  
  ## Extract/Define usage
  if (simulations == FALSE) {
    usage = usage(traps(capthist))
  }else{
    usage = rep(1, times = nrow(traps(capthist)))
  }
  
  ## Extract traps object
  traps = traps(capthist)
  
  ## distance from traps to each mesh point
  distmat = edist(traps,mesh) 
  
  ## Set up the parameter vector
  if(!is.null(startparams)){  
    params = startparams
  }else{
    autopars = autoini(capthist,mesh)
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
  capthist = t(capthist[,1,])   ###Remove a dimension as not using occasions
  
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
      
      lambda_x_mat <- t(apply(t(distmat), 2,lambda_x ,lambda0 = lambda0,sigma = sigma)%*%diag(as.numeric(usage)))
      
      lik_all_ind <- sapply(1:ncol(capthist),function(x) log(mean(exp(colSums(dpois(capthist[,x],lambda_x_mat,log = T)))*D*a)))
      
      lam <- colSums(lambda_x_mat) 
      prob_no_capt <- exp(-lam) ### Probability of 0 encounters
      prob_sing_capt <- lam*exp(-lam) ### Probability of 1 encounter 
      prob_capt <- 1 - prob_no_capt - prob_sing_capt  ### probability of atleast 2 encounters
      
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




scr2_lik_binom <- function(capthist, mesh, simulations = FALSE, model = NULL, startparams = NULL,hessian = F){
  
  ## Define the standard formula
  correct_model = list(D~1,g0~1,sigma~1)
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
  capthist <- remove_single_captures(capthist)
  
  ## Create design matrix
  desmat <- lapply(model,model.matrix,data = data.frame(cbind(D = 1,g0 = 1,sigma = 1, covariates(mesh))))
  
  ## Calculate number of parameters
  npars = lapply(desmat,ncol)
  
  ### extract number of individuals detected
  n <- nrow(capthist)  
  
  ### Area of each pixel
  a <- attr(mesh,'spacing')^2/100^2  
  
  
  ## Extract/Define usage
  if (simulations == FALSE) {
    usage = usage(traps(capthist))
  }else{
    usage = rep(1, times = nrow(traps(capthist)))
  }
  
  ## Extract traps object
  traps = traps(capthist)
  
  ## distance from traps to each mesh point
  distmat = edist(traps,mesh) 
  
  ## Set up the parameter vector
  if(!is.null(startparams)){  
    params = startparams
  }else{
    autopars = autoini(capthist,mesh)
    D = numeric(length = npars$D)
    D[1] = log(autopars$D)
    g0 = numeric(length = npars$g0)
    g0[1] = logit(autopars$g0)
    sigma = numeric(length = npars$sigma)
    sigma[1] = log(autopars$sigma)
    params = c(D,g0,sigma)
  }
  
  ## Vector of the index of each parameter
  args.index = split(1:length(params),rep(c('D','g0','sigma'),times = c(npars$D,npars$g0,npars$sigma)))
  
  capthist1 = apply(capthist,c(1,3),sum)
  
  ## Local function to optimise 
  optimiser <- local ({
    lik_opt <- function(params){
      D = as.vector(exp(desmat$D %*% params[args.index$D]))
      g0 = as.vector(invlogit(desmat$g0 %*% params[args.index$g0]))
      sigma = as.vector(exp(desmat$sigma %*% params[args.index$sigma]))
      
      g_x_mat <- t(apply(t(distmat), 2, g_hn ,g0 = g0,sigma = sigma))
      
      lik_all_ind <- sapply(1:nrow(capthist1),function(i)
        log(mean(exp(colSums(dbinom(capthist1[i,],ncol(capthist),g_x_mat,log = T)))*D*a)))
      
      haz = -log(1-g_x_mat)
      lam <- ncol(capthist)*colSums(haz) 
      prob_capt <- 1 - exp(-lam) - lam*exp(-lam)### Probability of 0 encounters
      
      ### the full likelihood
      loglik = - sum(D * a * prob_capt) - lfactorial(n) + sum(lik_all_ind)  
      
      return(-loglik)
    }
  })
  
  args = list(f = optimiser,p = params,hessian = hessian, print.level = 1)
  
  ## Run optimiser 
  mod = do.call(nlm,args)
  
  ## Output optimised parameters, g matrix and p.
  
  return(mod)
}



##### Simulate capture histories

simulate_scr_files <- function(spacing = NULL,        #spacing between mesh points
                               lambda0 = NULL,        
                               g0 = NULL,
                               sigma = 300, 
                               D = .1, 
                               d.cov = NULL, 
                               nx_traps = 5,          # No of  X x X traps 
                               trap_spacing = NULL,   # Spacing between traps in grid
                               spatial_smooth = 10){  # Spacial covariate smooth factor
  
  #Define spacing between traps as is recommended
  if(is.null(trap_spacing)) trap_spacing = 1.5*sigma
  
  #Create trapping grid
  raw_trap = expand.grid(x = seq( - ((nx_traps-1)*trap_spacing/2), ((nx_traps-1)*trap_spacing/2),length.out = nx_traps),
                         y = seq(- ((nx_traps-1)*trap_spacing/2),((nx_traps-1)*trap_spacing/2),length.out = nx_traps))
  
  ### Chnage detector type based on specified model
  if(is.null(lambda0) == is.null(g0)) stop('Provide either lambda0 or g0')
  
  if(is.null(lambda0) == F){
    detector = 'count'
    detectfn = 'hhn'
  }else{
    detector = 'proximity'
    detectfn = 'hn'
  }
  
  ## Create trap object
  traps = read.traps(data = raw_trap %>% data.frame(.,row.names = paste(1:nrow(raw_trap),'_trap',sep = '')), detector = detector)
  
  ## Create mesh object at recommended buffer of 4*sigma
  mesh = make.mask(traps, buffer = 4*sigma, spacing = spacing)
  
  if(is.null(spacing)) spacing = spacing(mesh)
  
  if(is.null(d.cov)){
    
    sim_pops = sim.popn(D,mesh, buffer = 0)
    
  }else{
    
    ## Create spatial covariate and smooth with a gaussian kernel
    distmat = as.matrix(dist(mesh))
    
    V = exp(-distmat/(spacing*spatial_smooth))
    
    run = runif(nrow(mesh),min = -1.5,max = 1.5)
    
    cov = t(chol(V))%*%run
    
    rm(distmat,V)
    
    covariates(mesh)$cov = cov
    set.seed(NULL)
    
    #simulate population
    sim_pops = sim.popn(D = exp(log(D) + d.cov*cov), model2D = 'IHP', core = mesh)
  }
  
  #simulate capture hisotry
  capthist = sim.capthist(traps,sim_pops,detectfn = detectfn,detectpar = list(sigma = sigma,lambda0 = lambda0),noccasions = 1,renumber = F)
  
  return(list(mesh = mesh, pop = sim_pops,capthist = capthist))
}


