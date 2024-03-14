library(secr)
library(tidyverse)
library(ggpubr)

source('Scripts/2Encounter_functions.R')


##### Fixed parameters #####
true_sigma = 300
true_d1 = .1
nsims = 100

##### Simulation Scenarios
lowD = .03
highD = .06
low_lambda0 = 1
high_lambda0 = 2

### Read simulated mesh objects
mesh = readRDS('Simulations/simulations_mesh.rds')

### Read simulated populations (high and low)
lp_pop = readRDS('Simulations/sim_pop_lp.rds')
hp_pop = readRDS('Simulations/sim_pop_hp.rds')

### Read simulated capture histories
lp_ld_capthist = readRDS('Simulations/sim_capthists_lp_ld.rds')
lp_hd_capthist = readRDS('Simulations/sim_capthists_lp_hd.rds')
hp_ld_capthist = readRDS('Simulations/sim_capthists_hp_ld.rds')
hp_hd_capthist = readRDS('Simulations/sim_capthists_hp_hd.rds')

### Read model fit objects
lp_ld_scr = readRDS('Simulations/simulations_scr_lp_ld.rds')
lp_ld_scr2 = readRDS('Simulations/simulations_scr2_lp_ld.rds')
lp_hd_scr = readRDS('Simulations/simulations_scr_lp_hd.rds')
lp_hd_scr2 = readRDS('Simulations/simulations_scr2_lp_hd.rds')
hp_ld_scr = readRDS('Simulations/simulations_scr_hp_ld.rds')
hp_ld_scr2 = readRDS('Simulations/simulations_scr2_hp_ld.rds')
hp_hd_scr = readRDS('Simulations/simulations_scr_hp_hd.rds')
hp_hd_scr2 = readRDS('Simulations/simulations_scr2_hp_hd.rds')


scr2_lik_hess_only <- function(capthist, mesh, detectfn = 'HHN', simulations = TRUE, model = NULL, startparams = NULL){
  
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
  
  args = list(par = params, fn = optimiser)
  
  ## Run optimiser 
  hessian = do.call(optimHess,args)
  
  ## Output optimised parameters, lambda matrix and p.
  
  return(hessian)
  
}

all_mesh = c(mesh,mesh,mesh,mesh)


#### Fit models

scr2_lik_hess_only(lp_ld_capthist[[1]],mesh[[1]],model = list(D~cov),startparams = lp_ld_scr2[[1]]$estimate)


lp_ld_min2_scr_hessian <- lapply(1:length(lp_ld_capthist), function(x) tryCatch(scr2_lik_hess_only(lp_ld_capthist[[x]],all_mesh[[x]],model = list(D~cov),startparams = lp_ld_scr2[[x]]$estimate),error = function(e) NULL))

saveRDS(lp_ld_min2_scr_hessian,'Simulations/simulations_scr2_lp_ld_hessian.rds')


lp_hd_min2_scr_hessian <- lapply(1:length(lp_hd_capthist), function(x) tryCatch(scr2_lik_hess_only(lp_hd_capthist[[x]],all_mesh[[x]],model = list(D~cov),startparams = lp_hd_scr2[[x]]$estimate),error = function(e) NULL))

saveRDS(lp_hd_min2_scr_hessian,'Simulations/simulations_scr2_lp_hd_hessian.rds')


hp_ld_min2_scr_hessian <- lapply(1:length(hp_ld_capthist), function(x) tryCatch(scr2_lik_hess_only(hp_ld_capthist[[x]],all_mesh[[x]],model = list(D~cov),startparams = hp_ld_scr2[[x]]$estimate),error = function(e) NULL))

saveRDS(hp_ld_min2_scr_hessian,'Simulations/simulations_scr2_hp_ld_hessian.rds')


hp_hd_min2_scr_hessian <- lapply(1:length(hp_hd_capthist), function(x) tryCatch(scr2_lik_hess_only(hp_hd_capthist[[x]],all_mesh[[x]],model = list(D~cov),startparams = hp_hd_scr2[[x]]$estimate),error = function(e) NULL))

saveRDS(hp_hd_min2_scr_hessian,'Simulations/simulations_scr2_hp_hd_hessian.rds') 

