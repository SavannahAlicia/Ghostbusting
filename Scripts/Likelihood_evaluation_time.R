source('scripts/2Encounter_functions.R')
capthists = read_rds('Simulations/sim_capthists_lp_hd.rds')[1:100]
mesh = read_rds('Simulations/simulations_mesh.rds')

set.seed(123)

#### Custom function changed to just evaluate the likelihood 

scr2_lik_only <- function(capthist, mesh, K = 2, detectfn = 'HHN', simulations = T, model = list(D~cov), startparams = NULL,hessian = F){
  
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
  
  if (K == 2){
    ## Remove single encounters
    capthist <- remove_single_captures(capthist)    
  }
  
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
  
  D = exp(desmat$D %*% params[args.index$D])
  lambda0 = exp(desmat$lambda0 %*% params[args.index$lambda0])
  sigma = exp(desmat$sigma %*% params[args.index$sigma])  
  
  lambda_x_mat <- t(apply(t(distmat), 2,lambda_x ,lambda0 = lambda0,sigma = sigma)%*%diag(as.numeric(usage)))
  
  lik_all_ind <- sapply(1:ncol(capthist),function(x) log(mean(exp(colSums(dpois(capthist[,x],lambda_x_mat,log = T)))*D*a)))
  
  lam <- colSums(lambda_x_mat) 
  prob_no_capt = exp(-lam) ### Probability of 0 encounters
  prob_capt = 1 - prob_no_capt   ### probability of atleast 2 encounters
  
  if (K == 2){
    prob_sing_capt = lam*exp(-lam) ### Probability of 1 encounter 
    prob_capt = prob_capt - prob_sing_capt
  }
  
  
  ### the full likelihood
  loglik = - sum(D * a * prob_capt) - lfactorial(n) + sum(lik_all_ind)  
  
  return(-loglik)
}

lik_eval_time1 = system.time(sapply(1:length(capthists), function(x) scr2_lik_only(capthists[[x]], mesh[[x]], K = 1)))
lik_eval_time2 = system.time(sapply(1:length(capthists), function(x) scr2_lik_only(capthists[[x]], mesh[[x]], K = 2)))


