rm(list = ls())

library(secr)
library(tidyverse)
library(raster)
library(ggpubr)
library(patchwork)
library(rlang)

source('2Encounter_functions.R')

site_names <- c("Munkhkhairkhan","Uyench","KhoridolSaridag","KharkhiraaTurgen",
                "SairKhatuu", "Darvi","KhasagtKhairkhan",   "Khuvch",
                "BaatarKhairkhan","BaruunShiluustei","Nemegt","Noyon","Tost","Zoolon",
                "Gurvainsaikhan", "JargalantKhairkhan", "Tsambagarav")   

mesh = read_rds('mesh.rds')
all_traps = read_rds('traps.rds')
caps = read.csv('captures.csv')

#### Areas of each site
data.frame(sites = site_names, area = sapply(mesh, function(x) spacing(x)^2*nrow(x)/1000^2))

#### Make capture history of data k  = 1
Mongolia1 <- make.capthist(captures = caps,traps = all_traps,bysession = T)

summary(Mongolia1,terse = T)

#### Make capture history of data k  = 2
Mongolia2 <- caps  %>% 
  group_by(ID) %>% 
  filter(n()>1) %>% 
  data.frame() %>% 
  make.capthist(captures = .,traps = all_traps,bysession = T) 

summary(Mongolia2,terse = T)

###### secr individual models
####homogeneous

setNumThreads(6)
hom_hex_ind_mod = lapply(1:length(Mongolia1),function(x)  tryCatch(secr.fit(Mongolia1[[x]],mesh[[x]],model = list(D~value), detectfn = 'HEX'),error = function(e) c(0,0,0,0)))

sapply(hom_hex_ind_mod[c(-3,-10)],function(x)region.N(x)[1,1])
sapply(hom_hex_ind_mod[c(-3,-10)],function(x)coef(x)[,1])

coef(hom_hex_ind_mod[[1]])

inhom_hex_shared_d1_mod = secr.fit(Mongolia1,model=list(D~session + value, lambda0~session, sigma~session), detectfn="HEX",mask=mesh) 

region.N(inhom_hex_shared_d1_mod)

saveRDS(hom_hex_ind_mod,'Mongolia_hom_individual_sites_hex.rds')

saveRDS(inhom_hex_shared_d1_mod,'Mongolia_shared_d1_hex.rds')
###Run models with 2 detection functions
###Exponential detection function has lower AIC


# hhn_mod <- secr.fit(Mongolia1,model=list(D~session + value, lambda0~1, sigma~1), detectfn="HHN",
#                      mask=mesh)
# 
# hex_mod <- secr.fit(Mongolia1,model=list(D~session + value, lambda0~1, sigma~1), detectfn="HEX",
#                      mask=mesh)
# 
# 
# saveRDS(list(hhn_mod,hex_mod), 'multisess_scrfit_hhn_hex.rds')

mods = readRDS('multisess_scrfit_hhn_hex.rds')

caphist = Mongolia2

##Number of detected individuals in each site
n = sapply(caphist, nrow)

##Area of each pixel
a = spacing(mesh[[1]])^2/100^2

## Extract traps object
traps = traps(caphist)

##Extract effort of each trap 
usage = sapply(traps, usage)

## distance from traps to each mesh point
distlist = lapply(1:length(traps),function(x) edist(traps[[x]],mesh[[x]]))

## Reduce dimensionality (remove occasions) for custom function
caphist = lapply(caphist, function(x) t(x[,1,]))   

## fix site with single individual
caphist$`3` <- t(caphist$`3`)

## Extract density covariate for each site
cov = lapply(mesh, function(x) covariates(x)$value)

## Function to get the likelihood for each site
site_lik = function(D1,lambda01,sigma1,caphist1,distmat1,usage1,n1){
  
  if(length(caphist1)== 0){
    return(0)
  }else{
    ## compute encounter matrix
    lambda_x_mat <- t(apply(t(distmat1), 2,lambda_hex ,lambda0 = lambda01,sigma = sigma1)%*%diag(as.numeric(usage1)))
    
    ## Compute the conditional likelihood for each individual
    lik_all_ind <- sapply(1:ncol(caphist1),function(x) log(mean(exp(colSums(dpois(caphist1[,x],lambda_x_mat,log = T)))*D1*a,na.rm = T)))
    
    ## Compute the probability of detecting an individual atleast twice given an activity centre
    lam <- colSums(lambda_x_mat) 
    prob_no_capt <- exp(-lam) ### Probability of 0 encounters
    prob_sing_capt <- lam*exp(-lam) ### Probability of 1 encounter 
    prob_capt <- 1 - prob_no_capt - prob_sing_capt  ### probability of atleast 2 encounters
    
    ### the full likelihood
    loglik = - sum(D1 * a * prob_capt,na.rm = T) - lfactorial(n1) + sum(lik_all_ind)  
    
    return(-loglik)
  }
}


## Optimiser for fitting the multisession model 
optimiser <- local ({
  lik_opt <- function(params){
    
    ## Density parameters. Covariate effect + intercept for each site.
    D = c(params[1], params[1] + params[2:(length(params)-3)])
    lambda0 = exp(params[(length(params)-1)])
    sigma = exp(params[length(params)])
    
    sum(sapply(c(1:(length(params)-3)), function(x) site_lik(D1 = exp(D[x] + params[(length(params)-2)]*cov[[x]]),
                                                             lambda01 = lambda0,
                                                             sigma1 = sigma,
                                                             caphist1 = caphist[[x]],
                                                             distmat1 = distlist[[x]],
                                                             usage1 = usage[[x]],
                                                             n1 = n[x])))
    
  }
})



hex_mod <- read_rds('multisess_scrfit_hhn_hex.rds')[[2]]

params = coef(hex_mod)[,1]

method = 'Nelder-Mead'
hessian = T
args = list(par = params,fn = optimiser,method = method,hessian = hessian, control = list(trace = 2,maxit = 5000))

## Run optimiser 

#mod = do.call(optim,args)
# saveRDS(mod,'Mongolia_min2Captures_contD_mod_hex1.rds')

mod <- readRDS('Mongolia_min2Captures_contD_mod_hex1.rds')

###Extract variance covariance matrix
mod$beta.vcv <- solve(mod$hessian)


###Covariate effect on Density and it's SE of SCR2 model
mod$par[18]
sqrt(mod$beta.vcv[18,18])

###Lambda0 estimate and SE of SCR2 model
exp(mod$par[19])

exp(mod$par[19])*sqrt(mod$beta.vcv[19,19])

###Sigma estimate and SE of SCR2 model
exp(mod$par[20])

exp(mod$par[20])*sqrt(mod$beta.vcv[20,20])


###Derive encounter rates for each survey from estimates parameters
lambda_x_list <- lapply(1:length(caphist), function(x) t(apply(t(distlist[[x]]), 2,lambda_hex ,lambda0 = exp(mod$par[length(caphist)+2]),sigma = exp(mod$par[length(caphist)+3]))%*%diag(as.numeric(usage[[x]]))))


cap_stats <- function(caphist){
  data.frame(caphist) %>% 
    group_by(ID) %>% 
    mutate(dets = n())  %>%
    ungroup() %>% 
    summarise(detections = n(), individuals = length(unique(ID)), single_detections = sum(dets == 1))
}

###capture history  summaries
captures_table <- data.frame(t(sapply(Mongolia1,cap_stats))) 


###Abundance of each site from SCR model 
scr_abund <- region.N(hex_mod)

abund_dat <- data.frame(t(sapply(scr_abund, function(x) x[1,1:2])))


### Abundance and SE for SCR2 model
scr2_nhpp_abund_f1 <- function(site,betas,bcov){
  site_ind <- c(1,rep(0,16))
  site_ind[site] <- 1
  b0 = as.numeric(site_ind%*%betas)
  sum(exp(b0 + bcov*covariates(mesh[[site]])$value)*attr(mesh[[site]],'area'),na.rm = T)
}

scr2_nhpp_abund <- sapply(1:length(caphist), function(x) scr2_nhpp_abund_f1(x,mod$par[1:17],mod$par[18]))

scr2_se <- vector(length = 17)

### Se of abundance using delta method
for(i in 1:17){
  grad2 <- nlme::fdHess(pars = c(mod$par[1:18]), function(x) scr2_nhpp_abund_f1(i,x[1:17],x[18]))$grad
  
  scr2_se[i] <- sqrt(t(grad2)%*%solve(mod$hessian[c(1:18),c(1:18)])%*%grad2)
}

###combine both abundance estimates from both methods 
cvdat <- data.frame(site_names,captures_table[,1:2],cbind(abund_dat,scr2_nhpp_abund,scr2_se)) %>% unnest()

###Compute CVs of abundance estimates and difference in abundance and CV between the methods
cvdat <- cvdat %>% 
  mutate(scr_cv = SE.estimate/estimate*100, 
         scr2_cv = scr2_se/scr2_nhpp_abund*100,
         deltaN = (scr2_nhpp_abund-estimate)/scr2_nhpp_abund*100,
         delatcv = scr2_cv - scr_cv)


###Clean up and order table
cvdat[,-1] <- apply(cvdat[,-1],MARGIN = 2, round,digits = 2)

cvdat <- cvdat[,c(1:5,8,6:7,9:11)]

colnames(cvdat) <- c('site','detections','individuals','scr_en','scr_se','scr_cv','scr2_en','scr2_se','scr2_cv','delta_en','delta_cv')

###Table 1 of manuscript
cvdat %>% arrange(delta_cv)

write.csv(cvdat,'outputs/min2_captures_multisession_results_CVtable.csv')


D0 = mod$par[1] + c(0,mod$par[2:17])
D1 = mod$par[18]
nhpp_exp_sing_f <- function(site,D0,D1){
  tryCatch({
    lam <- colSums(lambda_x_list[[site]])
    sum(exp(D0 + D1*covariates(mesh[[site]])$value)*attr(mesh[[site]],'area')*lam*exp(-lam),na.rm = T)
  },error = function(e) 0)
}


###Compute expected number of single detections using SCR2 parameters
exp_sing_dat <- sapply(1:length(caphist), function(x) nhpp_exp_sing_f(x,D0[x],D1))


###Run hypothesis test for each site
tab1 <- cbind(site = site_names,
              captures_table,
              data.frame(scr_D = as.numeric(abund_dat[,1]),
                         scr2_D = scr2_nhpp_abund,
                         exp_sing_cap = exp_sing_dat)) %>%   
  mutate(p_value = ppois(as.numeric(single_detections)-1,exp_sing_dat,lower.tail = F))


###Clean up and order table
tab1 <- tab1 %>% 
  dplyr::select(site,detections,individuals,scr_D,scr2_D,single_detections,exp_sing_cap,p_value) %>% 
  data.frame()

tab1[,-1] <- apply(tab1[,-1],MARGIN = 2, as.numeric)

tab1[,-1] <- apply(tab1[,-1],MARGIN = 2, round,digits = 2)

###Table 2 of manuscript
arrange(tab1,p_value)

write.csv(tab1,'outputs/min2_captures_multisession_results_hypothesisTest.csv')


##########Simulations
##########
##Simulating Tost 

set.seed(123)

mod_index = 13

###compute probability of detecting an individual once and 
###probability of detecting an individual at least twice
###given theta and an activity centre

pdot_p2dot <- function(site){
  lam <- colSums(lambda_x_list[[site]])
  pdot <- 1 - exp(-lam)
  p2dot = pdot - lam*exp(-lam)
  return(list(pdot = pdot, p2dot = p2dot))
}

tost_pdots <- pdot_p2dot(mod_index)

###Visualise pdot
tostp1 <- ggplot(mesh[[mod_index]],aes(x =x, y = y))+
  geom_tile(aes(fill = tost_pdots[[1]]))+
  geom_point(data = traps[[mod_index]], col = 'firebrick',shape = 3)+
  scale_fill_viridis_c(name = '',limits = c(0,1))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_void()+
  ggtitle('')+
  coord_equal()+
  theme(plot.title = element_text(hjust = .5),legend.key.height = unit(1,'inch'),plot.margin = unit(c(0,0,0,0),'cm'))

###Visualise p2dot
tostp2 <- ggplot(mesh[[mod_index]],aes(x =x, y = y))+
  geom_tile(aes(fill = tost_pdots[[2]]))+
  geom_point(data = traps[[mod_index]], col = 'firebrick',shape = 3)+
  scale_fill_viridis_c(name = '',limits = c(0,1))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_void()+
  ggtitle('')+
  coord_equal()+
  theme(plot.title = element_text(hjust = .5),legend.key.height = unit(1,'inch'),plot.margin = unit(c(0,0,0,0),'cm'))


library(gridExtra)

ggarrange(tostp1,tostp2,nrow = 2,common.legend = T, legend = 'right')

ggsave('outputs/tost_pdot.png',width = 8, units = 'in', dpi = 300)

###Extract model parameters for simulation
d1 = mod$par[18]
true_lambda0 = exp(mod$par[19])
true_sigma = exp(mod$par[20])


site_d0 = mod$par[1] + mod$par[mod_index]
site_mesh = mesh[[mod_index]]
site_traps = traps[[mod_index]]
site_usage = usage[[mod_index]]

site_EN <- sum(exp(site_d0 + d1*covariates(site_mesh)$value)*attr(site_mesh,'area'))

set.seed(213)

###Simulate 100 populations using model estimate of Density
site_pop <- replicate(100,sim.popn(exp(site_d0 + d1*covariates(site_mesh)$value),core = site_mesh,buffer = 0,model2D = 'IHP'),simplify = F)

###Simulate a capture history for each population using model observation parameters (theta)
site_caphist <- lapply(site_pop, function(x) sim.capthist(site_traps,x,detectfn = 'hex',detectpar = list(lambda0 = true_lambda0,sigma = true_sigma), noccasions = 1))

###Mean number of individuals detected
mean(sapply(site_caphist,nrow))
###Mean number of detections overall
mean(sapply(site_caphist,sum))

###Introduce ghosts into capture histories at 10%, 20% and 30% of individuals detected
ghost_capthist1 <- lapply(site_caphist, function(x) introduce_ghost(x,.1))
ghost_capthist2 <- lapply(site_caphist, function(x) introduce_ghost(x,.2))
ghost_capthist3 <- lapply(site_caphist, function(x) introduce_ghost(x,.3))
ghost_capthist <- c(site_caphist,ghost_capthist1,ghost_capthist2,ghost_capthist3)

setNumThreads(15)

###Fit standard SCR model
scr_fits <- lapply(ghost_capthist, function(x) tryCatch(secr.fit(x,site_mesh,model = list(D~value), detectfn = 'hex',start = c(D = site_d0,D.value = d1,lambda0 = log(true_lambda0), sigma = log(true_sigma))),error = function(e) NULL))

###Fir 
# min2_scr_fits <- lapply(ghost_capthist, function(x) tryCatch(scr2_lik(x,site_mesh,model = list(D~value), detectfn = 'HEX',startparams = c(D = site_d0,D.value = d1,lambda0 = log(true_lambda0), sigma = log(true_sigma))),error = function(e) NULL))
# 
# saveRDS(list(scr_fits,min2_scr_fits,ghost_capthist),'outputs/Tost_100sim.rds')

scr_tost_sim <- readRDS('outputs/Tost_100sim.rds')[[1]]
min2_scr_tost_sim <- readRDS('outputs/Tost_100sim.rds')[[2]]


scr_coefs <- data.frame(t(sapply(scr_tost_sim, function(x) coef(x)[,1]))) %>% mutate(mod = 'scr',trueD = site_d0,d1 = d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = 100))
scr2_coefs <- data.frame(t(sapply(min2_scr_tost_sim, function(x) x$optim$estimate))) %>% mutate(mod = 'scr2',trueD = site_d0,d1 = d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = 100))

colnames(scr_coefs)[1:4] = c('D','D.value','lambda0','sigma')

colnames(scr2_coefs) = colnames(scr_coefs)


tost_bias <- rbind(scr_coefs,scr2_coefs) %>% 
  rowwise() %>% 
  mutate(E_N = sum(exp(trueD + d1*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         mod_EN = sum(exp(D + D.value*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         EN_b = (mod_EN - E_N)/abs(E_N),
         D_b = (D- trueD)/abs(trueD),
         d1_b = (D.value - d1)/abs(d1),
         lambda0_b = (exp(lambda0) - true_lambda0)/abs(true_lambda0),
         sigma_b = (exp(sigma)- true_sigma)/abs(true_sigma)) %>% 
  group_by(mod,ghost_prop) %>% 
  summarise(EN_Bias = mean(EN_b), EN_lcl = EN_Bias - 1.96*sd(EN_b)/sqrt(n()),EN_ucl = EN_Bias + 1.96*sd(EN_b)/sqrt(n()),
            D_Bias = mean(D_b), D_lcl = D_Bias - 1.96*sd(D_b)/sqrt(n()),D_ucl = D_Bias + 1.96*sd(D_b)/sqrt(n()),
            d1_Bias = mean(d1_b), d1_lcl = d1_Bias - 1.96*sd(d1_b)/sqrt(n()),d1_ucl = d1_Bias + 1.96*sd(d1_b)/sqrt(n()),
            lambda0_Bias = mean(lambda0_b), lambda0_lcl = lambda0_Bias - 1.96*sd(lambda0_b)/sqrt(n()),lambda0_ucl = lambda0_Bias + 1.96*sd(lambda0_b)/sqrt(n()),
            sigma_Bias = mean(sigma_b), sigma_lcl = sigma_Bias - 1.96*sd(sigma_b)/sqrt(n()),sigma_ucl = sigma_Bias + 1.96*sd(sigma_b)/sqrt(n())) %>% 
  pivot_longer(3:17,names_sep = '_',names_to = c('par','estimate')) %>% 
  pivot_wider(names_from = 'estimate',values_from = value) %>% 
  filter(!par %in% c('D','d1')) 

tost_bias$par <- factor(tost_bias$par,labels = c(expression('E[N'~']'),bquote(lambda[0]),bquote(sigma)))


tost_bias_plot <- tost_bias %>% 
  ggplot(aes(x = ghost_prop, y = Bias, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = F)+
  geom_errorbar(aes(ymax = ucl,ymin = lcl),position = position_dodge(width = .05), width = .01,show.legend = F)+
  geom_hline(yintercept = 0)+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  #  geom_blank(data = zoolon_bias)+
  facet_grid(par~., scales = 'free',labeller = label_parsed)+
  ylab('')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))


tost_bias_plot

ggsave('outputs/Tost_100simulation.png',width = 10, units = 'in',dpi = 300)


mod_var <- rbind(scr_coefs,scr2_coefs) %>%
  rowwise() %>%
  mutate(E_N = sum(exp(trueD + d1*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         mod_EN = sum(exp(D + D.value*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         EN_v = (mod_EN - E_N)^2,
         D_v = (D- trueD)^2,
         d1_v = (D.value - d1)^2,
         lambda0_v = (exp(lambda0) - true_lambda0)^2,
         sigma_v = (exp(sigma) - true_sigma)^2) %>%
  group_by(mod,ghost_prop) %>%
  summarise(EN_variance = mean(EN_v),
            D_variance = mean(D_v),
            d1_variance = mean(d1_v),
            lambda0_variance = mean(lambda0_v),
            sigma_variance = mean(sigma_v)) %>%
  pivot_longer(3:7,names_sep = '_',names_to = c('par','estimate')) %>%
  pivot_wider(names_from = 'estimate',values_from = value)  %>%
  filter(!par %in% c('D','d1'))

mod_var %>% 
  ggplot(aes(x = ghost_prop, y = variance, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = F)+
  #  geom_errorbar(aes(ymax = ucl,ymin = lcl),position = position_dodge(width = .05), width = .01,show.legend = F)+
  geom_hline(yintercept = 0)+
  scale_x_continuous(labels = scales::percent)+
  #  scale_y_continuous(labels = scales::percent)+
  facet_grid(par~. ,scales = 'free',labeller = label_parsed)+
  ylab('MSE')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))

ggsave('outputs/Tost_100simulation_MSE.png',width = 10, units = 'in',dpi = 300)

##########
##Simulating zoolon

set.seed(123)
mod_index = 14



pdot_p2dot <- function(site){
  lam <- colSums(lambda_x_list[[site]])
  pdot <- 1 - exp(-lam)
  p2dot = pdot - lam*exp(-lam)
  return(list(pdot = pdot, p2dot = p2dot))
}

zoolon_pdots <- pdot_p2dot(mod_index)

zoolonp1 <- ggplot(mesh[[mod_index]],aes(x =x, y = y))+
  geom_tile(aes(fill = zoolon_pdots[[1]]))+
  geom_point(data = traps[[mod_index]], col = 'firebrick',shape = 3)+
  scale_fill_viridis_c(name = '',limits = c(0,1))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_void()+
  ggtitle('')+
  coord_equal()+
  theme(plot.title = element_text(hjust = .5),legend.key.height = unit(1,'inch'),plot.margin = unit(c(0,0,0,0),'cm'))

zoolonp2 <- ggplot(mesh[[mod_index]],aes(x =x, y = y))+
  geom_tile(aes(fill = zoolon_pdots[[2]]))+
  geom_point(data = traps[[mod_index]], col = 'firebrick',shape = 3)+
  scale_fill_viridis_c(name = '',limits = c(0,1))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_void()+
  ggtitle('')+
  coord_equal()+
  theme(plot.title = element_text(hjust = .5),legend.key.height = unit(1,'inch'),plot.margin = unit(c(0,0,0,0),'cm'))

ggarrange(zoolonp1,zoolonp2,nrow = 2,common.legend = T, legend = 'right')

ggsave('outputs/zoolon_pdot.png',width = 8, units = 'in', dpi = 300)



d1 = mod$par[18]
true_lambda0 = exp(mod$par[19])
true_sigma = exp(mod$par[20])


site_d0 = mod$par[1] + mod$par[mod_index]
site_mesh = mesh[[mod_index]]
site_traps = traps[[mod_index]]
site_usage = usage[[mod_index]]

site_EN <- sum(exp(site_d0 + d1*covariates(site_mesh)$value)*attr(site_mesh,'area'))

site_pop <- replicate(100,sim.popn(exp(site_d0 + d1*covariates(site_mesh)$value),core = site_mesh,buffer = 0,model2D = 'IHP'),simplify = F)

site_caphist <- lapply(site_pop, function(x) sim.capthist(site_traps,x,detectfn = 'hex',detectpar = list(lambda0 = true_lambda0,sigma = true_sigma), noccasions = 1))

mean(sapply(site_caphist,nrow))
mean(sapply(site_caphist,sum))

ghost_capthist1 <- lapply(site_caphist, function(x) introduce_ghost(x,.1))
ghost_capthist2 <- lapply(site_caphist, function(x) introduce_ghost(x,.2))
ghost_capthist3 <- lapply(site_caphist, function(x) introduce_ghost(x,.3))
ghost_capthist <- c(site_caphist,ghost_capthist1,ghost_capthist2,ghost_capthist3)

scr_fits <- lapply(ghost_capthist, function(x) tryCatch(secr.fit(x,site_mesh,model = list(D~value), detectfn = 'hex',start = c(D = site_d0,D.value = d1,lambda0 = log(true_lambda0), sigma = log(true_sigma))),error = function(e) NULL))

min2_scr_fits <- lapply(ghost_capthist, function(x) tryCatch(scr2_lik(x,site_mesh,model = list(D~value), detectfn = 'HEX',startparams = c(D = site_d0,D.value = d1,lambda0 = log(true_lambda0), sigma = log(true_sigma))),error = function(e) NULL))

saveRDS(list(scr_fits,min2_scr_fits,ghost_capthist),'outputs/Zoolon_100sim.rds')

scr_zoolon_sim <- readRDS('outputs/Zoolon_100sim.rds')[[1]]
min2_scr_zoolon_sim <- readRDS('outputs/Zoolon_100sim.rds')[[2]]


scr_coefs <- data.frame(t(sapply(scr_zoolon_sim, function(x) coef(x)[,1]))) %>% mutate(mod = 'scr',trueD = site_d0,d1 = d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = 100))
scr2_coefs <- data.frame(t(sapply(min2_scr_zoolon_sim, function(x) x$optim$estimate))) %>% mutate(mod = 'scr2',trueD = site_d0,d1 = d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = 100))

colnames(scr_coefs)[1:4] = c('D','D.value','lambda0','sigma')

colnames(scr2_coefs) = colnames(scr_coefs)


zoolon_bias <- rbind(scr_coefs,scr2_coefs) %>% 
  rowwise() %>% 
  mutate(E_N = sum(exp(trueD + d1*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         mod_EN = sum(exp(D + D.value*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         EN_b = (mod_EN - E_N)/abs(E_N),
         D_b = (D- trueD)/abs(trueD),
         d1_b = (D.value - d1)/abs(d1),
         lambda0_b = (exp(lambda0) - true_lambda0)/abs(true_lambda0),
         sigma_b = (exp(sigma)- true_sigma)/abs(true_sigma)) %>% 
  group_by(mod,ghost_prop) %>%
  filter(mod_EN < 100) %>% 
  summarise(EN_Bias = mean(EN_b), EN_lcl = EN_Bias - 1.96*sd(EN_b)/sqrt(n()),EN_ucl = EN_Bias + 1.96*sd(EN_b)/sqrt(n()),
            D_Bias = mean(D_b), D_lcl = D_Bias - 1.96*sd(D_b)/sqrt(n()),D_ucl = D_Bias + 1.96*sd(D_b)/sqrt(n()),
            d1_Bias = mean(d1_b), d1_lcl = d1_Bias - 1.96*sd(d1_b)/sqrt(n()),d1_ucl = d1_Bias + 1.96*sd(d1_b)/sqrt(n()),
            lambda0_Bias = mean(lambda0_b), lambda0_lcl = lambda0_Bias - 1.96*sd(lambda0_b)/sqrt(n()),lambda0_ucl = lambda0_Bias + 1.96*sd(lambda0_b)/sqrt(n()),
            sigma_Bias = mean(sigma_b), sigma_lcl = sigma_Bias - 1.96*sd(sigma_b)/sqrt(n()),sigma_ucl = sigma_Bias + 1.96*sd(sigma_b)/sqrt(n())) %>% 
  pivot_longer(3:17,names_sep = '_',names_to = c('par','estimate')) %>% 
  pivot_wider(names_from = 'estimate',values_from = value) %>% 
  filter(!par %in% c('D','d1'))  


zoolon_bias$par <- factor(zoolon_bias$par,labels = c(expression('E[N'~']'),bquote(lambda[0]),bquote(sigma)))

zoolon_bias_plot <- zoolon_bias %>% 
  ggplot(aes(x = ghost_prop, y = Bias, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = F)+
  geom_errorbar(aes(ymax = ucl,ymin = lcl),position = position_dodge(width = .05), width = .01,show.legend = F)+
  geom_hline(yintercept = 0)+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(par~., scales = 'free', labeller = label_parsed)+
  ylab('Relative Bias')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))


zoolon_bias_plot

ggsave('outputs/Zoolon_100simulation.png',width = 8, units = 'in',dpi = 300)

mod_var <- rbind(scr_coefs,scr2_coefs) %>%
  rowwise() %>%
  mutate(E_N = sum(exp(trueD + d1*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         mod_EN = sum(exp(D + D.value*covariates(site_mesh)$value)*attr(site_mesh,'area')),
         EN_v = (mod_EN - E_N)^2,
         D_v = (D- trueD)^2,
         d1_v = (D.value - d1)^2,
         lambda0_v = (exp(lambda0) - true_lambda0)^2,
         sigma_v = (exp(sigma) - true_sigma)^2) %>%
  group_by(mod,ghost_prop) %>%
  filter(mod_EN < 100) %>% 
  summarise(EN_variance = mean(EN_v),
            D_variance = mean(D_v),
            d1_variance = mean(d1_v),
            lambda0_variance = mean(lambda0_v),
            sigma_variance = mean(sigma_v)) %>%
  pivot_longer(3:7,names_sep = '_',names_to = c('par','estimate')) %>%
  pivot_wider(names_from = 'estimate',values_from = value)  %>%
  filter(!par %in% c('D','d1'))

mod_var %>% 
  ggplot(aes(x = ghost_prop, y = variance, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = F)+
  #  geom_errorbar(aes(ymax = ucl,ymin = lcl),position = position_dodge(width = .05), width = .01,show.legend = F)+
#  geom_hline(yintercept = 0)+
  scale_x_continuous(labels = scales::percent)+
  #  scale_y_continuous(labels = scales::percent)+
  facet_grid(par~. ,scales = 'free',labeller = label_parsed)+
  ylab('MSE')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))



figure <- ggarrange(ggarrange(zoolon_bias_plot,tost_bias_plot,ncol = 1,labels = c('A','B')),ggarrange(zoolonp1,zoolonp2,tostp1,tostp2,nrow = 4,common.legend = T, legend = 'right',labels = c('C(i)',' (ii)','D(i)',' (ii)')))

annotate_figure(figure,left = text_grob('Relative Bias',rot = 90,size = 20),bottom = text_grob('Proportion of Ghosts Introduced',size = 20,hjust = 1.5))

ggsave('outputs/sim_outputs.png',width = 15, dpi = 300, units = 'in',height = 10)


ggarrange(zoolonp1,tostp1,zoolonp2,tostp2,ncol= 2,nrow = 2,common.legend = T, legend = 'right', labels = c('Zoolon','Tost','',''),label.x = 0.5)

ggsave('outputs/probdet_outputs.png',width = 12,height = 5, dpi = 300, units = 'in')

ggarrange(zoolon_bias_plot,tost_bias_plot,labels = c('Zoolon','Tost'),label.x = 0.5)

ggsave('outputs/sim_outputs.png',width = 12,height = 5, dpi = 300, units = 'in')


for(i in 1:17){
  png(paste('outputs/trap_layouts/',i,'_',site_names[i],'.png',sep = ''),bg = 'transparent')
  plot(attr(mesh[[i]],'polygon'), col = 'white')
  plot(traps[[i]],add = TRUE)
  dev.off()
}


cap_stats <- function(caphist){
  data.frame(caphist) %>% 
    group_by(ID) %>% 
    mutate(dets = n())  %>%
    ungroup() %>% 
    summarise(detections = n(), individuals = length(unique(ID)), single_detections = sum(dets == 1))
}

covariates(mesh)

capthist = readRDS('outputs/Zoolon_100sim.rds')[[3]]

p_values_sim <- function(scr_fits,min2_scr_fits,mesh,capthist){

  nsims = 100
  
  scr_coefs <- data.frame(t(sapply(scr_fits, function(x) coef(x)[,1]))) %>% mutate(mod = 'scr',ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  scr2_coefs <- data.frame(t(sapply(min2_scr_fits, function(x) x$optim$estimate))) %>% mutate(mod = 'scr2',ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  colnames(scr_coefs)[1:4] = c('D','D.value','lambda0','sigma')
  
  colnames(scr2_coefs) = colnames(scr_coefs)
  
  
  distmat = edist(traps(capthist[[1]]),mesh)
  
  lambda_x_list_scr <- lapply(1:length(capthist), function(x) t(apply(t(distmat), 2,lambda_hhn ,lambda0 = exp(scr_coefs$lambda0[x]),sigma =  exp(scr_coefs$sigma[x]))))
  lambda_x_list_scr2 <- lapply(1:length(capthist), function(x) t(apply(t(distmat), 2,lambda_hhn ,lambda0 = exp(scr2_coefs$lambda0[x]),sigma =  exp(scr2_coefs$sigma[x]))))
  
  scr_exp_sing <- sapply(1:length(capthist), function(x){
    tryCatch({
      lam <- colSums(lambda_x_list_scr[[x]])
      sum(exp(scr_coefs$D[x] + scr_coefs$D.value[x]*covariates(mesh)$value)*attr(mesh,'area')*lam*exp(-lam),na.rm = T)
    },error = function(e) 0)
  })
  
  scr2_exp_sing <- sapply(1:length(capthist), function(x){
    tryCatch({
      lam <- colSums(lambda_x_list_scr2[[x]])
      sum(exp(scr2_coefs$D[x] + scr2_coefs$D.value[x]*covariates(mesh)$value)*attr(mesh,'area')*lam*exp(-lam),na.rm = T)
    },error = function(e) 0)
  })
  
  obs_sing <- t(sapply(capthist, function(x) cap_stats(x)))
  
  
  tab1 <- cbind(obs_sing,scr_exp_sing = scr_exp_sing,scr2_exp_sing = scr2_exp_sing) %>%
    data.frame() %>% 
    unnest() %>% 
    mutate(ghost_prop = rep(c(0,.1,.2,.3),each = nsims)) %>% 
    mutate(scr_p_value = ppois(as.numeric(single_detections)-1,scr_exp_sing,lower.tail = F),
           scr2_p_value = ppois(as.numeric(single_detections)-1,scr2_exp_sing,lower.tail = F)) 
  
  tab1 = tab1 %>% 
    group_by(ghost_prop) %>% 
    summarise(scr_hyp = sum(scr_p_value < .05),scr2_hyp = sum(scr2_p_value < .05)) 
  
  return(tab1)
}

capthist_tost = readRDS('outputs/Tost_100sim.rds')[[3]]

p_values_sim(scr_tost_sim,min2_scr_tost_sim,mesh[[13]],capthist_tost)

capthist_zoolon = readRDS('outputs/zoolon_100sim.rds')[[3]]

p_values_sim(scr_zoolon_sim,min2_scr_zoolon_sim,mesh[[14]],capthist_zoolon)
