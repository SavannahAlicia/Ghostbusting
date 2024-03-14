library(secr)
library(tidyverse)

source('Scripts/2Encounter_functions.R')

##### Fixed parameters #####
true_sigma = 300
true_d1 = .1
nsims = 100
set.seed(123)

##### Simulation Scenarios
lowD = .03
highD = .06
low_lambda0 = 1
high_lambda0 = 2

##### Simulate 100 landscapes

scr_objects = replicate(nsims,simulate_scr_files(D = lowD,
                                                 d.cov = true_d1, 
                                                 lambda0 = low_lambda0,
                                                 sigma = true_sigma,
                                                 nx_traps = 5),simplify = F)

mesh = lapply(scr_objects, function(x) x$mesh)

saveRDS(mesh,'Simulations/simulations_mesh.rds')

##### extract traps
traps = traps(scr_objects[[1]]$capthist)

##### Plot survey
ggplot(mesh[[1]], aes(x = x, y = y))+
  geom_tile(aes(fill = covariates(mesh[[1]])$cov))+
  geom_point(data = traps)



########### low pop
lp_pop = lapply(mesh,function(x) sim.popn(D = exp(log(lowD) + true_d1*covariates(x)$cov), model2D = 'IHP', core = x))

write_rds(lp_pop,'Simulations/sim_pop_lp.rds')
########### low detection data

lp_ld_capthists = lapply(lp_pop, function(x) sim.capthist(traps,x,detectfn = 'HHN',detectpar = list(lambda0 = low_lambda0,sigma = true_sigma),noccasions = 1))


lp_ld_dat = c(mean_pop = mean(sapply(lp_pop, nrow)),
              mean_inds = mean(sapply(lp_ld_capthists,nrow)),
              mean_dets = mean(sapply(lp_ld_capthists,sum)),
              mean_recaps = mean(sapply(lp_ld_capthists,max)),
              mean_sing_dets = mean(sapply(lp_ld_capthists, function(x) sum(rowSums(x)== 1))))


###########  high detection data


lp_hd_capthists = lapply(lp_pop, function(x) sim.capthist(traps,x,detectfn = 'HHN',detectpar = list(lambda0 = high_lambda0,sigma = true_sigma),noccasions = 1))

lp_hd_dat = c(mean_pop = mean(sapply(lp_pop, nrow)),
              mean_inds = mean(sapply(lp_hd_capthists,nrow)),
              mean_dets = mean(sapply(lp_hd_capthists,sum)),
              mean_recaps = mean(sapply(lp_hd_capthists,max)),
              mean_sing_dets = mean(sapply(lp_hd_capthists, function(x) sum(rowSums(x)== 1))))



########### high pop
hp_pop = lapply(mesh,function(x) sim.popn(D = exp(log(highD) + true_d1*covariates(x)$cov), model2D = 'IHP', core = x))

write_rds(hp_pop, 'Simulations/sim_pop_hp.rds')
########### low detection data


hp_ld_capthists = lapply(hp_pop, function(x) sim.capthist(traps,x,detectfn = 'HHN',detectpar = list(lambda0 = low_lambda0,sigma = true_sigma),noccasions = 1))

hp_ld_dat = c(mean_pop = mean(sapply(hp_pop, nrow)),
              mean_inds = mean(sapply(hp_ld_capthists,nrow)),
              mean_dets = mean(sapply(hp_ld_capthists,sum)),
              mean_recaps = mean(sapply(hp_ld_capthists,max)),
              mean_sing_dets = mean(sapply(hp_ld_capthists, function(x) sum(rowSums(x)== 1))))


###########  high detection data

hp_hd_capthists = lapply(hp_pop, function(x) sim.capthist(traps,x,detectfn = 'HHN',detectpar = list(lambda0 = high_lambda0,sigma = true_sigma),noccasions = 1))

hp_hd_dat = c(mean_pop = mean(sapply(hp_pop, nrow)),
              mean_inds = mean(sapply(hp_hd_capthists,nrow)),
              mean_dets = mean(sapply(hp_hd_capthists,sum)),
              mean_recaps = mean(sapply(hp_hd_capthists,max)),
              mean_sing_dets = mean(sapply(hp_hd_capthists, function(x) sum(rowSums(x)== 1))))

lp_ld_dat
lp_hd_dat
hp_ld_dat
hp_hd_dat

###### Create ghost individuals in each capture history

lp_ld_capthists_ghost = c(lp_ld_capthists,
                          lapply(lp_ld_capthists, function(x) introduce_ghost(x,.1)),
                          lapply(lp_ld_capthists, function(x) introduce_ghost(x,.2)),
                          lapply(lp_ld_capthists, function(x) introduce_ghost(x,.3)))

lp_hd_capthists_ghost = c(lp_hd_capthists,
                          lapply(lp_hd_capthists, function(x) introduce_ghost(x,.1)),
                          lapply(lp_hd_capthists, function(x) introduce_ghost(x,.2)),
                          lapply(lp_hd_capthists, function(x) introduce_ghost(x,.3)))


hp_ld_capthists_ghost = c(hp_ld_capthists,
                          lapply(hp_ld_capthists, function(x) introduce_ghost(x,.1)),
                          lapply(hp_ld_capthists, function(x) introduce_ghost(x,.2)),
                          lapply(hp_ld_capthists, function(x) introduce_ghost(x,.3)))

hp_hd_capthists_ghost = c(hp_hd_capthists,
                          lapply(hp_hd_capthists, function(x) introduce_ghost(x,.1)),
                          lapply(hp_hd_capthists, function(x) introduce_ghost(x,.2)),
                          lapply(hp_hd_capthists, function(x) introduce_ghost(x,.3)))



saveRDS(lp_ld_capthists_ghost,'Simulations/sim_capthists_lp_ld.rds')
saveRDS(lp_hd_capthists_ghost,'Simulations/sim_capthists_lp_hd.rds')
saveRDS(hp_ld_capthists_ghost,'Simulations/sim_capthists_hp_ld.rds')
saveRDS(hp_hd_capthists_ghost,'Simulations/sim_capthists_hp_hd.rds')


all_mesh = c(mesh,mesh,mesh,mesh)


#### Fit models

lp_ld_scr_fits <- lapply(1:length(lp_ld_capthists_ghost), function(x) tryCatch(secr.fit(lp_ld_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',start = c(D = log(lowD),D.value = true_d1,lambda0 = log(low_lambda0), sigma = log(true_sigma))),error = function(e) NULL))

saveRDS(lp_ld_scr_fits,'Simulations/simulations_scr_lp_ld.rds')

lp_ld_min2_scr_fits <- lapply(1:length(lp_ld_capthists_ghost), function(x) tryCatch(scr2_lik(lp_ld_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',startparams = c(D = log(lowD),D.value = true_d1,lambda0 = log(low_lambda0), sigma = log(true_sigma)),simulations = T, hessian = T),error = function(e) NULL))

saveRDS(lp_ld_min2_scr_fits,'Simulations/simulations_scr2_lp_ld.rds')

lp_hd_scr_fits <- lapply(1:length(lp_hd_capthists_ghost), function(x) tryCatch(secr.fit(lp_hd_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',start = c(D = log(lowD),D.value = true_d1,lambda0 = log(high_lambda0), sigma = log(true_sigma))),error = function(e) NULL))

saveRDS(lp_hd_scr_fits,'Simulations/simulations_scr_lp_hd.rds')

lp_hd_min2_scr_fits <- lapply(1:length(lp_hd_capthists_ghost), function(x) tryCatch(scr2_lik(lp_hd_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',startparams = c(D = log(lowD),D.value = true_d1,lambda0 = log(high_lambda0), sigma = log(true_sigma)),simulations = T, hessian = T),error = function(e) NULL))

saveRDS(lp_hd_min2_scr_fits,'Simulations/simulations_scr2_lp_hd.rds')

hp_ld_scr_fits <- lapply(1:length(hp_ld_capthists_ghost), function(x) tryCatch(secr.fit(hp_ld_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',start = c(D = log(highD),D.value = true_d1,lambda0 = log(low_lambda0), sigma = log(true_sigma))),error = function(e) NULL))

saveRDS(hp_ld_scr_fits,'Simulations/simulations_scr_hp_ld.rds')

hp_ld_min2_scr_fits <- lapply(1:length(hp_ld_capthists_ghost), function(x) tryCatch(scr2_lik(hp_ld_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',startparams = c(D = log(highD),D.value = true_d1,lambda0 = log(low_lambda0), sigma = log(true_sigma)),simulations = T, hessian = T),error = function(e) NULL))

saveRDS(hp_ld_min2_scr_fits,'Simulations/simulations_scr2_hp_ld.rds')

hp_hd_scr_fits <- lapply(1:length(hp_hd_capthists_ghost), function(x) tryCatch(secr.fit(hp_hd_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',start = c(D = log(highD),D.value = true_d1,lambda0 = log(high_lambda0), sigma = log(true_sigma))),error = function(e) NULL))

saveRDS(hp_hd_scr_fits,'Simulations/simulations_scr_hp_hd.rds')

hp_hd_min2_scr_fits <- lapply(1:length(hp_hd_capthists_ghost), function(x) tryCatch(scr2_lik(hp_hd_capthists_ghost[[x]],all_mesh[[x]],model = list(D~cov), detectfn = 'HHN',startparams = c(D = log(highD),D.value = true_d1,lambda0 = log(high_lambda0), sigma = log(true_sigma)),simulations = T, hessian = T),error = function(e) NULL))

saveRDS(hp_hd_min2_scr_fits,'Simulations/simulations_scr2_hp_hd.rds')