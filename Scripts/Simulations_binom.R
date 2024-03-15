library(secr)
library(tidyverse)

setwd("D:/ar337/Ghostbusting/")
source('Scripts/2Encounter_functions.R')

##### Fixed parameters #####
true_sigma = 300
true_d1 = .1
nsims = 100
set.seed(23)

##### Simulation Scenarios
lowD = .02
highD = .05
low_lambda0 = 1
high_lambda0 = 3
noccasions = 10

low_g0 = 1 - exp(-low_lambda0/noccasions)
high_g0 = 1 - exp(-high_lambda0/noccasions)

mesh = readRDS('Simulations/simulations_mesh.rds')

traps = traps(simulate_scr_files(D = lowD,
                           g0 = low_g0,
                           sigma = true_sigma,
                           nx_traps = 5)$capthist)

########### Low pop
lp_pop = readRDS('Simulations/sim_pop_lp.rds')

########### Low detection data
lp_ld_capthists = lapply(lp_pop, function(x) sim.capthist(traps, x, detectfn = 'HN', detectpar = list(g0 = low_g0,sigma = true_sigma), noccasions = noccasions))

lp_ld_dat = c(mean_pop = mean(sapply(lp_pop, nrow)),
              mean_inds = mean(sapply(lp_ld_capthists,nrow)),
              mean_dets = mean(sapply(lp_ld_capthists,sum)),
              mean_recaps = mean(sapply(lp_ld_capthists,max)),
              mean_sing_dets = mean(sapply(lp_ld_capthists, function(x) sum(rowSums(x)== 1))))


########### High detection data
lp_hd_capthists = lapply(lp_pop, function(x) sim.capthist(traps, x, detectfn = 'HN', detectpar = list(g0 = high_g0,sigma = true_sigma), noccasions = noccasions))

lp_hd_dat = c(mean_pop = mean(sapply(lp_pop, nrow)),
              mean_inds = mean(sapply(lp_hd_capthists,nrow)),
              mean_dets = mean(sapply(lp_hd_capthists,sum)),
              mean_recaps = mean(sapply(lp_hd_capthists,max)),
              mean_sing_dets = mean(sapply(lp_hd_capthists, function(x) sum(rowSums(x)== 1))))

########### High pop
hp_pop = readRDS('Simulations/sim_pop_hp.rds')

########### Low detection data
hp_ld_capthists = lapply(hp_pop, function(x) sim.capthist(traps, x, detectfn = 'HN', detectpar = list(g0 = low_g0, sigma = true_sigma), noccasions = noccasions))

hp_ld_dat = c(mean_pop = mean(sapply(hp_pop, nrow)),
              mean_inds = mean(sapply(hp_ld_capthists,nrow)),
              mean_dets = mean(sapply(hp_ld_capthists,sum)),
              mean_recaps = mean(sapply(hp_ld_capthists,max)),
              mean_sing_dets = mean(sapply(hp_ld_capthists, function(x) sum(rowSums(x)== 1))))

########### High detection data
hp_hd_capthists = lapply(hp_pop, function(x) sim.capthist(traps, x, detectfn = 'HN', detectpar = list(g0 = high_g0, sigma = true_sigma), noccasions = noccasions))

hp_hd_dat = c(mean_pop = mean(sapply(hp_pop, nrow)),
              mean_inds = mean(sapply(hp_hd_capthists, nrow)),
              mean_dets = mean(sapply(hp_hd_capthists, sum)),
              mean_recaps = mean(sapply(hp_hd_capthists, max)),
              mean_sing_dets = mean(sapply(hp_hd_capthists, function(x) sum(rowSums(x)== 1))))

lp_ld_dat
lp_hd_dat
hp_ld_dat
hp_hd_dat

lp_ld_capthists_ghost = c(lp_ld_capthists,
                          lapply(lp_ld_capthists, function(x) introduce_ghost(x, .1, noccasions = noccasions)),
                          lapply(lp_ld_capthists, function(x) introduce_ghost(x, .2, noccasions = noccasions)),
                          lapply(lp_ld_capthists, function(x) introduce_ghost(x, .3, noccasions = noccasions)))

lp_hd_capthists_ghost = c(lp_hd_capthists,
                          lapply(lp_hd_capthists, function(x) introduce_ghost(x, .1, noccasions = noccasions)),
                          lapply(lp_hd_capthists, function(x) introduce_ghost(x, .2, noccasions = noccasions)),
                          lapply(lp_hd_capthists, function(x) introduce_ghost(x, .3, noccasions = noccasions)))

hp_ld_capthists_ghost = c(hp_ld_capthists,
                          lapply(hp_ld_capthists, function(x) introduce_ghost(x, .1, noccasions = noccasions)),
                          lapply(hp_ld_capthists, function(x) introduce_ghost(x, .2, noccasions = noccasions)),
                          lapply(hp_ld_capthists, function(x) introduce_ghost(x, .3, noccasions = noccasions)))

hp_hd_capthists_ghost = c(hp_hd_capthists,
                          lapply(hp_hd_capthists, function(x) introduce_ghost(x, .1, noccasions = noccasions)),
                          lapply(hp_hd_capthists, function(x) introduce_ghost(x, .2, noccasions = noccasions)),
                          lapply(hp_hd_capthists, function(x) introduce_ghost(x, .3, noccasions = noccasions)))

saveRDS(lp_ld_capthists_ghost, 'Simulations/sim_capthists_lp_ld_prox.rds')
saveRDS(lp_hd_capthists_ghost, 'Simulations/sim_capthists_lp_hd_prox.rds')
saveRDS(hp_ld_capthists_ghost, 'Simulations/sim_capthists_hp_ld_prox.rds')
saveRDS(hp_hd_capthists_ghost, 'Simulations/sim_capthists_hp_hd_prox.rds')

all_mesh = c(mesh, mesh, mesh, mesh)

lp_ld_scr_fits <- lapply(1:length(lp_ld_capthists_ghost), function(x) tryCatch(secr.fit(lp_ld_capthists_ghost[[x]], all_mesh[[x]], model = list(D~cov), detectfn = 'HN', start = c(D = log(lowD), D.value = true_d1, g0 = logit(low_g0), sigma = log(true_sigma))), error = function(e) NULL))

saveRDS(lp_ld_scr_fits, 'Simulations/simulations_scr_lp_ld_g0.rds')

lp_ld_min2_scr_fits <- lapply(1:length(lp_ld_capthists_ghost), function(x) tryCatch(scr2_lik_binom(lp_ld_capthists_ghost[[x]], all_mesh[[x]], model = list(D~cov), startparams = c(D = log(lowD), D.value = true_d1, g0 = logit(low_g0), sigma = log(true_sigma)), simulations = T), error = function(e) NULL))

saveRDS(lp_ld_min2_scr_fits,'Simulations/simulations_scr2_lp_ld_g0.rds')

lp_hd_scr_fits <- lapply(1:length(lp_hd_capthists_ghost), function(x) tryCatch(secr.fit(lp_hd_capthists_ghost[[x]], all_mesh[[x]], model = list(D~cov), detectfn = 'HN', start = c(D = log(lowD), D.value = true_d1,g0 = logit(high_g0), sigma = log(true_sigma))), error = function(e) NULL))

saveRDS(lp_hd_scr_fits,'Simulations/simulations_scr_lp_hd_g0.rds')

lp_hd_min2_scr_fits <- lapply(1:length(lp_hd_capthists_ghost), function(x) tryCatch(scr2_lik_binom(lp_hd_capthists_ghost[[x]], all_mesh[[x]], model = list(D~cov), startparams = c(D = log(lowD), D.value = true_d1,g0 = logit(high_g0), sigma = log(true_sigma)), simulations = T), error = function(e) NULL))

saveRDS(lp_hd_min2_scr_fits,'Simulations/simulations_scr2_lp_hd_g0.rds')

hp_ld_scr_fits <- lapply(1:length(hp_ld_capthists_ghost), function(x) tryCatch(secr.fit(hp_ld_capthists_ghost[[x]],all_mesh[[x]], model = list(D~cov), detectfn = 'HN', start = c(D = log(highD), D.value = true_d1, g0 = logit(low_g0), sigma = log(true_sigma))), error = function(e) NULL))

saveRDS(hp_ld_scr_fits,'Simulations/simulations_scr_hp_ld_g0.rds')

hp_ld_min2_scr_fits <- lapply(1:length(hp_ld_capthists_ghost), function(x) tryCatch(scr2_lik_binom(hp_ld_capthists_ghost[[x]],all_mesh[[x]], model = list(D~cov), startparams = c(D = log(highD), D.value = true_d1, g0 = logit(low_g0), sigma = log(true_sigma)), simulations = T), error = function(e) NULL))

saveRDS(hp_ld_min2_scr_fits,'Simulations/simulations_scr2_hp_ld_g0.rds')

hp_hd_scr_fits <- lapply(1:length(hp_hd_capthists_ghost), function(x) tryCatch(secr.fit(hp_hd_capthists_ghost[[x]],all_mesh[[x]], model = list(D~cov), detectfn = 'HN', start = c(D = log(highD), D.value = true_d1, g0 = logit(high_g0), sigma = log(true_sigma))), error = function(e) NULL))

saveRDS(hp_hd_scr_fits,'Simulations/simulations_scr_hp_hd_g0.rds')

hp_hd_min2_scr_fits <- lapply(1:length(hp_hd_capthists_ghost), function(x) tryCatch(scr2_lik_binom(hp_hd_capthists_ghost[[x]],all_mesh[[x]], model = list(D~cov), startparams = c(D = log(highD), D.value = true_d1, g0 = logit(high_g0), sigma = log(true_sigma)), simulations = T), error = function(e) NULL))

saveRDS(hp_hd_min2_scr_fits,'Simulations/simulations_scr2_hp_hd_g0.rds')
