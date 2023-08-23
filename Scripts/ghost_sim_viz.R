rm(list = ls())

library(secr)
library(tidyverse)

setwd("D:/ar337/2encounter-SCR/ghosts_repo")
source('2Encounter_functions.R')


true_sigma = 300
true_d1 = .1
nsims = 100
lowD = .03
highD = .1
low_lambda0 = 1
high_lambda0 = 2

mesh = readRDS('simulations_mesh.rds')

lp_ld_capthist = readRDS('sim_capthists_lp_ld.rds')
lp_hd_capthist = readRDS('sim_capthists_lp_hd.rds')
hp_ld_capthist = readRDS('sim_capthists_hp_ld.rds')
hp_hd_capthist = readRDS('sim_capthists_hp_hd.rds')

lp_ld_scr = readRDS('simulations_scr_lp_ld.rds')
lp_ld_scr2 = readRDS('simulations_scr2_lp_ld.rds')

# scr2_unfit = which(lengths(lp_ld_scr2) == 0)
# 
# lp_ld_scr2[[scr2_unfit]] = lp_ld_scr2[[1]]
# 
# lp_ld_scr2[[scr2_unfit]]$code  = 5


lp_hd_scr = readRDS('simulations_scr_lp_hd.rds')
lp_hd_scr2 = readRDS('simulations_scr2_lp_hd.rds')
hp_ld_scr = readRDS('simulations_scr_hp_ld.rds')
hp_ld_scr2 = readRDS('simulations_scr2_hp_ld.rds')
hp_hd_scr = readRDS('simulations_scr_hp_hd.rds')
hp_hd_scr2 = readRDS('simulations_scr2_hp_hd.rds')

scr_fit = lp_ld_scr
scr2_fit = lp_ld_scr2
true_D = lowD
true_lambda0 = low_lambda0

bias_dat <- function(scr_fit,
                     scr2_fit,
                     mesh,
                     true_D,
                     true_lambda0, 
                     true_sigma = 300,
                     true_d1 = .1,
                     nsims = 100){
  
  scr_coefs <- data.frame(t(sapply(scr_fit, function(x) coef(x)[,1]))) %>% mutate(mod = 'scr',trueD = true_D,d1 = true_d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  scr2_coefs <- data.frame(t(sapply(scr2_fit, function(x) x$estimate))) %>% mutate(mod = 'scr2',trueD = true_D,d1 = true_d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  colnames(scr_coefs)[1:4] = c('D','D.value','lambda0','sigma')
  colnames(scr2_coefs) = colnames(scr_coefs)
  
  E_N = sapply(mesh, function(x) sum(exp(log(true_D) + true_d1*covariates(x)$cov)*attr(x,'area')))
  
  scr_mod_EN = sapply(1:nrow(scr_coefs), function(x) sum(exp(scr_coefs[x,]$D + scr_coefs[x,]$D.value*covariates(mesh[[((x-1) %% nsims) + 1]])$cov)*attr(mesh[[((x-1) %% nsims) + 1]],'area')))
  scr2_mod_EN = sapply(1:nrow(scr2_coefs), function(x) sum(exp(scr2_coefs[x,]$D + scr_coefs[x,]$D.value*covariates(mesh[[((x-1) %% nsims) + 1]])$cov)*attr(mesh[[((x-1) %% nsims) + 1]],'area')))
  
  scr_mod_conv = sapply(scr_fit, function(x) x$fit$code)
  scr2_mod_conv = sapply(scr2_fit, function(x) x$code)
  
  mod_bias <- rbind(scr_coefs,scr2_coefs) %>%
    mutate(E_N =  rep(E_N,8),
           mod_EN = c(scr_mod_EN,scr2_mod_EN),
           conv = c(scr_mod_conv,scr2_mod_conv)) %>%
    rowwise() %>%
    filter(!is.na(mod_EN)) %>% 
    filter(mod_EN < 5*E_N,conv < 3) %>% 
    mutate(EN_b = (mod_EN - E_N)/abs(E_N),
           D_b = (D- trueD)/abs(trueD),
           d1_b = (D.value - d1)/abs(d1),
           lambda0_b = (exp(lambda0) - true_lambda0)/abs(true_lambda0),
           sigma_b = (exp(sigma) - true_sigma)/abs(true_sigma)) %>%
    group_by(mod,ghost_prop) %>%
    summarise(EN_Bias = mean(EN_b), EN_lcl = EN_Bias - 1.96*sd(EN_b)/sqrt(n()),EN_ucl = EN_Bias + 1.96*sd(EN_b)/sqrt(n()),
              D_Bias = mean(D_b), D_lcl = D_Bias - 1.96*sd(D_b)/sqrt(n()),D_ucl = D_Bias + 1.96*sd(D_b)/sqrt(n()),
              d1_Bias = mean(d1_b), d1_lcl = d1_Bias - 1.96*sd(d1_b)/sqrt(n()),d1_ucl = d1_Bias + 1.96*sd(d1_b)/sqrt(n()),
              lambda0_Bias = mean(lambda0_b), lambda0_lcl = lambda0_Bias - 1.96*sd(lambda0_b)/sqrt(n()),lambda0_ucl = lambda0_Bias + 1.96*sd(lambda0_b)/sqrt(n()),
              sigma_Bias = mean(sigma_b), sigma_lcl = sigma_Bias - 1.96*sd(sigma_b)/sqrt(n()),sigma_ucl = sigma_Bias + 1.96*sd(sigma_b)/sqrt(n())) %>%
    pivot_longer(3:17,names_sep = '_',names_to = c('par','estimate')) %>%
    pivot_wider(names_from = 'estimate',values_from = value) %>%
    filter(!par %in% c('D','d1')) 
  
  mod_bias$par <- factor(mod_bias$par,labels = c(expression('E[N'~']'),bquote(lambda[0]),bquote(sigma)))
  
  return(mod_bias)
}



lp_ld_dat = bias_dat(scr_fit = lp_ld_scr,
                     scr2_fit = lp_ld_scr2,
                     mesh = mesh,
                     true_D = lowD,
                     true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'lp_ld')

lp_hd_dat = bias_dat(scr_fit = lp_hd_scr,
                     scr2_fit = lp_hd_scr2,
                     mesh = mesh,
                     true_D = lowD,
                     true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'lp_hd')

hp_ld_dat = bias_dat(scr_fit = hp_ld_scr,
                     scr2_fit = hp_ld_scr2,
                     mesh = mesh,
                     true_D = highD,
                     true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'hp_ld')

hp_hd_dat = bias_dat(scr_fit = hp_hd_scr,
                     scr2_fit = hp_hd_scr2,
                     mesh = mesh,
                     true_D = highD,
                     true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'hp_hd')

all_sims = bind_rows(lp_ld_dat,lp_hd_dat,hp_ld_dat,hp_hd_dat) 

all_sims$sim_name = factor(all_sims$sim_name, labels = c(expression('D = 0.15'~lambda[0]~'= 2'),
                                                         expression('D = 0.15'~lambda[0]~'= 1'),
                                                         expression('D = 0.03'~lambda[0]~'= 2'),
                                                         expression('D = 0.03'~lambda[0]~'= 1')))

all_sims %>% 
  ggplot(aes(x = ghost_prop, y = Bias, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  scale_shape_discrete(name = 'Model')+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = T)+
  geom_errorbar(aes(ymax = ucl,ymin = lcl),position = position_dodge(width = .05), width = .01,show.legend = F)+
  geom_hline(yintercept = 0)+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(par~sim_name, scales = 'free',labeller = label_parsed)+
  ylab('')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))



######Estimate variance

bias_var_dat <- function(scr_fit,
                         scr2_fit,
                         mesh,
                         true_D,
                         true_lambda0, 
                         true_sigma = 300,
                         true_d1 = .1,
                         nsims = 100){
  
  scr_coefs <- data.frame(t(sapply(scr_fit, function(x) coef(x)[,1]))) %>% mutate(mod = 'scr',trueD = true_D,d1 = true_d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  scr2_coefs <- data.frame(t(sapply(scr2_fit, function(x) x$estimate))) %>% mutate(mod = 'scr2',trueD = true_D,d1 = true_d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  colnames(scr_coefs)[1:4] = c('D','D.value','lambda0','sigma')
  colnames(scr2_coefs) = colnames(scr_coefs)
  
  E_N = sapply(mesh, function(x) sum(exp(log(true_D) + true_d1*covariates(x)$cov)*attr(x,'area')))
  
  scr_mod_EN = sapply(1:nrow(scr_coefs), function(x) sum(exp(scr_coefs[x,]$D + scr_coefs[x,]$D.value*covariates(mesh[[((x-1) %% nsims) + 1]])$cov)*attr(mesh[[((x-1) %% nsims) + 1]],'area')))
  scr2_mod_EN = sapply(1:nrow(scr2_coefs), function(x) sum(exp(scr2_coefs[x,]$D + scr_coefs[x,]$D.value*covariates(mesh[[((x-1) %% nsims) + 1]])$cov)*attr(mesh[[((x-1) %% nsims) + 1]],'area')))
  
  scr_mod_conv = sapply(scr_fit, function(x) x$fit$code)
  scr2_mod_conv = sapply(scr2_fit, function(x) x$code)
  
  mod_var <- rbind(scr_coefs,scr2_coefs) %>%
    mutate(E_N = rep(E_N,8),
           mod_EN = c(scr_mod_EN,scr2_mod_EN),
           conv = c(scr_mod_conv,scr2_mod_conv)) %>%
    rowwise() %>%
    filter(!is.na(mod_EN)) %>% 
    filter(mod_EN < 3*E_N,conv < 3) %>% 
    mutate(EN_v = (mod_EN - E_N)^2,
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
  
  mod_var$par <- factor(mod_var$par,labels = c(expression('E[N'~']'),bquote(lambda[0]),bquote(sigma)))
  
  return(mod_var)
}


lp_ld_dat = bias_var_dat(scr_fit = lp_ld_scr,
                         scr2_fit = lp_ld_scr2,
                         mesh = mesh,
                         true_D = lowD,
                         true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'lp_ld')

lp_hd_dat = bias_var_dat(scr_fit = lp_hd_scr,
                         scr2_fit = lp_hd_scr2,
                         mesh = mesh,
                         true_D = lowD,
                         true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'lp_hd')

hp_ld_dat = bias_var_dat(scr_fit = hp_ld_scr,
                         scr2_fit = hp_ld_scr2,
                         mesh = mesh,
                         true_D = highD,
                         true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'hp_ld')

hp_hd_dat = bias_var_dat(scr_fit = hp_hd_scr,
                         scr2_fit = hp_hd_scr2,
                         mesh = mesh,
                         true_D = highD,
                         true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'hp_hd')



all_sims = bind_rows(lp_ld_dat,lp_hd_dat,hp_ld_dat,hp_hd_dat) 

all_sims$sim_name = factor(all_sims$sim_name, labels = c(expression('D = 0.15'~lambda[0]~'= 2'),
                                                         expression('D = 0.15'~lambda[0]~'= 1'),
                                                         expression('D = 0.03'~lambda[0]~'= 2'),
                                                         expression('D = 0.03'~lambda[0]~'= 1')))

all_sims %>% 
  ggplot(aes(x = ghost_prop, y = variance, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  scale_shape_discrete(name = 'Model')+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = T)+
  geom_hline(yintercept = 0)+
  scale_x_continuous(labels = scales::percent)+
  facet_grid(par~sim_name, scales = 'free',labeller = label_parsed)+
  ylab('MSE')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))

all_sims %>% 
  ggplot(aes(x = ghost_prop, y = log(variance), group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  scale_shape_discrete(name = 'Model')+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = T)+
  scale_x_continuous(labels = scales::percent)+
  facet_grid(par~sim_name, scales = 'free',labeller = label_parsed)+
  ylab('log(MSE)')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))



cap_stats <- function(caphist){
  data.frame(caphist) %>% 
    group_by(ID) %>% 
    mutate(dets = n())  %>%
    ungroup() %>% 
    summarise(detections = n(), individuals = length(unique(ID)), single_detections = sum(dets == 1))
}



p_values_sim <- function(scr2_fit,
                         capthist,
                         mesh,
                         true_D,
                         true_lambda0, 
                         true_sigma = 300,
                         true_d1 = .1,
                         nsims = 100){
  

  scr2_coefs <- data.frame(t(sapply(scr2_fit, function(x) x$estimate))) %>% mutate(mod = 'scr2',trueD = true_D,d1 = true_d1,true_lambda0 = true_lambda0,true_sigma = true_sigma,ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  colnames(scr2_coefs)[1:4] = c('D','D.value','lambda0','sigma')

  distmat = edist(traps(capthist[[1]]),mesh[[1]])
  lambda_x_list_scr2 <- lapply(1:length(capthist), function(x) t(apply(t(distmat), 2,lambda_hhn ,lambda0 = exp(scr2_coefs$lambda0[x]),sigma =  exp(scr2_coefs$sigma[x]))))
  

  scr2_exp_sing <- sapply(1:length(capthist), function(x){
    tryCatch({
      lam <- colSums(lambda_x_list_scr2[[x]])
      sum(exp(scr2_coefs$D[x] + scr2_coefs$D.value[x]*covariates(mesh[[((x-1) %% nsims) + 1]])$cov)*attr(mesh[[((x-1) %% nsims) + 1]],'area')*lam*exp(-lam),na.rm = T)
    },error = function(e) 0)
  })
  
  obs_sing <- t(sapply(capthist, function(x) cap_stats(x)))
  
  
  tab1 <- cbind(obs_sing,scr2_exp_sing = scr2_exp_sing) %>%
    data.frame() %>% 
    unnest() %>% 
    mutate(ghost_prop = rep(c(0,.1,.2,.3),each = nsims)) %>% 
    mutate(scr2_p_value = ppois(as.numeric(single_detections)-1,scr2_exp_sing,lower.tail = F)) 
  
  tab1 = tab1 %>% 
    group_by(ghost_prop) %>% 
    summarise(scr2_hyp = sum(scr2_p_value <= .05)) 
  
  return(tab1)
}

lp_ld_hp = p_values_sim(scr2_fit = lp_ld_scr2,
                        capthist = lp_ld_capthist,
                        mesh = mesh,
                        true_D = lowD,
                        true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'lp_ld')

lp_hd_hp = p_values_sim(scr2_fit = lp_hd_scr2,
                        capthist = lp_hd_capthist,
                        mesh = mesh,
                        true_D = lowD,
                        true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'lp_hd')

hp_ld_hp = p_values_sim(scr2_fit = hp_ld_scr2,
                        capthist = hp_ld_capthist,
                        mesh = mesh,
                        true_D = highD,
                        true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'hp_ld')

hp_hd_hp = p_values_sim(scr2_fit = hp_hd_scr2,
                        capthist = hp_hd_capthist,
                        mesh = mesh,
                        true_D = highD,
                        true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'hp_hd')


all_hyp = rbind(lp_ld_hp,
                lp_hd_hp,
                hp_ld_hp,
                hp_hd_hp)

all_hyp$sim = factor(all_hyp$sim_name, labels = c(expression('D = 0.15'~lambda[0]~'= 2'),
                                             expression('D = 0.15'~lambda[0]~'= 1'),
                                             expression('D = 0.03'~lambda[0]~'= 2'),
                                             expression('D = 0.03'~lambda[0]~'= 1')))


all_hyp %>% 
  ggplot(aes(x = ghost_prop, y = scr2_hyp)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  geom_bar(position = position_dodge(width = .05), show.legend = F, stat = 'identity')+
  geom_hline(yintercept = 95)+
  geom_hline(yintercept = 5,linetype = 2)+
  scale_x_continuous(labels = scales::percent)+
  facet_grid(~sim, scales = 'free',labeller = label_parsed)+
  ylab('')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))


mesh_ex <- mesh[[1]]
trap_ex = traps(capthist[[1]])

ggplot(mesh_ex,aes(x = x, y = y))+
  geom_tile(aes(fill = covariates(mesh_ex)$cov))+
  geom_point(data = trap_ex, shape = 3, col = 'firebrick', stroke = 2)+
  scale_fill_viridis_c(name = 'cov')+
  coord_equal()
