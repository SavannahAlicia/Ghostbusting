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

### function to evaluate number of ghosts created
n_ghosts = function(capthist){
  inds = sapply(capthist,nrow)
  data.frame(split(inds,rep(1:4,each = 100))) %>% 
    pivot_longer(2:4) %>% 
    mutate(ghosts = value-X1) %>% 
    group_by(name) %>% 
    summarise(mean = mean(ghosts),sd = sd(ghosts))
}

n_ghosts(lp_ld_capthist)
n_ghosts(lp_hd_capthist)
n_ghosts(hp_ld_capthist)
n_ghosts(hp_hd_capthist)

### Read model fit objects
lp_ld_scr = readRDS('Simulations/simulations_scr_lp_ld.rds')
lp_ld_scr2 = readRDS('Simulations/simulations_scr2_lp_ld.rds')
lp_hd_scr = readRDS('Simulations/simulations_scr_lp_hd.rds')
lp_hd_scr2 = readRDS('Simulations/simulations_scr2_lp_hd.rds')
hp_ld_scr = readRDS('Simulations/simulations_scr_hp_ld.rds')
hp_ld_scr2 = readRDS('Simulations/simulations_scr2_hp_ld.rds')
hp_hd_scr = readRDS('Simulations/simulations_scr_hp_hd.rds')
hp_hd_scr2 = readRDS('Simulations/simulations_scr2_hp_hd.rds')


### Function to compute bias of estimates
bias_dat <- function(scr_fit,
                     scr2_fit,
                     mesh,
                     true_D,
                     true_lambda0, 
                     true_sigma = 300,
                     true_d1 = .1,
                     nsims = 100){
  
  ### Extract maximum like
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

all_sims$sim_name = factor(all_sims$sim_name, labels = c(expression('D = 0.06'~lambda[0]~'= 2'),
                                                         expression('D = 0.06'~lambda[0]~'= 1'),
                                                         expression('D = 0.03'~lambda[0]~'= 2'),
                                                         expression('D = 0.03'~lambda[0]~'= 1')))


################################
#### CI COverage


EN_f = function(D,mesh){
  sum(exp(D[1] + D[2]*covariates(mesh)$cov)*attr(mesh,'area'))
}


CI_coverage_dat <- function(scr_fit,
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
  
  
  scr_EN_dat = cbind(coef = 'EN',t(sapply(scr_fit, function(x) region.N(x)[1,]))[,c('estimate','lcl','ucl')])
  
  scr2_EN = sapply(1:nrow(scr2_coefs), function(x) EN_f(as.numeric(scr2_coefs[x,1:2]),mesh[[((x-1) %% nsims) + 1]]))  
  scr2_vcv = lapply(scr2_fit, function(x) solve(x$hessian))
  scr2_grad = t(sapply(1:nrow(scr2_coefs),function(x) nlme::fdHess(as.numeric(scr2_coefs[x,1:2]),EN_f, mesh = mesh[[((x-1) %% nsims) + 1]])$gradient))
  scr2_EN_se = sapply(1:nrow(scr2_coefs), function(x) sqrt(scr2_grad[x,]%*%scr2_vcv[[x]][1:2,1:2]%*%scr2_grad[x,]))
  scr2_EN_df = data.frame(estimate = scr2_EN,SE.estimate = scr2_EN_se)
  scr2_EN_dat = cbind(coef = 'EN',secr:::add.cl(df = scr2_EN_df, alpha = .05,loginterval = T)[,c('estimate','lcl','ucl')])
  
  scr_cis = lapply(1:length(scr_fit), function(x) cbind(rbind(scr_EN_dat[x,],predict(scr_fit[[x]])[2:3,c(2,4,5)] %>% rownames_to_column(var = 'coef')),true = c(E_N[((x-1) %% nsims) + 1],true_lambda0,true_sigma)))
  scr2_cis = lapply(1:nrow(scr2_coefs), function(x){
    se = sqrt(diag(scr2_vcv[[x]]))[3:4]
    cbind(
      rbind(scr2_EN_dat[x,],
            data.frame(coef = c('lambda0','sigma'),
                       estimate = exp(as.numeric(scr2_coefs[x,3:4])),
                       lcl = exp(as.numeric(scr2_coefs[x,3:4]) - 1.96*se),
                       ucl = exp(as.numeric(scr2_coefs[x,3:4]) + 1.96*se))),
      true = c(E_N[((x-1) %% nsims) + 1],true_lambda0,true_sigma))
  })
  
  scr_ci_coverage = data.frame(t(sapply(scr_cis,function(df){
    as.numeric((df$true > df$lcl) & (df$true < df$ucl))
  }))) %>%  mutate(mod = 'scr',ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  scr2_ci_coverage = data.frame(t(sapply(scr2_cis,function(df){
    as.numeric((df$true > df$lcl) & (df$true < df$ucl))
  }))) %>%  mutate(mod = 'scr2',ghost_prop = rep(c(0,.1,.2,.3),each = nsims))
  
  colnames(scr_ci_coverage)[1:3] = colnames(scr2_ci_coverage)[1:3] = c('EN','lambda0','sigma')
  
  scr_mod_conv = sapply(scr_fit, function(x) x$fit$code)
  scr2_mod_conv = sapply(scr2_fit, function(x) x$code)
  
  all_dat = bind_rows(scr_ci_coverage %>% 
              mutate(conv = scr_mod_conv) %>% 
              filter(conv < 3) %>% 
              pivot_longer(1:3,names_to = 'par') %>% 
              group_by(ghost_prop,par) %>% 
              summarise(cic = round(sum(value)/n()*100,0)) %>% 
              mutate(mod = 'scr'),
            
            scr2_ci_coverage %>% 
              mutate(conv = scr2_mod_conv) %>% 
              filter(conv < 3) %>% 
              pivot_longer(1:3,names_to = 'par') %>% 
              group_by(ghost_prop,par) %>% 
              summarise(cic = round(sum(value)/n()*100,0))%>% 
              mutate(mod = 'scr2'))
  all_dat$par <- factor(all_dat$par,labels = c(expression('E[N'~']'),bquote(lambda[0]),bquote(sigma)))
  
  return(all_dat)
}

lp_ld_ci = CI_coverage_dat(scr_fit = lp_ld_scr,
                         scr2_fit = lp_ld_scr2,
                         mesh = mesh,
                         true_D = lowD,
                         true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'lp_ld')

lp_hd_ci = CI_coverage_dat(scr_fit = lp_hd_scr,
                         scr2_fit = lp_hd_scr2,
                         mesh = mesh,
                         true_D = lowD,
                         true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'lp_hd')

hp_ld_ci = CI_coverage_dat(scr_fit = hp_ld_scr,
                         scr2_fit = hp_ld_scr2,
                         mesh = mesh,
                         true_D = highD,
                         true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'hp_ld')

hp_hd_ci = CI_coverage_dat(scr_fit = hp_hd_scr,
                         scr2_fit = hp_hd_scr2,
                         mesh = mesh,
                         true_D = highD,
                         true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'hp_hd')



all_ci = bind_rows(lp_ld_ci,lp_hd_ci,hp_ld_ci,hp_hd_ci) 

all_ci$sim_name = factor(all_ci$sim_name, labels = c(expression('D = 0.06'~lambda[0]~'= 2'),
                                                         expression('D = 0.06'~lambda[0]~'= 1'),
                                                         expression('D = 0.03'~lambda[0]~'= 2'),
                                                         expression('D = 0.03'~lambda[0]~'= 1')))

all_ci %>% 
  ggplot(aes(x = ghost_prop, y = cic/100, group = mod, fill = mod)) +
  scale_fill_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  geom_bar(position = position_dodge(width = .05),  show.legend = T, stat = 'identity', width = .04)+
  geom_hline(yintercept = .95, linetype = 2)+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(par~sim_name,labeller = label_parsed)+
  ylab('Confidence Interval Coverage')+
  xlab('')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))

ggsave('Simulations/cic_sim_results.png',width = 8, height = 5, dpi =300)

bias_plot = all_sims %>% 
  ggplot(aes(x = ghost_prop, y = Bias, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  scale_shape_discrete(name = 'Model')+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = T)+
  geom_errorbar(aes(ymax = ucl,ymin = lcl),position = position_dodge(width = .05), width = .01,show.legend = F)+
  geom_hline(yintercept = 0)+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(par~sim_name, scales = 'free',labeller = label_parsed)+
  ylab('Relative Bias')+
  xlab('')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))

######Estimate MSE

bias_mse_dat <- function(scr_fit,
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
  
  mod_mse <- rbind(scr_coefs,scr2_coefs) %>%
    mutate(E_N = rep(E_N,8),
           mod_EN = c(scr_mod_EN,scr2_mod_EN),
           conv = c(scr_mod_conv,scr2_mod_conv)) %>%
    rowwise() %>%
    filter(!is.na(mod_EN)) %>% 
    filter(mod_EN < 5*E_N,conv < 3) %>% 
    group_by(mod,ghost_prop) %>%
    summarise(EN_mse = mean((mod_EN - E_N)^2),
              D_mse = mean((D- trueD)^2),
              d1_mse = mean((D.value - d1)^2),
              lambda0_mse = mean((exp(lambda0) - true_lambda0)^2),
              sigma_mse = mean((exp(sigma) - true_sigma)^2),
              
              EN_var = var(mod_EN),
              D_var = var(D),
              d1_var = var(D.value),
              lambda0_var = var(exp(lambda0)),
              sigma_var = var(exp(sigma))) %>%
    pivot_longer(3:12,names_sep = '_',names_to = c('par','estimate')) %>%
    pivot_wider(names_from = 'estimate',values_from = value)  %>%
    filter(!par %in% c('D','d1'))
  
  mod_mse$par <- factor(mod_mse$par,labels = c(expression('E[N'~']'),bquote(lambda[0]),bquote(sigma)))
  
  return(mod_mse)
}


lp_ld_dat = bias_mse_dat(scr_fit = lp_ld_scr,
                         scr2_fit = lp_ld_scr2,
                         mesh = mesh,
                         true_D = lowD,
                         true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'lp_ld')

lp_hd_dat = bias_mse_dat(scr_fit = lp_hd_scr,
                         scr2_fit = lp_hd_scr2,
                         mesh = mesh,
                         true_D = lowD,
                         true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'lp_hd')

hp_ld_dat = bias_mse_dat(scr_fit = hp_ld_scr,
                         scr2_fit = hp_ld_scr2,
                         mesh = mesh,
                         true_D = highD,
                         true_lambda0 = low_lambda0) %>% 
  mutate(sim_name = 'hp_ld')

hp_hd_dat = bias_mse_dat(scr_fit = hp_hd_scr,
                         scr2_fit = hp_hd_scr2,
                         mesh = mesh,
                         true_D = highD,
                         true_lambda0 = high_lambda0) %>% 
  mutate(sim_name = 'hp_hd')



all_sims = bind_rows(lp_ld_dat,lp_hd_dat,hp_ld_dat,hp_hd_dat) 

all_sims$sim_name = factor(all_sims$sim_name, labels = c(expression('D = 0.06'~lambda[0]~'= 2'),
                                                         expression('D = 0.06'~lambda[0]~'= 1'),
                                                         expression('D = 0.03'~lambda[0]~'= 2'),
                                                         expression('D = 0.03'~lambda[0]~'= 1')))

all_sims %>% 
  ggplot(aes(x = ghost_prop, y = mse, group = mod, col = mod, shape = mod)) +
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

log_mse_plot = all_sims %>% 
  ggplot(aes(x = ghost_prop, y = log(mse), group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  scale_shape_discrete(name = 'Model')+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = T)+
  scale_x_continuous(labels = scales::percent)+
  facet_grid(par~sim_name, scales = 'free',labeller = label_parsed)+
  ylab('log(MSE)')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))

var_plot = all_sims %>% 
  ggplot(aes(x = ghost_prop, y = var, group = mod, col = mod, shape = mod)) +
  scale_color_discrete(name = 'Model',type = c('firebrick','navyblue'))+
  scale_shape_discrete(name = 'Model')+
  geom_point(position = position_dodge(width = .05), size = 3, show.legend = T)+
  scale_x_continuous(labels = scales::percent)+
  facet_grid(par~sim_name, scales = 'free',labeller = label_parsed)+
  ylab('Variance')+
  xlab('')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'))

ggarrange(bias_plot,log_mse_plot,nrow = 2, common.legend = T, legend =  'bottom', labels = c('A','B'))

ggsave('Simulations/sim_results.png',width = 8, height = 10, dpi =300)




#########################


sim_capt_summaries <- function(pop,capthist){
  pop_sim = sapply(pop,nrow)
  
  total_dets = sapply(capthist,sum)
  individuals_dets = sapply(capthist,nrow)
  single_dets = sapply(capthist,function(x) sum(rowSums(x) == 1))
  
  data.frame(total_dets = total_dets,individuals_dets = individuals_dets,single_dets = single_dets, pop = pop_sim, ghost_prop = rep(c(0,.1,.2,.3), each = 100)) %>% 
    pivot_longer(1:4) %>% 
    group_by(ghost_prop,name) %>% 
    summarise(mean = mean(value), sd = sd(value)) %>% 
    mutate(est = paste0(round(mean,1),' (',round(sd,1),')')) %>% 
    dplyr::select(ghost_prop,name,est) %>% 
    pivot_wider(names_from = name,values_from = est) %>% 
    dplyr::select(ghost_prop,pop,total_dets,individuals_dets,single_dets) %>% 
    data.frame()
}

tab1 = rbind(sim_capt_summaries(lp_pop, lp_ld_capthist) %>% mutate(scenario = 'lp_ld'),
      sim_capt_summaries(lp_pop, lp_hd_capthist) %>% mutate(scenario = 'lp_hd'),
      sim_capt_summaries(hp_pop, hp_ld_capthist) %>% mutate(scenario = 'hp_ld'),
      sim_capt_summaries(hp_pop, hp_hd_capthist) %>% mutate(scenario = 'hp_hd'))

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

all_hyp$sim = factor(all_hyp$sim_name, labels = c(expression('D = 0.06'~lambda[0]~'= 2'),
                                                  expression('D = 0.06'~lambda[0]~'= 1'),
                                                  expression('D = 0.03'~lambda[0]~'= 2'),
                                                  expression('D = 0.03'~lambda[0]~'= 1')))

all_hyp$type = ifelse(all_hyp$ghost_prop == 0, 'False','True')



hyp_test_lambda0 <- all_hyp %>% 
  ggplot(aes(x = ghost_prop, y = scr2_hyp)) +
  geom_bar(aes(fill  = type), position = position_dodge(width = .05), stat = 'identity',col = 'black')+
  scale_fill_viridis_d(option = 'H', direction = -1,name = 'Positives')+
  geom_hline(yintercept = 95,linetype = 2)+
  geom_hline(yintercept = 5,linetype = 2)+
  scale_x_continuous(labels = scales::percent)+
  facet_grid(~sim, scales = 'free',labeller = label_parsed)+
  ylab('# Rejected')+
  xlab('Proportion of Ghosts introduced')+
  theme_bw()+
  theme(text = element_text(size = 15),strip.text.y = element_text(angle = 0, face = 'bold'),plot.margin = unit(c(1,0,0,1),units = 'cm'),legend.position = 'bottom')

hyp_test_lambda0
ggsave('Simulations/hypothesis_test.png', width = 8, height = 5, dpi = 300, units = 'in')

mesh_ex <- mesh[[1]]
trap_ex = traps(lp_hd_capthist[[1]])

ggplot(mesh_ex,aes(x = x, y = y))+
  geom_tile(aes(fill = covariates(mesh_ex)$cov))+
  geom_point(data = trap_ex, shape = 3, col = 'firebrick', stroke = 1)+
  scale_fill_viridis_c(name = 'Spatial\nCovariate')+
  coord_equal()+
  theme_classic()

ggsave('Simulations/mesh_realisation.png',dpi = 300, width = 6, height = 3, units = 'in')

