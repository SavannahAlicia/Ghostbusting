library(secr)
library(tidyverse)
library(raster)
library(ggpubr)
library(patchwork)
library(rlang)
library(ggspatial)

source('Scripts/2Encounter_functions.R')

### Read data
mesh = read_rds('Data/mesh.rds')
all_traps = read_rds('Data/traps.rds')
caps = read.csv('Data/captures.csv')

site_names <- c("Munkhkhairkhan", "Tost")   

#### Areas of each site
data.frame(sites = site_names, area = sapply(mesh, function(x) spacing(x)^2 * nrow(x)/1000^2))

#### Make capture history of data k = 1
Mongolia1 <- make.capthist(captures = caps, traps = all_traps, bysession = T)

summary(Mongolia1, terse = T)

#### Make capture history of data k  = 2
Mongolia2 <- caps  %>% 
  group_by(ID) %>% 
  filter(n()>1) %>% 
  data.frame() %>% 
  make.capthist(captures = ., traps = all_traps, bysession = T) 

summary(Mongolia2, terse = T)

###### Secr models
scr_mods = lapply(1:length(site_names), function(x) secr.fit(capthist = Mongolia1[[x]], mask = mesh[[x]], model = list(D~value), detectfn = 'HHN'))

scr2_mods = lapply(1:length(site_names), function(x) scr2_lik(Mongolia1[[x]], mesh[[x]], model = list(D~value), hessian = T))

###### Extract estimates
scr2_abund_f = function(phi,mesh){
  b0 = phi[1]
  b1 = phi[2]
  sum(exp(b0 + b1 * covariates(mesh)$value) * attr(mesh,'area'), na.rm = T)  
}

scr2_abund_se_f = function(phi,mesh,vcov){
  grad2 <- nlme::fdHess(pars = phi, function(x) scr2_abund_f(x,mesh))$grad
  scr2_se <- sqrt(t(grad2) %*% vcov %*% grad2)
}

scr2_coefs_f = function(fit,mesh){
  EN = scr2_abund_f(fit$estimate[1:2],mesh)
  EN_se = scr2_abund_se_f(fit$estimate[1:2], mesh, solve(fit$hessian)[1:2,1:2])  
  theta = exp(fit$estimate[3:4])
  se = sqrt(diag(solve(fit$hessian)))[3:4] * theta
  return(data.frame(coef = c('EN', 'lambda0', 'sigma'), estimate = c(EN,theta), SE.estimate = c(EN_se,se)))
}

scr_coefs_f = function(fit){
  EN = region.N(fit)[1,1:2]
  theta = predict(scr_mods[[1]])[2:3,2:3]
  (rbind(EN,theta) %>% mutate(coef = c('EN', 'lambda0', 'sigma')) )[,c(3,1,2)]
}

scr_coefs = lapply(1:length(site_names), function(x) scr_coefs_f(scr_mods[[x]]) %>% mutate(site = site_names[x], mod = 'scr'))
scr2_coefs = lapply(1:length(site_names), function(x) scr2_coefs_f(scr2_mods[[x]],mesh = mesh[[x]]) %>% mutate(site = site_names[x], mod = 'scr2'))

est_dat = do.call(rbind,c(scr_coefs,scr2_coefs)) %>% 
  mutate(CV = SE.estimate/estimate)
rownames(est_dat) <- NULL

est_dat

cap_stats <- function(caphist){
  data.frame(caphist) %>% 
    group_by(ID) %>% 
    mutate(dets = n())  %>%
    ungroup() %>% 
    summarise(detections = n(), individuals = length(unique(ID)), single_detections = sum(dets == 1))
}

###capture history  summaries
captures_table <- data.frame(t(sapply(Mongolia1, cap_stats))) 

mk = ggplot(mesh[[1]], aes(x = x, y = y)) +
  geom_tile(aes(fill = covariates(mesh[[1]])$value)) +
  scale_fill_viridis_c(name = 'Occupancy\nprobability', limits = c(0,1)) +
  geom_point(data = traps(Mongolia1)[[1]], shape = 3, col = 'firebrick',stroke = 1.5) +
  coord_equal() +
  theme_void() +
  annotation_scale()

tost = ggplot(mesh[[2]], aes(x = y, y = x)) +
  geom_tile(aes(fill = covariates(mesh[[2]])$value)) +
  scale_fill_viridis_c(name = 'Occupancy\nprobability', limits = c(0,1)) +
  geom_point(data = traps(Mongolia1)[[2]], shape = 3, col = 'firebrick',stroke = 1.5) +
  coord_equal() +
  theme_void() +
  annotation_scale()

ggarrange(mk, tost, nrow = 1,labels = c('A','B'), common.legend = T, legend = 'right')

ggsave('CaseStudy/site_plots.png', dpi = 300, width = 9, height = 6)
