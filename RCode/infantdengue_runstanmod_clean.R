#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  
  # Model options
  include_ade=1 # ADE among children born to seropositive mothers?
  include_sevagefrailty=0 # Differece in severe disease odds by age
  include_infagefrailty=1 # Difference in infection hazard by age
  include_agereporting=0 # Difference in reporting odds by age
  include_maternalprotection=1 # Difference in infection hazard between infants born to seropositive and seronegative mothers?
  alphaconst=0 # Reporting odds constant over time
  alphalin=0 # Reporting odds log-linear with calendar year?
  filesuffix_extra=''
  foi_dir = './Data/02-foi/stanfit_v3/t51_rta_async_noType_neglectZika_final/'
  foi_type = 't51_rta_async_noType_neglectZika_final'
  numiter=1500
  lambda_uncertainty = 0 # Incorporate uncertainty in FOI and maternal seroprevalence?
  
} else {
  for (iarg in 1:length(args)) { 
    eval (parse (text = args[[iarg]] )) 
  }
}
  

require(dplyr)
require(tidyverse)
require(rstan)

all_den = read.csv('./Data/02-analysisData/DengueCasesByStateYear_2000to2024.csv')

# Restrict to 2000-2013
all_den = all_den %>% arrange(state,year_not) %>% filter(year_not<=2014)
  
statenums = unique(all_den$state)

# Get maternal seroprevalence and FOI
maternal_sp_dat = array(NA,dim=c(2014-2000+1,length(statenums),3),dimnames=list(2000:2014,
                                                                                statenums,
                                                                                c("SP","Smult","S1")))

foi_dat = array(NA,dim=c(2014-2000+1,length(statenums)),dimnames=list(2000:2014,
                                                                      statenums))

for (s in statenums) {
    
  if (file.exists(paste0(foi_dir,'S_mother/',s,'.csv'))) {
    
    t = read.csv(paste0(foi_dir,'S_mother/',s,'.csv'))
    
    t = t %>% filter(year>=2000 & year<=2014) %>% mutate(SP = Smult+S1) %>%
      group_by(year) %>% dplyr::summarise(Smult=mean(Smult),SP=mean(SP),S1=mean(S1))
    
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("SP")] = t$SP
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("Smult")] = t$Smult
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("S1")] = t$S1
    
  }
  
  if (file.exists(paste0(foi_dir,'lambda_t/',s,'.csv'))) {
    
    t = read.csv(paste0(foi_dir,'lambda_t/',s,'.csv'))
    t = t %>% filter(year>=2000 & year<=2014) %>%
      group_by(year) %>% dplyr::summarise(FOI = mean(val))
    
    foi_dat[,dimnames(foi_dat)[[2]]==s] = t$FOI
    
  }
  
}

all_den$maternal_sp = sapply(1:nrow(all_den),function(x) (maternal_sp_dat[dimnames(maternal_sp_dat)[[1]]==all_den$year_not[x],
                                                                          dimnames(maternal_sp_dat)[[2]]==all_den$state[x],
                                                                          c("SP")]))
all_den$maternal_s1 = sapply(1:nrow(all_den),function(x) (maternal_sp_dat[dimnames(maternal_sp_dat)[[1]]==all_den$year_not[x],
                                                                          dimnames(maternal_sp_dat)[[2]]==all_den$state[x],
                                                                          c("S1")]))
all_den$maternal_smult = sapply(1:nrow(all_den),function(x) (maternal_sp_dat[dimnames(maternal_sp_dat)[[1]]==all_den$year_not[x],
                                                                             dimnames(maternal_sp_dat)[[2]]==all_den$state[x],
                                                                             c("Smult")]))
all_den$foi = sapply(1:nrow(all_den),function(x) (foi_dat[dimnames(foi_dat)[[1]]==all_den$year_not[x],
                                                          dimnames(foi_dat)[[2]]==all_den$state[x]]))

# Parameters to output
namepars = c("mat_protection_decay",
             "baseline_hazmatprot_sn",
             "baseline_hazmatprot_sp",
             "baseline_hazmult_sn",
             "baseline_hazmult_sp",
             "baseline_replogor",
             "oneyear_probsev",
             "baseline_sevlogor",
             "adepeak_logor",
             "ademode",
             "adekappa")

if (lambda_uncertainty==1) {

  model_code = stanc(file="./RCode/fit_infantdengue_lambdauncertainty_modeloption_clean.stan")
  pars_to_output = c(namepars,'alphast_param','lambda_st','msp_st')
  
} else {
  
  model_code = stanc(file="./RCode/fit_infantdengue_nolambdauncertainty_modeloption_clean.stan")
  pars_to_output = c(namepars,'alphast_param')
  
}
compiled_model_working = stan_model(stanc_ret = model_code)

options(mc.cores = 4) # checks number of cores without having later to specify the cores argument
rstan_options(auto_write = TRUE) # extended packages to use stan

runmodel = function(include_ade,include_sevagefrailty,include_infagefrailty,include_agereporting,include_maternalprotection,
                    lambda_uncertainty,
                    alphaconst,alphalin,numiter,filesuffix_extra) {
  
  filesuffix = paste0('_ade_',include_ade,
                      '_sevage_',include_sevagefrailty,
                      '_infage_',include_infagefrailty,
                      '_repage_',include_agereporting,
                      '_matprot_',include_maternalprotection,
                      '_alphaconst_',alphaconst,
                      '_alphalin_',alphalin,
                      '_foitype_',foi_type,
                      '_lambdauncertainty_',lambda_uncertainty,
                      filesuffix_extra)
  
  if (lambda_uncertainty==1) {
    
    lambda_prior = read.csv(paste0(foi_dir,'FOI_Prior2MixNorm.csv'))
    lambda_prior_bounds = read.csv(paste0(foi_dir,'FOI_PriorBounds.csv'))
    msp_prior = read.csv(paste0(foi_dir,'MSP_Prior2MixNorm.csv'))
    msp_prior_bounds = read.csv(paste0(foi_dir,'MSP_PriorBounds.csv'))
    
    data_input = list(
      nmonth = 12,
      nst = nrow(all_den),
      ns = length(unique(all_den$state)),
      ny = 2014-2000+1,
      y_st_bymonth=as.matrix(all_den[,paste0('infantcases_',0:11,'month')]),
      n_st_bymonth=as.matrix(all_den[,paste0('infantbirths_',1:12,'month')]),
      sev_st_bymonth=as.matrix(all_den[,paste0('infantseverecases_',0:11,'month')]),
      lambda_prior=lambda_prior,
      lambda_prior_lbounds=lambda_prior_bounds[,1],
      lambda_prior_ubounds=lambda_prior_bounds[,2],
      msp_prior=msp_prior,
      msp_prior_lbounds=msp_prior_bounds[,1],
      msp_prior_ubounds=msp_prior_bounds[,2],
      alpha_constant_t = alphaconst,
      alpha_linear = alphalin,
      include_ade = include_ade,
      include_sevagefrailty = include_sevagefrailty,
      include_infagefrailty = include_infagefrailty,
      include_agereporting = include_agereporting,
      include_maternalprotection = include_maternalprotection
    )
  } else {
    data_input = list(
      nmonth = 12,
      nst = nrow(all_den),
      ns = length(unique(all_den$state)),
      ny = 2014-2000+1,
      y_st_bymonth=as.matrix(all_den[,paste0('infantcases_',0:11,'month')]),
      n_st_bymonth=as.matrix(all_den[,paste0('infantbirths_',1:12,'month')]),
      sev_st_bymonth=as.matrix(all_den[,paste0('infantseverecases_',0:11,'month')]),
      lambda_st=all_den$foi,
      msp1_st=all_den$maternal_s1,
      mspmult_st=all_den$maternal_smult,
      alpha_constant_t = alphaconst,
      alpha_linear = alphalin,
      include_ade = include_ade,
      include_sevagefrailty = include_sevagefrailty,
      include_infagefrailty = include_infagefrailty,
      include_agereporting = include_agereporting,
      include_maternalprotection = include_maternalprotection
    )
  }
  
  if (file.exists(paste0('./infantmodel',filesuffix,'.rds'))) {
    infant_model = readRDS(paste0('./infantmodel',filesuffix,'.rds'))
  } else {
    
    system.time(
      infant_model <-
        sampling(
          compiled_model_working,
          data = data_input,
          chains = 4,
          iter = numiter,
          # warmup = ,
          thin = 5,
          control = list(adapt_delta = 0.80),
          include=TRUE,
          pars=pars_to_output
        )
    )
    infant_model@stanmodel@dso <- new("cxxdso")
    saveRDS(infant_model,paste0('./infantmodel',filesuffix,'.rds'))
  }
  
}

runmodel(include_ade,include_sevagefrailty,include_infagefrailty,include_agereporting,include_maternalprotection,
         lambda_uncertainty,
         alphaconst,alphalin,numiter,filesuffix_extra)
