setwd('C:/Users/mhitchings/Dropbox (UFL)/DengueInfants/')

require(dplyr)
require(tidyverse)
require(loo)
require(grid)
require(gridExtra)
require(arrow)
require(rstan)
require(ggplot2); theme_set(theme_classic())

include_1day = 0

include_maternalprotection=1
alphaconst=0
alphalin=0
include_reportingdiffsev=0
include_vtrans=0
alphaprior=0
model_type=0
zeroday_missing=0
negbin=0
foi_dir = './Data/02-foi/stanfit_v3/31_rta_async_sevShift_long/'
foi_type='31_rta_async_sevShift_long'
filesuffix_extra=''
lambda_uncertainty = 0

pred_cases = function(par,
                      nmonth,
                      nst,
                      ns,
                      ny,
                      y_st_bymonth,
                      n_st_bymonth,
                      sev_st_bymonth,
                      y_0d_st,
                      sev_0d_st,
                      lambda_st,
                      msp_st,
                      msp1_st,
                      mspmult_st,
                      mpinfsn_st,
                      mpinfs1_st,
                      mpinfmult_st,
                      model_type,
                      ade_dist,
                      include_vtrans,
                      vtrans_prob,
                      prop_missage,
                      fixedpar) {
  
  expit = function(x) {1/(1+exp(-x))}
  logit = function(x) {log(x/(1-x))}
  
  if (length(intersect(names(par),names(fixedpar)))>0) {stop("Repetition of parameters in par and fixedpar")}
  
  allpars = c(par,fixedpar)
  
  mat_protection_decay = allpars[['mat_protection_decay']]
  baseline_hazmatprot_sn = allpars[['baseline_hazmatprot_sn']]
  baseline_hazmatprot_sp = allpars[['baseline_hazmatprot_sp']]
  baseline_hazinfmult = allpars[['baseline_hazinfmult']]
  baseline_hazinfmult = allpars[['baseline_hazinfmult']]
  ageinfection_decay = allpars[['ageinfection_decay']]
  oneyear_probsev = allpars[['oneyear_probsev']]
  titersev_mult = allpars[['titersev_mult']]
  titersev_meanlog = allpars[['titersev_meanlog']]
  titersev_sdlog = allpars[['titersev_sdlog']]
  
  baseline_replogor = allpars[['baseline_replogor']]
  agereporting_decay = allpars[['agereporting_decay']]
  baseline_sevlogor = allpars[['baseline_sevlogor']]
  ageseverity_decay = allpars[['ageseverity_decay']]
  reportingsev_logor = allpars[['reportingsev_logor']]
  
  titersevmulti_meanlog = titersev_meanlog
  titersevmulti_sdlog=titersev_sdlog
  
  alphast = rep(NA,nst)
  for (i in 1:nst) {
    alphast[i] = allpars[[paste0('alpha',i)]]
  }
  
  lambda_mat = matrix(lambda_st,nrow=nst,ncol=nmonth,byrow=F)
  msp_mat = matrix(msp_st,nrow=nst,ncol=nmonth,byrow=F)
  msp1_mat = matrix(msp1_st,nrow=nst,ncol=nmonth,byrow=F)
  mspmult_mat = matrix(mspmult_st,nrow=nst,ncol=nmonth,byrow=F)
  mpinfsn_mat=matrix(mpinfsn_st,nrow=nst,ncol=nmonth,byrow=F)
  mpinfs1_mat=matrix(mpinfs1_st,nrow=nst,ncol=nmonth,byrow=F)
  mpinfmult_mat=matrix(mpinfmult_st,nrow=nst,ncol=nmonth,byrow=F)
  
  lambdaprevyear_mat = lambda_mat
  lambdaprevyear_mat[2:nst,] = lambda_mat[1:(nst-1),]
  lambdaprevyear_mat[seq(1,nst,by=ny),] = lambda_mat[seq(1,nst,by=ny),]
  
  n_st_bymonth_prevyear = n_st_bymonth
  n_st_bymonth_prevyear[2:nst,] = n_st_bymonth[1:(nst-1),]
  n_st_bymonth_prevyear[seq(1,nst,by=ny),] = n_st_bymonth[seq(1,nst,by=ny),]
  
  cumhaz_sp=matrix(NA,nrow=nst,ncol=nmonth)
  cumhaz_sn=matrix(NA,nrow=nst,ncol=nmonth)
  
  denom_st_bymonth=matrix(NA,nrow=nst,ncol=nmonth)
  
  Ms = 0:(nmonth-1)
  
  hr_maternalprotection_sn = matrix(1 - baseline_hazmatprot_sn * exp(-mat_protection_decay * Ms),nrow=nst,ncol=nmonth,byrow=T)
  hr_maternalprotection_sp = matrix(1 - baseline_hazmatprot_sp * exp(-mat_protection_decay * Ms),nrow=nst,ncol=nmonth,byrow=T)
  hr_agefrailty = matrix(1 + baseline_hazinfmult * exp(-ageinfection_decay * Ms),nrow=nst,ncol=nmonth,byrow=T)
  
  p_sev_sn = matrix(expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*Ms)),nrow=nst,ncol=nmonth,byrow=T)
  
  # Note that if model_type is 0 or 1, p_sev_s1 = p_sev_sm
  if (ade_dist==0) {
    p_sev_s1 = matrix(expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dlnorm(Ms,titersev_meanlog, titersev_sdlog)),
                      nrow=nst,ncol=nmonth,byrow=T)
    
    p_sev_sm = matrix(expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dlnorm(Ms,titersevmulti_meanlog, titersevmulti_sdlog)),
                      nrow=nst,ncol=nmonth,byrow=T)
  } else {
    p_sev_s1 = matrix(expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dnorm(Ms,titersev_meanlog, titersev_sdlog)),
                      nrow=nst,ncol=nmonth,byrow=T)
    
    p_sev_sm = matrix(expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dnorm(Ms,titersevmulti_meanlog, titersevmulti_sdlog)),
                      nrow=nst,ncol=nmonth,byrow=T)
  }
  
  p_rep_mild = expit(logit(matrix(alphast,nrow=nst,ncol=nmonth,byrow=F)) + 
                       matrix(baseline_replogor * exp(-agereporting_decay * Ms),nrow=nst,ncol=nmonth,byrow=T))
  p_rep_sev = expit(logit(matrix(alphast,nrow=nst,ncol=nmonth,byrow=F)) + 
                      matrix(baseline_replogor * exp(-agereporting_decay * Ms),nrow=nst,ncol=nmonth,byrow=T) + reportingsev_logor)
  
  # Calculate probability of vertical transmission in seronegative and seropositive
  p_sn_vt = vtrans_prob * mpinfsn_mat
  p_sp_vt = vtrans_prob * (mpinfs1_mat * mspmult_mat + mpinfmult_mat * mspmult_mat) / mspmult_mat
  p_vt = p_sn_vt * (1 - mspmult_mat) + p_sp_vt * mspmult_mat
  
  # We will need the total children born in each calendar year
  n_byyear = apply(n_st_bymonth,1,sum)
  
  for (M in 0:(nmonth-1)) {
    
    # Build up cumulative hazard by state-year and month of age
    cumhaz_sp[,M+1] = 0
    cumhaz_sn[,M+1] = 0
    
    # Need to calculate the denominator (who is at risk of having case at age M in year t)
    # For this, we add up all children who could have been M months old during year t
    # i.e. children who were born in months 12-M+1 to 12 of the previous year plus children
    # born in months 1 to 12-M of the current year
    denom_st_bymonth[,M+1] = 0
    
    # Need to calculate cumulative hazard up to age M for a child who could have experienced
    # dengue infection at age M in year t
    # For this, need to loop over children born in calendar month 12-M+1 of the previous year to 
    # calendar month 12-M of the current year
    for (c in 1:12) {
      if (c>12-M) {
        
        # Children born in months 12-M+1 to 12 of the previous year were of age M during year t
        denom_st_bymonth[,M+1] = denom_st_bymonth[,M+1] + n_st_bymonth_prevyear[,c]
        
      } else {
        # Children who were born in months 1 to 12-M were of age M during year t
        denom_st_bymonth[,M+1] = denom_st_bymonth[,M+1] + n_st_bymonth[,c]
        
      }
      
      cumhaz_sn_stmcm = matrix(0,nrow=nst,ncol=nmonth)
      cumhaz_sp_stmcm = matrix(0,nrow=nst,ncol=nmonth)
      
      if (c>12-M) {
        
        cumhaz_sn_stmcm[,c:12] = lambdaprevyear_mat[,c:12]/12 * hr_agefrailty[,1:(12-c+1)] * hr_maternalprotection_sn[,1:(12-c+1)]
        cumhaz_sn_stmcm[,1:(M+c-12)] = lambda_mat[,1:(M+c-12)]/12 * hr_agefrailty[,(12-c+2):(M+1)] * hr_maternalprotection_sn[,(12-c+2):(M+1)]
        
        cumhaz_sp_stmcm[,c:12] = lambdaprevyear_mat[,c:12]/12 * hr_agefrailty[,1:(12-c+1)] * hr_maternalprotection_sp[,1:(12-c+1)]
        cumhaz_sp_stmcm[,1:(M+c-12)] = lambda_mat[,1:(M+c-12)]/12 * hr_agefrailty[,(12-c+2):(M+1)] * hr_maternalprotection_sp[,(12-c+2):(M+1)]
        
      } else {
        
        cumhaz_sn_stmcm[,1:(M+1)] = lambda_mat[,1:(M+1)]/12 * hr_agefrailty[,1:(M+1)] * hr_maternalprotection_sn[,1:(M+1)]
        cumhaz_sp_stmcm[,1:(M+1)] = lambda_mat[,1:(M+1)]/12 * hr_agefrailty[,1:(M+1)] * hr_maternalprotection_sp[,1:(M+1)]
        
      }
      
      # Finally, add this cumulative hazard weighted by the relative number of children born in that month
      cumhaz_sn[,M+1] = cumhaz_sn[,M+1] + rowSums(cumhaz_sn_stmcm) * n_st_bymonth[,c]/n_byyear
      cumhaz_sp[,M+1] = cumhaz_sp[,M+1] + rowSums(cumhaz_sp_stmcm) * n_st_bymonth[,c]/n_byyear
      
    }
    
  }
  
  # Now if we are including vertical transmission, treat the first month differently
  # Do the 0-2 week cases (i.e. including vertical transmission risk)
  if (include_vtrans==1) {
    
    p_inf_sp = (1 - exp(-lambda_mat/12 * hr_maternalprotection * hr_agefrailty))
    p_inf_sn = (1 - exp(-lambda_mat/12 * hr_agefrailty))
    
    p_inf_sp[,1] = p_inf_sp[,1] + p_sp_vt
    p_inf_sn[,1] = p_inf_sn[,1] + p_sn_vt
    
  } else {
    
    p_inf_sp = (1 - p_sp_vt) * 
      exp(-cumhaz_sp) * (1 - exp(-lambda_mat/12 * hr_maternalprotection_sp * hr_agefrailty))
    p_inf_sn = (1 - p_sn_vt) * 
      exp(-cumhaz_sp) * (1 - exp(-lambda_mat/12 * hr_maternalprotection_sn * hr_agefrailty))
    
  }
  
  p_inf = msp_mat * p_inf_sp + (1 - msp_mat) * p_inf_sn
  
  p_sev = (mspmult_mat*p_sev_sm + msp1_mat*p_sev_s1) * p_inf_sp + (1-msp_mat)*p_sev_sn * p_inf_sn
  
  pred_y_st_month = denom_st_bymonth * p_inf * (p_rep_mild + p_sev * (p_rep_sev - p_rep_mild)) * (1 - prop_missage)
  pred_sev_st_month = denom_st_bymonth * p_sev * p_rep_sev * (1 - prop_missage)
  
  log_lik_case = dpois(y_st_bymonth,pred_y_st_month,log=T)
  dim(log_lik_case)=dim(pred_y_st_month)
  
  log_lik_sev = dpois(sev_st_bymonth,pred_sev_st_month,log=T)
  dim(log_lik_sev)=dim(pred_sev_st_month)
  
  
  return(list(pred_y_st_month,pred_sev_st_month,log_lik_case,log_lik_sev))
}

ageprofile = function(par,
                      nmonth,
                      nst,
                      ns,
                      ny,
                      y_st_bymonth,
                      n_st_bymonth,
                      sev_st_bymonth,
                      y_0d_st,
                      sev_0d_st,
                      lambda_st,
                      msp_st,
                      msp1_st,
                      mspmult_st,
                      mpinfsn_st,
                      mpinfs1_st,
                      mpinfmult_st,
                      model_type,
                      ade_dist,
                      include_vtrans,
                      vtrans_prob,
                      prop_missage,
                      fixedpar) {
  
  expit = function(x) {1/(1+exp(-x))}
  logit = function(x) {log(x/(1-x))}
  
  if (length(intersect(names(par),names(fixedpar)))>0) {stop("Repetition of parameters in par and fixedpar")}
  
  allpars = c(par,fixedpar)
  
  mat_protection_decay = allpars[['mat_protection_decay']]
  baseline_hazmatprot_sn = allpars[['baseline_hazmatprot_sn']]
  baseline_hazmatprot_sp = allpars[['baseline_hazmatprot_sp']]
  baseline_hazinfmult = allpars[['baseline_hazinfmult']]
  baseline_hazinfmult = allpars[['baseline_hazinfmult']]
  ageinfection_decay = allpars[['ageinfection_decay']]
  oneyear_probsev = allpars[['oneyear_probsev']]
  titersev_mult = allpars[['titersev_mult']]
  titersev_meanlog = allpars[['titersev_meanlog']]
  titersev_sdlog = allpars[['titersev_sdlog']]
  
  baseline_replogor = allpars[['baseline_replogor']]
  agereporting_decay = allpars[['agereporting_decay']]
  baseline_sevlogor = allpars[['baseline_sevlogor']]
  ageseverity_decay = allpars[['ageseverity_decay']]
  reportingsev_logor = allpars[['reportingsev_logor']]
  
  titersevmulti_meanlog = titersev_meanlog
  titersevmulti_sdlog=titersev_sdlog
  
  alphast = rep(NA,nst)
  for (i in 1:nst) {
    alphast[i] = allpars[[paste0('alpha',i)]]
  }
  
  lambda_mat = matrix(lambda_st,nrow=nst,ncol=nmonth,byrow=F)
  msp_mat = matrix(msp_st,nrow=nst,ncol=nmonth,byrow=F)
  msp1_mat = matrix(msp1_st,nrow=nst,ncol=nmonth,byrow=F)
  mspmult_mat = matrix(mspmult_st,nrow=nst,ncol=nmonth,byrow=F)
  mpinfsn_mat=matrix(mpinfsn_st,nrow=nst,ncol=nmonth,byrow=F)
  mpinfs1_mat=matrix(mpinfs1_st,nrow=nst,ncol=nmonth,byrow=F)
  mpinfmult_mat=matrix(mpinfmult_st,nrow=nst,ncol=nmonth,byrow=F)
  
  lambdaprevyear_mat = lambda_mat
  lambdaprevyear_mat[2:nst,] = lambda_mat[1:(nst-1),]
  lambdaprevyear_mat[seq(1,nst,by=ny),] = lambda_mat[seq(1,nst,by=ny),]
  
  n_st_bymonth_prevyear = n_st_bymonth
  n_st_bymonth_prevyear[2:nst,] = n_st_bymonth[1:(nst-1),]
  n_st_bymonth_prevyear[seq(1,nst,by=ny),] = n_st_bymonth[seq(1,nst,by=ny),]
  
  cumhaz_sp=matrix(NA,nrow=nst,ncol=nmonth)
  cumhaz_sn=matrix(NA,nrow=nst,ncol=nmonth)
  
  denom_st_bymonth=matrix(NA,nrow=nst,ncol=nmonth)
  
  Ms = 0:(nmonth-1)
  
  hr_maternalprotection_sn = 1 - baseline_hazmatprot_sn * exp(-mat_protection_decay * Ms)
  hr_maternalprotection_sp = 1 - baseline_hazmatprot_sp * exp(-mat_protection_decay * Ms)
  hr_agefrailty = 1 + baseline_hazinfmult * exp(-ageinfection_decay * Ms)
  
  hr_sn = hr_maternalprotection_sn * hr_agefrailty
  hr_sp = hr_maternalprotection_sp * hr_agefrailty
  
  p_rep_mild = exp(baseline_replogor * exp(-agereporting_decay * Ms))
  
  rep_sn = hr_sn * p_rep_mild
  rep_sp = hr_sp * p_rep_mild
  
  p_sev_sn = exp(baseline_sevlogor * exp(-ageseverity_decay*Ms))
  
  # Note that if model_type is 0 or 1, p_sev_s1 = p_sev_sm
  if (ade_dist==0) {
    p_sev_s1 = exp(baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dlnorm(Ms,titersev_meanlog, titersev_sdlog))
    
    p_sev_sm = exp(baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dlnorm(Ms,titersevmulti_meanlog, titersevmulti_sdlog))
  } else {
    p_sev_s1 = exp(baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dnorm(Ms,titersev_meanlog, titersev_sdlog))
    
    p_sev_sm = exp(baseline_sevlogor * exp(-ageseverity_decay*Ms) + 
                              titersev_mult * dnorm(Ms,titersevmulti_meanlog, titersevmulti_sdlog))
  }
  
  prepsev_sn = rep_sn * p_sev_sn
  prepsev_sp = rep_sp * p_sev_sm
  
  return(list(hr_sn,hr_sp,rep_sn,rep_sp,p_sev_sn,p_sev_sm,prepsev_sn,prepsev_sp))
}

expit = function(x) {1/(1+exp(-x))}
logit = function(x) {log(x/(1-x))}


### Generate FOI and maternal seroprevalence across states for 100 years
### Use birth and maternal age data to generate hazards for infants in each state-year
### Use hazards to generate number of infant cases and severe cases by month

state_codes = read.table('C:/Users/mhitchings/Dropbox (UFL)/DengueInfants/Data/01-processedData/mapping/state.txt')
colnames(state_codes) = c("letter","number")

if (include_1day==1) {
  all_den = read.csv('C:/Users/mhitchings/Dropbox (UFL)/DengueInfants/Data/02-analysisData/DengueCasesByStateYear_0517.csv')
} else {
  all_den = read.csv('C:/Users/mhitchings/Dropbox (UFL)/DengueInfants/Data/02-analysisData/DengueCasesByStateYear_0517_No1Day.csv')
}

all_den = all_den %>% arrange(state,year_not)

statenums = unique(all_den$state)

maternal_sp_dat = array(NA,dim=c(2014-2000+1,length(statenums),3),dimnames=list(2000:2014,
                                                                                statenums,
                                                                                c("SP","Smult","S1")))

foi_dat = array(NA,dim=c(2014-2000+1,length(statenums)),dimnames=list(2000:2014,
                                                                      statenums))

maternal_prev_dat = array(0,dim=c(2014-2000+1,length(statenums),3),dimnames=list(2000:2014,
                                                                                 statenums,
                                                                                 c("MP_SN","MP_S1","MP_Smult")))

nst = nrow(all_den)
ns = length(unique(all_den$state))
ny = 2014-2000+1
nmonth=12

include_vtrans = 0
vtrans_prob=0
prop_missage=0

y_st_bymonth=as.matrix(all_den[,paste0('infantcases_',0:11,'month')])
sev_st_bymonth=as.matrix(all_den[,paste0('infantseverecases_',0:11,'month')])
y_2weeks_st=all_den$infantcases_2weeks
sev_2weeks_st=all_den$infantseverecases_2weeks

ade_dist = 1
include_ade = 1
include_sevagefrailty = 1
include_infagefrailty = 1
include_agereporting = 1

namepars = c("mat_protection_decay",
   "baseline_hazmatprot_sn",
   "baseline_hazmatprot_sp",
   "baseline_hazinfmult",
   "ageinfection_decay",
   "baseline_replogor",
   "agereporting_decay",
   "oneyear_probsev",
   "baseline_sevlogor",
   "ageseverity_decay",
   "titersev_mult",
   "titersev_meanlog",
   "titersev_sdlog",
   "reportingsev_logor")

namepars_fortable = c("mat_protection_decay",
            "baseline_hazmatprot_sn",
            "baseline_hazmatprot_sp",
            "baseline_hazinfmult",
            "ageinfection_decay",
            "baseline_replogor",
            "agereporting_decay",
            "oneyear_probsev",
            "baseline_sevlogor",
            "ageseverity_decay",
            "titersev_mult",
            "titersev_mu",
            "titersev_sd",
            "reportingsev_logor")

filesuffix = paste0('_ade_',include_ade,
          '_sevage_',include_sevagefrailty,
          '_infage_',include_infagefrailty,
          '_repage_',include_agereporting,
          '_matprot_',include_maternalprotection,
          '_alphaconst_',alphaconst,
          '_alphalin_',alphalin,
          '_foitype_',foi_type,
          '_lambdauncertainty_',lambda_uncertainty,
          '_adedist_',c('LN','N')[ade_dist+1],
          filesuffix_extra)

p = readRDS(paste0('./stan_fits/NewFOIModels/infantmodel',filesuffix,'.rds'))

post_pars = as.data.frame(as.matrix(p))
post_pars = post_pars[,c(namepars,grep('alphast_param',colnames(post_pars),value=T))]

if (alphaconst==1) {
  post_pars_forpred = cbind(post_pars[,1:14],post_pars[,rep(15:41,each=ny)])
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
} else if (alphalin==1) {
  alphaparams = post_pars[,grepl('alphast_param',colnames(post_pars))]
  allalphas = matrix(NA,nrow=nrow(post_pars),ncol=nst)

  for (t in 1:ny) {
    for (s in 1:ns) {
      allalphas[,(s - 1)*ny + t] = expit(logit(alphaparams[,s]) + (t-1) * alphaparams[,s+ns]);
    }
  }

  post_pars_forpred = cbind(post_pars[,1:14],allalphas)
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
} else {
  post_pars_forpred = post_pars
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
}

meanfitpars = apply(post_pars_forpred,2,mean)
meancifitpars = apply(post_pars_forpred,2,function(x) quantile(x,c(0.5,0.025,0.975)))

### HR at birth for seropositive vs. seronegative infants
rel_hr_atbirth = (1-post_pars_forpred$baseline_hazmatprot_sp)/(1-post_pars_forpred$baseline_hazmatprot_sn)
quantile(rel_hr_atbirth,c(0.5,0.025,0.975))

meancifitpars = meancifitpars[,colnames(meancifitpars) %in% namepars]
meancifitpars[,'baseline_hazmatprot_sn'] = 1 - meancifitpars[,'baseline_hazmatprot_sn']
meancifitpars[,'baseline_hazmatprot_sp'] = 1 - meancifitpars[,'baseline_hazmatprot_sp']
meancifitpars[,'baseline_replogor'] = exp(meancifitpars[,'baseline_replogor'])
meancifitpars[,'baseline_sevlogor'] = exp(meancifitpars[,'baseline_sevlogor'])
meancifitpars[,'oneyear_probsev'] = 100*expit(meancifitpars[,'oneyear_probsev'])
partable = data.frame('Parameter'=namepars_fortable,
                'MedianCrI'=paste0(round(meancifitpars[1,],2),' (',round(meancifitpars[2,],2),',',round(meancifitpars[3,],2),')'))
write.csv(partable,'~/UF/Research/DengueInfants/Manuscripts/BestModParTable_0512.csv',row.names=F)

reprate_forplot = post_pars_forpred[,grepl('alpha',colnames(post_pars_forpred))]
reprate_forplot_summ = as.data.frame(t(apply(reprate_forplot,2,function(x) quantile(x,c(0.5,0.025,0.975)))))
colnames(reprate_forplot_summ) = c('Median','LCI','UCI')
reprate_forplot_summ$State = all_den$statecode
reprate_forplot_summ$Year = all_den$year_not

ggsave(paste0('~/UF/Research/DengueInfants/Manuscripts/Figure_S14_ReportingRate.png'),
       reprate_forplot_summ %>% ggplot(aes(x=Year,y=Median)) + geom_line() +
         geom_ribbon(aes(x=Year,ymin=LCI,ymax=UCI),col='grey',alpha=0.5)+
         facet_wrap(vars(State),nrow=6,ncol=6)+
         scale_x_continuous(breaks=c(2000,2005,2010,2015),labels=c('00','05','10','15'))+ylab('Reporting rate')+scale_y_continuous(limits=c(0,1)),
       height=5,width=5,units='in',device='png')

maternal_sp_dat = array(NA,dim=c(2014-2000+1,length(statenums),3),dimnames=list(2000:2014,
                                                                          statenums,
                                                                          c("SP","Smult","S1")))

foi_dat = array(NA,dim=c(2014-2000+1,length(statenums)),dimnames=list(2000:2014,
                                                                statenums))

maternal_prev_dat = array(0,dim=c(2014-2000+1,length(statenums),3),dimnames=list(2000:2014,
                                                                           statenums,
                                                                           c("MP_SN","MP_S1","MP_Smult")))

if (foi_dir=='./Data/02-foi/stanfit_diff_reduced/11_rta_init_long/state/all/') {

  for (s in statenums) {
  
    if (file.exists(paste0(foi_dir,s,'/S_infantAb.csv'))) {
    
      t = read.csv(paste0(foi_dir,s,'/S_infantAb.csv'))
      
      
      t = t %>% filter(year>=1999) %>% mutate(Smult = S2+S3+S4,
                                              SP = Smult+S1) %>%
        group_by(year) %>% summarise(Smult=mean(Smult),SP=mean(SP),
                                     S1=mean(S1),S1_I = mean(S1_I),
                                     S2=mean(S2),S2_I = mean(S2_I),
                                     S3=mean(S3),S3_I=mean(S3_I),
                                     S4=mean(S4),S4_I=mean(S4_I))
      t = t %>% mutate(MP_SN = S1_I * S1,
                       MP_S1 = S2_I * S2,
                       MP_Smult = S3_I * S3 + S4_I * S4)
      
      maternal_prev_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("MP_SN")] = (head(t$MP_SN,-1)+tail(t$MP_SN,-1))/2
      maternal_prev_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("MP_S1")] = (head(t$MP_S1,-1)+tail(t$MP_S1,-1))/2
      maternal_prev_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("MP_Smult")] = (head(t$MP_Smult,-1)+tail(t$MP_Smult,-1))/2
      
      maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("SP")] = (head(t$SP,-1)+tail(t$SP,-1))/2
      maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("Smult")] = (head(t$Smult,-1)+tail(t$Smult,-1))/2
      maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("S1")] = (head(t$S1,-1)+tail(t$S1,-1))/2
    
    }
    
    if (file.exists(paste0(foi_dir,s,'/hazard_infant.csv'))) {
    
      t = read.csv(paste0(foi_dir,s,'/hazard_infant.csv'))
      t = t %>% filter(year>=2000) %>%
        group_by(year) %>% summarise(FOI = mean(hazard))
      
      foi_dat[,dimnames(foi_dat)[[2]]==s] = t$FOI
    
    }

  }

} else if (grepl('stanfit_v3',foi_dir)) {

  for (s in statenums) {
  
    if (file.exists(paste0(foi_dir,'S_mother/',s,'.csv'))) {
    
    t = read.csv(paste0(foi_dir,'S_mother/',s,'.csv'))
    
    t = t %>% filter(year>=2000 & year<=2014) %>% mutate(SP = Smult+S1) %>%
      group_by(year) %>% summarise(Smult=mean(Smult),SP=mean(SP),S1=mean(S1))
    
    if (s == 43 & foi_dir == "./Data/02-foi/stanfit_v3/31_rta_async_sevShift_long/") {
      t = t %>% bind_rows(data.frame('year'=2013:2014,
                                     'Smult'=rep(t$Smult[13],2),
                                     'SP'=rep(t$SP[13],2),
                                     'S1'=rep(t$S1[13],2)
      )
      )
    }
    
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("SP")] = t$SP
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("Smult")] = t$Smult
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("S1")] = t$S1
    
    }
    
  if (file.exists(paste0(foi_dir,'lambda_t/',s,'.csv'))) {
  
    t = read.csv(paste0(foi_dir,'lambda_t/',s,'.csv'))
    t = t %>% filter(year>=2000 & year<=2014) %>%
      group_by(year) %>% summarise(FOI = mean(val))
    
    foi_dat[,dimnames(foi_dat)[[2]]==s] = t$FOI
    
    }
  
  }

} else {

  for (s in statenums) {
  
    if (file.exists(paste0(foi_dir,'S_infantAb/',s,'.csv'))) {
    
    t = read.csv(paste0(foi_dir,'S_infantAb/',s,'.csv'))
    
    t = t %>% filter(year>=2000 & year<=2014) %>% mutate(SP = Smult+S1) %>%
      group_by(year) %>% summarise(Smult=mean(Smult),SP=mean(SP),S1=mean(S1))
    
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("SP")] = t$SP
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("Smult")] = t$Smult
    maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("S1")] = t$S1
    
    }
    
    if (file.exists(paste0(foi_dir,'hazard_infant/',s,'.csv'))) {
    
    t = read.csv(paste0(foi_dir,'hazard_infant/',s,'.csv'))
    t = t %>% filter(year>=2000 & year<=2014) %>%
      group_by(year) %>% summarise(FOI = mean(hazard))
    
    foi_dat[,dimnames(foi_dat)[[2]]==s] = t$FOI
    
    }
    
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

all_den$maternal_prev_sn = sapply(1:nrow(all_den),function(x) (maternal_prev_dat[dimnames(maternal_prev_dat)[[1]]==all_den$year_not[x],
                                                                           dimnames(maternal_prev_dat)[[2]]==all_den$state[x],
                                                                           c("MP_SN")]))
all_den$maternal_prev_s1 = sapply(1:nrow(all_den),function(x) (maternal_prev_dat[dimnames(maternal_prev_dat)[[1]]==all_den$year_not[x],
                                                                           dimnames(maternal_prev_dat)[[2]]==all_den$state[x],
                                                                           c("MP_S1")]))
all_den$maternal_prev_smult = sapply(1:nrow(all_den),function(x) (maternal_prev_dat[dimnames(maternal_prev_dat)[[1]]==all_den$year_not[x],
                                                                              dimnames(maternal_prev_dat)[[2]]==all_den$state[x],
                                                                              c("MP_Smult")]))

lambda_st=all_den$foi
msp1_st=all_den$maternal_s1
mspmult_st=all_den$maternal_smult
msp_st = msp1_st+mspmult_st
y_0d_st=all_den$ncase_noage
sev_0d_st=all_den$nsevcase_noage
n_st_bymonth = as.matrix(all_den[,paste0('infantbirths_',1:12,'month')])
mpinfsn_st=all_den$maternal_prev_sn
mpinfs1_st=all_den$maternal_prev_s1
mpinfmult_st=all_den$maternal_prev_smult

pred = apply(post_pars_forpred,1,function(x)
pred_cases(x,
       nmonth,
       nst,
       ns,
       ny,
       y_st_bymonth,
       n_st_bymonth,
       sev_st_bymonth,
       y_0d_st,
       sev_0d_st,
       lambda_st,
       msp_st,
       msp1_st,
       mspmult_st,
       mpinfsn_st,
       mpinfs1_st,
       mpinfmult_st,
       model_type,
       ade_dist,
       include_vtrans,
       vtrans_prob,
       prop_missage,
       c())
)

mean_matrix_predictedcases = Reduce("+", lapply(pred,function(x) x[[1]]))/ length(pred)
mean_matrix_predictedsevcases = Reduce("+", lapply(pred,function(x) x[[2]]))/ length(pred)

cc = c('lightslategrey','seagreen','darkorchid2','lightsalmon','navy')

### Predicted and observed cases by state, year, and month of age
pred_vs_obs_cases = all_den[,c("year_not","statecode","region",paste0("infantcases_",0:11,"month"))] %>% 
  mutate(region = factor(region,levels=c("North","Northeast","Central-West","Southeast","South")))
colnames(pred_vs_obs_cases) = c("year_not","statecode","region",0:11)

pred_vs_obs_cases = pred_vs_obs_cases %>% pivot_longer(cols = '0':'11',names_to = 'Month', values_to = 'Observed') %>%
  mutate(Fit = c(t(mean_matrix_predictedcases)),
         Month = as.numeric(Month),
         y5 = floor(year_not/5),
         y5_label = case_when(y5 == 400 ~ "2000-2004",
                              y5 == 401 ~ "2005-2009",
                              y5 == 402 ~ "2010-2014",
                              y5 == 403 ~ "2015-2019",
                              y5 == 404 ~ "2020-2024")
  )

pred_vs_obs_cases_byregion_5year = pred_vs_obs_cases %>% 
  group_by(region,y5_label,Month) %>% summarise(Observed = sum(Observed),
                                                Fit = sum(Fit))

### Predicted and observed severe cases by state, year, and month of age
pred_vs_obs_sevcases = all_den[,c("year_not","statecode","region",paste0("infantseverecases_",0:11,"month"))] %>% 
  mutate(region = factor(region,levels=c("North","Northeast","Central-West","Southeast","South")))
colnames(pred_vs_obs_sevcases) = c("year_not","statecode","region",0:11)

pred_vs_obs_sevcases = pred_vs_obs_sevcases %>% pivot_longer(cols = '0':'11',names_to = 'Month', values_to = 'Observed') %>%
  mutate(Fit = c(t(mean_matrix_predictedsevcases)),
         Month = as.numeric(Month),
         y5 = floor(year_not/5),
         y5_label = case_when(y5 == 400 ~ "2000-2004",
                              y5 == 401 ~ "2005-2009",
                              y5 == 402 ~ "2010-2014",
                              y5 == 403 ~ "2015-2019",
                              y5 == 404 ~ "2020-2024"))

pred_vs_obs_sevcases_byregion_5year = pred_vs_obs_sevcases %>% 
  group_by(region,y5_label,Month) %>% summarise(Observed = sum(Observed),
                                                Fit = sum(Fit))

pred_vs_obs_cases_brazil_5year = pred_vs_obs_cases %>% 
  group_by(y5_label,Month) %>% summarise(Observed=sum(Observed),
                                         Fit=sum(Fit))
pred_vs_obs_sevcases_brazil_5year = pred_vs_obs_sevcases %>% 
  group_by(y5_label,Month) %>% summarise(Observed=sum(Observed),
                                         Fit=sum(Fit))

### Predicted and observed cases across Brazil over time
pred_vs_obs_cases_byiter_vec = unlist(lapply(pred,function(x) c(t(x[[1]]))))
pred_vs_obs_cases_byiter = data.frame('iter'=rep(1:length(pred),each=405*12),
                                      'statecode'=rep(state_codes$letter,each=15*12),
                                      'year_not'=rep(2000:2014,each=12),
                                      'Month'=0:11,
                                      'Pred'=pred_vs_obs_cases_byiter_vec)
pred_vs_obs_cases_byiter = pred_vs_obs_cases_byiter %>% left_join(pred_vs_obs_cases %>% select(year_not,statecode,region,Month,Observed,y5,y5_label),by=c('year_not','statecode','Month'))

pred_vs_obs_sevcases_byiter_vec = unlist(lapply(pred,function(x) c(t(x[[2]]))))
pred_vs_obs_sevcases_byiter = data.frame('iter'=rep(1:length(pred),each=405*12),
                                      'statecode'=rep(state_codes$letter,each=15*12),
                                      'year_not'=rep(2000:2014,each=12),
                                      'Month'=0:11,
                                      'Pred'=pred_vs_obs_sevcases_byiter_vec)
pred_vs_obs_sevcases_byiter = pred_vs_obs_sevcases_byiter %>% left_join(pred_vs_obs_sevcases %>% select(year_not,statecode,region,Month,Observed,y5,y5_label),by=c('year_not','statecode','Month'))

pred_vs_obs_cases_brazil_5year = pred_vs_obs_cases_byiter %>% 
  group_by(y5_label,Month,iter) %>% summarise(Observed=sum(Observed),
                                         Fit=sum(Pred))

ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS_CaseAgeDistribution_5Year_Brazil.png',
  ggplot() + 
  geom_ribbon(data=pred_vs_obs_cases_brazil_5year %>% 
                group_by(y5_label,Month) %>% 
                summarise(Observed=mean(Observed),
                          MeanFit=mean(Fit),
                          FitUCI=quantile(Fit,0.975),
                          FitLCI=quantile(Fit,0.025)),
              aes(x=Month,ymin=FitLCI,ymax=FitUCI),alpha=0.2,colour='grey') + 
  geom_line(data=pred_vs_obs_cases_brazil_5year %>% 
  group_by(y5_label,Month) %>% 
  summarise(Observed=mean(Observed),
            MeanFit=mean(Fit),
            FitUCI=quantile(Fit,0.975),
            FitLCI=quantile(Fit,0.025)) %>% 
  pivot_longer(cols = 'Observed':'MeanFit'),aes(x=Month,y=value,linetype=name)) +
  facet_wrap(vars(y5_label)) +
  scale_x_continuous(name='Age',breaks=c(0,6,12))+
  ylab('Count')+
  scale_linetype_manual(name='',values=c(1,2),labels=c('Fit','Observed')),
  height=3,width=7,units='in',device='png')

pred_vs_obs_sevcases_brazil_5year = pred_vs_obs_sevcases_byiter %>% 
  group_by(y5_label,Month,iter) %>% summarise(Observed=sum(Observed),
                                              Fit=sum(Pred))

ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS_SevCaseAgeDistribution_5Year_Brazil.png',
       ggplot() + 
         geom_ribbon(data=pred_vs_obs_sevcases_brazil_5year %>% 
                       group_by(y5_label,Month) %>% 
                       summarise(Observed=mean(Observed),
                                 MeanFit=mean(Fit),
                                 FitUCI=quantile(Fit,0.975),
                                 FitLCI=quantile(Fit,0.025)),
                     aes(x=Month,ymin=FitLCI,ymax=FitUCI),alpha=0.2,colour='grey') + 
         geom_line(data=pred_vs_obs_sevcases_brazil_5year %>% 
                     group_by(y5_label,Month) %>% 
                     summarise(Observed=mean(Observed),
                               MeanFit=mean(Fit),
                               FitUCI=quantile(Fit,0.975),
                               FitLCI=quantile(Fit,0.025)) %>% 
                     pivot_longer(cols = 'Observed':'MeanFit'),aes(x=Month,y=value,linetype=name)) +
         facet_wrap(vars(y5_label)) +
         scale_x_continuous(name='Age',breaks=c(0,6,12))+
         ylab('Count')+
         scale_linetype_manual(name='',values=c(1,2),labels=c('Fit','Observed')),
       height=3,width=7,units='in',device='png')


pred_vs_obs_cases_brazil_5year %>% 
  pivot_longer(cols = 'Observed':'Fit') %>% ggplot() + #geom_point(aes(x=Month,y=value,shape=name)) + 
  geom_line(aes(x=Month,y=value,linetype=name)) + 
  facet_wrap(vars(y5_label))

pred_vs_obs_sevcases_brazil_5year %>% 
  pivot_longer(cols = 'Observed':'Fit') %>% ggplot() + #geom_point(aes(x=Month,y=value,shape=name)) + 
  geom_line(aes(x=Month,y=value,linetype=name)) + 
  facet_wrap(vars(y5_label))


# Predicted and observed proportion of cases among neonates
predprop_byregion_5year = pred_vs_obs_cases_byregion_5year %>% 
  mutate(PredCase=Fit*as.numeric(Month==0),
         ObsCase=Observed*as.numeric(Month==0)) %>% 
  group_by(region,y5_label) %>% 
  summarise(PropObs=sum(ObsCase/sum(Observed)),
            PropPred=sum(PredCase)/sum(Fit))

predsevprop_byregion_5year = pred_vs_obs_sevcases_byregion_5year %>% 
  mutate(PredCase=Fit*as.numeric(Month>=5),
         ObsCase=Observed*as.numeric(Month>=5)) %>% 
  group_by(region,y5_label) %>% 
  summarise(PropObs=sum(ObsCase/sum(Observed)),
            PropPred=sum(PredCase)/sum(Fit))

allpredprop_byregion_5year = 
  bind_rows(predprop_byregion_5year %>% mutate(Case = "All cases"),
            predsevprop_byregion_5year %>% mutate(Case = "Severe cases"))


ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS_11_FitProportions.png',
       allpredprop_byregion_5year %>% 
         pivot_longer(cols = 'PropObs':'PropPred') %>% 
         mutate(y5 = case_when(y5_label == '2000-2004' ~ 1,
                               y5_label == '2005-2009' ~ 2,
                               y5_label == '2010-2014' ~ 3
         ),
         fitlabel = case_when(name=="PropObs" ~ "Observed",
                              name=="PropPred" ~ "Fit")) %>%
         ggplot() + 
         geom_point(aes(x=y5,y=value,colour=region)) + 
         geom_line(aes(x=y5,y=value,colour=region)) + 
         scale_colour_manual(name='Region',values=cc) + 
         facet_grid(Case ~ fitlabel) + 
         ylab('Proportion') + 
         scale_x_continuous(name='Year',breaks=seq(1:3),labels=c('00-04','05-09','10-14'))+
         theme(legend.position='right',
               axis.text = element_text(size=6)),
       height=5,width=6,units='in',device='png')


predcase_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[1]]))),nrow=length(pred),ncol=12,byrow=T)
predsev_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[2]]))),nrow=length(pred),ncol=12,byrow=T)

data_forplot = data.frame('Month'=rep(1:12,times=6),
                    'Severe'=c(rep('Cases',36),rep('Severe',36)),
                    'Type'=rep(c(rep('Observed',12),rep('Fit',24)),2),
                    'Mean'=c(colSums(y_st_bymonth),apply(predcase_age_dist,2,mean),rep(NA,12),
                             colSums(sev_st_bymonth),apply(predsev_age_dist,2,mean),rep(NA,12)),
                    'LCI'=c(rep(NA,24),apply(predcase_age_dist,2,function(x) quantile(x,0.025)),
                            rep(NA,24),apply(predsev_age_dist,2,function(x) quantile(x,0.025))),
                    'UCI'=c(rep(NA,24),apply(predcase_age_dist,2,function(x) quantile(x,0.975)),
                            rep(NA,24),apply(predsev_age_dist,2,function(x) quantile(x,0.975)))
)

# Plot by state, sum across years
state_indices = matrix(1:405,nrow=ns,ncol=ny,byrow=T)

predcase_age_dist_bystate = lapply(pred,function(p) t(sapply(1:ns,function(x) colSums(p[[1]][state_indices[x,],]))))

predy_bystate_total = data.frame('Month' = rep(1:12,ns))
predy_bystate_total$State = paste0(rep(unique(all_den$state),each=12),"-",
                             rep(unique(all_den$statecode),each=12))
trow = 1
for (s in 1:ns) {

  for (m in 1:12) {

    predy_bystate_total[trow,c("Median","LCI","UCI")] = quantile(unlist(lapply(predcase_age_dist_bystate,function(x) x[s,m])),c(0.5,0.025,0.975))

    trow = trow + 1

  }

}

predy_bystate_total$Group = 'Fit'

obsy_bystate_total = data.frame('Median'=c(t(sapply(0:11,function(x)
sapply(unique(all_den$state),
   function(s) sum(all_den[all_den$state==s,paste0("infantcases_",x,"month")]))))),
'LCI'=NA,
'UCI'=NA,
'Month'=rep(1:12,ns),
'State'=paste0(rep(unique(all_den$state),each=12),"-",
           rep(unique(all_den$statecode),each=12)),
'Group'='Observed')
predy_bystate_total = bind_rows(predy_bystate_total,obsy_bystate_total)

## Severe cases
predsev_age_dist_bystate = lapply(pred,function(p) t(sapply(1:ns,function(x) colSums(p[[2]][state_indices[x,],]))))

predsev_bystate_total = data.frame('Month' = rep(1:12,ns))
predsev_bystate_total$State = paste0(rep(unique(all_den$state),each=12),"-",
                               rep(unique(all_den$statecode),each=12))
trow = 1
for (s in 1:ns) {

  for (m in 1:12) {

    predsev_bystate_total[trow,c("Median","LCI","UCI")] = quantile(unlist(lapply(predsev_age_dist_bystate,function(x) x[s,m])),c(0.5,0.025,0.975))

    trow = trow + 1

  }

}

predsev_bystate_total$Group = 'Fit'

obssev_bystate_total = data.frame('Median'=c(t(sapply(0:11,function(x)
  sapply(unique(all_den$state),
         function(s) sum(all_den[all_den$state==s,paste0("infantseverecases_",x,"month")]))))),
  'LCI'=NA,
  'UCI'=NA,
  'Month'=rep(1:12,ns),
  'State'=paste0(rep(unique(all_den$state),each=12),"-",
             rep(unique(all_den$statecode),each=12)),
  'Group'='Observed')
predsev_bystate_total = bind_rows(predsev_bystate_total,obssev_bystate_total)

profiles = apply(post_pars_forpred,1,function(x)
  ageprofile(x,
             nmonth,
             nst,
             ns,
             ny,
             y_st_bymonth,
             n_st_bymonth,
             sev_st_bymonth,
             y_0d_st,
             sev_0d_st,
             lambda_st,
             msp_st,
             msp1_st,
             mspmult_st,
             mpinfsn_st,
             mpinfs1_st,
             mpinfmult_st,
             model_type,
             ade_dist,
             include_vtrans,
             vtrans_prob,
             prop_missage,
             c())
)

hr_sn = matrix(unlist(lapply(profiles,function(x) x[[1]])),nrow=length(profiles),ncol=12,byrow=T)
hr_sp = matrix(unlist(lapply(profiles,function(x) x[[2]])),nrow=length(profiles),ncol=12,byrow=T)
rep_sn = matrix(unlist(lapply(profiles,function(x) x[[3]])),nrow=length(profiles),ncol=12,byrow=T)
rep_sp = matrix(unlist(lapply(profiles,function(x) x[[4]])),nrow=length(profiles),ncol=12,byrow=T)
p_sev_sn = matrix(unlist(lapply(profiles,function(x) x[[5]])),nrow=length(profiles),ncol=12,byrow=T)
p_sev_sm = matrix(unlist(lapply(profiles,function(x) x[[6]])),nrow=length(profiles),ncol=12,byrow=T)
prepsev_sn = matrix(unlist(lapply(profiles,function(x) x[[7]])),nrow=length(profiles),ncol=12,byrow=T)
prepsev_sm = matrix(unlist(lapply(profiles,function(x) x[[8]])),nrow=length(profiles),ncol=12,byrow=T)

data_forplot = data.frame('Month'=rep(1:12,times=8),
                          'EventType'=c(rep('Infection (HR)',24),
                                        rep('Reported Case (OR)',24),
                                        rep('Severe Disease (OR)',24),
                                        rep('Reported Sev. Case',24)),
                          'Serostatus'=rep(c(rep('Seronegative',12),
                                             rep('Seropositive',12)),4),
                          'Mean'=c(apply(hr_sn,2,mean),
                                   apply(hr_sp,2,mean),
                                   apply(rep_sn,2,mean),
                                   apply(rep_sp,2,mean),
                                   apply(p_sev_sn,2,mean),
                                   apply(p_sev_sm,2,mean),
                                   apply(prepsev_sn,2,mean),
                                   apply(prepsev_sm,2,mean)),
                          'LCI'=c(apply(hr_sn,2,function(x) quantile(x,0.025)),
                                  apply(hr_sp,2,function(x) quantile(x,0.025)),
                                  apply(rep_sn,2,function(x) quantile(x,0.025)),
                                  apply(rep_sp,2,function(x) quantile(x,0.025)),
                                  apply(p_sev_sn,2,function(x) quantile(x,0.025)),
                                  apply(p_sev_sm,2,function(x) quantile(x,0.025)),
                                  apply(prepsev_sn,2,function(x) quantile(x,0.025)),
                                  apply(prepsev_sm,2,function(x) quantile(x,0.025))),
                          'UCI'=c(apply(hr_sn,2,function(x) quantile(x,0.975)),
                                  apply(hr_sp,2,function(x) quantile(x,0.975)),
                                  apply(rep_sn,2,function(x) quantile(x,0.975)),
                                  apply(rep_sp,2,function(x) quantile(x,0.975)),
                                  apply(p_sev_sn,2,function(x) quantile(x,0.975)),
                                  apply(p_sev_sm,2,function(x) quantile(x,0.975)),
                                  apply(prepsev_sn,2,function(x) quantile(x,0.975)),
                                  apply(prepsev_sm,2,function(x) quantile(x,0.975)))
)
data_forplot$EventType = factor(data_forplot$EventType,levels=c('Infection (HR)',
                                                                'Severe Disease (OR)',
                                                                'Reported Case (OR)',
                                                                'Reported Sev. Case'))

textsize=8
### Profile plots for model diagram
prof_cc = scales::seq_gradient_pal("#C03830", "#317EC2", "Lab")(seq(0,1,length.out=3))

sp = c(0.5)
data_forplot_comb = data.frame('Month'=1:12,
                               'EventType'="Reported Case (OR)",
                               'Seroprevalence'=rep(sp,each=12),
                               'Serostatus'=rep(c('Medium'),each=12)
)
data_forplot_comb$Mean = data_forplot_comb$Seroprevalence * 
  data_forplot$Mean[data_forplot$EventType=="Reported Case (OR)" & data_forplot$Serostatus=="Seropositive"] +
  (1-data_forplot_comb$Seroprevalence) * 
  data_forplot$Mean[data_forplot$EventType=="Reported Case (OR)" & data_forplot$Serostatus=="Seronegative"]

data_forplot_comb = bind_rows(data_forplot,data_forplot_comb)

profileplot_hr = ggplot() +
  geom_line(data=data_forplot_comb %>% filter(EventType=="Reported Case (OR)"),
            aes(x=as.numeric(Month),y=Mean,
                colour=factor(Serostatus,levels=c('Seropositive','Medium','Seronegative')),
                alpha=factor(as.numeric(Serostatus %in% c('Medium')))),show.legend=F) +
  scale_x_continuous(name='Age',breaks=c(0,6,12))+
  scale_colour_manual(name=NULL,labels=c('Seropositive','Medium','Seronegative'),
                      values=prof_cc) + 
  scale_alpha_manual(name=NULL,values=c(0.5,1)) + 
  ylab('Case hazard\n(relative to 12-month-old)')+
  expand_limits(y=c(0))+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize))
ggsave('~/UF/Research/DengueInfants/Manuscripts/ProfileHR_0520.png',
       profileplot_hr,
       device='png',width=2,height=2,units='in')

profileplot_repsev = ggplot() +
  geom_line(data=data_forplot %>% filter(EventType=="Severe Disease (OR)"),
            aes(x=as.numeric(Month),y=Mean,colour=Serostatus),alpha=0.5,show.legend=F) +
  scale_x_continuous(name='Age',breaks=c(0,6,12))+
  scale_colour_manual(name=NULL,values=c('#317EC2','#C03830')) + 
  ylab('Severe disease odds\n(relative to 12-month-old)')+
  expand_limits(y=c(0,3))+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize))
ggsave('~/UF/Research/DengueInfants/Manuscripts/ProfileSevOR.png',
       profileplot_repsev,
       device='png',width=2,height=2,units='in')

sp = c(0.5)
data_forplot_comb = data.frame('Month'=1:12,
                               'EventType'="Reported Sev. Case",
                               'Seroprevalence'=rep(sp,each=12),
                               'Serostatus'=rep(c('Medium'),each=12)
)
data_forplot_comb$Mean = data_forplot_comb$Seroprevalence * 
  data_forplot$Mean[data_forplot$EventType=="Reported Sev. Case" & data_forplot$Serostatus=="Seropositive"] +
  (1-data_forplot_comb$Seroprevalence) * 
  data_forplot$Mean[data_forplot$EventType=="Reported Sev. Case" & data_forplot$Serostatus=="Seronegative"]

data_forplot_comb = bind_rows(data_forplot,data_forplot_comb)

profileplot_sevcase = ggplot() +
  geom_line(data=data_forplot_comb %>% filter(EventType=="Reported Sev. Case"),
            aes(x=as.numeric(Month),y=Mean,
                colour=factor(Serostatus,levels=c('Seropositive','Medium','Seronegative')),
                alpha=factor(as.numeric(Serostatus %in% c('Medium')))),show.legend=F) +
  scale_x_continuous(name='Age',breaks=c(0,6,12))+
  scale_colour_manual(name=NULL,labels=c('Seropositive','Medium','Seronegative'),
                      values=prof_cc) + 
  scale_alpha_manual(name=NULL,values=c(0.5,1)) + 
  ylab('Severe case hazard\n(relative to 12-month-old)')+
  expand_limits(y=c(0))+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize))
ggsave('~/UF/Research/DengueInfants/Manuscripts/ProfileSevCase_0520.png',
       profileplot_sevcase,
       device='png',width=2,height=2,units='in')

sp = c(0.5)
data_forplot_comb = data.frame('Month'=1:12,
                               'EventType'="Reported Sev. Case",
                               'Seroprevalence'=rep(sp,each=12),
                               'Serostatus'=rep(c('Seroprevalence 50%'),each=12)
)
data_forplot_comb$Mean = data_forplot_comb$Seroprevalence * 
  data_forplot$Mean[data_forplot$EventType=="Reported Sev. Case" & data_forplot$Serostatus=="Seropositive"] +
  (1-data_forplot_comb$Seroprevalence) * 
  data_forplot$Mean[data_forplot$EventType=="Reported Sev. Case" & data_forplot$Serostatus=="Seronegative"]

data_forplot_comb = bind_rows(data_forplot,data_forplot_comb)

profileplot_sevcase = ggplot() +
  geom_line(data=data_forplot_comb %>% filter(EventType=="Reported Sev. Case"),
            aes(x=as.numeric(Month),y=Mean,
                colour=factor(Serostatus,levels=c('Seropositive','Seroprevalence 50%','Seronegative')))) +
  scale_x_continuous(name='Age',breaks=c(0,6,12))+
  scale_colour_manual(name=NULL,labels=c('Seropositive','Seroprevalence 50%','Seronegative'),
                      values=prof_cc) + 
  ylab('Severe case hazard\n(relative to 12-month-old)')+
  expand_limits(y=c(0))+
  theme(legend.position='bottom',
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize))
ggsave('~/UF/Research/DengueInfants/Manuscripts/ProfileSevCase_WithLegend_0520.png',
       profileplot_sevcase,
       device='png',width=6,height=2,units='in')


profileplot_repsev = ggplot() +
  geom_line(data=data_forplot %>% filter(EventType=="Severe Disease (OR)"),
            aes(x=as.numeric(Month),y=Mean,colour=Serostatus),linewidth=1,show.legend=F) +
  scale_x_continuous(name='Age',breaks=c(0,6,12))+
  scale_colour_manual(name=NULL,values=c('#317EC2','#C03830')) + 
  ylab('Severe disease odds\n(relative to 12-month-old)')+
  expand_limits(y=c(0,3))+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize))


profileplot = ggplot() +
  geom_ribbon(data=data_forplot,aes(x=Month,ymin=LCI,ymax=UCI,linetype=Serostatus),alpha=0.2,fill='grey',col=NA,show.legend = F)+
  geom_line(data=data_forplot,aes(x=as.numeric(Month),y=Mean,linetype=Serostatus)) +
  scale_x_continuous(name='Month',breaks=c(0,3,6,9,12))+
  scale_linetype_manual(name=NULL,values=c(1,2)) + 
  ylab('Ratio')+
  expand_limits(y=c(0))+
  facet_wrap(vars(EventType),ncol=2,scales = 'free')+
  guides(
    linetype = guide_legend(position='inside')
  )+
  theme(legend.position.inside = c(0.3,0.9),
        legend.background = element_blank(),
        strip.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize))
ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS10_Profiles.png',
       profileplot,
       device='png',width=6,height=6,units='in')

predcase_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[1]]))),nrow=length(pred),ncol=12,byrow=T)
predsev_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[2]]))),nrow=length(pred),ncol=12,byrow=T)

data_forplot = data.frame('Month'=rep(1:12,times=6),
                          'Severe'=c(rep('Cases',36),rep('Severe',36)),
                          'Type'=rep(c(rep('Observed',12),rep('Fit',24)),2),
                          'Mean'=c(colSums(y_st_bymonth),apply(predcase_age_dist,2,mean),rep(NA,12),
                                   colSums(sev_st_bymonth),apply(predsev_age_dist,2,mean),rep(NA,12)),
                          'LCI'=c(rep(NA,24),apply(predcase_age_dist,2,function(x) quantile(x,0.025)),
                                  rep(NA,24),apply(predsev_age_dist,2,function(x) quantile(x,0.025))),
                          'UCI'=c(rep(NA,24),apply(predcase_age_dist,2,function(x) quantile(x,0.975)),
                                  rep(NA,24),apply(predsev_age_dist,2,function(x) quantile(x,0.975)))
)

predplot = ggplot() +
  geom_ribbon(data=data_forplot %>% filter(is.na(Mean)),aes(x=Month,ymin=LCI,ymax=UCI),alpha=0.2,fill='grey',col=NA)+
  geom_line(data=data_forplot %>% filter(!is.na(Mean)),aes(x=as.numeric(Month),y=Mean,linetype=Type)) +
  scale_x_continuous(name='Month',breaks=c(0,3,6,9,12))+
  ylab('Cases')+
  scale_linetype_manual(name=NULL,values=c(1,2)) + 
  facet_wrap(vars(Severe),ncol=2,scales = 'free') + 
  guides(
    linetype = guide_legend(position='inside')
  )+
  theme(legend.position.inside = c(0.3,0.8),
        legend.background = element_blank(),
        strip.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize))
ggsave('~/UF/Research/DengueInfants/Manuscripts/Figure2.png',
       predplot,device='png',width=4,height=2,units='in')

