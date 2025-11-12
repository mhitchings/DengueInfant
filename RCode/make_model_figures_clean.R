require(dplyr)
require(tidyverse)
require(loo)
require(grid)
require(gridExtra)
require(arrow)
require(rstan)
require(ggplot2); theme_set(theme_classic())

include_maternalprotection=1
alphaconst=0
alphalin=0
foi_dir = './Data/02-foi/stanfit_v3/t51_rta_async_noType_neglectZika_final/'
foi_type='t51_rta_final'
filesuffix_extra=''
lambda_uncertainty = 0

# Function to describe rise and fall in severe disease odds in children of seropositive mothers
ade_func = function(y,mode,kappa,mult,a,c) {
  
  res = rep(0,length(y))
  
  res[y>=a & y<=c] = mult * ((y[y>=a & y<=c] - a) / (mode * (c - a))) ^ (mode * kappa)  * ((c - y[y>=a & y<=c]) / ((1 - mode) * (c - a))) ^ ((1 - mode) * kappa)
  
  return(res)
}

pred_cases = function(par,
                      nmonth,
                      nst,
                      ns,
                      ny,
                      y_st_bymonth,
                      n_st_bymonth,
                      sev_st_bymonth,
                      lambda_st,
                      msp_st,
                      msp1_st,
                      mspmult_st,
                      fixedpar) {
  
  expit = function(x) {1/(1+exp(-x))}
  logit = function(x) {log(x/(1-x))}
  
  if (length(intersect(names(par),names(fixedpar)))>0) {stop("Repetition of parameters in par and fixedpar")}
  
  allpars = c(par,fixedpar)
  
  mat_protection_decay = allpars[['mat_protection_decay']]
  baseline_hazmatprot_sn = allpars[['baseline_hazmatprot_sn']]
  baseline_hazmatprot_sp = allpars[['baseline_hazmatprot_sp']]
  baseline_hazmult_sn = allpars[['baseline_hazmult_sn']]
  baseline_hazmult_sp = allpars[['baseline_hazmult_sp']]
  oneyear_probsev = allpars[['oneyear_probsev']]
  adepeak_logor = allpars[['adepeak_logor']]
  ademode = allpars[['ademode']]
  adekappa = allpars[['adekappa']]
  
  baseline_replogor = allpars[['baseline_replogor']]
  baseline_sevlogor = allpars[['baseline_sevlogor']]
  
  alphast = rep(NA,nst)
  for (i in 1:nst) {
    alphast[i] = allpars[[paste0('alpha',i)]]
  }
  
  lambda_mat = matrix(lambda_st,nrow=nst,ncol=nmonth,byrow=F)
  msp_mat = matrix(msp_st,nrow=nst,ncol=nmonth,byrow=F)
  msp1_mat = matrix(msp1_st,nrow=nst,ncol=nmonth,byrow=F)
  mspmult_mat = matrix(mspmult_st,nrow=nst,ncol=nmonth,byrow=F)
  
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
  
  # Hazard ratio of infection in for infants born to seronegative and seropositive mothers
  hr_sp = matrix(c(1 + baseline_hazmult_sp,1 - baseline_hazmatprot_sp * exp(-mat_protection_decay * 0:(nmonth-2))),nrow=nst,ncol=nmonth,byrow=T)
  hr_sn = matrix(c(1 + baseline_hazmult_sn,1 - baseline_hazmatprot_sn * exp(-mat_protection_decay * 0:(nmonth-2))),nrow=nst,ncol=nmonth,byrow=T)
  
  # Odds of severe disease upon infection in infants born to seronegative and seropositive mothers
  p_sev_sn = matrix(expit(oneyear_probsev + baseline_sevlogor * as.numeric(Ms==0)),nrow=nst,ncol=nmonth,byrow=T)
  
  p_sev_sp = matrix(expit(oneyear_probsev + baseline_sevlogor * as.numeric(Ms==0) + 
                            ade_func(Ms,ademode,adekappa,adepeak_logor,1,nmonth-1)),
                    nrow=nst,ncol=nmonth,byrow=T)
  
  # Reporting odds
  p_rep = expit(logit(matrix(alphast,nrow=nst,ncol=nmonth,byrow=F)) + 
                  matrix(baseline_replogor * as.numeric(Ms==0),nrow=nst,ncol=nmonth,byrow=T))
  
  # We will need the total children born in each calendar year
  n_byyear = apply(n_st_bymonth,1,sum)
  
  cumhaz_sn[,1] = 0
  cumhaz_sp[,1] = 0
  
  denom_st_bymonth[,1] = rowSums(n_st_bymonth)
  
  for (M in 1:(nmonth-1)) {
    
    # Build up cumulative hazard by state-year and month of age
    cumhaz_sp[,M+1] = 0
    cumhaz_sn[,M+1] = 0
    
    # Need to calculate the denominator (who is at risk of having case at age M in year t)
    # For this, we add up all children who could have been M months old during year t
    # i.e. children who were born in months 12-M+1 to 12 of the previous year plus children
    # born in months 1 to 12-M of the current year
    denom_st_bymonth[,M+1] = 0
    
    # Need to calculate cumulative hazard up to age M for a child who could have experienced
    # dengue infection at age m in year t
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
        
        cumhaz_sn_stmcm[,c:12] = lambdaprevyear_mat[,c:12]/12 * hr_sn[,1:(12-c+1)]
        cumhaz_sn_stmcm[,1:(M+c-12)] = lambda_mat[,1:(M+c-12)]/12 * hr_sn[,(12-c+2):(M+1)]
        
        cumhaz_sp_stmcm[,c:12] = lambdaprevyear_mat[,c:12]/12 * hr_sp[,1:(12-c+1)]
        cumhaz_sp_stmcm[,1:(M+c-12)] = lambda_mat[,1:(M+c-12)]/12 * hr_sp[,(12-c+2):(M+1)]
        
      } else {
        
        cumhaz_sn_stmcm[,1:M] = lambda_mat[,1:M]/12 * hr_sn[,1:M]
        
        cumhaz_sp_stmcm[,1:M] = lambda_mat[,1:M]/12 * hr_sp[,1:M]
        
      }
      
      # Finally, add this cumulative hazard weighted by the relative number of children born in that month
      cumhaz_sn[,M+1] = cumhaz_sn[,M+1] + rowSums(cumhaz_sn_stmcm) * n_st_bymonth[,c]/n_byyear
      cumhaz_sp[,M+1] = cumhaz_sp[,M+1] + rowSums(cumhaz_sp_stmcm) * n_st_bymonth[,c]/n_byyear
      
    }
    
  }
  
  p_inf_sp = exp(-cumhaz_sp) * (1 - exp(-lambda_mat/12 * hr_sp))
  p_inf_sn = exp(-cumhaz_sn) * (1 - exp(-lambda_mat/12 * hr_sn))
  
  p_inf = msp_mat * p_inf_sp + (1 - msp_mat) * p_inf_sn
  
  p_sev = (msp_mat*p_sev_sp) * p_inf_sp + (1-msp_mat)*p_sev_sn * p_inf_sn
  
  pred_y_st_month = denom_st_bymonth * p_inf * p_rep
  pred_sev_st_month = denom_st_bymonth * p_sev * p_rep
  
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
                      lambda_st,
                      msp_st,
                      msp1_st,
                      mspmult_st,
                      fixedpar) {
  
  expit = function(x) {1/(1+exp(-x))}
  logit = function(x) {log(x/(1-x))}
  
  if (length(intersect(names(par),names(fixedpar)))>0) {stop("Repetition of parameters in par and fixedpar")}
  
  allpars = c(par,fixedpar)
  
  mat_protection_decay = allpars[['mat_protection_decay']]
  baseline_hazmatprot_sn = allpars[['baseline_hazmatprot_sn']]
  baseline_hazmatprot_sp = allpars[['baseline_hazmatprot_sp']]
  baseline_hazmult_sn = allpars[['baseline_hazmult_sn']]
  baseline_hazmult_sp = allpars[['baseline_hazmult_sp']]
  oneyear_probsev = allpars[['oneyear_probsev']]
  adepeak_logor = allpars[['adepeak_logor']]
  ademode = allpars[['ademode']]
  adekappa = allpars[['adekappa']]
  
  baseline_replogor = allpars[['baseline_replogor']]
  baseline_sevlogor = allpars[['baseline_sevlogor']]
  
  alphast = rep(NA,nst)
  for (i in 1:nst) {
    alphast[i] = allpars[[paste0('alpha',i)]]
  }
  
  lambda_mat = matrix(lambda_st,nrow=nst,ncol=nmonth,byrow=F)
  msp_mat = matrix(msp_st,nrow=nst,ncol=nmonth,byrow=F)
  msp1_mat = matrix(msp1_st,nrow=nst,ncol=nmonth,byrow=F)
  mspmult_mat = matrix(mspmult_st,nrow=nst,ncol=nmonth,byrow=F)
  
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
  
  hr_sn = c(1 + baseline_hazmult_sn,1 - baseline_hazmatprot_sn * exp(-mat_protection_decay * 0:(nmonth-2)))
  hr_sp = c(1 + baseline_hazmult_sp,1 - baseline_hazmatprot_sp * exp(-mat_protection_decay * 0:(nmonth-2)))
  
  or_rep = exp(baseline_replogor * as.numeric(Ms==0))
  
  rep_sn = hr_sn * or_rep
  rep_sp = hr_sp * or_rep
  
  or_sev_sn = exp(baseline_sevlogor * as.numeric(Ms==0))
  
  or_sev_sp = exp(baseline_sevlogor * as.numeric(Ms==0) + 
                    ade_func(Ms,ademode,adekappa,adepeak_logor,1,nmonth-1))
  
  orrepsev_sn = rep_sn * or_sev_sn
  orrepsev_sp = rep_sp * or_sev_sp
  
  return(list(hr_sn,hr_sp,rep_sn,rep_sp,or_sev_sn,or_sev_sp,orrepsev_sn,orrepsev_sp))
}

expit = function(x) {1/(1+exp(-x))}
logit = function(x) {log(x/(1-x))}


### Generate FOI and maternal seroprevalence across states for 100 years
### Use birth and maternal age data to generate hazards for infants in each state-year
### Use hazards to generate number of infant cases and severe cases by month

state_codes = read.table('C:/Users/mhitchings/Dropbox (UFL)/DengueInfants/Data/01-processedData/mapping/state.txt')
colnames(state_codes) = c("letter","number")

all_den = read.csv('./Data/02-analysisData/DengueCasesByStateYear_2000to2024.csv')

all_den = all_den %>% arrange(state,year_not) %>% filter(year_not <= 2014)

statenums = unique(all_den$state)

nst = nrow(all_den)
ns = length(unique(all_den$state))
ny = 2014-2000+1
nmonth=12


y_st_bymonth=as.matrix(all_den[,paste0('infantcases_',0:11,'month')])
sev_st_bymonth=as.matrix(all_den[,paste0('infantseverecases_',0:11,'month')])

include_ade = 1
include_sevagefrailty = 0
include_infagefrailty = 1
include_agereporting = 0

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

namepars_fortable = c("mat_protection_decay",
                      "onemonth_protection_sn",
                      "onemonth_protection_sp",
                      "neonate_hr_sn",
                      "neonate_hr_sp",
                      "baseline_replogor",
                      "oneyear_probsev",
                      "baseline_sevlogor",
                      "adepeak_logor",
                      "ademode",
                      "adekappa")

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

# Read model
p = readRDS(paste0('./stan_fits/infantmodel',filesuffix,'.rds'))

# Model parameters into a matrix and apply some transformations to interpretable model parameters
post_pars = as.data.frame(as.matrix(p))
post_pars = post_pars[,c(namepars,grep('alphast_param',colnames(post_pars),value=T))]

if (alphaconst==1) {
  post_pars_forpred = cbind(post_pars[,1:length(namepars)],
                            post_pars[,rep((length(namepars)+1):(length(namepars)+ns),each=ny)])
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
} else if (alphalin==1) {
  alphaparams = post_pars[,grepl('alphast_param',colnames(post_pars))]
  allalphas = matrix(NA,nrow=nrow(post_pars),ncol=nst)
  
  for (t in 1:ny) {
    for (s in 1:ns) {
      allalphas[,(s - 1)*ny + t] = expit(logit(alphaparams[,s]) + (t-1) * alphaparams[,s+ns]);
    }
  }
  
  post_pars_forpred = cbind(post_pars[,1:length(namepars)],allalphas)
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
} else {
  post_pars_forpred = post_pars
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
}

meancifitpars = apply(post_pars_forpred,2,function(x) quantile(x,c(0.5,0.025,0.975)))

meancifitpars = meancifitpars[,colnames(meancifitpars) %in% namepars]
meancifitpars[,'baseline_hazmult_sn'] = 1 + meancifitpars[,'baseline_hazmult_sn']
meancifitpars[,'baseline_hazmult_sp'] = 1 + meancifitpars[,'baseline_hazmult_sp']
meancifitpars[,'baseline_hazmatprot_sn'] = 1 - meancifitpars[,'baseline_hazmatprot_sn']
meancifitpars[,'baseline_hazmatprot_sp'] = 1 - meancifitpars[,'baseline_hazmatprot_sp']
meancifitpars[,'baseline_replogor'] = exp(meancifitpars[,'baseline_replogor'])
meancifitpars[,'baseline_sevlogor'] = exp(meancifitpars[,'baseline_sevlogor'])
meancifitpars[,'oneyear_probsev'] = 100*expit(meancifitpars[,'oneyear_probsev'])
meancifitpars[,'adepeak_logor'] = exp(meancifitpars[,'adepeak_logor'])
meancifitpars[,'ademode'] = 2 + 10 * meancifitpars[,'ademode']

partable = data.frame('Parameter'=namepars_fortable,
                      'MedianCrI'=paste0(sprintf('%.2f',round(meancifitpars[1,],2)),
                                         ' (',sprintf('%.2f',round(meancifitpars[2,],2)),
                                         ',',sprintf('%.2f',round(meancifitpars[3,],2)),')'))

write.csv(partable,'./Table_S2.csv',row.names=F)

### Figure S13: trace plot
png('./Figure_S13.png',width=8,height=9,units='in',res=180)
plot(stan_trace(p,c(setdiff(namepars,c('baseline_replogor','baseline_sevlogor')),'lp__'),ncol=3))
dev.off()

### Figure S14: pairs plot
png('./Figure_S14.png',width=15,height=15,units='in',res=180)
pairs(p,pars=namepars)
dev.off()

### HR at birth for seropositive vs. seronegative infants
rel_hr_atbirth = (1+post_pars_forpred$baseline_hazmult_sp)/(1+post_pars_forpred$baseline_hazmult_sn)
quantile(rel_hr_atbirth,c(0.5,0.025,0.975))

rel_hr_at1month = (1-post_pars_forpred$baseline_hazmatprot_sp)/(1-post_pars_forpred$baseline_hazmatprot_sn)
quantile(rel_hr_at1month,c(0.5,0.025,0.975))

hr_atbirth_sp = (1+post_pars_forpred$baseline_hazmult_sp)
quantile(hr_atbirth_sp,c(0.5,0.025,0.975))
hr_atbirth_sn = (1+post_pars_forpred$baseline_hazmult_sn)
quantile(hr_atbirth_sn,c(0.5,0.025,0.975))

hr_onemonth_sp = (1-post_pars_forpred$baseline_hazmatprot_sp)
quantile(hr_onemonth_sp,c(0.5,0.025,0.975))
hr_onemonth_sn = (1-post_pars_forpred$baseline_hazmatprot_sn)
quantile(hr_onemonth_sn,c(0.5,0.025,0.975))



reprate_forplot = post_pars_forpred[,grepl('alpha',colnames(post_pars_forpred))]
reprate_forplot_summ = as.data.frame(t(apply(reprate_forplot,2,function(x) quantile(x,c(0.5,0.025,0.975)))))
colnames(reprate_forplot_summ) = c('Median','LCI','UCI')
reprate_forplot_summ$State = all_den$statecode
reprate_forplot_summ$Year = all_den$year_not

ggsave('./Figure_S20.png',
       reprate_forplot_summ %>% ggplot(aes(x=Year,y=Median)) + 
         geom_ribbon(aes(x=Year,ymin=LCI,ymax=UCI),col='grey',alpha=0.5)+
         geom_line() +
         facet_wrap(vars(State),nrow=6,ncol=6)+
         scale_x_continuous(breaks=c(2000,2005,2010,2015),
                            labels=c('00','05','10','15'))+
         ylab('Reporting Rate')+scale_y_continuous(limits=c(0,1)),
       height=5,width=5,units='in',device='png')

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
    
    if (s == 16 & foi_dir == "./Data/02-foi/stanfit_v3/t46_rta_async_noType/") {
      t = t %>% bind_rows(data.frame('year'=2000,
                                     'Smult'=rep(t$Smult[1],1),
                                     'SP'=rep(t$SP[1],1),
                                     'S1'=rep(t$S1[1],1)
      )
      ) %>% arrange(year)
    }
    
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


lambda_st=all_den$foi
msp1_st=all_den$maternal_s1
mspmult_st=all_den$maternal_smult
msp_st = msp1_st+mspmult_st
n_st_bymonth = as.matrix(all_den[,paste0('infantbirths_',1:12,'month')])

pred = apply(post_pars_forpred,1,function(x)
  pred_cases(x,
             nmonth,
             nst,
             ns,
             ny,
             y_st_bymonth,
             n_st_bymonth,
             sev_st_bymonth,
             lambda_st,
             msp_st,
             msp1_st,
             mspmult_st,
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
  scale_linetype_manual(name='',values=c(1,2),labels=c('Fit','Observed'))

pred_vs_obs_sevcases_brazil_5year = pred_vs_obs_sevcases_byiter %>% 
  group_by(y5_label,Month,iter) %>% summarise(Observed=sum(Observed),
                                              Fit=sum(Pred))

ggsave('./Figure_S23.png',
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
       height=2,width=6,units='in',device='png')

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
  bind_rows(predprop_byregion_5year %>% mutate(Case = "All Cases"),
            predsevprop_byregion_5year %>% mutate(Case = "Severe Cases"))

ggsave('./Figure_S17.png',
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
       height=5,width=5,units='in',device='png')


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
predy_bystate_total$State = rep(unique(all_den$statecode),each=12)
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
'State'=rep(unique(all_den$statecode),each=12),
'Group'='Observed')
predy_bystate_total = bind_rows(predy_bystate_total,obsy_bystate_total)

ggsave('./Figure_S21.png',
  ggplot(predy_bystate_total) + geom_line(aes(x=Month,y=Median,linetype=Group)) +
    geom_ribbon(aes(x=Month,ymin=LCI,ymax=UCI,fill=Group),show.legend=F,alpha=0.5)+
    ylim(c(0,max(c(predy_bystate_total$Median,predy_bystate_total$UCI),na.rm=T)))+
    ylab('Reported Cases')+
    facet_wrap(vars(State),nrow=5,ncol=6)+
    scale_x_continuous(breaks=seq(0,12,3))+
    scale_linetype_manual(name='',values=c(1,2),labels=c('Fit','Observed'))+
    scale_fill_manual(name='',values=c('grey','grey'),labels=c('Fit','Observed'))+
    theme(legend.position=c(0.9,0.1)),
  height=6,width=6,units='in',device='png')

## Severe cases
predsev_age_dist_bystate = lapply(pred,function(p) t(sapply(1:ns,function(x) colSums(p[[2]][state_indices[x,],]))))

predsev_bystate_total = data.frame('Month' = rep(1:12,ns))
predsev_bystate_total$State = rep(unique(all_den$statecode),each=12)
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
  'State'=rep(unique(all_den$statecode),each=12),
  'Group'='Observed')
predsev_bystate_total = bind_rows(predsev_bystate_total,obssev_bystate_total)

ggsave('./Figure_S22.png',
  ggplot(predsev_bystate_total) + geom_line(aes(x=Month,y=Median,linetype=Group)) +
    geom_ribbon(aes(x=Month,ymin=LCI,ymax=UCI,fill=Group),show.legend=F,alpha=0.5)+
    ylim(c(0,max(c(predsev_bystate_total$Median,predsev_bystate_total$UCI),na.rm=T)))+
    ylab('Reported Severe Cases')+
    facet_wrap(vars(State),nrow=5,ncol=6)+
    scale_x_continuous(breaks=seq(0,12,3))+
    scale_linetype_manual(name='',values=c(1,2),labels=c('Fit','Observed'))+
    scale_fill_manual(name='',values=c('grey','grey'),labels=c('Fit','Observed'))+
    theme(legend.position=c(0.9,0.1)),
  height=6,width=6,units='in',device='png')

## Age profiles for Figure 4
profiles = apply(post_pars_forpred,1,function(x)
  ageprofile(x,
             nmonth,
             nst,
             ns,
             ny,
             y_st_bymonth,
             n_st_bymonth,
             sev_st_bymonth,
             lambda_st,
             msp_st,
             msp1_st,
             mspmult_st,
             c())
)

hr_sn = matrix(unlist(lapply(profiles,function(x) x[[1]])),nrow=length(profiles),ncol=12,byrow=T)
hr_sp = matrix(unlist(lapply(profiles,function(x) x[[2]])),nrow=length(profiles),ncol=12,byrow=T)
rep_sn = matrix(unlist(lapply(profiles,function(x) x[[3]])),nrow=length(profiles),ncol=12,byrow=T)
rep_sp = matrix(unlist(lapply(profiles,function(x) x[[4]])),nrow=length(profiles),ncol=12,byrow=T)
or_sev_sn = matrix(unlist(lapply(profiles,function(x) x[[5]])),nrow=length(profiles),ncol=12,byrow=T)
or_sev_sp = matrix(unlist(lapply(profiles,function(x) x[[6]])),nrow=length(profiles),ncol=12,byrow=T)
orrepsev_sn = matrix(unlist(lapply(profiles,function(x) x[[7]])),nrow=length(profiles),ncol=12,byrow=T)
orrepsev_sp = matrix(unlist(lapply(profiles,function(x) x[[8]])),nrow=length(profiles),ncol=12,byrow=T)

data_forplot = data.frame('Month'=rep(1:12,times=8),
                          'EventType'=c(rep('Infection Hazard Ratio',24),
                                        rep('Reported Case (OR)',24),
                                        rep('Severe Disease Odds Ratio',24),
                                        rep('Reported Sev. Case',24)),
                          'Serostatus'=rep(c(rep('Seronegative',12),
                                             rep('Seropositive',12)),4),
                          'Mean'=c(apply(hr_sn,2,mean),
                                   apply(hr_sp,2,mean),
                                   apply(rep_sn,2,mean),
                                   apply(rep_sp,2,mean),
                                   apply(or_sev_sn,2,mean),
                                   apply(or_sev_sp,2,mean),
                                   apply(orrepsev_sn,2,mean),
                                   apply(orrepsev_sp,2,mean)),
                          'LCI'=c(apply(hr_sn,2,function(x) quantile(x,0.025)),
                                  apply(hr_sp,2,function(x) quantile(x,0.025)),
                                  apply(rep_sn,2,function(x) quantile(x,0.025)),
                                  apply(rep_sp,2,function(x) quantile(x,0.025)),
                                  apply(or_sev_sn,2,function(x) quantile(x,0.025)),
                                  apply(or_sev_sp,2,function(x) quantile(x,0.025)),
                                  apply(orrepsev_sn,2,function(x) quantile(x,0.025)),
                                  apply(orrepsev_sp,2,function(x) quantile(x,0.025))),
                          'UCI'=c(apply(hr_sn,2,function(x) quantile(x,0.975)),
                                  apply(hr_sp,2,function(x) quantile(x,0.975)),
                                  apply(rep_sn,2,function(x) quantile(x,0.975)),
                                  apply(rep_sp,2,function(x) quantile(x,0.975)),
                                  apply(or_sev_sn,2,function(x) quantile(x,0.975)),
                                  apply(or_sev_sp,2,function(x) quantile(x,0.975)),
                                  apply(orrepsev_sn,2,function(x) quantile(x,0.975)),
                                  apply(orrepsev_sp,2,function(x) quantile(x,0.975)))
)
data_forplot$EventType = factor(data_forplot$EventType,levels=c('Infection Hazard Ratio',
                                                                'Severe Disease Odds Ratio',
                                                                'Reported Case (OR)',
                                                                'Reported Sev. Case'))

textsize=8


profileplot = ggplot() +
  geom_ribbon(data=data_forplot %>% 
                filter(EventType %in% c('Infection Hazard Ratio',
                                        'Severe Disease Odds Ratio')),
              aes(x=Month,ymin=LCI,ymax=UCI,linetype=Serostatus),alpha=0.2,fill='grey',col=NA,show.legend = F)+
  geom_line(data=data_forplot  %>% 
              filter(EventType %in% c('Infection Hazard Ratio',
                                      'Severe Disease Odds Ratio')),
            aes(x=as.numeric(Month),y=Mean,linetype=Serostatus),linewidth=0.25) +
  scale_x_continuous(name='Month',breaks=c(1,3,6,9,12),limits=c(1,12))+
  scale_linetype_manual(name=NULL,values=c(3,4)) + 
  ylab('Ratio')+
  expand_limits(y=c(0))+
  facet_wrap(vars(EventType),ncol=2,scales = 'free')+
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

predcase_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[1]]))),nrow=length(pred),ncol=12,byrow=T)
predsev_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[2]]))),nrow=length(pred),ncol=12,byrow=T)

data_forplot_pred = data.frame('Month'=rep(1:12,times=6),
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
  geom_ribbon(data=data_forplot_pred %>% filter(is.na(Mean)),aes(x=Month,ymin=LCI,ymax=UCI),alpha=0.2,fill='grey',col=NA)+
  geom_line(data=data_forplot_pred %>% filter(!is.na(Mean)),aes(x=as.numeric(Month),y=Mean,linetype=Type),linewidth=0.25) +
  scale_x_continuous(name='Month',breaks=c(1,3,6,9,12),limits=c(1,12))+
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

ggsave('./Figure_4.png',
  grid.arrange(predplot,profileplot,nrow=2),
  height=4,width=4,units='in',device='png')

### Figure S18 and S19: model fit and age profiles incorporating uncertainty in FOI
# maternal seroprevalence
lambda_uncertainty = 1

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

# Read model
p = readRDS(paste0('./stan_fits/infantmodel',filesuffix,'.rds'))

# Model parameters into a matrix and apply some transformations to interpretable model parameters
post_pars = as.data.frame(as.matrix(p))
post_pars = post_pars[,c(namepars,grep('alphast_param',colnames(post_pars),value=T))]

if (alphaconst==1) {
  post_pars_forpred = cbind(post_pars[,1:length(namepars)],
                            post_pars[,rep((length(namepars)+1):(length(namepars)+ns),each=ny)])
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
} else if (alphalin==1) {
  alphaparams = post_pars[,grepl('alphast_param',colnames(post_pars))]
  allalphas = matrix(NA,nrow=nrow(post_pars),ncol=nst)
  
  for (t in 1:ny) {
    for (s in 1:ns) {
      allalphas[,(s - 1)*ny + t] = expit(logit(alphaparams[,s]) + (t-1) * alphaparams[,s+ns]);
    }
  }
  
  post_pars_forpred = cbind(post_pars[,1:length(namepars)],allalphas)
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
} else {
  post_pars_forpred = post_pars
  colnames(post_pars_forpred) = c(namepars,paste0('alpha',1:nst))
}

pred = apply(post_pars_forpred,1,function(x)
  pred_cases(x,
             nmonth,
             nst,
             ns,
             ny,
             y_st_bymonth,
             n_st_bymonth,
             sev_st_bymonth,
             lambda_st,
             msp_st,
             msp1_st,
             mspmult_st,
             c())
)

profiles = apply(post_pars_forpred,1,function(x)
  ageprofile(x,
             nmonth,
             nst,
             ns,
             ny,
             y_st_bymonth,
             n_st_bymonth,
             sev_st_bymonth,
             lambda_st,
             msp_st,
             msp1_st,
             mspmult_st,
             c())
)

hr_sn = matrix(unlist(lapply(profiles,function(x) x[[1]])),nrow=length(profiles),ncol=12,byrow=T)
hr_sp = matrix(unlist(lapply(profiles,function(x) x[[2]])),nrow=length(profiles),ncol=12,byrow=T)
rep_sn = matrix(unlist(lapply(profiles,function(x) x[[3]])),nrow=length(profiles),ncol=12,byrow=T)
rep_sp = matrix(unlist(lapply(profiles,function(x) x[[4]])),nrow=length(profiles),ncol=12,byrow=T)
or_sev_sn = matrix(unlist(lapply(profiles,function(x) x[[5]])),nrow=length(profiles),ncol=12,byrow=T)
or_sev_sp = matrix(unlist(lapply(profiles,function(x) x[[6]])),nrow=length(profiles),ncol=12,byrow=T)
orrepsev_sn = matrix(unlist(lapply(profiles,function(x) x[[7]])),nrow=length(profiles),ncol=12,byrow=T)
orrepsev_sp = matrix(unlist(lapply(profiles,function(x) x[[8]])),nrow=length(profiles),ncol=12,byrow=T)

data_forplot = data.frame('Month'=rep(1:12,times=8),
                          'EventType'=c(rep('Infection Hazard Ratio',24),
                                        rep('Reported Case Odds Ratio',24),
                                        rep('Severe Disease Odds Ratio',24),
                                        rep('Reported Sev. Case',24)),
                          'Serostatus'=rep(c(rep('Seronegative',12),
                                             rep('Seropositive',12)),4),
                          'Mean'=c(apply(hr_sn,2,mean),
                                   apply(hr_sp,2,mean),
                                   apply(rep_sn,2,mean),
                                   apply(rep_sp,2,mean),
                                   apply(or_sev_sn,2,mean),
                                   apply(or_sev_sp,2,mean),
                                   apply(orrepsev_sn,2,mean),
                                   apply(orrepsev_sp,2,mean)),
                          'LCI'=c(apply(hr_sn,2,function(x) quantile(x,0.025)),
                                  apply(hr_sp,2,function(x) quantile(x,0.025)),
                                  apply(rep_sn,2,function(x) quantile(x,0.025)),
                                  apply(rep_sp,2,function(x) quantile(x,0.025)),
                                  apply(or_sev_sn,2,function(x) quantile(x,0.025)),
                                  apply(or_sev_sp,2,function(x) quantile(x,0.025)),
                                  apply(orrepsev_sn,2,function(x) quantile(x,0.025)),
                                  apply(orrepsev_sp,2,function(x) quantile(x,0.025))),
                          'UCI'=c(apply(hr_sn,2,function(x) quantile(x,0.975)),
                                  apply(hr_sp,2,function(x) quantile(x,0.975)),
                                  apply(rep_sn,2,function(x) quantile(x,0.975)),
                                  apply(rep_sp,2,function(x) quantile(x,0.975)),
                                  apply(or_sev_sn,2,function(x) quantile(x,0.975)),
                                  apply(or_sev_sp,2,function(x) quantile(x,0.975)),
                                  apply(orrepsev_sn,2,function(x) quantile(x,0.975)),
                                  apply(orrepsev_sp,2,function(x) quantile(x,0.975)))
)
data_forplot$EventType = factor(data_forplot$EventType,levels=c('Infection Hazard Ratio',
                                                                'Severe Disease Odds Ratio',
                                                                'Reported Case Odds Ratio',
                                                                'Reported Sev. Case'))

textsize=8


profileplot = ggplot() +
  geom_ribbon(data=data_forplot,
              aes(x=Month,ymin=LCI,ymax=UCI,linetype=Serostatus),alpha=0.2,fill='grey',col=NA,show.legend = F)+
  geom_line(data=data_forplot,
            aes(x=as.numeric(Month),y=Mean,linetype=Serostatus),linewidth=0.25) +
  scale_x_continuous(name='Month',breaks=c(1,3,6,9,12),limits=c(1,12))+
  scale_linetype_manual(name=NULL,values=c(3,4)) + 
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

predcase_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[1]]))),nrow=length(pred),ncol=12,byrow=T)
predsev_age_dist = matrix(unlist(lapply(pred,function(x) colSums(x[[2]]))),nrow=length(pred),ncol=12,byrow=T)

data_forplot_pred = data.frame('Month'=rep(1:12,times=6),
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
  geom_ribbon(data=data_forplot_pred %>% filter(is.na(Mean)),aes(x=Month,ymin=LCI,ymax=UCI),alpha=0.2,fill='grey',col=NA)+
  geom_line(data=data_forplot_pred %>% filter(!is.na(Mean)),aes(x=as.numeric(Month),y=Mean,linetype=Type),linewidth=0.25) +
  scale_x_continuous(name='Month',breaks=c(1,3,6,9,12),limits=c(1,12))+
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

ggsave('./Figure_S18.png',
       predplot,
       height=2,width=4,units='in',device='png')
ggsave('./Figure_S19.png',
       profileplot,
       height=4,width=4,units='in',device='png')

### Figure S24: change in incidence and mean infant age with maternal seroprevalence
## according to model

## Change pred function to have fixed input reporting rate
pred_cases_fixedreprate = function(par,
                      nmonth,
                      nst,
                      ns,
                      ny,
                      y_st_bymonth,
                      n_st_bymonth,
                      sev_st_bymonth,
                      lambda_st,
                      msp_st,
                      msp1_st,
                      mspmult_st,
                      reporting_rate,
                      fixedpar) {
  
  expit = function(x) {1/(1+exp(-x))}
  logit = function(x) {log(x/(1-x))}
  
  if (length(intersect(names(par),names(fixedpar)))>0) {stop("Repetition of parameters in par and fixedpar")}
  
  allpars = c(par,fixedpar)
  
  mat_protection_decay = allpars[['mat_protection_decay']]
  baseline_hazmatprot_sn = allpars[['baseline_hazmatprot_sn']]
  baseline_hazmatprot_sp = allpars[['baseline_hazmatprot_sp']]
  baseline_hazmult_sn = allpars[['baseline_hazmult_sn']]
  baseline_hazmult_sp = allpars[['baseline_hazmult_sp']]
  oneyear_probsev = allpars[['oneyear_probsev']]
  adepeak_logor = allpars[['adepeak_logor']]
  ademode = allpars[['ademode']]
  adekappa = allpars[['adekappa']]
  
  baseline_replogor = allpars[['baseline_replogor']]
  baseline_sevlogor = allpars[['baseline_sevlogor']]
  
  alphast = reporting_rate
  
  lambda_mat = matrix(lambda_st,nrow=nst,ncol=nmonth,byrow=F)
  msp_mat = matrix(msp_st,nrow=nst,ncol=nmonth,byrow=F)
  msp1_mat = matrix(msp1_st,nrow=nst,ncol=nmonth,byrow=F)
  mspmult_mat = matrix(mspmult_st,nrow=nst,ncol=nmonth,byrow=F)
  
  lambdaprevyear_mat = lambda_mat
  
  n_st_bymonth_prevyear = n_st_bymonth
  
  cumhaz_sp=matrix(NA,nrow=nst,ncol=nmonth)
  cumhaz_sn=matrix(NA,nrow=nst,ncol=nmonth)
  
  denom_st_bymonth=matrix(NA,nrow=nst,ncol=nmonth)
  
  Ms = 0:(nmonth-1)
  
  # Hazard ratio of infection in for infants born to seronegative and seropositive mothers
  hr_sp = matrix(c(1 + baseline_hazmult_sp,1 - baseline_hazmatprot_sp * exp(-mat_protection_decay * 0:(nmonth-2))),nrow=nst,ncol=nmonth,byrow=T)
  hr_sn = matrix(c(1 + baseline_hazmult_sn,1 - baseline_hazmatprot_sn * exp(-mat_protection_decay * 0:(nmonth-2))),nrow=nst,ncol=nmonth,byrow=T)
  
  # Odds of severe disease upon infection in infants born to seronegative and seropositive mothers
  p_sev_sn = matrix(expit(oneyear_probsev + baseline_sevlogor * as.numeric(Ms==0)),nrow=nst,ncol=nmonth,byrow=T)
  
  p_sev_sp = matrix(expit(oneyear_probsev + baseline_sevlogor * as.numeric(Ms==0) + 
                            ade_func(Ms,ademode,adekappa,adepeak_logor,1,nmonth-1)),
                    nrow=nst,ncol=nmonth,byrow=T)
  
  # Reporting odds
  p_rep = expit(logit(matrix(alphast,nrow=nst,ncol=nmonth,byrow=F)) + 
                  matrix(baseline_replogor * as.numeric(Ms==0),nrow=nst,ncol=nmonth,byrow=T))
  
  # We will need the total children born in each calendar year
  n_byyear = apply(n_st_bymonth,1,sum)
  
  cumhaz_sn[,1] = 0
  cumhaz_sp[,1] = 0
  
  denom_st_bymonth[,1] = rowSums(n_st_bymonth)
  
  for (M in 1:(nmonth-1)) {
    
    # Build up cumulative hazard by state-year and month of age
    cumhaz_sp[,M+1] = 0
    cumhaz_sn[,M+1] = 0
    
    # Need to calculate the denominator (who is at risk of having case at age M in year t)
    # For this, we add up all children who could have been M months old during year t
    # i.e. children who were born in months 12-M+1 to 12 of the previous year plus children
    # born in months 1 to 12-M of the current year
    denom_st_bymonth[,M+1] = 0
    
    # Need to calculate cumulative hazard up to age M for a child who could have experienced
    # dengue infection at age m in year t
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
        
        cumhaz_sn_stmcm[,c:12] = lambdaprevyear_mat[,c:12]/12 * hr_sn[,1:(12-c+1)]
        cumhaz_sn_stmcm[,1:(M+c-12)] = lambda_mat[,1:(M+c-12)]/12 * hr_sn[,(12-c+2):(M+1)]
        
        cumhaz_sp_stmcm[,c:12] = lambdaprevyear_mat[,c:12]/12 * hr_sp[,1:(12-c+1)]
        cumhaz_sp_stmcm[,1:(M+c-12)] = lambda_mat[,1:(M+c-12)]/12 * hr_sp[,(12-c+2):(M+1)]
        
      } else {
        
        cumhaz_sn_stmcm[,1:M] = lambda_mat[,1:M]/12 * hr_sn[,1:M]
        
        cumhaz_sp_stmcm[,1:M] = lambda_mat[,1:M]/12 * hr_sp[,1:M]
        
      }
      
      # Finally, add this cumulative hazard weighted by the relative number of children born in that month
      cumhaz_sn[,M+1] = cumhaz_sn[,M+1] + rowSums(cumhaz_sn_stmcm) * n_st_bymonth[,c]/n_byyear
      cumhaz_sp[,M+1] = cumhaz_sp[,M+1] + rowSums(cumhaz_sp_stmcm) * n_st_bymonth[,c]/n_byyear
      
    }
    
  }
  
  p_inf_sp = exp(-cumhaz_sp) * (1 - exp(-lambda_mat/12 * hr_sp))
  p_inf_sn = exp(-cumhaz_sn) * (1 - exp(-lambda_mat/12 * hr_sn))
  
  p_inf = msp_mat * p_inf_sp + (1 - msp_mat) * p_inf_sn
  
  p_sev = (msp_mat*p_sev_sp) * p_inf_sp + (1-msp_mat)*p_sev_sn * p_inf_sn
  
  pred_y_st_month = denom_st_bymonth * p_inf * p_rep
  pred_sev_st_month = denom_st_bymonth * p_sev * p_rep
  
  return(list(pred_y_st_month,pred_sev_st_month))
}

lambda_uncertainty=0
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

p = readRDS(paste0('./stan_fits/infantmodel',filesuffix,'.rds'))

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

meanfitpars = apply(post_pars_forpred,2,median)

nst=1
ns=1
ny=1
n_st_bymonth = matrix(10000,nrow=1,ncol=nmonth)
nmonth=12

lambda_st = 0.05
msp_st = 0
msp1_st = 0
mspmult_st = 0

reporting_rate = 0.2
msps = seq(0,1,by=0.05)

for (i in 1:length(msps)) {
  msp = msps[i]
  
  pred = pred_cases_fixedreprate(meanfitpars,
             nmonth,
             nst,
             ns,
             ny,
             y_st_bymonth,
             n_st_bymonth,
             sev_st_bymonth,
             lambda_st,
             msp,
             msp/2,
             msp/2,
             reporting_rate,
             c())
  
  if (i == 1) {
    
    df_casenums = data.frame('MaternalSP'=msp,
                             'Month'=rep(0:11,2),
                             'Type'=rep(c('All Cases','Severe Cases'),each=12),
                             'Cases'=c(pred[[1]],pred[[2]]))
    
  } else {
    
    df_casenums = rbind(df_casenums,
                        data.frame('MaternalSP'=msp,
                                   'Month'=rep(0:11,2),
                                   'Type'=rep(c('All Cases','Severe Cases'),each=12),
                                   'Cases'=c(pred[[1]],pred[[2]]))
    )
    
    
  }
  
  
  
}

df_ir = df_casenums %>% group_by(MaternalSP,Type) %>% 
  summarise(IR = sum(Cases)/sum(n_st_bymonth))

p_ir = df_ir %>% group_by(Type) %>% mutate(IRR = IR/IR[1]) %>% 
  ggplot() + geom_line(aes(x=MaternalSP,y=IRR)) + 
  facet_wrap(vars(Type)) + 
  ylim(c(0,1.5)) + 
  xlab('') + 
  ylab('Incidence Rate Ratio Relative\n To No Seroprevalence')

p_age = df_casenums %>% group_by(MaternalSP,Type) %>% 
  summarise(MeanAge = sum(Cases*Month)/sum(Cases)) %>%
  ggplot() + geom_line(aes(x=MaternalSP,y=MeanAge)) + 
  facet_wrap(vars(Type),scales='free') + 
  xlab('Maternal Seroprevalence') + 
  ylab('Mean Age Of Infant Case')

ggsave('./Figure_S24.png',
       grid.arrange(p_ir,p_age,nrow=2),
       height=6,width=6,units='in',device='png')
