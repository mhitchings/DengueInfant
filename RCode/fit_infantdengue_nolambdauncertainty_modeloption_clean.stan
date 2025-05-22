//
functions {
  
  real expit(real x) {
    
    return 1/(1+exp(-x));
  }
  
}

data {
  int nmonth; // Number of months of infant dengue cases
  int nst; // Number of state-years of data
  int ns; // Number of states
  int ny; // Number of years
  int y_st_bymonth[nst,nmonth]; // The number of infant dengue cases by month, and by state-year
  real n_st_bymonth[nst,nmonth]; // The number of infants born in each state-year
  int sev_st_bymonth[nst,nmonth]; // The number of infant severe dengue cases by month, and by state-year
  real lambda_st[nst]; // Dengue FOI in each state-year
  real msp1_st[nst]; // Maternal monotypic seroprevalence in each state-year
  real mspmult_st[nst]; // Maternal multitypic seroprevalence in each state-year
  int<lower=0,upper=1> alpha_constant_t; //Whether reporting rate is constant over time
  int<lower=0,upper=1> alpha_linear; //Whether reporting rate is linear over time
  int<lower=0,upper=1> include_ade; //Whether or not to include ADE of severe disease
  int<lower=0,upper=1> include_sevagefrailty; //Whether or not to include increasing risk of severe disease with younger age
  int<lower=0,upper=1> include_infagefrailty; //Whether or not to include increasing risk of infection with younger age
  int<lower=0,upper=1> include_agereporting; //Whether or not to include increasing risk of reporting with younger age
  int<lower=0,upper=1> include_maternalprotection; //Whether or not to include maternal protection against infection
  int<lower=0,upper=1> ade_dist; // Distribution for ADE odds ratio (0=log-normal, 1=normal)
}

parameters {
  real<lower=0,upper=1> alphast_param[(alpha_linear + alpha_constant_t) ? ns*(2 * alpha_linear + alpha_constant_t):nst];
  
  // Maternal protection against infection
  real<lower=0,upper=8> mat_protection_decay_param[include_maternalprotection ? 1:0];
  real<lower=0,upper=1> baseline_hazmatprot_sp_param[include_maternalprotection ? 1:0];
  real<upper=1> baseline_hazmatprot_sn_param[include_maternalprotection ? 1:0];
  
  // Decaying risk of infection with age
  real<lower=0> baseline_hazinfmult_param[include_infagefrailty ? 1:0];
  real<lower=0,upper=8> ageinfection_decay_param[include_infagefrailty ? 1:0];
  
  // Decaying risk of reporting with age
  real<lower=0,upper=5> baseline_replogor_param[include_agereporting ? 1:0];
  real<lower=0,upper=8> agereporting_decay_param[include_agereporting ? 1:0];
  
  // Changing risk of severe infection by age - elevated at birth and exponential decrease
  real<lower=0,upper=8> ageseverity_decay_param[include_sevagefrailty ? 1:0];
  real<lower=-3,upper=3> baseline_sevlogor_param[include_sevagefrailty ? 1:0];
  real<lower=-10,upper=0> oneyear_probsev;
  
  // Elevated risk of severe infection in children born to seropositive mothers - lognormal
  real<lower=0,upper=30> titersev_mult_param[include_ade ? 1:0];
  real<lower=0,upper=12> titersev_meanlog_param[include_ade ? 1:0];
  real<lower=0,upper=3> titersev_sdlog_param[include_ade ? 1:0];
  
}

transformed parameters {
  
  // Maternal protection against infection
  real<lower=0,upper=1> alphast[nst];
  real<lower=0,upper=8> mat_protection_decay;
  real<lower=0,upper=1> baseline_hazmatprot_sp;
  real<upper=1> baseline_hazmatprot_sn;
  real<lower=0> baseline_hazinfmult;
  real<lower=0,upper=8> ageinfection_decay;
  real<lower=0,upper=5> baseline_replogor;
  real<lower=0,upper=8> agereporting_decay;
  real<lower=-3,upper=3> baseline_sevlogor;
  real<lower=0,upper=8> ageseverity_decay;
  real<lower=0,upper=30> titersev_mult;
  real<lower=0,upper=12> titersev_meanlog;
  real<lower=0,upper=3> titersev_sdlog;
  real<lower=0,upper=3> titersevmulti_meanlogdiff;
  real<lower=0,upper=3> titersevmulti_sdlog;
  
  real<lower=0> phi_case;
  real<lower=0> phi_sev;
  
  if (alpha_constant_t) {
    for (t in 1:ny) {
      for (s in 1:ns) {
        alphast[(s - 1)*ny + t] = alphast_param[s];
      }
    }
  } else if (alpha_linear) {
    for (t in 1:ny) {
      for (s in 1:ns) {
        alphast[(s - 1)*ny + t] = expit(logit(alphast_param[s]) + (t-1) * alphast_param[s+ns]);
      }
    }
  } else {
    alphast = alphast_param;
  }
  
  if (include_maternalprotection) {
    mat_protection_decay = mat_protection_decay_param[1];
    baseline_hazmatprot_sp = baseline_hazmatprot_sp_param[1];
    baseline_hazmatprot_sn = baseline_hazmatprot_sn_param[1];
  } else {
    mat_protection_decay = 1;
    baseline_hazmatprot_sp = 0;
    baseline_hazmatprot_sn = 0;
    
  }
  
  if (include_infagefrailty) {
    baseline_hazinfmult = baseline_hazinfmult_param[1];
    ageinfection_decay = ageinfection_decay_param[1];
  } else {
    baseline_hazinfmult = 0;
    ageinfection_decay = 1;
  }
  
  if (include_agereporting) {
    baseline_replogor = baseline_replogor_param[1];
    agereporting_decay = agereporting_decay_param[1];
  } else {
    baseline_replogor = 0;
    agereporting_decay = 1;
  }
  
  if (include_sevagefrailty) {
    baseline_sevlogor = baseline_sevlogor_param[1];
    ageseverity_decay = ageseverity_decay_param[1];
  } else {
    baseline_sevlogor = 0;
    ageseverity_decay = 1;
  }
  
  if (include_ade) {
    titersev_mult = titersev_mult_param[1];
    titersev_meanlog = titersev_meanlog_param[1];
    titersev_sdlog = titersev_sdlog_param[1];
    
  } else {
    titersev_mult = 0;
    titersev_meanlog = log(100);
    titersev_sdlog = 0.1;
  }
  
}

model {
  
  real msp_st[nst];
  
  real hr_maternalprotection_sp[nmonth];
  real hr_maternalprotection_sn[nmonth];
  real hr_agefrailty[nmonth];
  real cumhaz_sp;
  real cumhaz_sn;
  real p_inf_sp;
  real p_inf_sn;
  real p_sev_sp[nmonth];
  real p_sev_sn[nmonth];
  real p_sev;
  real p_inf;
  real p_rep_mild;
  real p_rep_sev;
  
  // Need to keep track of the children that are at risk of having a case at age 0 in calendar year t
  real denom_st_bymonth;
  
  // Some temporary variables to help with cumulative hazard calculation
  real n_byyear;
  real cumhaz_sn_stmc;
  real cumhaz_sp_stmc;
  
  for (M in 0:(nmonth-1)) {
    
    hr_maternalprotection_sp[M+1] = 1 - baseline_hazmatprot_sp * exp(-mat_protection_decay * M);
    hr_maternalprotection_sn[M+1] = 1 - baseline_hazmatprot_sn * exp(-mat_protection_decay * M);
    hr_agefrailty[M+1] = 1 + baseline_hazinfmult * exp(-ageinfection_decay * M);
    
    p_sev_sn[M+1] = expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*M));
      
    // probability of severe disease by month
    if (ade_dist==0) {
      p_sev_sp[M+1] = expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*M) + titersev_mult * exp(lognormal_lpdf(M | titersev_meanlog, titersev_sdlog)));
    } else {
      p_sev_sp[M+1] = expit(oneyear_probsev + baseline_sevlogor * exp(-ageseverity_decay*M) + titersev_mult * exp(normal_lpdf(M | titersev_meanlog, titersev_sdlog)));
    }
  }
  
  for (st in 1:nst) {
    
    msp_st[st] = msp1_st[st] + mspmult_st[st];
    
    // We will need the total children born in each calendar year
    n_byyear = 0;
      
    for (c in 1:12) {
      n_byyear = n_byyear + n_st_bymonth[st,c];
    }
    
    for (M in 0:(nmonth-1)) {

      // Build up cumulative hazard by state-year and month of age
      cumhaz_sp = 0;
      cumhaz_sn = 0;
      
      // Need to calculate cumulative hazard up to age M for a child who could have experienced
      // dengue infection at age M in year t
      // For this, need to loop over children born in calendar month 12-M+1 of the previous year to 
      // calendar month 12-M of the current year
      for (c in 1:12) {
        cumhaz_sn_stmc = 0;
        cumhaz_sp_stmc = 0;
        
        // If c > 12-M, this means they were born in the previous year, so some months 
        // of a child's life were experienced in the previous year
        if (c > 12-M) {
        
          for (m in 0:M) {
          
            // If their being m months old occurred in the previous calendar year (i.e. if m<=12-c) use last year's
            // population FOI
            if (m<=12-c) {
              
              // Special case as I don't have FOI from 1999, so just 2000
              if ((st-1) % ny ==0) {
                
                cumhaz_sn_stmc += lambda_st[st]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sn[m+1];
                cumhaz_sp_stmc += lambda_st[st]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sp[m+1];
                
              } else {
                
                // Data is ordered by state first, then by year, so st-1 is the index of the same state's previous year
                cumhaz_sn_stmc += lambda_st[st-1]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sn[m+1];
                cumhaz_sp_stmc += lambda_st[st-1]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sp[m+1];
                
              }
              
            } else {
              // If their being m months old occurred in the same calendar year (i.e. if m>12-c) use this year's population FOI
              
              cumhaz_sn_stmc += lambda_st[st]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sn[m+1];
              cumhaz_sp_stmc += lambda_st[st]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sp[m+1];
              
            }
          
          }
          
        } else {
          
          for (m in 0:M) {
            
            cumhaz_sn_stmc += lambda_st[st]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sn[m+1];
            cumhaz_sp_stmc += lambda_st[st]/12 * hr_agefrailty[m+1] * hr_maternalprotection_sp[m+1];
            
          }
          
        }
          
        // Finally, add this cumulative hazard weighted by the relative number of children born in that month
        cumhaz_sn += cumhaz_sn_stmc * n_st_bymonth[st,c]/n_byyear;
        cumhaz_sp += cumhaz_sp_stmc * n_st_bymonth[st,c]/n_byyear;
        
      }
      
      p_inf_sp = exp(-cumhaz_sp) * (1 - exp(-lambda_st[st]/12 * hr_maternalprotection_sp[M+1] * hr_agefrailty[M+1]));
      p_inf_sn = exp(-cumhaz_sn) * (1 - exp(-lambda_st[st]/12 * hr_maternalprotection_sn[M+1] * hr_agefrailty[M+1]));
        
      p_inf = msp_st[st] * p_inf_sp + (1 - mspmult_st[st] - msp1_st[st]) * p_inf_sn;
      
      p_rep_mild = expit(logit(alphast[st]) + baseline_replogor * exp(-agereporting_decay * M));
      
      p_sev = (msp_st[st]*p_sev_sp[M+1]) * p_inf_sp + (1-msp_st[st])*p_sev_sn[M+1] * p_inf_sn;
      
      p_rep_sev = expit(logit(alphast[st]) + baseline_replogor * exp(-agereporting_decay * M));
      
      // Need to calculate the denominator (who is at risk of having case at age M in year t)
      // For this, we add up all children who could have been M months old during year t
      // i.e. children who were born in months 12-M+1 to 12 of the previous year plus children
      // born in months 1 to 12-M of the current year
      denom_st_bymonth = 0;
      for (c in 1:12) {
        
        if (c>12-M) {
          
          // Children born in months 12-M+1 to 12 of the previous year were of age M during year t
            
          // Special case as I don't have n from 1999, so just 2000
          if ((st-1) % ny ==0) {
        
            denom_st_bymonth += n_st_bymonth[st,c];
          
          } else {
              
            denom_st_bymonth += n_st_bymonth[st-1,c];
              
          }
            
        } else {
          // Children who were born in months 1 to 12-M were of age M during year t
          denom_st_bymonth += n_st_bymonth[st,c];
            
        }
        
      }
      
      y_st_bymonth[st,M+1] ~ poisson(denom_st_bymonth * p_inf * (p_rep_mild + p_sev * (p_rep_sev - p_rep_mild)));
      sev_st_bymonth[st,M+1] ~ poisson(denom_st_bymonth * p_sev * p_rep_sev);
        
    }
    
  }
  
}
