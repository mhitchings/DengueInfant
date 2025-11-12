//
functions {
  
  real expit(real x) {
    
    return 1/(1+exp(-x));
  }
  
  real ade_func(real y, real mode, real kappa, real logOR, real a, real c) {

    return(logOR * ((y - a) / (mode * (c - a))) ^ (mode * kappa)  * ((c - y) / ((1 - mode) * (c - a))) ^ ((1 - mode) * kappa));
  }
    
  real bvnormmix(real x, real lambda1, real lambda2, real mu1, real mu2, real sigma1, real sigma2) {
    
    return log(lambda1 * exp(normal_lpdf(x | mu1, sigma1)) + lambda2 * exp(normal_lpdf(x | mu2, sigma2)));
    
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
  real lambda_prior[nst,6]; //Parameters to define bivariate normal mixture prior for FOI
  real lambda_prior_lbounds[nst]; 
  real lambda_prior_ubounds[nst]; 
  real msp_prior[nst,6]; //Parameters to define bivariate normal mixture prior for maternal seropositivity
  real msp_prior_lbounds[nst]; 
  real msp_prior_ubounds[nst]; 
  int<lower=0,upper=1> alpha_constant_t; //Whether reporting rate is constant over time
  int<lower=0,upper=1> alpha_linear; //Whether reporting rate is linear over time
  int<lower=0,upper=1> include_ade; //Whether or not to include ADE of severe disease
  int<lower=0,upper=1> include_sevagefrailty; //Whether or not to include increasing risk of severe disease with younger age
  int<lower=0,upper=1> include_infagefrailty; //Whether or not to include increasing risk of infection with younger age
  int<lower=0,upper=1> include_agereporting; //Whether or not to include increasing risk of reporting with younger age
  int<lower=0,upper=1> include_maternalprotection; //Whether or not to include maternal protection against infection
}

parameters {
  real<lower=0,upper=1> alphast_param[(alpha_linear + alpha_constant_t) ? ns*(2 * alpha_linear + alpha_constant_t):nst];
  
  real<lower=lambda_prior_lbounds,upper=lambda_prior_ubounds> lambda_st[nst]; // Dengue FOI in each state-year
  real<lower=msp_prior_lbounds,upper=msp_prior_ubounds> msp_st[nst]; // Maternal monotypic seroprevalence in each state-year
  
  // Profile of infant risk among those born to seronegative and seropositive mothers
  // For both groups, we estimate the hazard ratio in month 0 and month 1, and from month
  // two onwards the hazard ratio decays towards 1 at an exponential rate
  real<lower=0.1,upper=1> mat_protection_decay_param[include_maternalprotection ? 1:0];
  
  // Decreased risk in infants 1 month vs. general population
  real<lower=0,upper=1> baseline_hazmatprot_sp_param[include_maternalprotection ? 1:0];
  real<lower=0,upper=1> baseline_hazmatprot_sn_param[include_maternalprotection ? 1:0];
  
  // Increase risk in neonate vs. general population
  real<lower=0> baseline_hazmult_sp_param[include_infagefrailty ? 1:0];
  real<lower=0> baseline_hazmult_sn_param[include_infagefrailty ? 1:0];
  
  // Increased odds of reporting in neonates
  real<lower=0,upper=5> baseline_replogor_param[include_agereporting ? 1:0];
  
  // Changing risk of severe infection in nenoates
  real<lower=-3,upper=3> baseline_sevlogor_param[include_sevagefrailty ? 1:0];
  
  // Log odds of severe disease in a one-year old
  real<lower=-10,upper=0> oneyear_probsev;
  
  // Elevated risk of severe infection in children born to seropositive mothers - has the shape of a beta distribution, with a mode, peak and "concentration"
  real<lower=0,upper=3> adepeak_logor_param[include_ade ? 1:0];
  real<lower=0,upper=1> ademode_param[include_ade ? 1:0];
  real<lower=0,upper=10> adekappa_param[include_ade ? 1:0];
  
}

transformed parameters {
  
  // Maternal protection against infection
  real<lower=0,upper=1> alphast[nst];
  real<lower=0.1,upper=1> mat_protection_decay;
  real<lower=0,upper=1> baseline_hazmatprot_sp;
  real<upper=1> baseline_hazmatprot_sn;
  real<lower=0> baseline_hazmult_sp;
  real<lower=0> baseline_hazmult_sn;
  real<lower=0,upper=5> baseline_replogor;
  real<lower=-3,upper=3> baseline_sevlogor;
  real<lower=0,upper=3> adepeak_logor;
  real<lower=0,upper=1> ademode;
  real<lower=0,upper=10> adekappa;
  
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
    baseline_hazmult_sp = baseline_hazmult_sp_param[1];
    baseline_hazmult_sn = baseline_hazmult_sn_param[1];
  } else {
    baseline_hazmult_sp = 0;
    baseline_hazmult_sn = 0;
  }
  
  if (include_agereporting) {
    baseline_replogor = baseline_replogor_param[1];
  } else {
    baseline_replogor = 0;
  }
  
  if (include_sevagefrailty) {
    baseline_sevlogor = baseline_sevlogor_param[1];
  } else {
    baseline_sevlogor = 0;
  }
  
  if (include_ade) {
    adepeak_logor = adepeak_logor_param[1];
    ademode = ademode_param[1];
    adekappa = adekappa_param[1];
    
  } else {
    adepeak_logor = 0;
    ademode = 0.5;
    adekappa = 2;
  }
  
}

model {
  
  real hr_sp[nmonth];
  real hr_sn[nmonth];
  real cumhaz_sp;
  real cumhaz_sn;
  real p_inf_sp;
  real p_inf_sn;
  real p_sev_sp[nmonth];
  real p_sev_sn[nmonth];
  real p_sev;
  real p_inf;
  real p_rep;
  
  // Need to keep track of the children that are at risk of having a case at age 0 in calendar year t
  real denom_st_bymonth;
  
  // Some temporary variables to help with cumulative hazard calculation
  real n_byyear;
  real cumhaz_sn_stmc;
  real cumhaz_sp_stmc;
  
    // Priors
  for (st in 1:nst) {
    target += bvnormmix(lambda_st[st],lambda_prior[st,1],lambda_prior[st,2],lambda_prior[st,3],lambda_prior[st,4],lambda_prior[st,5],lambda_prior[st,6]) ; 
    target += bvnormmix(msp_st[st],msp_prior[st,1],msp_prior[st,2],msp_prior[st,3],msp_prior[st,4],msp_prior[st,5],msp_prior[st,6]) ; 
  }

  for (M in 0:(nmonth-1)) {
    
    if (M==0) {
    
      hr_sp[M+1] = (1 + baseline_hazmult_sp);
      hr_sn[M+1] = (1 + baseline_hazmult_sn);
      
      // probability of severe disease by month
      p_sev_sn[M+1] = expit(oneyear_probsev + baseline_sevlogor);
      p_sev_sp[M+1] = expit(oneyear_probsev + baseline_sevlogor);
      
    } else {
      
      hr_sp[M+1] = (1 - baseline_hazmatprot_sp * exp(-mat_protection_decay * (M-1)));
      hr_sn[M+1] = (1 - baseline_hazmatprot_sn * exp(-mat_protection_decay * (M-1)));
      
      // probability of severe disease by month
      p_sev_sn[M+1] = expit(oneyear_probsev);
      p_sev_sp[M+1] = expit(oneyear_probsev + ade_func(M,ademode,adekappa,adepeak_logor,1,nmonth-1));
      
    }

  }
  
  for (st in 1:nst) {
    
    // We will need the total children born in each calendar year
    n_byyear = 0;
      
    for (c in 1:12) {
      n_byyear = n_byyear + n_st_bymonth[st,c];
    }
    
    for (M in 0:(nmonth-1)) {

      // Build up cumulative hazard by state-year and month of age
      cumhaz_sn = 0;
      cumhaz_sp = 0;
      
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
        
          for (m in 0:(M-1)) {
          
            // If their being m months old occurred in the previous calendar year (i.e. if m<=12-c) use last year's
            // population FOI
            if (m<=12-c) {
              
              // Special case as I don't have FOI from 1999, so just 2000
              if ((st-1) % ny ==0) {
                
                cumhaz_sn_stmc += lambda_st[st]/12 * hr_sn[m+1];
                cumhaz_sp_stmc += lambda_st[st]/12 * hr_sp[m+1];
                
              } else {
                
                // Data is ordered by state first, then by year, so st-1 is the index of the same state's previous year
                cumhaz_sn_stmc += lambda_st[st-1]/12 * hr_sn[m+1];
                cumhaz_sp_stmc += lambda_st[st-1]/12 * hr_sp[m+1];
                
              }
              
            } else {
              // If their being m months old occurred in the same calendar year (i.e. if m>12-c) use this year's population FOI
              
              cumhaz_sn_stmc += lambda_st[st]/12 * hr_sn[m+1];
              cumhaz_sp_stmc += lambda_st[st]/12 * hr_sp[m+1];
              
            }
          
          }
          
        } else {
          
          for (m in 0:(M-1)) {
            
            cumhaz_sn_stmc += lambda_st[st]/12 * hr_sn[m+1];
            cumhaz_sp_stmc += lambda_st[st]/12 * hr_sp[m+1];
            
          }
          
        }
          
        // Finally, add this cumulative hazard weighted by the relative number of children born in that month
        cumhaz_sn += cumhaz_sn_stmc * n_st_bymonth[st,c]/n_byyear;
        cumhaz_sp += cumhaz_sp_stmc * n_st_bymonth[st,c]/n_byyear;
        
      }
      
      p_inf_sn = exp(-cumhaz_sn) * (1 - exp(-lambda_st[st]/12 * hr_sn[M+1]));
      p_inf_sp = exp(-cumhaz_sp) * (1 - exp(-lambda_st[st]/12 * hr_sp[M+1]));
        
      p_inf = msp_st[st] * p_inf_sp + (1 - msp_st[st]) * p_inf_sn;
      
      if (M==0) {
        p_rep = expit(logit(alphast[st]) + baseline_replogor);
      } else {
        p_rep = expit(logit(alphast[st]));
      }
      
      p_sev = (msp_st[st]*p_sev_sp[M+1]) * p_inf_sp + (1-msp_st[st])*p_sev_sn[M+1] * p_inf_sn;
      
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
      
      y_st_bymonth[st,M+1] ~ poisson(denom_st_bymonth * p_inf * p_rep);
      sev_st_bymonth[st,M+1] ~ poisson(denom_st_bymonth * p_sev * p_rep);
        
    }
    
  }
  
}

