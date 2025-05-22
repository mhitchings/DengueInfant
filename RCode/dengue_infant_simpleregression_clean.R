require(dplyr)
require(tidyverse)
require(ggplot2)
require(grid)
require(gridExtra)
require(lubridate)
require(broom)
require(broom.helpers)
require(broom.mixed)
require(lme4)
require(rms)
require(rstan)
require(arrow)
require(mgcv)
require(lme4)

state_codes = read.table('./Data/01-processedData/mapping/state.txt')
colnames(state_codes) = c("letter","number")

all_den = read.csv('./Data/02-analysisData/DengueCasesByStateYear_2000to2024_No1Day.csv')

all_den = all_den %>% arrange(state,year_not)

statenums = unique(all_den$state)

all_den = all_den %>% mutate(agediff = mean6pluscaseage-meanpopage,
                             logir = log(ncase/pop*100000),
                             loginfantir = log(infantcases/infantpop*100000),
                             w=sqrt(infantcases)/sum(sqrt(all_den$infantcases)),
                             loginfantsevir = log(infantseverecases/infantpop*100000),
                             wsev = sqrt(infantseverecases)/sum(sqrt(all_den$infantseverecases)),
                             w1to5 = sqrt(cases1to5)/sum(sqrt(all_den$cases1to5)))

### Linear regression of infant case age in each state-year against state-year-level variables

m_meanage_linear = lm(meaninfantcaseage ~ year_not + mean6pluscaseage + meanpopage + logir + meanmaternalage + proppos + statecode,data=all_den,
                      weights = w)
tab = broom.helpers::tidy_and_attach(m_meanage_linear)
tab = tab %>% mutate(estconfi = paste0(sprintf('%.2f',estimate),' (',sprintf('%.2f',conf.low),',',sprintf('%.2f',conf.high),')'),
                     pval = sprintf('%.2f',p.value)) %>%
  dplyr::select(term,estconfi,pval)

# Infant incidence rate (log) as outcome
m_infantir_linear = lm(loginfantir ~ year_not + mean6pluscaseage + meanpopage + logir + meanmaternalage + proppos + statecode,data=all_den %>% filter(loginfantir>-Inf),weights = w)
tab = broom.helpers::tidy_and_attach(m_infantir_linear)
tab = tab %>% mutate(estconfi = paste0(sprintf('%.2f',estimate),' (',sprintf('%.2f',conf.low),',',sprintf('%.2f',conf.high),')'),
                     pval = sprintf('%.3f',p.value)) %>%
  dplyr::select(term,estconfi,pval)

# Now for age of severe infant cases
m_meansevage = lm(meaninfantsevcaseage ~ year_not + mean6pluscaseage + meanpopage + logir + meanmaternalage + proppos + statecode,data=all_den,weights = wsev)
tab = broom.helpers::tidy_and_attach(m_meansevage)
tab = tab %>% mutate(estconfi = paste0(sprintf('%.2f',estimate),' (',sprintf('%.2f',conf.low),',',sprintf('%.2f',conf.high),')'),
                     pval = sprintf('%.2f',p.value)) %>%
  dplyr::select(term,estconfi,pval)

# Same but for children 1 to 5 (should be no effect of maternal SP)
m_meanage1to5_linear = glm(mean1to5caseage ~ year_not + mean6pluscaseage + meanpopage + logir + meanmaternalage + proppos + statecode,data=all_den,weights = w1to5)
tab = broom.helpers::tidy_and_attach(m_meanage1to5_linear)
tab = tab %>% mutate(estconfi = paste0(sprintf('%.2f',estimate),' (',sprintf('%.2f',conf.low),',',sprintf('%.2f',conf.high),')'),
                     pval = sprintf('%.2f',p.value)) %>%
  dplyr::select(term,estconfi,pval)


