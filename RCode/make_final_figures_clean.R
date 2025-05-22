require(dplyr)
require(tidyverse)
require(ggplot2); theme_set(theme_classic())
require(grid)
require(gridExtra)
require(lubridate)
require(broom)
require(broom.helpers)
require(lme4)
require(rms)
require(rstan)
require(arrow)
require(geodata)
require(sf)
require(scales)
require(ggh4x)
require(cowplot)

setwd("C:/Users/mhitchings/Dropbox (UFL)/DengueInfants")
state_codes = read.table('./Data/01-processedData/mapping/state.txt')
colnames(state_codes) = c("letter","number")

den_cases = readRDS('./Data/02-analysisData/DengueLineList_Cases.rds')

all_den = read.csv('./Data/02-analysisData/DengueCasesByStateYear_2000to2024_No1Day.csv')

all_den = all_den %>% arrange(state,year_not)

statenums = unique(all_den$state)

ncaseunder18 = den_cases %>% 
  filter(!is.na(state) & age_years<=18 & year_not>=2000 & year_not<=2024) %>% 
  group_by(year_not,state) %>% 
  summarise(ncaseunder18 = n())
all_den = all_den %>% left_join(ncaseunder18,
                                by=c('year_not','state'))

nseverecaseunder18 = den_cases %>% 
  filter(!is.na(state) & age_years<=18 & year_not>=2000 & year_not<=2024 & severe=="Severe") %>% 
  group_by(year_not,state) %>% 
  summarise(nseverecaseunder18 = n())
all_den = all_den %>% left_join(nseverecaseunder18,
                                by=c('year_not','state'))


mean6pluscaseage = den_cases %>% 
  filter(!is.na(state) & age_years>=6 & year_not>=2000 & year_not<=2024) %>% 
  group_by(year_not,state) %>% 
  summarise(mean6pluscaseage = mean(age_years))
all_den = all_den %>% left_join(mean6pluscaseage,
                                by=c('year_not','state'))

all_den = all_den %>% mutate(agediff = meanage-meanpopage,
                             logir = log(ir),
                             loginfantir = log(infantcases/infantpop*100000))

### Heat map
state_codes = state_codes %>% mutate(region=case_when(
  is.na(letter) ~ "None",
  letter %in% c("AC","AM","AP","PA","RO","RR","TO") ~ "North",
  letter %in% c("AL","BA","CE","MA","PB","PE","PI","RN","SE") ~ "Northeast",
  letter %in% c("GO","MT","MS","DF") ~ "Central-West",
  letter %in% c("ES","MG","RJ","SP") ~ "Southeast",
  letter %in% c("PR","RS","SC") ~ "South"
  
)) %>% mutate(region=factor(region,levels=c('South','Southeast','Central-West','Northeast','North')),
              region_num = case_when(region == "North" ~ 1,
                                     region == "Northeast" ~ 2,
                                     region == "Central-West" ~ 3,
                                     region == "Southeast" ~ 4,
                                     region== "South" ~ 5)) %>% 
  arrange(region)
state_codes$staterank=1:nrow(state_codes)

all_den = all_den %>% left_join(state_codes %>% dplyr::select(letter,staterank,region_num) %>% dplyr::rename(statecode=letter),by='statecode')

lastyear = 2024

### Plots by region
all_den_byregion = all_den %>% group_by(year_not,region_num) %>%
  summarise(regioninfantcases = sum(infantcases),
            regioninfantpop = sum(infantpop),
            regioninfantseverecases = sum(infantseverecases),
            regionpropinfcaseneonate = sum(infantcases_0month)/sum(infantcases),
            regionpropinfsevcasegt4 = 1 - sum(infantseverecases_0month+infantseverecases_1month+infantseverecases_2month+
                                                infantseverecases_3month+infantseverecases_4month)/sum(infantseverecases),
            regionmeaninfantcaseage = sum(meaninfantcaseage*infantcases,na.rm=T)/sum(infantcases),
            regionmeaninfantsevcaseage = sum(meaninfantsevcaseage*infantseverecases,na.rm=T)/sum(infantseverecases))
all_den_byregion_y = all_den %>% mutate(y_label = case_when(year_not <=2009 ~ "2000-2010",
                                                             year_not>2009 & year_not<=2019 ~ "2010-2019",
                                                             year_not>=2020 ~ "2020-2024"
                                       )) %>% group_by(y_label,region_num) %>%
  summarise(regioninfantcases = sum(infantcases),
            regioninfantpop = sum(infantpop),
            regioninfantseverecases = sum(infantseverecases),
            regionpropinfcaseneonate = sum(infantcases_0month)/sum(infantcases),
            regionpropinfsevcasegt4 = 1 - sum(infantseverecases_0month+infantseverecases_1month+infantseverecases_2month+
                                                infantseverecases_3month+infantseverecases_4month)/sum(infantseverecases),
            regionmedianpropinfcaseneonate = median(1-propinfcasegt1,na.rm=T),
            regionmedianpropinfsevcasegt4 = median(propinfsevcasegt4,na.rm=T)
  ) %>% 
  mutate(ir = regioninfantcases/regioninfantpop,
         sevir = regioninfantseverecases/regioninfantpop)

all_den_brazil = all_den %>% group_by(year_not) %>%
  summarise(brazilinfantcases = sum(infantcases),
            brazilinfantpop = sum(infantpop),
            brazilinfantseverecases = sum(infantseverecases),
            propinfcaseneonate = sum(infantcases_0month)/sum(infantcases),
            propinfsevcasegt4 = 1 - sum(infantseverecases_0month+infantseverecases_1month+infantseverecases_2month+
                                                infantseverecases_3month+infantseverecases_4month)/sum(infantseverecases),
            brazilmeaninfantcaseage = sum(meaninfantcaseage*infantcases,na.rm=T)/sum(infantcases),
            brazilmeaninfantsevcaseage = sum(meaninfantsevcaseage*infantseverecases,na.rm=T)/sum(infantseverecases))


textsize = 10

all_den_brazil_5y = all_den %>% mutate(y5 = floor(year_not/5),
                                       y5_label = case_when(y5 == 400 ~ "2000-2004",
                                                            y5 == 401 ~ "2005-2009",
                                                            y5 == 402 ~ "2010-2014",
                                                            y5 == 403 ~ "2015-2019",
                                                            y5 == 404 ~ "2020-2024"
                                       )) %>% group_by(y5) %>%
  summarise(brazilinfantcases = sum(infantcases),
            brazilinfantpop = sum(infantpop),
            brazilinfantseverecases = sum(infantseverecases),
            brazilpropinfcaseneonate = sum(infantcases_0month)/sum(infantcases),
            brazilpropinfsevcasegt4 = 1 - sum(infantseverecases_0month+infantseverecases_1month+infantseverecases_2month+
                                                infantseverecases_3month+infantseverecases_4month)/sum(infantseverecases)
  ) %>% 
  mutate(ir = brazilinfantcases/brazilinfantpop,
         sevir = brazilinfantseverecases/brazilinfantpop)


all_den_brazil %>% filter(year_not %in% c(2000,2024)) %>% dplyr::select(year_not,brazilinfantcases,brazilinfantpop,brazilinfantseverecases) %>% 
  pivot_wider(names_from='year_not',values_from='brazilinfantcases':'brazilinfantseverecases') %>% 
  mutate(ir2000=brazilinfantcases_2000/brazilinfantpop_2000,
         ir2024=brazilinfantcases_2024/brazilinfantpop_2024,
         irsev2000=brazilinfantseverecases_2000/brazilinfantpop_2000,
         irsev2024=brazilinfantseverecases_2024/brazilinfantpop_2024,
         foldburden=ir2024/ir2000,
         foldsevburden=irsev2024/irsev2000)

all_den_byy = all_den %>% 
  mutate(y_label = case_when(year_not <=2009 ~ "2000-2010",
                             year_not>2009 & year_not<=2019 ~ "2010-2019",
                             year_not>=2020 ~ "2020-2024"
         )) %>% group_by(y_label,statecode) %>%
  summarise(infantcases = sum(infantcases),
            infantpop = sum(infantpop),
            infantseverecases = sum(infantseverecases),
            ypropinfcaseneonate = sum(infantcases_0month)/sum(infantcases),
            ypropinfsevcasegt4 = 1 - sum(infantseverecases_0month+infantseverecases_1month+infantseverecases_2month+
                                                infantseverecases_3month+infantseverecases_4month)/sum(infantseverecases),
            propmedianinfcaseneonate = median(1-propinfcasegt1,na.rm=T),
            propmedianinfsevcasegt4 = median(propinfsevcasegt4,na.rm=T),
            region_num=mean(region_num))
all_den_byy = all_den_byy %>% 
  mutate(ir = infantcases/infantpop,
         sevir = infantseverecases/infantpop)

# Figure 1 will be a line graph of incidence rate, line graph of case age, and histogram
# An inset map will give the legend

# Histogram of case age

### Age distribution by region
dailycountsbyregion = matrix(NA,nrow=5,ncol=48)
dailysevcountsbyregion = matrix(NA,nrow=5,ncol=48)

regions = c("North","Northeast","Central-West","Southeast","South")

for (row in 1:5) {
  
  x=hist((den_cases %>% filter(!is.na(state) & age_years<1 & age_months<12 & region==regions[row]))$age_days,
         breaks=c(seq(1,29,7),seq(59,360,30)))
  counts = x$counts
  dailycounts = c(counts[1:4],
                  rep(counts[5:length(counts)]/4,each=4)
  )
  dailycountsbyregion[row,] = dailycounts
  
  x=hist((den_cases %>% filter(!is.na(state) & age_years<1 & age_months<12 & region==regions[row] & severe=="Severe"))$age_days,
         breaks=c(seq(1,29,7),seq(59,360,30)))
  counts = x$counts
  dailysevcounts = c(counts[1:4],
                     rep(counts[5:length(counts)]/4,each=4)
  )
  dailysevcountsbyregion[row,]=dailysevcounts
}
df_forplot = data.frame('Age'=rep(seq(1,52,length.out=48),2*5),
                        'Region'=rep(rep(regions,each=48),2),
                        'Counts'=c(c(t(dailycountsbyregion)),c(t(dailysevcountsbyregion))),
                        'Type'=c(rep('All Cases',length(dailycounts)*5),
                                 rep('Severe Cases',length(dailycounts)*5))
)

cc = c('lightslategrey','seagreen','darkorchid2','lightsalmon','navy')

map_legend = ggplot() + 
  geom_sf(data = brazil_regions,aes(fill=factor(region_num)),show.legend=F) + 
  geom_sf(data = brazil,linetype=1,fill=NA,linewidth=0.5,colour='black')+
  scale_fill_manual(name='Region',values=cc) + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA,linewidth=1),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

# Figure 1A: annual incidence rate in infants by region, with panels for all cases and severe cases
figure1_a =  
  ggplot() + 
  geom_line(data = all_den_byregion %>% mutate(infantsevir = regioninfantseverecases/regioninfantpop*100000,
                                               infantir = regioninfantcases/regioninfantpop*100000) %>% 
              pivot_longer(cols=c('infantir','infantsevir')) %>% 
              mutate(label = case_when(name=="infantir" ~ "All cases",
                                       name=="infantsevir" ~ "Severe cases")),
            aes(x=year_not,y=value,color=factor(region_num)),alpha=0.5,show.legend = F) + 
  geom_line(data = all_den_brazil %>% mutate(infantsevir = brazilinfantseverecases/brazilinfantpop*100000,
                                             infantir = brazilinfantcases/brazilinfantpop*100000) %>% 
              pivot_longer(cols=c('infantir','infantsevir')) %>% 
              mutate(label = case_when(name=="infantir" ~ "All cases",
                                       name=="infantsevir" ~ "Severe cases")),
            aes(x=year_not,y=value),color='black',linewidth=0.75,linetype=1,show.legend = F) + 
  facet_wrap(vars(label),nrow=2,scales='free_y') + 
  scale_x_continuous(name='Year',breaks=c(2000,2010,2020)) + 
  ylab('Infant IR per 100,000') + 
  ggtitle('A')+
  scale_color_manual(name='Region',values=cc)+
  theme(legend.position='bottom',
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        strip.text = element_text(size=textsize,hjust=0))

figure1_a_with_inset = 
  ggdraw() + 
  draw_plot(figure1_a) + 
  draw_plot(map_legend,x=0.3,y=0.6,width=0.3,height=0.3)

### Figure 1B is empirical age distribution of all cases and severe cases in Brazil,
### aggregated over all states and over 2000-2024, as well as the age distribution
### of severe cases from O'Driscoll et al (Thailand)
thailand_sevcases = read.csv('~/UF/Research/DengueInfants/Manuscripts/InfantCases.csv')
thailand_sevcases$dens = thailand_sevcases$cases/sum(thailand_sevcases$cases)
thailand_sevcases$age_months = thailand_sevcases$age_months * 30.75
thailand_sevcases$label='Severe Cases'

df_forplot_thailand = data.frame('Age'=seq(1,52,length.out=52),
                                 'Counts'=rep(thailand_sevcases$dens/4,each=4),
                                 'label'=rep('Severe Cases',52))

## Add in the number of cases in 1-year-olds as a separate point
ncase_1to2 = den_cases %>% 
  filter(!is.na(state) & age_years==1) %>% count()
nsevcase_1to2 = den_cases %>% 
  filter(!is.na(state) & severe=="Severe" & age_years==1) %>% count()

df_forplot_brazil = data.frame('Age'=rep(c(seq(1,52,length.out=48),65),2),
                               'Counts'=c(colSums(dailycountsbyregion),ncase_1to2$n/52,
                                          colSums(dailysevcountsbyregion),nsevcase_1to2$n/52),
                               'Type'=c(rep('All Cases',length(dailycounts)+1),
                                        rep('Severe Cases',length(dailycounts)+1))
)

# Normalize counts to get empirical density
df_forplot_brazil = df_forplot_brazil %>% group_by(Type) %>%
  mutate(Dens=Counts/sum(Counts))

figure1_b = ggplot() + 
  geom_col(data=df_forplot_thailand,aes(x=Age,y=Counts),fill='grey',colour='grey') + 
  geom_point(data=df_forplot_brazil %>% filter(Age==65) %>% rename(label=Type),
             aes(x=Age,y=Dens),colour='black') + 
  geom_segment(data=df_forplot_brazil %>% filter(Age==65) %>% rename(label=Type),
               aes(x=Age,xend=Age,y=Dens,yend=0),colour='black') + 
  geom_area(data = df_forplot_brazil %>% filter(Age<65) %>% 
              rename(label=Type),aes(x=Age,y=Dens),colour='black',
            fill='grey',alpha=0.5,
            linewidth=0.5,linetype=1) + 
  ylab('Proportion')+
  scale_x_continuous(name='Age (months)',
                     breaks=c(c(0,3,6,9,12)*4*1.085,65),
                     labels=c(0,3,6,9,12,'12-24'),
                     limits=c(0,70)) +
  ggtitle('B')+
  facet_wrap(vars(label),nrow=2,
             ,scales='free_y'
  )+
  theme(legend.position='bottom',
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        strip.text = element_text(size=textsize,hjust=0))



## Figure 1C is the proportion of cases and severe cases that fall in the upper "mode"
## aggregated across three time periods (2000-2009, 2010-2019, and 2020-2024). States
## will be as transparent points, regions will be as solid lines
all_den_combined_byy = bind_rows(all_den_byy  %>% 
  pivot_longer(cols=c('propmedianinfcaseneonate','propmedianinfsevcasegt4')) %>% 
  mutate(label = case_when(name=="propmedianinfcaseneonate" ~ "Cases among neonates",
                           name=="propmedianinfsevcasegt4" ~ "Sev. Cases among 5-12 mo.")),
  all_den_byregion_y  %>% mutate(statecode='All') %>%
    pivot_longer(cols=c('regionmedianpropinfcaseneonate','regionmedianpropinfsevcasegt4')) %>% 
    mutate(label = case_when(name=="regionmedianpropinfcaseneonate" ~ "Cases among neonates",
                             name=="regionmedianpropinfsevcasegt4" ~ "Sev. Cases among 5-12 mo."))
) %>% mutate(y_numeric = case_when(y_label=="2000-2010" ~ 1,
                                   y_label=="2010-2019" ~ 2,
                                   y_label=="2020-2024" ~ 3))

figure1_c <- ggplot() + 
  geom_point(data = all_den_combined_byy %>% filter(statecode != "All"),
             aes(x=y_numeric,y=value,
                 fill=factor(region_num),
                 colour=factor(region_num)),
             alpha=0.2,
             show.legend=F)+
  geom_line(data = all_den_combined_byy %>% filter(statecode=="All"),
            aes(x=y_numeric,y=value,
                color=factor(region_num)),
            show.legend=F)  + 
  scale_x_continuous(name='Period',breaks=c(1:3),
                     labels=c('00-09','10-19','20-24'),
                     limits=c(0.75,3.25)) + 
  ylab('Proportion') +
  ggtitle('C')+
  scale_color_manual(name='Region',values=cc)+
  scale_fill_manual(name='Region',values=cc)+
  facet_wrap(vars(label),nrow=2)+
  theme(legend.position='bottom',
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        strip.text = element_text(size=textsize,hjust=0))

ggsave('~/UF/Research/DengueInfants/Manuscripts/Figure1_Landscape_0507.png',
       grid.arrange(figure1_a_with_inset,figure1_b,figure1_c,ncol=3),
       height=4,width=7.5,units='in',device='png')



######################################################################################
######################################################################################
######################################################################################
######################################################################################

### Figure S1: annual incidence rate by state in the whole population
ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS1_0512.png',
       all_den %>% 
         mutate(statecode = factor(statecode,levels=rev(state_codes$letter)),
                region = factor(region,levels=unique(rev(state_codes$region)))) %>% 
         ggplot() + 
         geom_line(aes(x=year_not,y=ir,colour=region),linewidth=1) + 
         facet_wrap(vars(statecode)) + 
         scale_color_manual(name='Region',values=cc)+
         ylab('Annual incidence rate per 100,000') +
         scale_x_continuous(name='Year',breaks=seq(2000,2025,5),labels=c('00','05','10','15','20','25')),
       height=6,width=8,units='in',device='png')

# Fit exponential growth in IR to each state
df = data.frame('State'=unique(all_den$statecode),
                'Slope'=rep(NA,27),
                'p'=rep(NA,27),
                'R2'=rep(NA,27))
row = 1
for (sc in unique(all_den$statecode)) {
  m_exp_growth = lm(loginfantir ~ year_not,data=all_den %>% filter(statecode==sc & loginfantir>-Inf))
  df$Slope[row] = exp(coef(m_exp_growth)[2])
  df$p[row] = summary(m_exp_growth)$coef[2,4]
  df$R2[row] = summary(m_exp_growth)$r.squared
  row = row+1
  
}

df %>% left_join(state_codes %>% dplyr::select(letter,staterank,region_num,region) %>% dplyr::rename(State=letter),by='State') %>%
  arrange(Slope)

brazil_states = st_as_sf(readRDS('./Data/SpatialData/gadm/gadm41_BRA_1_pk.rds')) %>% 
  mutate(letter = substr(HASC_1,4,5)) %>% 
  left_join(state_codes,by='letter')
brazil_regions = brazil_states %>% group_by(region_num) %>% summarise()
brazil = brazil_regions %>% group_by() %>% summarise()

## Plot of exponential growth by state
figs2 = ggplot() + 
  geom_sf(data = brazil_states %>% rename(State=letter) %>% left_join(df,by='State'),aes(fill=Slope)) + 
  scale_fill_continuous(name='Annual fold increase') + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6))

ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS2_ExpGrowth.png',
       figs2,
       height=4,width=4,units='in',device='png')

# Proportion of cases and severe cases that are in infants by region and time
all_den_byregion = all_den %>% group_by(year_not,region_num) %>%
  summarise(regioninfantcases = sum(infantcases,na.rm=T),
            regioninfantpop = sum(infantpop,na.rm=T),
            regioninfantseverecases = sum(infantseverecases,na.rm=T),
            regiontotalcasesunder18 = sum(ncaseunder18,na.rm=T),
            regiontotalseverecasesunder18 = sum(nseverecaseunder18,na.rm=T),
            regiontotalcases = sum(ncase),
            regiontotalseverecases = sum(nsevere))

all_den_brazil = all_den %>% group_by(year_not) %>%
  summarise(brazilinfantcases = sum(infantcases,na.rm=T),
            brazilinfantpop = sum(infantpop,na.rm=T),
            brazilinfantseverecases = sum(infantseverecases,na.rm=T),
            braziltotalcasesunder18 = sum(ncaseunder18,na.rm=T),
            braziltotalseverecasesunder18 = sum(nseverecaseunder18,na.rm=T),
            braziltotalcases = sum(ncase),
            braziltotalseverecases = sum(nsevere))

all_den_combined = bind_rows(all_den_byregion,
                             all_den_brazil %>% mutate(region_num=6) %>% 
                               rename(regioninfantcases=brazilinfantcases,
                                      regioninfantseverecases=brazilinfantseverecases,
                                      regiontotalcasesunder18=braziltotalcasesunder18,
                                      regiontotalseverecasesunder18=braziltotalseverecasesunder18))  %>% 
  mutate(caseratio = regioninfantcases/regiontotalcasesunder18,
         severecaseratio = regioninfantseverecases/regiontotalseverecasesunder18) %>%
  pivot_longer(cols=c('caseratio','severecaseratio')) %>% 
  mutate(label = case_when(name=="caseratio" ~ "All Cases",
                           name=="severecaseratio" ~ "Severe Cases"))


figure_propsev =  
  ggplot() + 
  geom_line(data = all_den_combined,
            aes(x=year_not,y=value,color=factor(region_num),alpha=factor(as.numeric(region_num==6))),linewidth=1) + 
  facet_wrap(vars(label),ncol=2) + 
  scale_x_continuous(name='Year',breaks=seq(2000,2020,5)) + 
  ylab('Proportion of child cases in infants') + 
  scale_alpha_manual(values=c(0.5,1))+
  scale_color_manual(name='Region',values=c(cc,'black'),labels=c('North','Northeast','Central-West','Southeast','South','Brazil'))+
  guides(alpha='none')+
  theme(legend.position='bottom',
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        strip.text = element_text(size=textsize,hjust=0))

ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS3_PropChildCases_0512.png',
       figure_propsev,
       height=4,width=6.6,units='in',device='png')


#### Supplementary figure of age distribution when we have date of birth
### Age distribution by region
x=hist((den_cases %>% filter(!is.na(dob) & !is.na(state) & age_days_calc_symp<=749 & age_days_calc_symp>0))$age_days_calc_symp,
       breaks=c(seq(1,29,7),seq(59,749,30)))
counts = x$counts
dailycounts = c(counts[1:4],
                rep(counts[5:length(counts)]/4,each=4)
)

x=hist((den_cases %>% filter(!is.na(dob) & !is.na(state) & age_days_calc_symp<=749 & age_days_calc_symp>0 & severe=="Severe"))$age_days_calc_symp,
       breaks=c(seq(1,29,7),seq(59,749,30)))
counts = x$counts
dailysevcounts = c(counts[1:4],
                   rep(counts[5:length(counts)]/4,each=4)
)


df_forplot = data.frame('Age'=rep(seq(1,104,length.out=100),2),
                        'Count'=c(dailycounts,dailysevcounts),
                        'Type'=c(rep('All Cases',length(dailycounts)),
                                 rep('Severe Cases',length(dailycounts)))
)
df_forplot = df_forplot %>% group_by(Type) %>%
  mutate(Dens = Count/sum(Count))

figures_agedist = ggplot() + 
  geom_area(data = df_forplot %>% 
              rename(label=Type),aes(x=Age,y=Dens),colour='black',
            fill='grey',alpha=0.5,
            linewidth=0.5,linetype=1) + 
  ylab('Proportion')+
  scale_x_continuous(name='Age (months)',
                     breaks=c(c(0,6,12,18,24)*4*1.085),
                     labels=c(0,6,12,18,24)) +
  facet_wrap(vars(label),nrow=2,
             ,scales='free_y'
  )+
  theme(legend.position='bottom',
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        strip.text = element_text(size=textsize,hjust=0))
ggsave('~/UF/Research/DengueInfants/Manuscripts/FigureS4_AgeDistribution_WithDOB.png',
       figures_agedist,
       height=4,width=2.5,units='in',device='png')

# Figure S5-9: Age distribution of severe cases by state and 5-year period in SE, NE, and CW
ggsave('~/UF/Research/DengueInfants/Manuscripts/Figure_S5_0512.png',
       den_cases  %>% 
         filter(!is.na(state) & age_years<1 & age_months<12 & severe=="Severe" & region=="North") %>% 
         mutate(y5 = floor(year_not/5),
                y5_label = case_when(y5 == 400 ~ "2000-2004",
                                     y5 == 401 ~ "2005-2009",
                                     y5 == 402 ~ "2010-2014",
                                     y5 == 403 ~ "2015-2019",
                                     y5 == 404 ~ "2020-2024"
                ),
                region = factor(region,levels = rev(unique(state_codes$region)))) %>% 
         ggplot(aes(x=age_months)) + geom_histogram(breaks=(-1):12) + 
         facet_grid(statecode~y5_label) + 
         xlab("Age (months)") + ylab("Case count"),
       height=8,width=6,units='in',device='png')

ggsave('~/UF/Research/DengueInfants/Manuscripts/Figure_S6_0512.png',
       den_cases  %>% 
         filter(!is.na(state) & age_years<1 & age_months<12 & severe=="Severe" & region=="Northeast") %>% 
         mutate(y5 = floor(year_not/5),
                y5_label = case_when(y5 == 400 ~ "2000-2004",
                                     y5 == 401 ~ "2005-2009",
                                     y5 == 402 ~ "2010-2014",
                                     y5 == 403 ~ "2015-2019",
                                     y5 == 404 ~ "2020-2024"
                ),
                region = factor(region,levels = rev(unique(state_codes$region)))) %>% 
         ggplot(aes(x=age_months)) + geom_histogram(breaks=(-1):12) + 
         facet_grid(statecode~y5_label) + 
         xlab("Age (months)") + ylab("Case count"),
       height=8,width=6,units='in',device='png')

ggsave('~/UF/Research/DengueInfants/Manuscripts/Figure_S7_0512.png',
       den_cases  %>% 
         filter(!is.na(state) & age_years<1 & age_months<12 & severe=="Severe" & region=="Central-West") %>% 
         mutate(y5 = floor(year_not/5),
                y5_label = case_when(y5 == 400 ~ "2000-2004",
                                     y5 == 401 ~ "2005-2009",
                                     y5 == 402 ~ "2010-2014",
                                     y5 == 403 ~ "2015-2019",
                                     y5 == 404 ~ "2020-2024"
                ),
                region = factor(region,levels = rev(unique(state_codes$region)))) %>% 
         ggplot(aes(x=age_months)) + geom_histogram(breaks=(-1):12) + 
         facet_grid(statecode~y5_label) + 
         xlab("Age (months)") + ylab("Case count"),
       height=8,width=6,units='in',device='png')

ggsave('~/UF/Research/DengueInfants/Manuscripts/Figure_S8_0512.png',
       den_cases  %>% 
         filter(!is.na(state) & age_years<1 & age_months<12 & severe=="Severe" & region=="Southeast") %>% 
         mutate(y5 = floor(year_not/5),
                y5_label = case_when(y5 == 400 ~ "2000-2004",
                                     y5 == 401 ~ "2005-2009",
                                     y5 == 402 ~ "2010-2014",
                                     y5 == 403 ~ "2015-2019",
                                     y5 == 404 ~ "2020-2024"
                ),
                region = factor(region,levels = rev(unique(state_codes$region)))) %>% 
         ggplot(aes(x=age_months)) + geom_histogram(breaks=(-1):12) + 
         facet_grid(statecode~y5_label) + 
         xlab("Age (months)") + ylab("Case count"),
       height=8,width=6,units='in',device='png')

ggsave('~/UF/Research/DengueInfants/Manuscripts/Figure_S9_0512.png',
       den_cases  %>% 
         filter(!is.na(state) & age_years<1 & age_months<12 & severe=="Severe" & region=="South") %>% 
         mutate(y5 = floor(year_not/5),
                y5_label = case_when(y5 == 400 ~ "2000-2004",
                                     y5 == 401 ~ "2005-2009",
                                     y5 == 402 ~ "2010-2014",
                                     y5 == 403 ~ "2015-2019",
                                     y5 == 404 ~ "2020-2024"
                ),
                region = factor(region,levels = rev(unique(state_codes$region)))) %>% 
         ggplot(aes(x=age_months)) + geom_histogram(breaks=(-1):12) + 
         facet_grid(statecode~y5_label) + 
         xlab("Age (months)") + ylab("Case count"),
       height=8,width=6,units='in',device='png')


## Heat map of maternal seroprevalence
all_den = all_den %>% filter(year_not<=2014)

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
        group_by(year) %>% summarise(Smult=median(Smult),SP=median(SP),
                                     S1=median(S1),S1_I = median(S1_I),
                                     S2=median(S2),S2_I = median(S2_I),
                                     S3=median(S3),S3_I=median(S3_I),
                                     S4=median(S4),S4_I=median(S4_I))
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
        group_by(year) %>% summarise(FOI = median(hazard))
      
      foi_dat[,dimnames(foi_dat)[[2]]==s] = t$FOI
      
    }
    
  }
  
} else if (grepl('stanfit_v3',foi_dir)) {
  
  for (s in statenums) {
    
    if (file.exists(paste0(foi_dir,'S_mother/',s,'.csv'))) {
      
      t = read.csv(paste0(foi_dir,'S_mother/',s,'.csv'))
      
      t = t %>% filter(year>=2000 & year<=2014) %>% mutate(SP = Smult+S1) %>%
        group_by(year) %>% summarise(Smult=median(Smult),SP=median(SP),S1=median(S1))
      
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
        group_by(year) %>% summarise(Smult=median(Smult),SP=median(SP),S1=median(S1))
      
      maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("SP")] = t$SP
      maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("Smult")] = t$Smult
      maternal_sp_dat[,dimnames(maternal_sp_dat)[[2]]==s,c("S1")] = t$S1
      
    }
    
    if (file.exists(paste0(foi_dir,'hazard_infant/',s,'.csv'))) {
      
      t = read.csv(paste0(foi_dir,'hazard_infant/',s,'.csv'))
      t = t %>% filter(year>=2000 & year<=2014) %>%
        group_by(year) %>% summarise(FOI = median(hazard))
      
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

lastyear=2014
ggsave('~/UF/Research/DengueInfants/Manuscripts/HeatMap_MaternalSeroprevalence.png',
       all_den %>% ggplot() + 
         geom_tile(aes(x=year_not,y=staterank,fill=maternal_sp)) + 
         geom_hline(yintercept = 3.5)+
         geom_hline(yintercept = 7.5)+
         geom_hline(yintercept = 11.5)+
         geom_hline(yintercept = 20.5)+
         annotate(geom="text",x=lastyear+1,y=2,label="South",angle=270)+
         annotate(geom="text",x=lastyear+1,y=5.5,label="Southeast",angle=270)+
         annotate(geom="text",x=lastyear+1,y=9.5,label="Central-West",angle=270)+
         annotate(geom="text",x=lastyear+1,y=16,label="Northeast",angle=270)+
         annotate(geom="text",x=lastyear+1,y=24,label="North",angle=270)+
         scale_y_continuous(name='State',
                            labels=state_codes$letter,
                            breaks=state_codes$staterank) + 
         scale_x_continuous(name='Year',
                            labels=seq(2000,lastyear,by=5),
                            breaks=seq(2000,lastyear,by=5)) + 
         scale_fill_gradient(name='Maternal seroprevalence',
                             low="white", high="blue") + 
         theme(legend.position = 'bottom'),
       device='png',height=8,width=6,units='in')

