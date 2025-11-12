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

state_codes = read.table('./Data/01-processedData/mapping/state.txt')
colnames(state_codes) = c("letter","number")

# Note that this file is not made publicly available as it contains date of birth
den_cases = try(readRDS('./Data/02-analysisData/DengueLineList_Cases.rds'))

all_den = read.csv('./Data/02-analysisData/DengueCasesByStateYear_2000to2024.csv')

all_den = all_den %>% arrange(state,year_not)

statenums = unique(all_den$state)

### Heat map
statenames = read.csv('./Data/SpatialData/municipality_match.csv')
statenames = statenames %>% distinct(State_iso,State) %>% rename(letter=State_iso,StateName=State)

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
  arrange(region) %>% left_join(statenames)
state_codes$staterank=1:nrow(state_codes)


all_den = all_den %>% left_join(state_codes %>% dplyr::select(letter,staterank,region_num,StateName) %>% dplyr::rename(statecode=letter),by='statecode')

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

brazil_states = st_as_sf(readRDS('./Data/SpatialData/gadm/gadm41_BRA_1_pk.rds')) %>% 
  mutate(letter = substr(HASC_1,4,5)) %>% 
  left_join(state_codes,by='letter')
brazil_regions = brazil_states %>% group_by(region_num) %>% summarise()
brazil = brazil_regions %>% group_by() %>% summarise()

# Figure 1 will be a line graph of incidence rate, line graph of case age, and histogram
# An inset map will give the legend

# Histogram of case age
### NOTE: If the line list of cases has not been read in, this code will use the case count
# data provided on github

regions = c("North","Northeast","Central-West","Southeast","South")

### Age distribution by region
if (!inherits(den_cases,'try-error')) {
  dailycountsbyregion = matrix(NA,nrow=5,ncol=48)
  dailysevcountsbyregion = matrix(NA,nrow=5,ncol=48)
  
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
} else {
  dailycountsbyregion = read.csv('./Data/02-analysisData/DailyCountsByRegion_ForFig1.csv')
  dailysevcountsbyregion = read.csv('./Data/02-analysisData/DailySevCountsByRegion_ForFig1.csv')
}
df_forplot = data.frame('Age'=rep(seq(1,52,length.out=48),2*5),
                        'Region'=rep(rep(regions,each=48),2),
                        'Counts'=c(c(t(dailycountsbyregion)),c(t(dailysevcountsbyregion))),
                        'Type'=c(rep('All Cases',length(t(dailycountsbyregion))),
                                 rep('Severe Cases',length(t(dailysevcountsbyregion))))
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
  ylab('Infant Incidence Rate per 100,000') + 
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
thailand_sevcases = read.csv('./Data/ThailandInfantCases.csv')
thailand_sevcases$dens = thailand_sevcases$cases/sum(thailand_sevcases$cases)
thailand_sevcases$age_months = thailand_sevcases$age_months * 30.75
thailand_sevcases$label='Severe Cases'

df_forplot_thailand = data.frame('Age'=seq(1,52,length.out=52),
                                 'Counts'=rep(thailand_sevcases$dens/4,each=4),
                                 'label'=rep('Severe Cases',52))

df_forplot_brazil = data.frame('Age'=rep(seq(1,52,length.out=48),2),
                               'Counts'=c(colSums(dailycountsbyregion),
                                          colSums(dailysevcountsbyregion)),
                               'Type'=c(rep('All Cases',48),
                                        rep('Severe Cases',48))
)

# Normalize counts to get empirical density
df_forplot_brazil = df_forplot_brazil %>% group_by(Type) %>%
  mutate(Dens=Counts/sum(Counts))

figure1_b = ggplot() + 
  geom_col(data=df_forplot_thailand,aes(x=Age,y=Counts),fill='grey',colour='grey') + 
  geom_area(data = df_forplot_brazil %>% filter(Age<65) %>% 
              rename(label=Type),aes(x=Age,y=Dens),colour='black',
            fill='grey',alpha=0.5,
            linewidth=0.5,linetype=1) + 
  ylab('Proportion')+
  scale_x_continuous(name='Age (Months)',
                     breaks=c(c(0,3,6,9,12)*4*1.085),
                     labels=c(0,3,6,9,12),
                     limits=c(0,12*4*1.085)) +
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
  mutate(label = case_when(name=="propmedianinfcaseneonate" ~ "Cases Among Neonates",
                           name=="propmedianinfsevcasegt4" ~ "Sev. Cases Among 5-12 Mo.")),
  all_den_byregion_y  %>% mutate(statecode='All') %>%
    pivot_longer(cols=c('regionmedianpropinfcaseneonate','regionmedianpropinfsevcasegt4')) %>% 
    mutate(label = case_when(name=="regionmedianpropinfcaseneonate" ~ "Cases Among Neonates",
                             name=="regionmedianpropinfsevcasegt4" ~ "Sev. Cases Among 5-12 Mo."))
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

ggsave('./Figure_1.png',
       grid.arrange(figure1_a_with_inset,figure1_b,figure1_c,ncol=3),
       height=4,width=8,units='in',device='png')

#### Figure 2: Estimated maternal seroprevalence
all_den_2014 = all_den %>% filter(year_not<=2014)
foi_dir = "./Data/02-foi/stanfit_v3/t51_rta_async_noType_neglectZika_final/"

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



all_den_2014$maternal_sp = sapply(1:nrow(all_den_2014),function(x) (maternal_sp_dat[dimnames(maternal_sp_dat)[[1]]==all_den_2014$year_not[x],
                                                                                    dimnames(maternal_sp_dat)[[2]]==all_den_2014$state[x],
                                                                                    c("SP")]))
all_den_2014$maternal_s1 = sapply(1:nrow(all_den_2014),function(x) (maternal_sp_dat[dimnames(maternal_sp_dat)[[1]]==all_den_2014$year_not[x],
                                                                                    dimnames(maternal_sp_dat)[[2]]==all_den_2014$state[x],
                                                                                    c("S1")]))
all_den_2014$maternal_smult = sapply(1:nrow(all_den_2014),function(x) (maternal_sp_dat[dimnames(maternal_sp_dat)[[1]]==all_den_2014$year_not[x],
                                                                                       dimnames(maternal_sp_dat)[[2]]==all_den_2014$state[x],
                                                                                       c("Smult")]))
all_den_2014$foi = sapply(1:nrow(all_den_2014),function(x) (foi_dat[dimnames(foi_dat)[[1]]==all_den_2014$year_not[x],
                                                                    dimnames(foi_dat)[[2]]==all_den_2014$state[x]]))


lastyear=2014

msp_2014 = brazil_states %>% left_join(all_den_2014 %>% filter(year_not==2014) %>% dplyr::select(StateName,maternal_sp))

msp_difference = brazil_states %>% left_join(all_den_2014 %>% filter(year_not %in% c(2000,2014)) %>% 
                                               group_by(StateName) %>% 
                                               mutate(mspdiff = maternal_sp - lag(maternal_sp)) %>% 
                                               filter(!is.na(mspdiff)) %>% 
                                               select(StateName,mspdiff))



fig2a = msp_2014 %>% ggplot() + geom_sf(aes(fill=maternal_sp)) + 
  ggtitle('A') +
  scale_fill_gradient(name='Maternal Seroprevalence\n In 2014',
                      low="white", high="blue") + 
  theme(legend.text = element_text(size=textsize-2),
        legend.title = element_text(size=textsize-2),
        legend.position = 'bottom')
fig2b = msp_difference %>% ggplot() + geom_sf(aes(fill=mspdiff)) + 
  ggtitle('B') +
  scale_fill_gradient(name='Change In Maternal\nSeroprevalence 2000-2014',
                      low="white", high="red") + 
  theme(legend.text = element_text(size=textsize-2),
        legend.title = element_text(size=textsize-2),
        legend.position = 'bottom')
ggsave('./Figure_2.png',
       grid.arrange(fig2a,fig2b,ncol=2),  
       device='png',height=4,width=7,units='in')

#### Figure 3 is schematic of model

#### Figure 4 is estimated risk profiles and fit to overall data: see "make_model_figures_clean.R"



######################################################################################
################### SUPPLEMENTARY FIGURES ############################################
######################################################################################
######################################################################################

### Figure S1: annual incidence rate by state in the whole population
ggsave('./Figure_S1.png',
       all_den %>% 
         mutate(StateName = factor(StateName,levels=rev(state_codes$StateName)),
                region = factor(region,levels=unique(rev(state_codes$region)))) %>% 
         ggplot() + 
         geom_line(aes(x=year_not,y=ir,colour=region),linewidth=1) + 
         facet_wrap(vars(StateName),ncol=5) + 
         scale_color_manual(name='Region',values=cc)+
         ylab('Annual Incidence Rate Per 100,000') +
         scale_x_continuous(name='Year',breaks=seq(2000,2025,5),labels=c('00','05','10','15','20','25'))+
         theme(legend.position = 'bottom'),
       height=6,width=9,units='in',device='png')

### Figure S2: Proportion of cases and severe cases that are in infants by region and time
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
  ylab('Proportion Of Child Cases In Infants') + 
  scale_alpha_manual(values=c(0.5,1))+
  scale_color_manual(name='Region',values=c(alpha(cc,0.5),'black'),labels=c('North','Northeast','Central-West','Southeast','South','Brazil'))+
  guides(alpha='none')+
  theme(legend.position='bottom',
        axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize),
        legend.text = element_text(size=textsize),
        legend.title = element_text(size=textsize),
        strip.text = element_text(size=textsize,hjust=0))

ggsave('./Figure_S2.png',
       figure_propsev,
       height=4,width=6.6,units='in',device='png')

### Figure S3: Age distribution of cases and severe cases up to age 2 when we have date of birth
### Figure S4-8: Age distribution of severe cases by state and 5-year period in SE, NE, and CW
### NOTE: can only run if line list of cases is available

### Age distribution by region
if (!inherits(den_cases,'try-error')) {
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
  ggsave('./FigureS3.png',
         figures_agedist,
         height=4,width=2.5,units='in',device='png')
  
  ### Figure S4-8: Age distribution of severe cases by state and 5-year period in SE, NE, and CW
  ggsave('./Figure_S4.png',
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
           facet_grid(StateName~y5_label) + 
           xlab("Age (Months)") + ylab("Case Count"),
         height=8,width=6,units='in',device='png')
  
  ggsave('./Figure_S5.png',
         den_cases  %>% 
           filter(!is.na(state) & age_years<1 & age_months<12 & severe=="Severe" & region=="Northeast") %>% 
           mutate(y5 = floor(year_not/5),
                  y5_label = case_when(y5 == 400 ~ "2000-2004",
                                       y5 == 401 ~ "2005-2009",
                                       y5 == 402 ~ "2010-2014",
                                       y5 == 403 ~ "2015-2019",
                                       y5 == 404 ~ "2020-2024"
                  ),
                  region = factor(region,levels = rev(unique(state_codes$region))),
                  StateName = gsub("Rio Grande do Norte","Rio Grande\ndo Norte",StateName)
           ) %>% 
           ggplot(aes(x=age_months)) + geom_histogram(breaks=(-1):12) + 
           facet_grid(StateName~y5_label) + 
           xlab("Age (Months)") + ylab("Case Count"),
         height=8,width=6,units='in',device='png')
  
  ggsave('./Figure_S6.png',
         den_cases  %>%
           filter(!is.na(state) & age_years<1 & age_months<12 & severe=="Severe" & region=="Central-West") %>% 
           mutate(y5 = floor(year_not/5),
                  y5_label = case_when(y5 == 400 ~ "2000-2004",
                                       y5 == 401 ~ "2005-2009",
                                       y5 == 402 ~ "2010-2014",
                                       y5 == 403 ~ "2015-2019",
                                       y5 == 404 ~ "2020-2024"
                  ),
                  region = factor(region,levels = rev(unique(state_codes$region)))
           ) %>% 
           ggplot(aes(x=age_months)) + geom_histogram(breaks=(-1):12) + 
           facet_grid(StateName~y5_label) + 
           xlab("Age (Months)") + ylab("Case Count"),
         height=8,width=6,units='in',device='png')
  
  ggsave('./Figure_S7.png',
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
           facet_grid(StateName~y5_label) + 
           xlab("Age (Months)") + ylab("Case Count"),
         height=8,width=6,units='in',device='png')
  
  ggsave('./Figure_S8.png',
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
           facet_grid(StateName~y5_label) + 
           xlab("Age (Months)") + ylab("Case Count"),
         height=8,width=6,units='in',device='png')
}

### Figure S9: Plots of mean case age, mean maternal age, and population IR against infant case age

p_modecaseage_infcaseage = all_den %>% 
  filter(mode6pluscaseage<100) %>% 
  filter(!is.na(meaninfantcaseage)) %>% mutate(agediff = meanage-meanpopage) %>%
  ggplot() + geom_point(aes(x=mode6pluscaseage,y=meaninfantcaseage))+
  geom_smooth(aes(x=mode6pluscaseage,y=meaninfantcaseage))+
  xlab('Peak Age Of Cases')+ylab('Mean Age Of Infant Case')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))
p_popir_infcaseage = all_den %>% mutate(agediff = meanage-meanpopage,
                                        log10ir = log10(ir)) %>% 
  filter(log10ir>-Inf & !is.na(log10ir) & !is.na(meaninfantcaseage)) %>%
  ggplot() + geom_point(aes(x=log10ir,y=meaninfantcaseage))+
  geom_smooth(aes(x=log10ir,y=meaninfantcaseage))+
  scale_x_continuous(name='Population Dengue Incidence Rate Per 100,000',
                     breaks = 0:4,
                     labels = 10^(0:4))+ylab('Mean Age Of Infant Case')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))
p_matage_infcaseage = all_den %>% 
  filter(!is.na(meaninfantcaseage)) %>% 
  mutate(agediff = meanage-meanpopage) %>% 
  ggplot() + geom_point(aes(x=meanmaternalage,y=meaninfantcaseage))+
  geom_smooth(aes(x=meanmaternalage,y=meaninfantcaseage))+
  xlab('Mean Maternal Age')+ylab('Mean Age Of Infant Case')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))
p_year_infcaseage = all_den %>% 
  mutate(agediff = meanage-meanpopage,
         log10ir = log10(ir)) %>% 
  filter(!is.na(meaninfantcaseage)) %>%
  ggplot() + geom_point(aes(x=year_not,y=meaninfantcaseage))+
  geom_smooth(aes(x=year_not,y=meaninfantcaseage))+
  xlab('Year')+
  ylab('Mean Age Of Infant Case')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))

p_popir_modecaseage = all_den %>% 
  mutate(agediff = meanage-meanpopage,
         log10ir = log10(ir)) %>% 
  filter(log10ir>-Inf) %>% 
  filter(mode6pluscaseage<100) %>% 
  ggplot() + geom_point(aes(x=log10ir,y=mode6pluscaseage))+
  geom_smooth(aes(x=log10ir,y=mode6pluscaseage))+
  scale_x_continuous(name='Population Dengue Incidence Rate Per 100,000',
                     breaks = 0:4,
                     labels = 10^(0:4))+
  ylab('Peak Age Of Cases')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))
p_matage_meancaseage = all_den %>% 
  filter(mode6pluscaseage<100) %>% 
  mutate(agediff = meanage-meanpopage) %>%
  ggplot() + geom_point(aes(x=meanmaternalage,y=mode6pluscaseage))+
  geom_smooth(aes(x=meanmaternalage,y=mode6pluscaseage))+
  xlab('Mean Maternal Age')+
  ylab('Peak Age Of Cases')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))
p_year_modecaseage = all_den %>% 
  filter(mode6pluscaseage<100) %>% 
  mutate(agediff = meanage-meanpopage,
         log10ir = log10(ir)) %>% 
  ggplot() + geom_point(aes(x=year_not,y=mode6pluscaseage))+
  geom_smooth(aes(x=year_not,y=mode6pluscaseage))+
  xlab('Year')+
  ylab('Peak Age Of Cases')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))

p_matage_popir = all_den %>% 
  mutate(agediff = meanage-meanpopage,
         log10ir = log10(ir)) %>% 
  filter(log10ir>-Inf) %>%
  ggplot() + geom_point(aes(x=meanmaternalage,y=log10ir))+
  geom_smooth(aes(x=meanmaternalage,y=log10ir))+
  xlab('Mean Maternal Age')+
  scale_y_continuous(name='Population Dengue Incidence Rate Per 100,000',
                     breaks = 0:4,
                     labels = 10^(0:4))+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))
p_year_popir = all_den %>% 
  mutate(agediff = meanage-meanpopage,
         log10ir = log10(ir)) %>% 
  filter(log10ir>-Inf) %>%
  ggplot() + geom_point(aes(x=year_not,y=log10ir))+
  geom_smooth(aes(x=year_not,y=log10ir))+
  xlab('Year')+
  scale_y_continuous(name='Population Dengue Incidence Rate Per 100,000',
                     breaks = 0:4,
                     labels = 10^(0:4))+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))

p_year_matage = all_den %>% 
  mutate(agediff = meanage-meanpopage,
         log10ir = log10(ir)) %>% 
  ggplot() + geom_point(aes(x=year_not,y=meanmaternalage))+
  geom_smooth(aes(x=year_not,y=meanmaternalage))+
  xlab('Year')+
  ylab('Mean Maternal Age')+
  theme(axis.title = element_text(size=textsize),
        axis.text = element_text(size=textsize/2),
        title = element_text(size=textsize))

blank <- grid.rect(gp=gpar(col="white"))

panelplot = grid.arrange(p_modecaseage_infcaseage,p_popir_infcaseage,p_year_infcaseage,
                         blank,p_popir_modecaseage,p_year_modecaseage,
                         blank,blank,p_year_popir,
                         ncol=3)

ggsave('./Figure_S9.png',
       panelplot,
       height=10,width=10,units='in',device='png')

### Figure S10: Plot of FOI against reported incidence among infants

## Summary for text
all_den_2014 %>% group_by(year_not) %>% summarise(MSP=mean(maternal_sp),
                                                  FOI=mean(foi))
all_den_2014 %>% group_by(region,StateName) %>% summarise(MSP=mean(maternal_sp),
                                                          FOI=mean(foi))
print(all_den_2014 %>% group_by(region,StateName) %>% summarise(MSP=mean(maternal_sp),
                                                                FOI=mean(foi)) %>% arrange(MSP),n=27)
print(all_den_2014 %>% group_by(region,StateName) %>% summarise(MSP=mean(maternal_sp),
                                                                FOI=mean(foi)),n=27)

print(all_den_2014 %>% filter(year_not %in% c(2000,2014)) %>% 
        group_by(statecode) %>% 
        mutate(mspdiff = maternal_sp - lag(maternal_sp)) %>% 
        filter(!is.na(mspdiff)) %>% 
        select(statecode,mspdiff) %>% 
        arrange(mspdiff),n=27)

ggsave('./Figure_S10.png',
       all_den_2014 %>% filter(infantir>0) %>% mutate(loginfantir = log10(infantir)) %>% ggplot() + 
         geom_point(aes(x=loginfantir,y=foi)) + 
         scale_y_continuous(name='Force Of Infection') + 
         scale_x_continuous(name='Infant Incidence Rate Per 100,000',
                            breaks = (-1):3,labels = 10^((-1):3)),
       device='png',height=5,width=5,units='in')

## Figure S11 is a summary of the seroprevalence review data

### Figure S12: Plot of maternal seroprevalence against reported incidence among infants
all_den_2014 %>% group_by(year_not) %>% summarise(MSP=mean(maternal_sp))
all_den_2014 %>% group_by(region) %>% summarise(MSP=mean(maternal_sp))

ggsave('./Figure_S12.png',
       all_den_2014  %>% mutate(agediff = meanage-meanpopage,log10ir = log10(ir),
                                foi_group = factor(ntile(foi,4),levels=1:4,
                                                   labels=c('Force Of Infection < 2.1%',
                                                            '2.2% < Force Of Infection < 4.1%',
                                                            '4.2% < Force Of Infection < 7.1%',
                                                            'Force Of Infection > 7.2%')
                                )) %>% 
         filter(log10ir>-Inf) %>%
         ggplot() + geom_point(aes(x=maternal_sp,y=log10ir))+
         geom_smooth(aes(x=maternal_sp,y=log10ir))+
         xlab('Maternal Seroprevalence')+
         scale_y_continuous(name='Infant Dengue Incidence Rate Per 100,000',
                            breaks = 0:4,
                            labels = 10^(0:4))+
         facet_wrap(vars(foi_group)) +
         theme(axis.title = element_text(size=textsize),
               axis.text = element_text(size=textsize),
               title = element_text(size=textsize)),
       height=6,width=6,units='in',device='png')

### Figure S13, 14, and 17-24 are all based on the infant model fits, see "make_model_figures_clean.R"


### Figure S25 and S26: Assess missing data in infant cases.
### Figure S25 is a flow chart

### The line list data is not included in this public repository and must be obtained
### from SINAN e.g. using the microdatasus package

### Line list of all cases including non-eligible
den_cases = try(read_parquet('~/UF/Research/DengueInfants/Data/denguecases_20002024.parquet'))

if (!inherits(den_cases,'try-error')) {
  den_cases = den_cases %>% mutate(age_months = ifelse(!is.na(age_crude_month),age_crude_month,
                                                       ifelse(!is.na(age_crude_days),0,NA)))
  
  den_cases = den_cases %>% filter(age_days<365 & age_years<1 & age_months<12 & severe != "Discarded") # 283,910
  
  den_cases %>% group_by(severe) %>% summarise(meanAge=mean(age_days))
  
  # First remove discarded, missing, or inconclusive
  den_cases %>% filter(severe=="Inconclusive") %>% count() # 225
  den_cases %>% filter(severe=="Missing") %>% count() # 7,511
  
  
  den_cases$birthyear = ifelse(!is.na(den_cases$dob),year(den_cases$dob),year(den_cases$symptomonsetdate-den_cases$age_days))
  den_cases$symptomyear = year(den_cases$symptomonsetdate)
  den_cases = den_cases %>% mutate(diff = abs(symptomyear-year_not),
                                   ms = month(symptomonsetdate))
  
  # Cases with difference between year of symptoms and year of notification >1 year, OR with difference of 1 year and month of symptoms in March-October
  den_cases %>% 
    filter(!(severe %in% c("Discarded","Inconclusive","Missing"))) %>% 
    filter(diff>1 | (diff==1 & !(ms %in% c(1,2,11,12)))) %>% count()  # 23,102
  
  # Missing data on state
  den_cases %>% 
    filter(!(severe %in% c("Discarded","Inconclusive","Missing"))) %>% 
    filter(diff==0 | (diff==1 & (ms %in% c(1,2,11,12)))) %>%
    filter(is.na(state)) %>% count() # 0
  
  den_cases = den_cases %>% mutate(age_days_calc_not = case_when(!is.na(dob) ~ as.numeric(notificationdate-dob),
                                                                 NA ~ age_days),
                                   age_days_calc_symp = case_when(!is.na(dob) ~ as.numeric(symptomonsetdate-dob),
                                                                  NA ~ age_days))
  
  den_cases %>% 
    filter(!(severe %in% c("Discarded","Inconclusive","Missing"))) %>% 
    filter(diff==0 | (diff==1 & (ms %in% c(1,2,11,12)))) %>%
    filter(!is.na(state)) %>% 
    filter(abs(age_days_calc_not-age_days)>90) %>% count()
  
  # Exclude those with age_crude_days>30 or age_crude_month>12
  den_cases %>% 
    filter(!(severe %in% c("Discarded","Inconclusive","Missing"))) %>% 
    filter(diff==0 | (diff==1 & (ms %in% c(1,2,11,12)))) %>%
    filter(!is.na(state)) %>% 
    filter(age_crude_days>30 | age_crude_month>12) %>% count() # 673
  
  # Exclude those with age_crude_month=0 or age_crude_years=0 (this appears to indicate
  # an erroneous entry. If age_crude_month=0, they should have age in days recorded, and similarly for age_crude_year=0.
  # Also, in this group, those with a DOB have ages recorded as much higher)
  den_cases %>% 
    filter(!(severe %in% c("Discarded","Inconclusive","Missing"))) %>% 
    filter(diff==0 | (diff==1 & (ms %in% c(1,2,11,12)))) %>%
    filter(!is.na(state)) %>% 
    filter(is.na(age_crude_days) | age_crude_days<=30) %>% 
    filter(is.na(age_crude_month) | age_crude_month<=12) %>% 
    filter(age_crude_month==0 | age_crude_years==0) %>% count() # 38
  
  den_cases %>% 
    filter(!(severe %in% c("Discarded","Inconclusive","Missing"))) %>% 
    filter(diff==0 | (diff==1 & (ms %in% c(1,2,11,12)))) %>%
    filter(!is.na(state)) %>% 
    filter(is.na(age_crude_days) | age_crude_days<=30) %>% 
    filter(is.na(age_crude_month) | age_crude_month<=12) %>% 
    filter(is.na(age_crude_month) | age_crude_month>0) %>% 
    filter(is.na(age_crude_years) | age_crude_years>0) %>%
    filter(age_days<=1) %>% count() # 65,523
  
  state_codes = read.table('./Data/01-processedData/mapping/state.txt')
  colnames(state_codes) = c("letter","number")
  
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
  
  den_cases = den_cases %>% left_join(state_codes %>% dplyr::select(letter,number,staterank,region_num) %>% dplyr::rename(state=number),by='state')
  
  den_cases = den_cases %>% 
    filter(severe != "Discarded") %>% mutate(
      missing =   1-as.numeric(!(severe %in% c("Inconclusive","Missing")) & 
                                 (diff==0 | (diff==1 & (ms %in% c(1,2,11,12)))) & 
                                 (!is.na(state)) & 
                                 (is.na(age_crude_days) | age_crude_days<=30) &
                                 (is.na(age_crude_month) | age_crude_month<=12) &
                                 (is.na(age_crude_month) | age_crude_month>0) &
                                 (is.na(age_crude_years) | age_crude_years>0) &
                                 (age_days>1)
      )
    )
  
  den_cases %>% filter(age_days>1) %>% ggplot() + geom_boxplot(aes(x=factor(missing),y=age_days))
  den_cases %>% mutate(Cat = case_when(severe %in% c("Inconclusive","Missing") ~ "Excluded",
                                       TRUE ~ "Included")) %>% 
    group_by(Cat) %>% summarise(meanAge=mean(age_days))
  
  lastyear=2024
  regions = c("North","Northeast","Central-West","Southeast","South")
  ggsave('./Figure_S26.png',
         den_cases %>% 
           filter(!is.na(state)) %>% 
           group_by(year_not,region_num) %>% summarise(PropMissing=mean(missing)) %>%
           ggplot() +
           geom_line(aes(x=year_not,y=PropMissing,colour=factor(region_num))) +
           ylab('Proportion excluded') + 
           scale_x_continuous(name='Year',
                              labels=seq(2000,lastyear+1,by=5),
                              breaks=seq(2000,lastyear+1,by=5)) + 
           scale_color_manual(name='Region',values=cc,labels=regions)+
           theme(legend.position = 'right',
                 legend.text=element_text(size=6),
                 legend.title=element_text(size=6),
                 axis.title=element_text(size=6),
                 axis.text=element_text(size=6)),
         height=3,width=4.5,units='in',device='png')
}