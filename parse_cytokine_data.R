library(tidyverse)
ff = list.files(path='~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/Data/Cytokines/',pattern = '*.csv',full.names = T)
dat = list()
for(i in 1:length(ff)){
  dat[[i]]=readr::read_csv(file = ff[i], skip = 1)
  dat[[i]]$plate=paste('plate',i)
}
dat_all = dplyr::bind_rows(dat)

dat_all$day = NA
dat_all$day[grep('D3',x = dat_all$Sample)]='D3'
dat_all$day[grep('D0',x = dat_all$Sample)]='D0'

cytokines = unique(dat_all$Assay)
dat_all = dat_all %>% filter(Sample!='Control sample')

library(tidyverse)
library(gridExtra)
pdf('cytokine_data.pdf')

for(cytk in cytokines){
  cytk_dat = dat_all %>% filter(Assay==cytk, Concentration>0 | is.na(Concentration))
  my_lims_x = range(cytk_dat$Signal)
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  h1 =
    cytk_dat %>% filter(is.na(Concentration)) %>% ggplot(aes(x=Signal, colour = day))+
    geom_histogram() + scale_x_log10(limits = my_lims_x)+theme_minimal()
  
  p1=cytk_dat %>% ggplot(aes(y=Concentration, x = Signal, colour = plate))+
    scale_x_log10()+
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,)+
    geom_point()+theme_minimal()+ylab('pg/ml')+
    # theme(legend.position = 'none')+
    geom_smooth(aes(group = NA))+
    annotation_logticks(sides='l')+ggtitle(cytk)
  print(grid.arrange(h1,p1))
  cytk_dat %>% filter(is.na(Concentration)) %>% ggplot(aes(x=Signal))+
    geom_histogram() + scale_x_log10(limits=c(1,10^6))
}
dev.off()

dat_all = dat_all %>% filter(is.na(Concentration))
dat_all$ID=NA
for(i in 1:nrow(dat_all)){
  xx = strsplit(x = dat_all$Sample[i],split = '-',fixed = T)[[1]]
  dat_all$ID[i] = paste(xx[1],xx[2],xx[4],sep = '-')
}
cytokines = c("IFN-g","IL-1RA","IP-10")

pcr_dat = read_csv('Analysis_Data/interim_all_analysis.csv')
pcr_dat = pcr_dat %>% filter(ID %in% dat_all$ID, Time<5)

library(brms)
mod_naive = brm(log10_viral_load | cens(censor) ~ Time+CT_RNaseP+(1+Time|ID),
                data = pcr_dat,
                prior = c(prior(lkj(2), class = "cor"),
                          prior(normal(5, 5), class = Intercept),
                          prior(normal(0,1), class='b'),
                          prior(constant(7), class = 'nu')),
                family = 'student', cores = 4, chains = 4)

summary(mod_naive)
pcr_dat$preds = colMeans(posterior_epred(mod_naive,newdata = data.frame(Time=pcr_dat$Time, CT_RNaseP=mean(pcr_dat$CT_RNaseP), ID=pcr_dat$ID)))

pcr_dat_summary = pcr_dat %>% group_by(ID) %>%
  mutate(
    baseline_VL = mean(log10_viral_load[Timepoint_ID==0]),
    slope = coef(lm(preds~Time))['Time']) %>%
  distinct(ID, .keep_all = T)

hist(pcr_dat_summary$slope, breaks = 15)

dat_all_2 = merge(dat_all, pcr_dat_summary[, c('ID','slope','baseline_VL')], by = 'ID')
dat_all_2 = dat_all_2 %>% filter(Assay %in% cytokines)
dat_all_2 = dat_all_2 %>% group_by(Assay, day) %>%
  mutate(Conc_standardised = scale(log10(`Calc. Concentration pg/ml`)))

dat_all_2 %>% ggplot(aes(x=`Calc. Concentration pg/ml`, y = slope))+
  geom_point() + geom_smooth(method = lm) + scale_x_log10()+
  facet_wrap(vars(Assay, day), nrow = 4, ncol = 3)

p_all_vc=dat_all_2 %>% ggplot(aes(x=Conc_standardised, y = slope))+
  geom_point() + geom_smooth(method = lm) + theme_minimal() +
  xlab('Standardised concentration')+ylab('Viral clearance (log10 change per day)')+
  facet_wrap(vars(Assay, day), nrow = 4, ncol = 3) 

p_all_vl=dat_all_2 %>% filter(day=='D0') %>%
  ggplot(aes(x=Conc_standardised, y = baseline_VL))+
  geom_point() + geom_smooth(method = lm) + theme_minimal() +
  xlab('Standardised concentration')+ylab('Baseline log10 viral load')+
  facet_wrap(vars(Assay), nrow = 4, ncol = 3)

dat_all_2 %>% filter(day=='D3') %>%
  ggplot(aes(x=Conc_standardised, y = baseline_VL))+
  geom_point() + geom_smooth(method = lm) + theme_minimal() +
  xlab('Standardised concentration')+ylab('Baseline log10 viral load')+
  facet_wrap(vars(Assay), nrow = 4, ncol = 3)

ggsave(filename =  '~/Downloads/first_100_baseline_VL.pdf', plot = p_all_vl)
ggsave(filename =  '~/Downloads/first_100_viral_clearance.pdf', plot = p_all_vc)


summary(lm(baseline_VL ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IFN-g')))
summary(lm(baseline_VL ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IL-1RA')))
summary(lm(baseline_VL ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IP-10')))

summary(lm(slope ~ log10(`Calc. Concentration pg/ml`), data = dat_all_2%>%filter(day=='D0',Assay=='IFN-g')))

summary(lm(slope ~ Conc_standardised, data = dat_all_2%>%filter(day=='D3',Assay=='IL-1RA')))
summary(lm(slope ~ log10(`Calc. Concentration pg/ml`), data = dat_all_2%>%filter(day=='D3',Assay=='IL-1RA')))

summary(lm(slope ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IP-10')))
summary(lm(slope ~ log10(`Calc. Concentration pg/ml`), data = dat_all_2%>%filter(day=='D0',Assay=='IP-10')))


xx=dat_all_2 %>% filter(day=='D0') %>% 
  pivot_wider(names_from = Assay, 
              values_from = Conc_standardised, 
              id_cols = c(ID, slope, baseline_VL))
mod1 = lm(slope ~ . , data = xx[, !colnames(xx) %in% c('ID','baseline_VL')])
summary(mod1)

mod2 = lm(baseline_VL ~ . , data = xx[, !colnames(xx) %in% c('ID','slope')])
summary(mod2)
