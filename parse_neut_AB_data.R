library(tidyverse)
neut_sero = readr::read_csv('~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/Data/Serology/Neutralising_Ab_csv_files/sVNT reults (PLATCOV)_200 samples_Feb2025.csv')
neut_sero = neut_sero %>% arrange(`Subject code`) #%>% filter(Day==0)

colnames(neut_sero)
key_vars = c('Ancestral',"Delta","BA.1","BA.2",'XBB',"SAR-CoV-1")
par(las=1)
pairs(neut_sero[, key_vars],col=neut_sero$Day+1,xlim=c(0,100),ylim=c(0,100))

pp=GGally::ggpairs(data = neut_sero,columns = key_vars,
                   upper = list(continuous = "cor"),
                   lower = list(continuous = "smooth"))+
  theme_minimal()
pp
ggplot2::ggsave(pp,filename = '~/Downloads/neutralising_Ab.pdf')


df_long <- neut_sero %>% 
  #mutate(across(-c(`Subject code`,day), as.character)) %>%
  pivot_longer(-c(`Subject code`,Day), names_to = "VoC", values_to = "Inhibition") %>%
  arrange(`Subject code`, Day) %>%
  mutate(ID = `Subject code`)
pcr_dat = readr::read_csv('Analysis_Data/interim_all_analysis.csv')
pcr_dat = pcr_dat %>% filter(ID %in% df_long$ID, Time<5)
baseline_dat = pcr_dat %>% filter(Timepoint_ID==0) %>% distinct(ID, .keep_all = T) %>%
  select(ID, Rand_date,Trt,N_dose,Age,Sex,Variant2)

df_long = merge(df_long, baseline_dat, by='ID')
df_long %>% filter(Day==0, !VoC %in% c('Mu')) %>% ggplot(aes(x=Rand_date, y = Inhibition))+
  geom_point() + geom_smooth() + xlab('Date of randomisation')+
  ylab('Inhibition (%)')+
  facet_wrap(vars(VoC))+theme_minimal()



pcr_dat = readr::read_csv('Analysis_Data/interim_all_analysis.csv')
pcr_dat = pcr_dat %>% filter(ID %in% df_long$ID, Time<5)

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
  distinct(ID, .keep_all = T) %>% select(ID, baseline_VL, slope)


df_long = merge(df_long, pcr_dat_summary, by = 'ID')

df_long %>% filter(Day==0, VoC %in% key_vars) %>% 
  ggplot(aes(x=Inhibition, y = slope))+
  geom_point() + geom_smooth(method=lm) + ylab('Viral clearance slope')+
  xlab('Inhibition (%)')+
  facet_wrap(vars(VoC))+theme_minimal()

