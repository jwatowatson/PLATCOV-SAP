library(tidyverse)
library(brms)
library(ggforce)
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

RUN_MODEL=F

platcov_dat = readr::read_csv('Analysis_Data/interim_all_analysis.csv')
length(unique(platcov_dat$ID))
xx=platcov_dat %>% filter(Timepoint_ID==0) %>% 
  group_by(ID) %>%
  mutate(Baseline_VL = mean(log10_viral_load)) %>%
  ungroup() %>% mutate(n = as.numeric(as.factor(ID))) %>%
  distinct(ID, .keep_all = T) 
pp=xx%>%
  ggplot(aes(x=Rand_date, y=Baseline_VL))+
  geom_point() + theme_minimal() + geom_smooth(aes(group=NA))+
  ylab('Baseline viral load (genomes/ml)')+
  xlab('Date of randomisation') + ggtitle(paste('n =',max(xx$n)))+
  geom_hline(yintercept = log10(250))
ggsave(filename = 'Plots/Baseline_VL_time.pdf',plot = pp)

unique(platcov_dat$Trt)
trts_zero_effect = c('No study drug',"Ivermectin","Favipiravir","Hydroxychloroquine")
trts_blinded = "Nitazoxanide"
platcov_dat = platcov_dat %>% filter(#Trt %in% trts_zero_effect,
  !is.na(CT_RNaseP),
  !Trt %in% trts_blinded,
  Timepoint_ID < 6)


range(platcov_dat$Rand_date)
length(unique(platcov_dat$ID))

platcov_dat = platcov_dat %>% 
  mutate(Epoch = cut(Rand_date,breaks = 10,include.lowest = T),
         CT_RNaseP_scaled = scale(CT_RNaseP)) %>% 
  group_by(ID) %>%
  mutate(Baseline_VL = mean(log10_viral_load[Timepoint_ID==0]),
         Any_FUP = any(Timepoint_ID>0),
         n_pcrs = n()) %>%
  ungroup() %>%
  filter(Baseline_VL>log10(250), Any_FUP, n_pcrs>4)

baseline_dat = platcov_dat %>% distinct(ID, .keep_all = T)

baseline_dat %>% 
  ggplot(aes(x=Rand_date, colour = Epoch)) + 
  geom_histogram(binwidth = 7) + 
  xlab('Date of randomisation')+
  theme_minimal()
# set the reference level to be early 2023
cc=levels(platcov_dat$Epoch)
cc = cc[cc != "2023-04-24"]
platcov_dat$Epoch = factor(platcov_dat$Epoch, levels = c("2023-04-24",cc))
trts = unique(platcov_dat$Trt); trts = trts[trts!='No study drug']
platcov_dat$Trt = factor(platcov_dat$Trt, levels = c('No study drug',trts))

#f_path = 'Rout/brms_fit_NSD.RData'
f_path = 'Rout/brms_fit_all.RData'
if(RUN_MODEL){
  ff = as.formula(log10_viral_load | cens(censor) ~ Time*Trt - Trt + Time*Epoch + Site + Symptom_onset + CT_RNaseP_scaled + (1+Time|ID))
  mod1 = brm(formula = ff,
             data = platcov_dat,
             prior = c(prior(lkj(2), class = "cor"),
                       prior(normal(5, 5), class = Intercept),
                       prior(normal(0,2), class='b'),
                       prior(constant(7), class = 'nu')),
             family = 'student', cores = 4, chains = 4, 
             iter = 1500, warmup = 500)
  
  mod_HCQ = brm(log10_viral_load | cens(censor) ~ Time*Trt - Trt + Time*Epoch + Site + Symptom_onset + CT_RNaseP_scaled + (1+Time|ID),
             data = platcov_dat,
             prior = c(prior(lkj(2), class = "cor"),
                       prior(normal(5, 5), class = Intercept),
                       prior(normal(0,2), class='b'),
                       prior(constant(7), class = 'nu')),
             family = 'student', cores = 4, chains = 4, 
             iter = 1500, warmup = 500)
  
  save(mod1, file = f_path)
} else {
  load(f_path)
}
brms::mcmc_plot(mod1,variable ='Time:Epoch',regex = T)+theme_minimal()
brms::mcmc_plot(mod1,variable ='Trt',regex = T)+theme_minimal()

brms::mcmc_plot(mod1,variable ='CT_RNaseP',regex = T)
dat_pred = platcov_dat %>% mutate(CT_RNaseP_scaled=0)
platcov_dat$pred_VL =
  colMeans(posterior_epred(mod1,newdata = dat_pred))

platcov_summary = platcov_dat %>% 
  group_by(ID) %>%
  mutate(
    baseline_VL = mean(log10_viral_load[Timepoint_ID==0]),
    slope = coef(lm(pred_VL~Time))['Time'],
    half_life = ifelse(slope< -0.05, 24*log10(1/2)/slope, NA)
  ) %>%
  distinct(ID, .keep_all = T) %>%
  select(ID, baseline_VL, slope, half_life,Rand_date, Trt, 
         Age,Weight,Time_since_last_dose,Sex, Symptom_onset)

platcov_summary %>% ggplot(aes(x=Rand_date, y = slope, colour = Trt))+geom_point()+
  # geom_smooth(aes(group=NA))+
  theme_minimal()+ xlab('Date of randomisation')+
  ylab('Slope: log10 change per day')
24*log10(1/2)/(-0.4)
24*log10(1/2)/(-0.8)

platcov_summary %>% ggplot(aes(x=Rand_date, y = half_life,colour = Trt))+
  geom_point()+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  geom_quantile(aes(group=NA),method = "rqss", lambda = 100, quantiles = 0.5)+
  theme_minimal()+ xlab('Date of randomisation')+
  ylab('Clearance half life (hours)')+
  ggtitle('No study drug or ineffective antiviral')


platcov_summary %>% ggplot(aes(x=Age, y = slope,colour = Trt))+
  geom_point()+
  geom_smooth(aes(group=NA),method = lm)+
  theme_minimal()+ xlab('Age (years)')+
  ylab('Clearance half life (hours)')+
  ggtitle('No study drug or ineffective antiviral')

mod_slope = mgcv::gam(slope ~ Age + Sex + Symptom_onset+Weight+Time_since_last_dose,
                      family = gaussian(), data = platcov_summary)
summary(mod_slope)
platcov_summary %>% ggplot(aes(x=baseline_VL, y=slope))+
  geom_point()+theme_minimal()+xlab('Baseline viral load')+
  ylab('Clearance slope')

pdf('Plots/Individual_fits_NSD.pdf')
p = platcov_dat %>%
  #  filter(ID %in% unique(platcov_dat$ID)[1:40]) %>%
  ggplot(aes(x=Time, y=log10_viral_load,colour=censor))+
  geom_point()+ylim(range(platcov_dat$log10_viral_load))+
  xlim(0,6)+ ylab('Viral genomes per ml (log10)')+ 
  xlab('Days from randomisation')+
  geom_line(aes(x=Time, y = pred_VL, colour = 'Model'))+
  facet_wrap_paginate(vars(ID), nrow = 5, ncol = 5)+theme_minimal()
k = n_pages(p)

for(i in 1:k){
  p = platcov_dat %>%
    #  filter(ID %in% unique(platcov_dat$ID)[1:40]) %>%
    ggplot(aes(x=Time, y=log10_viral_load,colour=censor))+
    geom_point()+ylim(range(platcov_dat$log10_viral_load))+
    xlim(0,6)+ ylab('Viral genomes per ml (log10)')+ 
    xlab('Days from randomisation')+
    geom_line(aes(x=Time, y = pred_VL, colour = 'Model'))+
    facet_wrap_paginate(vars(ID), nrow = 5, ncol = 5,page = i)+theme_minimal()
  print(p)
}
dev.off()
