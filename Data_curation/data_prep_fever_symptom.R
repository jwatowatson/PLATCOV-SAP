library(lubridate)
library(readxl)
library(ggplot2)
library(tidyverse)

##******************** Clinical database *******************
##*********************************************************
##*********************************************************
##*

clin_data = haven::read_dta('~/Dropbox/MORU/Adaptive Trials/PLATCOV/Data/InterimEnrolment.dta')
temp_data = haven::read_dta('~/Dropbox/MORU/Adaptive Trials/PLATCOV/Data/InterimFUTemp.dta')
vita_data = haven::read_dta('~/Dropbox/MORU/Adaptive Trials/PLATCOV/Data/InterimVitalSigns.dta')
# x0 = vita_data[, c('Site', 'Label', 'visit', 'vd_dat', "vs_tim", "vs_temp")]
x1 = temp_data[, c('Site', 'Label', 'visit', 'fut_amdat', "fut_amtim", "fut_amtemp")]
x2 = temp_data[, c('Site', 'Label', 'visit', 'fut_pmdat', "fut_pmtim", "fut_pmtemp")]
 colnames(x1) = colnames(x2) = 
  c('Site', 'Label', 'visit', 'fut_dat', "fut_tim", "fut_temp")
temp_data = rbind(x1, x2)
temp_data = temp_data[(!is.na(temp_data$fut_dat) & 
                         !is.na(temp_data$fut_tim) &
                         !is.na(temp_data$fut_temp)), ]
# summary(lm(fut_temp ~ Site, data = temp_data[temp_data$visit=='D4',]))
# ggplot(temp_data, aes(x=visit, y=fut_temp, fill=Site)) + 
# geom_boxplot()

ind_low = temp_data$fut_temp<36
sum(ind_low)
temp_data$fut_temp[ind_low]=36

ind_high = temp_data$fut_temp>42
sum(ind_high)
temp_data$fut_temp[ind_high]=42

temp_data$temp_time = 
  as.character(as.POSIXct(apply(temp_data[, c('fut_dat','fut_tim')],
                                1, paste, collapse =' '),
                          format = '%F %H:%M:%S'))
temp_data = merge(temp_data, 
                  clin_data, 
                  by = c('Label','Site'),
                  all.x=T)

temp_data = temp_data %>% filter(!is.na(scrpassed), 
                                 scrpassed==1, 
                                 !is.na(randat),
                                 rantim != '') %>% 
  ungroup() %>%
  mutate(rangrp = sjlabelled::as_character(rangrp),
         Rand_date_time = as.character(as.POSIXct(paste(randat,
                                                        rantim, sep=' '))),
         Time = difftime(temp_time,Rand_date_time, units = 'days'),
         Timepoint_ID = gsub(x = visit,replacement = '',pattern='D'),
         Timepoint_ID = gsub(x = Timepoint_ID, replacement = '', pattern = 'H1'),
         Timepoint_ID = as.numeric(Timepoint_ID),
         Time = as.numeric(Time)) %>%
  group_by(Label) %>%
  mutate(Fever_Baseline = ifelse(any(fut_temp>37.4 & Time<1), 1, 0))

write.csv(x = temp_data, file = '../Analysis_Data/temperature_data.csv', row.names = F, quote = F)

xx=sort(table(temp_data$Label))
sort(names(xx[xx<=20]))



######################### ******* Symptoms *********** ##################################
symp=haven::read_dta('~/Dropbox/MORU/Adaptive Trials/PLATCOV/Data/InterimSymptoms.dta')
symp = symp[!is.na(symp$sq_yn),]
sort(table(symp$Label),decreasing = T)
ids_symp_data = names(which(table(symp$Label)>4)) ## need at least 4 records to be included

symp = symp %>% filter(Label %in% ids_symp_data, visit != 'D0H1') %>%
  mutate(Timepoint_ID = gsub(x=visit, pattern='D',replacement=''))



HR_data = vita_data[, c('Label','visit','vs_hr')]
HR_data$Timepoint_ID = gsub(pattern = 'D',replacement = '',x = HR_data$visit)
HR_data$Timepoint_ID = as.numeric(gsub(pattern = 'H1',replacement = '',x = HR_data$Timepoint_ID))
HR_data = HR_data[!is.na(HR_data$vs_hr), ]
HR_data$heart_rate = HR_data$vs_hr

HR_data = aggregate(heart_rate ~ Label + Timepoint_ID, data = HR_data, mean)

plot(aggregate(heart_rate ~ Timepoint_ID, data = HR_data, mean))

symp_data = merge(symp, HR_data[, c('Label', 'Timepoint_ID','heart_rate')],
                  all=T, by = c('Label', "Timepoint_ID"))
symp_data$Timepoint_ID = as.numeric(symp_data$Timepoint_ID)
symp_data = symp_data %>% arrange(Label, Timepoint_ID)

write.csv(x = symp_data, file = '../Analysis_Data/symptom_data.csv', row.names = F, quote = F)

rm(list=ls())

