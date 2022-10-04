library(lubridate)
library(readxl)

##******************** Clinical database *******************
##*********************************************************
##*********************************************************
##*

clin_data = haven::read_dta('../Data/InterimEnrolment.dta')
temp_data = haven::read_dta('../Data/InterimFUTemp.dta')
vita_data = haven::read_dta('../Data/InterimVitalSigns.dta')
x0 = vita_data[, c('Label', 'visit', 'vd_dat', "vs_tim", "vs_temp")]
x1 = temp_data[, c('Label', 'visit', 'fut_amdat', "fut_amtim", "fut_amtemp")]
x2 = temp_data[, c('Label', 'visit', 'fut_pmdat', "fut_pmtim", "fut_pmtemp")]
colnames(x0) = colnames(x1) = colnames(x2) = 
  c('Label', 'visit', 'fut_dat', "fut_tim", "fut_temp")
temp_data = rbind(x0, x1, x2)
temp_data = temp_data[(!is.na(temp_data$fut_dat) & 
                         !is.na(temp_data$fut_tim) &
                         !is.na(temp_data$fut_temp)), ]

ind_low = temp_data$fut_temp<36
sum(ind_low)
temp_data$fut_temp[ind_low]=36

ind_high = temp_data$fut_temp>42
sum(ind_high)
temp_data$fut_temp[ind_high]=42

temp_data$temp_time = 
  as.character(as.POSIXct(apply(temp_data[, c('fut_dat','fut_tim')], 1, paste, collapse =' '),
                          format = '%F %H:%M:%S'))

ind = !is.na(clin_data$scrpassed) & clin_data$scrpassed==0
clin_data = clin_data[!ind, ]
sort(unique(clin_data$Label))
clin_data$rangrp = sjlabelled::as_character(clin_data$rangrp)

clin_data$Rand_date_time = NA
for(i in 1:nrow(clin_data)){
  if(!is.na(clin_data$randat[i]) & !is.na(clin_data$rantim[i])){
    clin_data$Rand_date_time[i] = as.character(as.POSIXct(paste(clin_data$randat[i],
                                                                clin_data$rantim[i], sep=' ')))
  }
}
temp_data = merge(temp_data, 
                  clin_data[, c('Label','rangrp','Rand_date_time')], 
                  by = 'Label',all.x=T)



# Cross check with online randomisation app data
# Sites TH57 and TH58 did not use app so cannot cross check
rand_app_data = rbind(read.csv("~/Dropbox/PLATCOV/data-TH1.csv"),
                      read.csv("~/Dropbox/PLATCOV/data-BR3.csv"))
rand_app_data$ID = paste0('PLT-', rand_app_data$site,'-',
                          stringr::str_pad(rand_app_data$randomizationID, 3, pad = "0"))
rand_app_data$Site = plyr::mapvalues(x = rand_app_data$site, from=c('TH1','BR3'),to=c('th001','br003'))

writeLines('The following randomisation database IDs are not in clinical database:\n')
print(rand_app_data$ID[!rand_app_data$ID %in% clin_data$Label])

rand_app_data$Rand_Time = as.POSIXct(rand_app_data$Date, format = '%a %b %d %H:%M:%S %Y',tz = 'GMT')
rand_app_data$tzone = plyr::mapvalues(x = rand_app_data$site, 
                                      from = c('TH1','BR3'),
                                      to = c('Asia/Bangkok','America/Sao_Paulo'))
rand_app_data$Rand_Time_TZ=NA
for(i in 1:nrow(rand_app_data)){
  rand_app_data$Rand_Time_TZ[i] = as.character(with_tz(rand_app_data$Rand_Time[i], tzone = rand_app_data$tzone[i]))
}

temp_data = merge(temp_data, rand_app_data, 
                  by.x = c('Label'),
                  by.y = c("ID"), all.x = T)

ind = is.na(temp_data$Rand_date_time) & !is.na(temp_data$Rand_Time_TZ)
sum(ind)
# View(temp_data[ind, ])
temp_data$Rand_date_time[ind] = temp_data$Rand_Time_TZ[ind]

ind = !is.na(temp_data$Rand_date_time) & is.na(temp_data$Rand_Time_TZ)
sum(ind)
temp_data$Rand_Time_TZ[ind] = temp_data$Rand_date_time[ind]
Rand_diffs= apply(temp_data[,c('Rand_date_time','Rand_Time_TZ')],1,
                  function(x) difftime(x[1], x[2], units='mins'))

if(any(abs(Rand_diffs)>5)){
  writeLines(sprintf('More than 5 min difference in rand time for %s', 
                     temp_data$Label[which(abs(Rand_diffs)>5)]))
}
print(temp_data[which(abs(Rand_diffs)>15 & !duplicated(temp_data$Label)), c('Label','Rand_Time_TZ','Rand_date_time') ])


temp_data$Time = apply(temp_data[,c('temp_time','Rand_date_time')],1,
                          function(x) difftime(x[1], x[2], units='days'))
temp_data$Timepoint_ID = (gsub(pattern = 'D',replacement = '',x = temp_data$visit))
temp_data$Timepoint_ID = as.numeric(gsub(pattern = 'H1',replacement = '',x = temp_data$Timepoint_ID))

plot(temp_data$Timepoint_ID, temp_data$Timepoint_ID-temp_data$Time)
ind = which(temp_data$Timepoint_ID-temp_data$Time < -10)
View(temp_data[ind, ])
write.csv(x = temp_data[ind, ], file = '~/Downloads/data_problems.csv')
temp_data$Time[ind] = temp_data$Timepoint_ID[ind]
temp_data$temperature_ax = temp_data$fut_temp
hist(temp_data$temperature_ax)

temp_data$ID = temp_data$Label


cols= c('ID', 'Timepoint_ID', 'Time', 'temperature_ax')
temp_data = temp_data[, cols]
temp_data = dplyr::arrange(temp_data, ID, Time)
temp_data$Time = round(temp_data$Time, 4)
temp_data$ind_dup = duplicated(temp_data[, 1:4])
table(temp_data$ind_dup)
temp_data = temp_data[!temp_data$ind_dup, 1:4]

write.csv(x = temp_data, file = 'temperature_data.csv', row.names = F, quote = F)



symp=haven::read_dta('../Data/InterimSymptoms.dta')
symp = symp[!is.na(symp$sq_yn),]
ids_symp_data = names(which(table(symp$Label)>3))
symp = symp[symp$Label %in% ids_symp_data, ]
symp = symp[symp$visit != 'D0H1', ]


View(symp[, c(3,4,grep(pattern = 'yn', x = colnames(symp)))])
table(symp$visit)

symp$Timepoint_ID = as.numeric(gsub(pattern = 'D',replacement = '',x = symp$visit))
symp$ID = symp$Label

cols= c('ID', 'Timepoint_ID', 'sq_yn')
symp_data = symp[, cols]
symp_data = dplyr::arrange(symp_data, ID, Timepoint_ID)


HR_data = vita_data[, c('Label','visit','vs_hr')]
HR_data$Timepoint_ID = gsub(pattern = 'D',replacement = '',x = HR_data$visit)
HR_data$Timepoint_ID = as.numeric(gsub(pattern = 'H1',replacement = '',x = HR_data$Timepoint_ID))
HR_data = HR_data[!is.na(HR_data$vs_hr), ]
HR_data$heart_rate = HR_data$vs_hr
HR_data$ID = HR_data$Label

HR_data = aggregate(heart_rate ~ ID + Timepoint_ID, data = HR_data, mean)

symp_data = merge(symp_data, HR_data[, c('ID', 'Timepoint_ID','heart_rate')],
                  all=T, by = c('ID', "Timepoint_ID"))

write.csv(x = symp_data, file = 'symptom_data.csv', row.names = F, quote = F)

rm(list=ls())

