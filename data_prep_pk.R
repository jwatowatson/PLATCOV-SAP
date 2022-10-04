dosing_data = haven::read_dta('../Data/InterimStudyDrugIvermectin.dta')
dosing_data = dosing_data[!is.na(dosing_data$da_dat),]
clin_data = haven::read_dta('../Data/InterimEnrolment.dta')

dosing_data$Date_Time = NA
for(i in 1:nrow(dosing_data)){
  dosing_data$Date_Time[i] = as.character(as.POSIXct(paste(dosing_data$da_dat[i],
                                                           dosing_data$da_sttim[i], 
                                                           sep=' ')))
}
dosing_data$Time = NA
for(id in unique(dosing_data$Label)){
  ind = which(dosing_data$Label==id)
  t_0 = min(dosing_data$Date_Time[ind])
  for(j in ind){
    dosing_data$Time[j] = difftime(dosing_data$Date_Time[j],
                                   t_0, units = 'hours')
  }
}
dosing_data = merge(dosing_data, clin_data[, c('Label','weight','age_yr','sex')],by='Label')

weight_based_doses = dosing_data$da_tabs[!duplicated(dosing_data$Label)]*6/dosing_data$weight[!duplicated(dosing_data$Label)]

par(las=1, cex.lab=1.3, cex.axis=1.3, family='serif')
hist(weight_based_doses,ylab='Number of patients',
     main = '', xlab='Daily Ivermectin dose (ug/kg)',
     breaks = seq(0.48,0.6, by=.01))
mean(weight_based_doses)
median(weight_based_doses)
range(weight_based_doses)
plot(dosing_data$weight[!duplicated(dosing_data$Label)], weight_based_doses,xlab='Weight', ylab = 'Dose (ug/kg)',panel.first=grid())

pk_dat = readr::read_csv('../Data/PK/IVM_levels.csv')
pk_dat$Date_Time = NA
for(i in 1:nrow(pk_dat)){
  pk_dat$Date_Time[i] = as.character(as.POSIXct(paste(pk_dat$`Collection date`[i],
                                                      pk_dat$`Collection time`[i], 
                                                      sep=' '),
                                                format='%d-%b-%y %H:%M:%S'))
}

pk_dat$Time = NA
for(id in unique(pk_dat$`Patient code`)){
  ind = which(pk_dat$`Patient code`==id)
  t_0 = min(dosing_data$Date_Time[dosing_data$Label==id])
  for(j in ind){
    pk_dat$Time[j] = difftime(pk_dat$Date_Time[j],
                              t_0, units = 'hours')
  }
}
hist(pk_dat$Time, breaks = seq(-100,200,by=10))
pk_dat$Label = pk_dat$`Patient code`
pk_dat = merge(pk_dat, clin_data[, c('Label','weight','age_yr','sex')],by='Label')

IVM_nonmem1 = data.frame(ID = NA,
                         TIME = dosing_data$Time,
                         AMT = dosing_data$da_tabs*6,
                         ODV = NA,
                         LNDV = NA,
                         MDV = 1,
                         BLQ = 0,
                         CMT = 1,
                         EVID = 1,
                         AGE = dosing_data$age_yr,
                         SEX = as.numeric(dosing_data$sex),
                         BW = dosing_data$weight,
                         ID_true = dosing_data$Label)

IVM_nonmem2 = data.frame(ID=NA,
                         TIME = pk_dat$Time,
                         AMT = NA,
                         ODV = pk_dat$`Ivermectin (IVM) (ng/mL)`,
                         LNDV = log(pk_dat$`Ivermectin (IVM) (ng/mL)`),
                         MDV = 0,
                         BLQ = 0,
                         CMT = 2,
                         EVID = 0,
                         AGE = pk_dat$age_yr,
                         SEX = as.numeric(pk_dat$sex),
                         BW = pk_dat$weight,
                         ID_true = pk_dat$`Patient code`)


IVM_nonmem = rbind(IVM_nonmem1, IVM_nonmem2)
IVM_nonmem$ID = as.numeric(as.factor(IVM_nonmem$ID_true))
IVM_nonmem = dplyr::arrange(IVM_nonmem, ID, TIME)


colnames(IVM_nonmem)[1] = '#ID'
IVM_nonmem = rbind(c('#','hours','mg','ng/ml','log','','','','','years',
                     'M=1','kg',''),
                   IVM_nonmem)


plot(diff(as.numeric(IVM_nonmem$TIME[-1])), IVM_nonmem$LNDV[-(1:2)],
     xlim = c(0,100))

readr::write_csv(x = IVM_nonmem, file = '../Data/PK/IVM_NONMEM_data.csv',na = '.')


### Read NONMEM output

id_map = IVM_nonmem[!duplicated(IVM_nonmem$`#ID`), c('#ID', 'ID_true'),][-1,]
pk_summaries = readr::read_csv('Ivermectin/NONMEM_out.csv')[,c('ID','F1','AUC_72','CMAX_72','AUC_168','CMAX_168')]
pk_summaries$ID = plyr::mapvalues(pk_summaries$ID, from = id_map$`#ID`, to=id_map$ID_true)
readr::write_csv(x = pk_summaries, file = 'Ivermectin/PK_summaries.csv')
