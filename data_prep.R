##******************** Clinical database *******************
##*********************************************************
##*********************************************************
##*

clin_data = haven::read_dta('../Data/InterimEnrolment.dta')

# check data for missing values
ind = !is.na(clin_data$cov_symphr) & is.na(clin_data$cov_sympday)
if(sum(ind)>0) clin_data$cov_sympday[ind] = clin_data$cov_symphr[ind]/24
table(clin_data$cov_sympday, useNA = 'ifany')
clin_data$cov_sympday[is.na(clin_data$cov_sympday)] = 2

table(clin_data$cov_test, useNA = 'ifany')
table(is.na(clin_data$age_yr))
clin_data$age_yr[is.na(clin_data$age_yr)] = 25
clin_data$BMI = clin_data$weight/(clin_data$height/100)^2
table(is.na(clin_data$BMI))


##******************** Vaccine database *******************
##*********************************************************
##*********************************************************

vacc_data = haven::read_dta('../Data/InterimVaccine.dta')

## Some cleaning
vacc_data$vc_name = tolower(vacc_data$vc_name)
writeLines('\nVaccine names before cleaning:')
print(table(vacc_data$vc_name))


vacc_data$vc_name[grep(pattern = 'astrazeneca',x = vacc_data$vc_name)] = 'AZ'
vacc_data$vc_name[grep(pattern = 'sinopharm',x = vacc_data$vc_name)] = 'SPH'
vacc_data$vc_name[grep(pattern = 'sinovac',x = vacc_data$vc_name)] = 'SV'
vacc_data$vc_name[grep(pattern = 'moderna',x = vacc_data$vc_name)] = 'MD'
vacc_data$vc_name[grep(pattern = 'pfizer',x = vacc_data$vc_name)] = 'PZ'
vacc_data$vc_name[vacc_data$vc_name==''] = NA

writeLines('\nVaccine names after cleaning:')
print(table(vacc_data$vc_name, useNA = 'ifany'))

clin_data$Any_dose = NA
clin_data$N_dose = NA
vacc_date_cols = grep('dos', colnames(vacc_data))
for(i in 1:nrow(clin_data)){
  id = clin_data$Label[i]
  ind = which(vacc_data$Label==id)
  if(length(ind)==0){
    clin_data$Any_dose[i] = 'No'
    clin_data$N_dose[i] = 0
  } else {
    clin_data$Any_dose[i] = 'Yes'
    clin_data$N_dose[i] = sum(!vacc_data[ind, vacc_date_cols]=='')
  }
}


clin_data$sex = as.numeric(sjlabelled::as_character(clin_data$sex)=='Male')
clin_data$rangrp = sjlabelled::as_character(clin_data$rangrp)
clin_data$cov_test = sjlabelled::as_character(clin_data$cov_test)
clin_data$cov_test[clin_data$cov_test=='Not done']=NA
clin_data$cov_test = as.numeric(clin_data$cov_test=='Positive')


##******************** Log database ***********************
##*********************************************************
##*********************************************************

log_data = haven::read_dta('../Data/InterimSampleLog.dta')
log_data$sl_barc[log_data$sl_barc=='']=NA
log_data = log_data[!is.na(log_data$sl_barc), ]


##******************** qPCR database **********************
##*********************************************************
##*********************************************************


fnames = list.files('../Data/CSV files/')
for(i in 1:length(fnames)){
  if(i==1){
    Res = readr::read_csv(paste0('../Data/CSV files/', fnames[i],sep=''))
  } else {
    temp = readr::read_csv(paste0('../Data/CSV files/', fnames[i],sep=''))
    Res = rbind(Res,temp)
  }
}

Res = Res[!is.na(Res$`Sample ID`), ]
Res$`Lot no.` = tolower(Res$`Lot no.`)
Res$Plate = as.numeric(as.factor(Res$`Lot no.`))


writeLines('\nShowing the number of plates and the number of samples per plate:')
print(table(Res$Plate))

writeLines('\nAre there 96 samples per plate?')
print(table(table(Res$Plate)==96))

## Extract standard curve data by plate
ind = grep('std', Res$`Sample ID`)
SC = Res[ind, c('Sample ID','N/S Gene','Target conc. c/mL','Plate')]
SC$`N/S Gene`[SC$`N/S Gene`=='Undetermined'] = 40
SC$CT_NS = as.numeric(SC$`N/S Gene`)
SC$log10_true_density = log10(as.numeric(SC$`Target conc. c/mL`))
SC$ID = apply(SC[, c('Sample ID','Plate')], 1, function(x) paste(x[1],as.numeric(x[2]), sep='_'))

# Select columns
cols = c('ID','Plate','CT_NS','log10_true_density')
SC = SC[,cols]

## Extract sample data
Res = Res[!is.na(Res$BARCODE), ]

# ids_missing = sort(unique(Res$`SUBJECT ID`[!Res$`SUBJECT ID` %in% clin_data$Label]))
# clin_data_missing = data.frame(array(dim = c(length(ids_missing), ncol(clin_data))))
# colnames(clin_data_missing)=colnames(clin_data)
# clin_data_missing$Label=ids_missing
# data_TH1 <- readr::read_csv("~/Dropbox/PLATCOV/data-TH1.csv")
# data_TH1 = data_TH1[paste0('PLT-TH1-',data_TH1$randomizationID)%in%ids_missing,] 
# clin_data_missing$rangrp = data_TH1$Treatment
# xx=strsplit(as.character(as.POSIXct(data_TH1$Date,format="%a %b %d %H:%M:%S %Y")+7*60*60),split = ' ')
# 
# clin_data_missing$randat = sapply(xx, function(x) x[1])
# clin_data_missing$rantim = sapply(xx, function(x) x[2])
# clin_data_missing$Site = 'th001'
# clin_data_missing$Any_dose = 'Yes'
# clin_data_missing$N_dose = 2
# clin_data_missing$cov_test = 1
# clin_data_missing$age_yr = 20
# clin_data_missing$BMI = 20
# clin_data_missing$sex = 1
# clin_data_missing$cov_sympday = 2
# 
# clin_data = rbind(clin_data, clin_data_missing)


Res = Res[Res$`SUBJECT ID` %in% clin_data$Label, ]

## get missing data from log file
table(Res$BARCODE %in% log_data$sl_barc)
# xx=Res[!Res$BARCODE %in% log_data$sl_barc, c('BARCODE', 'SUBJECT ID')]
# write.csv(x = xx, file = '~/Downloads/missing_barcodes.csv')

writeLines('Number of samples per patient:')
range(table(Res$`SUBJECT ID`))

Res$`N/S Gene`[Res$`N/S Gene`=='Undetermined'] = 40
Res$RNaseP[Res$RNaseP=='Undetermined'] = 40

Res$`N/S Gene` = as.numeric(Res$`N/S Gene`)
Res$RNaseP = as.numeric(Res$RNaseP)


Res$CT_NS = Res$`N/S Gene`
Res$`Target conc. c/mL`[Res$`Target conc. c/mL`==0]=10
Res$log10_viral_load = log10(Res$`Target conc. c/mL`+1)
Res$CT_RNaseP = Res$RNaseP

Res$ID = Res$`SUBJECT ID`
Res$Swab_ID = Res$Location
Res$Swab_ID = gsub(Res$Swab_ID, pattern = '1',replacement = '')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = '2',replacement = '')


Res$Timepoint_ID = Res$`TIME-POINT`
table(Res$Timepoint_ID)
Res$Timepoint_ID[Res$Timepoint_ID=='D0H0']='0'
Res$Timepoint_ID[Res$Timepoint_ID=='D0PRE']='0'
Res$Timepoint_ID = gsub(Res$Timepoint_ID,pattern = 'D',replacement = '',fixed = T)
Res$Timepoint_ID = gsub(Res$Timepoint_ID,pattern = 'Tx',replacement = '',fixed = T)
Res$Timepoint_ID = as.numeric(Res$Timepoint_ID)
table(Res$Timepoint_ID, useNA = 'ifany')

Res$Site = NA
Res$Time = NA
Res$Rand_date = NA
Res$Trt = NA
Res$Any_dose = NA
Res$N_dose = NA
Res$Antibody_test = NA
Res$Age = NA
Res$Sex = NA
Res$BMI = NA
Res$Symptom_onset = NA

for(i in 1:nrow(Res)){
  id = Res$ID[i]
  ind = which(clin_data$Label==id)

  ## Calculate time since randomisation for the visit
  rand_time = as.POSIXct(paste(clin_data$randat[ind],
                               clin_data$rantim[ind], sep=' '))
  Res$Rand_date[i] = as.character(rand_time)

  sample_time = as.POSIXct(paste(Res$`COLLECTION DATE`[i],
                                 Res$`Time Collected`[i], sep=' '),
                           format = '%d-%b-%y %X')
  Res$Time[i] = difftime(sample_time,rand_time,units = 'days')

  ## Fill in clinical data
  if(length(ind)==0) print('error')
  Res$Site[i] = clin_data$Site[ind]
  Res$Trt[i] = clin_data$rangrp[ind]
  Res$Any_dose[i] = clin_data$Any_dose[ind]
  Res$N_dose[i] = clin_data$N_dose[ind]
  Res$Antibody_test[i] = clin_data$cov_test[ind]
  Res$Age[i] = clin_data$age_yr[ind]
  Res$BMI[i] = clin_data$BMI[ind]
  Res$Sex[i] = clin_data$sex[ind]
  Res$Symptom_onset[i] = clin_data$cov_sympday[ind]
}

##### Make up data - placeholder for when data are available
Res$Variant=NA
for(id in unique(Res$ID)){
  ind = Res$ID==id
  Res$Variant[ind] = sample(x = c('Alpha','Delta','Omicron'),1)
}


##***********************************************
cols = c('ID','Time','Trt','Site','Timepoint_ID',
         'BARCODE','Swab_ID','Plate','Rand_date',
         'Any_dose','N_dose','Antibody_test',
         'Age', 'Sex', 'Symptom_onset','Variant',
         'CT_NS','CT_RNaseP','log10_viral_load')
writeLines('\n column names:')
print(cols)
Res = Res[, cols]

Res = dplyr::arrange(Res, Site, ID, Time)
SC = dplyr::arrange(SC, Plate, ID)

###### Write csv files

write.csv(x = SC, file = 'interim_control_dat.csv', row.names = F)
write.csv(x = Res, file = 'interim_dat.csv', row.names = F)

rm(list=ls())
