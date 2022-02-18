## Clinical database
clin_data = haven::read_dta('../Data/InterimEnrolment.dta')
clin_data$rangrp = clin_data$rangrp[sample(nrow(clin_data),nrow(clin_data),replace = F)]
log_data = haven::read_dta('../Data/InterimSampleLog.dta')
log_data = log_data[log_data$sl_barc != "", ]

## PCR data
fnames = list.files('../Data/CSV files/')
for(i in 1:length(fnames)){
  if(i==1){
    Res = readr::read_csv(paste0('../Data/CSV files/', fnames[i],sep=''),skip = 1)
  } else {
    temp = readr::read_csv(paste0('../Data/CSV files/', fnames[i],sep=''))
    Res = rbind(Res,temp)
  }
}
Res = Res[!is.na(Res$`Sample ID (Mmoloec)`), ]
Res$Plate = as.numeric(as.factor(Res$`Lot no.`))
table(Res$Plate)

table(table(Res$Plate)==96)

## Extract standard curve data by plate
ind = grep('std', Res$`Sample ID (Mmoloec)`)
SC = Res[ind, c('Sample ID (Mmoloec)','N/S Gene','Target conc. c/mL','Plate')]
SC$CT_NS = as.numeric(SC$`N/S Gene`)
SC$log10_true_density = log10(as.numeric(SC$`Target conc. c/mL`))
SC$ID = apply(SC[, c('Sample ID (Mmoloec)','Plate')], 1, function(x) paste(x[1],as.numeric(x[2]), sep='_'))

# Write csv file
cols = c('ID','Plate','CT_NS','log10_true_density')
SC = SC[,cols]

SC = dplyr::arrange(SC, Plate, ID)

write.csv(x = SC, file = 'interim_control_dat.csv', row.names = F)


## Extract sample data
Res = Res[!is.na(Res$BARCODE), ]
all(table(Res$`SUBJECT ID`)==20)

Res$`N/S Gene`[Res$`N/S Gene`=='Undetermined'] = 40
Res$RNaseP[Res$RNaseP=='Undetermined'] = 40

Res$`N/S Gene` = as.numeric(Res$`N/S Gene`)
Res$RNaseP = as.numeric(Res$RNaseP)


Res$CT_NS = Res$`N/S Gene`
Res$`Target conc. c/mL`[Res$`Target conc. c/mL`==0]=1
Res$log10_viral_load = log10(Res$`Target conc. c/mL`)
Res$CT_RNaseP = Res$RNaseP

Res$ID = Res$`SUBJECT ID`
Res$Swab_ID = Res$Location
Res$Swab_ID = gsub(Res$Swab_ID, pattern = '1',replacement = '')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = '2',replacement = '')

Res$Site = NA
Res$Time = NA
Res$Trt = NA
Res$Vaccinated = NA
Res$Vaccine_type = NA
Res$Antibody_test = NA
Res$Age = NA
Res$Sex = NA
Res$Symptom_onset = NA


## Some cleaning
clin_data$vc_name1 = tolower(clin_data$vc_name1)
table(clin_data$vc_name1)
clin_data$vc_name1[grep(pattern = 'ast',x = clin_data$vc_name1)] = 'AZ'
clin_data$vc_name1[grep(pattern = 'azt',x = clin_data$vc_name1)] = 'AZ'
clin_data$vc_name1[grep(pattern = 'sinop',x = clin_data$vc_name1)] = 'SPH'
clin_data$vc_name1[grep(pattern = 'sinov',x = clin_data$vc_name1)] = 'SV'
clin_data$vc_name1[grep(pattern = 'sinva',x = clin_data$vc_name1)] = 'SV'
clin_data$vc_name1[grep(pattern = 'mod',x = clin_data$vc_name1)] = 'MD'
clin_data$vc_name1[grep(pattern = 'pf',x = clin_data$vc_name1)] = 'PZ'
table(clin_data$vc_name1)
clin_data$vc_name1[clin_data$vc_name1=='']='None'

for(i in 1:nrow(Res)){
  id = Res$ID[i]
  ind = which(clin_data$Label==id)

  ## Calculate time since randomisation for the visit
  rand_time = as.POSIXct(paste(clin_data$randat[ind],
                               clin_data$rantim[ind], sep=' '))
  sample_time = as.POSIXct(paste(Res$`COLLECTION DATE`[i],
                                 Res$`Time Collected`[i], sep=' '),
                           format = '%d-%b-%Y %X')
  Res$Time[i] = difftime(sample_time,rand_time,units = 'days')

  ## Fill in clinical data
  if(length(ind)==0) print('error')
  Res$Site[i] = clin_data$Site[ind]
  Res$Trt[i] = sjlabelled::as_character(clin_data$rangrp[ind])
  Res$Vaccinated[i] = clin_data$vc_statyn1[ind]
  Res$Vaccine_type[i] = clin_data$vc_name1[ind]
  Res$Antibody_test[i] = as.numeric(clin_data$cov_test[ind]==1)
  Res$Age[i] = clin_data$age_yr[ind]
  Res$Sex[i] = 1 ## NEEDS UPDATING
  Res$Symptom_onset[i] = clin_data$cov_sympday[ind]
}



##### Make up data
Res$Variant=NA
for(id in unique(Res$ID)){
  ind = Res$ID==id
  Res$Variant[ind] = sample(x = c('Alpha','Delta','Omicron'),1)
}
Res$Vaccinated[is.na(Res$Vaccinated)]=0
Res$Antibody_test[is.na(Res$Antibody_test)]=0

Res = dplyr::arrange(Res, Site, Plate, ID, Time)



###### Write csv file
cols = c('ID','Time','Trt','Site',
         'BARCODE','Swab_ID','Plate',
         'Vaccinated','Vaccine_type','Antibody_test',
         'Age', 'Sex', 'Symptom_onset','Variant',
         'CT_NS','CT_RNaseP','log10_viral_load')

Res = Res[, cols]
write.csv(x = Res, file = 'interim_dat.csv', row.names = F)

rm(list=ls())
