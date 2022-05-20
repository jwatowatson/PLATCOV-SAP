##******************** Clinical database *******************
##*********************************************************
##*********************************************************
##*

clin_data = haven::read_dta('../Data/InterimEnrolment.dta')


# extract screening failure data
ind = !is.na(clin_data$scrpassed) & clin_data$scrpassed==0
screen_failure = clin_data[ind, 
                           c("Trial","Site","scrid","scrdat",
                             "scrpassed","reason_failure","scrnote")]

screen_failure$reason_failure = sjlabelled::as_character(screen_failure$reason_failure)

clin_data = clin_data[!ind, ]
sort(unique(clin_data$Label))
# check data for missing values
ind = !is.na(clin_data$cov_symphr) & is.na(clin_data$cov_sympday)
if(sum(ind)>0) clin_data$cov_sympday[ind] = clin_data$cov_symphr[ind]/24

ind = is.na(clin_data$cov_sympday)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing time since symptom onset', clin_data$Label[ind]))
}

ind = is.na(clin_data$cov_test)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing serology test', clin_data$Label[ind]))
}

clin_data$age_yr[clin_data$Label=='PLT-TH57-002'] = 21
ind = is.na(clin_data$age_yr)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing age', clin_data$Label[ind]))
}
clin_data$BMI = clin_data$weight/(clin_data$height/100)^2

ind = is.na(clin_data$BMI)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing BMI', clin_data$Label[ind]))
}


trt_distcont_data = haven::read_dta('../Data/InterimChangeTreatment.dta')


## AB rapid test data
ab_rdt = readxl::read_excel('../Data/Rapid Ab tests.xlsx')
ab_rdt$scrid = apply(ab_rdt, 1, function(x) paste0('PLT-TH1-SCR10',x[1]))
ab_rdt = ab_rdt[ab_rdt$scrid %in% clin_data$scrid &
                  !is.na(ab_rdt$C), c(2:4, 8)]
table(ab_rdt$C)
table(ab_rdt$`IgM (-, +)`)
table(ab_rdt$`IgG (-, +, ++, +++)`)

my_from = c('-','+','++','+++')
my_to = 0:3
ab_rdt$IgM = plyr::mapvalues(x = ab_rdt$`IgM (-, +)`, 
                             from=my_from,to=my_to)
ab_rdt$IgG = plyr::mapvalues(x = ab_rdt$`IgG (-, +, ++, +++)`, 
                             from=my_from,to=my_to)
ab_rdt = ab_rdt[, c('scrid','IgG','IgM')]
clin_data = merge(clin_data, ab_rdt, all = T)


### Variant data
var_data = read.csv('../Data/Variant_csv_files/variant genotyping Run 1.csv')
for(i in 1:nrow(var_data)){ 
  # extract patient ID
  var_data$Sample.ID[i]=unlist(strsplit(var_data$Sample.ID[i],split = '_'))[1]
  # simplify genotyping
  var_data$Variant.genotyping[i]=unlist(strsplit(var_data$Variant.genotyping[i],split = '_'))[1]
}

##******************** Vaccine database *******************
##*********************************************************
##*********************************************************

vacc_data = haven::read_dta('../Data/InterimVaccine.dta')

writeLines(sprintf('No vaccine data for %s',
                   clin_data$Label[!clin_data$Label %in% vacc_data$Label]))

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
  clin_data$N_dose[i] = sum(!vacc_data[ind, vacc_date_cols]=='')
  clin_data$Any_dose[i] = c('No','Yes')[1+as.numeric(clin_data$N_dose[i]>0)]
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

writeLines('Clinical data from the following patients not in database:\n')
print(unique(Res$`SUBJECT ID`[!Res$`SUBJECT ID` %in%  clin_data$Label]))
Res = Res[Res$`SUBJECT ID` %in% clin_data$Label, ]

## get missing data from log file
table(Res$BARCODE %in% log_data$sl_barc)
xx=Res[!Res$BARCODE %in% log_data$sl_barc, c('BARCODE', 'SUBJECT ID')]
write.csv(x = xx, file = '~/Downloads/missing_barcodes.csv')

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
Res$Weight = NA
Res$Symptom_onset = NA
Res$Per_protocol_all = NA
Res$Per_protocol_sample = NA
Res$IgG = NA
Res$IgM = NA

sampling_time_conflicts = c()
for(i in 1:nrow(Res)){
  
  id = Res$ID[i]
  ind = which(clin_data$Label==id)
  ## Calculate time since randomisation for the visit
  rand_time = as.POSIXct(paste(clin_data$randat[ind],
                               clin_data$rantim[ind], sep=' '))
  Res$Rand_date[i] = as.character(rand_time)
  
  # Mallika's data
  sample_time = 
    as.POSIXct(paste(Res$`COLLECTION DATE`[i],
                     Res$`Time Collected`[i], sep=' '),
               format = '%d-%b-%y %X')
  
  # get time of sampling from log file
  barcode_i = Res$BARCODE[i]
  ind_log = which(log_data$sl_barc==barcode_i)
  if(length(ind_log)>1) {
    writeLines(sprintf('duplicated entries for barcode %s', barcode_i))
    ind_log = ind_log[1]
  }
  # get the sample time from the log file
  if(length(ind_log)>0){
    sample_time_log = 
      as.POSIXct(paste(log_data$sl_sampdat[ind_log],
                       log_data$sl_samptim[ind_log], sep=' '))
  } else {
    sample_time_log = NA # NA if missing
  }
  
  s_times = c(sample_time, sample_time_log)
  if(all(is.na(s_times))){
    writeLines(sprintf('No sample time for patient %s at timepoint %s',id, Res$`TIME-POINT`[i]))
  } else {
    if(all(!is.na(s_times))){
      if(sample_time_log != sample_time){
        writeLines(sprintf('Conflicting sample times for barcode %s',barcode_i))
        sampling_time_conflicts=c(sampling_time_conflicts,barcode_i)
        my_sample_time = s_times[2] # trust the log as updated by Padd
      } else {
        my_sample_time = s_times[1]
      }
    } else {
      my_sample_time = s_times[which(!is.na(s_times))]
    }
  }
 
  Res$Time[i] = difftime(sample_time,rand_time,units = 'days')
  if(is.na(Res$Time[i])) {
    writeLines(sprintf('Missing sample time for patient %s',id))
  } else {
    if(Res$Time[i] < 0){
      writeLines(sprintf('Negative sample time for patient %s',id))
    }
  }
  ## Fill in clinical data
  if(length(ind)==0) print('error')
  Res$Site[i] = clin_data$Site[ind]
  Res$Trt[i] = clin_data$rangrp[ind]
  Res$Any_dose[i] = clin_data$Any_dose[ind]
  Res$N_dose[i] = clin_data$N_dose[ind]
  Res$Antibody_test[i] = clin_data$cov_test[ind]
  Res$Age[i] = clin_data$age_yr[ind]
  Res$BMI[i] = clin_data$BMI[ind]
  Res$Weight[i] = clin_data$weight[ind]
  Res$Sex[i] = clin_data$sex[ind]
  Res$Symptom_onset[i] = clin_data$cov_sympday[ind]
  Res$IgG[i] = clin_data$IgG[ind]
  Res$IgM[i] = clin_data$IgM[ind]
  
  if(id %in% trt_distcont_data$Label){
    ind_discont=which(trt_distcont_data$Label==id)
    if(sample_time > as.POSIXct(trt_distcont_data$cm_stdat[ind_discont],
                                format = '%d/%m/%Y')){
      Res$Per_protocol_sample[i]=0
    } else {
      Res$Per_protocol_sample[i]=1
    }
  } else {
    Res$Per_protocol_sample[i]=1
  }
}


## Add genotyping data
Res$Variant=NA
Res$Variant_Imputed=1 #1: imputed; 0: genotyped

# Imputation based on date
ind_Delta = Res$Rand_date < as.POSIXct('2021-12-17')
ind_BA1 = Res$Rand_date >= as.POSIXct('2021-12-17') &
  Res$Rand_date < as.POSIXct('2022-02-13')
ind_BA2 = Res$Rand_date >= as.POSIXct('2022-02-13')
Res$Variant[ind_Delta] = 'Delta'
Res$Variant[ind_BA1] = 'BA.1'
Res$Variant[ind_BA2] = 'BA.2'

for(id in unique(Res$ID)){
  ind = Res$ID==id
  if(id %in% var_data$Sample.ID){
    k = which(var_data$Sample.ID==id)
    Res$Variant[ind] = var_data$Variant.genotyping[k]
    Res$Variant_Imputed[ind] = 0
  } 
}



##***********************************************
cols = c('ID','Time','Trt','Site','Timepoint_ID',
         'BARCODE','Swab_ID','Plate','Rand_date',
         'Any_dose','N_dose','Antibody_test','Weight','BMI',
         'Age', 'Sex', 'Symptom_onset','Variant','Variant_Imputed',
         'CT_NS','CT_RNaseP',
         'Per_protocol_sample','IgG','IgM')
writeLines('\n column names:')
print(cols)
Res = Res[, cols]

Res = dplyr::arrange(Res, Site, ID, Time)
SC = dplyr::arrange(SC, Plate, ID)


###### Transformation to RNA copies per mL
# We fit a mixed effects model to the control data to estimate standard curves for each qPCR plate. 
# This includes random slopes and random intercepts for each assay.
control_dat = dplyr::arrange(SC, CT_NS, Plate)
control_dat$CT = control_dat$CT_NS
control_dat$CT[control_dat$CT_NS==40]=NA
control_dat$batch = as.factor(control_dat$Plate)

library(lme4)
conv_mod = lmer(log10_true_density ~ 1 + CT + (1+CT|batch), 
                data = control_dat,
                control = lmerControl(optimizer ="Nelder_Mead"))

preds = predict(conv_mod)
plot(control_dat$CT, control_dat$log10_true_density, xlim=c(20,40))
for(bb in levels(control_dat$batch)){
  ind = control_dat$batch==bb
  lines(control_dat$CT[ind], preds[ind])
}

preds_all = 
  predict(conv_mod,
          newdata = data.frame(CT=Res$CT_NS,
                               batch=as.factor(Res$Plate)))
preds_cens = 
  predict(conv_mod, 
          newdata = data.frame(CT=rep(40,nrow(Res)),
                               batch=as.factor(Res$Plate)))
Res$log10_viral_load = preds_all
Res$log10_cens_vl = preds_cens

Res$log10_viral_load[Res$CT_NS==40]=
  Res$log10_cens_vl[Res$CT_NS==40]
table(Res$CT_NS==40)
table(Res$log10_viral_load == Res$log10_cens_vl)



###### Write csv files
# Overall data files
write.csv(x = SC, file = 'interim_control_dat.csv', row.names = F)
write.csv(x = Res, file = 'interim_dat.csv', row.names = F)
write.csv(x = screen_failure, file = 'interim_screening_dat.csv', row.names = F)


# Specific analysis data files
IDs_Ivermectin = unique(Res$ID[Res$Trt %in% c('Ivermectin',"No study drug")])
IDs_pos_control = unique(Res$ID[Res$Rand_date < "2021-12-17 00:00:00" & 
                                  Res$Trt == 'Regeneron' & 
                                  Res$Site == 'th001'])
Res_ivermectin = Res[Res$ID %in% c(IDs_Ivermectin, IDs_pos_control), ]
write.csv(x = Res_ivermectin, file = 'Ivermectin_analysis.csv', row.names = F)

writeLines('Barcodes of conflicting sample times:')
print(sampling_time_conflicts)
rm(list=ls())
