library(lubridate)
library(readxl)

##******************** Clinical database *******************
##*********************************************************
##*********************************************************
##*

clin_data = haven::read_dta('../Data/InterimEnrolment.dta')
final_status = haven::read_dta('../Data/InterimFinalStatus.dta')
AE_data = haven::read_dta('../Data/InterimAE.dta')

final_status = final_status[!is.na(final_status$fs_compyn), ]
# extract screening failure data
ind = !is.na(clin_data$scrpassed) & clin_data$scrpassed==0
screen_failure = 
  clin_data[ind, 
            c("Trial","Site","scrid","scrdat",
              "scrpassed","reason_failure","scrnote")]

screen_failure$reason_failure = sjlabelled::as_character(screen_failure$reason_failure)

clin_data = clin_data[!ind, ]
sort(unique(clin_data$Label))

clin_data$rangrp = sjlabelled::as_character(clin_data$rangrp)
clin_data$cov_test = sjlabelled::as_character(clin_data$cov_test)
clin_data$cov_test[clin_data$cov_test=='Not done']=NA
clin_data$cov_test = as.numeric(clin_data$cov_test=='Positive')

AE_data = merge(AE_data, clin_data[, c('Label','rangrp')])

# check data for missing values
ind = !is.na(clin_data$cov_symphr) & is.na(clin_data$cov_sympday)
if(sum(ind)>0) clin_data$cov_sympday[ind] = clin_data$cov_symphr[ind]/24

ind = is.na(clin_data$cov_sympday)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing time since symptom onset', clin_data$Label[ind]))
  clin_data$cov_sympday[is.na(clin_data$cov_sympday)]=2
}

ind = is.na(clin_data$cov_test)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing serology test', clin_data$Label[ind]))
  clin_data$cov_test[is.na(clin_data$cov_test)]=1
}

clin_data$age_yr[clin_data$Label=='PLT-TH57-002'] = 21
ind = is.na(clin_data$age_yr) & !is.na(clin_data$dob_my)

# calculate age at randomisation
for(i in which(ind)){
  clin_data$age_yr[i] = trunc((as.POSIXct(clin_data$dob_my[i], format='%d/%m/%Y') %--% as.POSIXct(clin_data$randat[i])) / years(1))
}

ind = is.na(clin_data$age_yr)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing age', clin_data$Label[ind]))
}
clin_data$BMI = clin_data$weight/(clin_data$height/100)^2

ind = is.na(clin_data$BMI)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing BMI', clin_data$Label[ind]))
}


# Cross check with online randomisation app data
# Sites TH57 and TH58 did not use app so cannot cross check
rand_app_data = rbind(read.csv("~/Dropbox/PLATCOV/data-TH1.csv"),
                      read.csv("~/Dropbox/PLATCOV/data-BR3.csv"))
rand_app_data$ID = paste0('PLT-', rand_app_data$site,'-',
                          stringr::str_pad(rand_app_data$randomizationID, 3, pad = "0"))
rand_app_data$Site = plyr::mapvalues(x = rand_app_data$site, from=c('TH1','BR3'),to=c('th001','br003'))

writeLines('The following randomisation database IDs are not in clinical database:\n')
print(rand_app_data$ID[!rand_app_data$ID %in% clin_data$Label])

# rand_app_data = rand_app_data[rand_app_data$ID%in% clin_data$Label, ]
rand_app_data$Rand_Time = as.POSIXct(rand_app_data$Date, format = '%a %b %d %H:%M:%S %Y',tz = 'GMT')
rand_app_data$tzone = plyr::mapvalues(x = rand_app_data$site, 
                                      from = c('TH1','BR3'),
                                      to = c('Asia/Bangkok','America/Sao_Paulo'))
rand_app_data$Rand_Time_TZ=NA
for(i in 1:nrow(rand_app_data)){
  rand_app_data$Rand_Time_TZ[i] = as.character(with_tz(rand_app_data$Rand_Time[i], tzone = rand_app_data$tzone[i]))
}

clin_data = merge(clin_data, rand_app_data, 
                  by.x = c('Label','Site'),
                  by.y = c("ID",'Site'), all = T)

# discrepancies between randomisation database and CRFs?
clin_data$sex_char = sjlabelled::as_character(clin_data$sex.x)
clin_data$Rand_date_time = NA
for(i in 1:nrow(clin_data)){
  if(!is.na(clin_data$randat[i]) & !is.na(clin_data$rantim[i])){
    clin_data$Rand_date_time[i] = as.character(as.POSIXct(paste(clin_data$randat[i],
                                                                clin_data$rantim[i], sep=' ')))
  }
}

# variables we need are:
# * age
# * sex
# * time since symptom onset
# * number of vaccine doses
# * Serology rapid test

####### Age #########
table(clin_data$rangrp)
clin_data$rangrp[clin_data$rangrp=='Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir'
clin_data$rangrp[clin_data$rangrp=='Nitazoxanida']='Nitazoxanide'
ind = is.na(clin_data$Treatment) & !is.na(clin_data$rangrp)
clin_data$Treatment[ind] = clin_data$rangrp[ind]

ind = !is.na(clin_data$Treatment) & is.na(clin_data$rangrp)
clin_data$rangrp[ind] = clin_data$Treatment[ind]


if(any(! clin_data$Treatment == clin_data$rangrp)) {
  writeLines(sprintf('Randomisation inconsistent for %s', clin_data$Label[clin_data$Treatment != clin_data$rangrp]))
}

####### Age #########
ind = is.na(clin_data$age_yr) & !is.na(clin_data$age)
clin_data$age_yr[ind] = clin_data$age[ind]
ind = !is.na(clin_data$age_yr) & is.na(clin_data$age)
clin_data$age[ind] = clin_data$age_yr[ind]
if(any(is.na(clin_data$age))){
  writeLines(sprintf('Age missing for %s', clin_data$Label[is.na(clin_data$age)]))
}

if(any(!clin_data$age == clin_data$age_yr)) {
  writeLines(sprintf('Age inconsistent for %s', clin_data$Label[clin_data$age != clin_data$age_yr]))
}
print(clin_data[clin_data$age != clin_data$age_yr, c('Label','age','age_yr')])

####### Sex #########
ind = is.na(clin_data$sex_char) & !is.na(clin_data$sex.y)
clin_data$sex_char[ind] = clin_data$sex.y[ind]
ind = !is.na(clin_data$sex_char) & is.na(clin_data$sex.y)
clin_data$sex.y[ind] = clin_data$sex_char[ind]
if(any(is.na(clin_data$sex_char))){
  writeLines(sprintf('Sex missing for %s', clin_data$Label[is.na(clin_data$sex_char)]))
}

if(any(!clin_data$sex_char == clin_data$sex.y)) {
  writeLines(sprintf('Sex inconsistent for %s', clin_data$Label[clin_data$sex_char != clin_data$sex.y]))
}
clin_data$sex = as.numeric(clin_data$sex_char=='Male')

####### Rand Time #########
ind = is.na(clin_data$Rand_date_time) & !is.na(clin_data$Rand_Time_TZ)
clin_data$Rand_date_time[ind] = clin_data$Rand_Time_TZ[ind]
ind = !is.na(clin_data$Rand_date_time) & is.na(clin_data$Rand_Time_TZ)
clin_data$Rand_Time_TZ[ind] = clin_data$Rand_date_time[ind]
Rand_diffs= apply(clin_data[,c('Rand_date_time','Rand_Time_TZ')],1,
                  function(x) difftime(x[1], x[2], units='mins'))
if(any(abs(Rand_diffs)>5)){
  writeLines(sprintf('More than 5 min difference in rand time for %s', 
                     clin_data$Label[which(abs(Rand_diffs)>5)]))
}
print(clin_data[which(abs(Rand_diffs)>15), c('Label','Rand_Time_TZ','Rand_date_time') ])


clin_data = clin_data[, c('Label','rangrp','cov_test','sex','age_yr','randat',
                          "Rand_date_time",'BMI','weight','cov_sympday','Site')]

## Per protocol for treatment data
trt_distcont_data = haven::read_dta('../Data/InterimChangeTreatment.dta')



### Variant data
## PCR variant data (up until July 2022)
fnames_var = list.files('../Data/PCR genotyping/',pattern = 'variant genotyping',full.names = T)
dat = lapply(fnames_var,read_excel)
for(i in 1:length(dat)) dat[[i]] = dat[[i]][, c("SUBJECT ID","Summary")]
dat = do.call(what = rbind, dat)
dat = dat[!is.na(dat$Summary),]
var_data = dat[!duplicated(dat$`SUBJECT ID`), ]
colnames(var_data)=c('ID', 'Summary_PCR')

var_data$Summary_PCR = plyr::mapvalues(x = var_data$Summary_PCR, 
                                       from = c('Delta_B1.617.2','Omicron BA.2_B.1.1.529',
                                                'Omicron_B.1.1.529','undetermined'),
                                       to = c('Delta','BA.2','BA.1',NA))
writeLines(sprintf('We have PCR variant genotyping for %s patients',nrow(var_data)))

## Nanopore data
# fnames_var = list.files('../Data/Nanopore sequencing/',full.names = T)
res_cols = c('ID', 'Lineage', "Scorpio call", 'QC')
res_no_issue=read_excel(path = '../Data/Nanopore sequencing/PLATCoV_classification.xlsx',sheet='No issue')
res_no_issue$ID = sapply(res_no_issue$`Sequence name`, FUN = function(x) strsplit(x,split = '_')[[1]][1])
res_no_issue$QC = 'no issue'
if(any(duplicated(res_no_issue$ID))) print('warning - some duplicates in the nanopore output - will be removed')
res_no_issue = res_no_issue[!duplicated(res_no_issue$ID), ]


res_flagged=read_excel(path = '../Data/Nanopore sequencing/PLATCoV_classification.xlsx',sheet='Flagged issues')
res_flagged$ID = sapply(res_flagged$`Sequence name`, FUN = function(x) strsplit(x,split = '_')[[1]][1])
res_flagged$QC = 'flagged issue'
if(any(duplicated(res_flagged$ID))) print('warning - some duplicates in the nanopore output - will be removed')
res_flagged = res_flagged[!duplicated(res_flagged$ID), ]

res_low_cov=read_excel(path = '../Data/Nanopore sequencing/PLATCoV_classification.xlsx',sheet='Low coverage')
res_low_cov$ID = sapply(res_low_cov$`Sequence name`, FUN = function(x) strsplit(x,split = '_')[[1]][1])
res_low_cov$QC = 'low coverage'
if(any(duplicated(res_low_cov$ID))) print('warning - some duplicates in the nanopore output - will be removed')
res_low_cov = res_low_cov[!duplicated(res_low_cov$ID), ]


res = rbind(res_no_issue[, res_cols], res_flagged[, res_cols], res_low_cov[ res_cols])
all(!duplicated(res$ID))
res$Type = 'Nanopore'

table(res$Lineage)
table(res$`Scorpio call`)
res$Summary_nanopore = plyr::mapvalues(x = res$`Scorpio call`, 
                                       from = c('Delta (B.1.617.2-like)',
                                                'Omicron (BA.1-like)',
                                                'Omicron (BA.2-like)',
                                                'Omicron (BA.4-like)',
                                                'Omicron (BA.5-like)',
                                                'Omicron (Unassigned)'),
                                       to = c('Delta',
                                              'BA.1',
                                              'BA.2',
                                              'BA.4',
                                              'BA.5',
                                              'Omicron'))

table(res$Summary_nanopore)

# dat = lapply(fnames_var, read.csv)
# for(i in 1:length(dat)) {
#   dat[[i]]$Clade[grep(pattern = 'Omicron', x = dat[[i]]$Clade,ignore.case = T)]='Omicron'
#   dat[[i]]$Summary = apply(dat[[i]][, c('Clade', 'Lineage')], 1, paste, collapse='_')
#   dat[[i]] = dat[[i]][, c("SUBJECT.ID","Summary")]
# }
# dat = do.call(what = rbind, dat)
# dat = dat[!is.na(dat$Summary),]
# dat = dat[!duplicated(dat$SUBJECT.ID), ]
writeLines(sprintf('We have nanopore variant typing for %s patients',nrow(res)))

variant_data = merge(var_data, res, by='ID', all = T)
variant_data = dplyr::arrange(variant_data, ID)
variant_data$Summary = variant_data$Summary_nanopore

ind = is.na(variant_data$Lineage) & !is.na(variant_data$Summary_PCR)
sum(ind)
variant_data$Summary[ind] = variant_data$Summary_PCR[ind]

ind = which(variant_data$Summary=='Omicron')
table(variant_data$Lineage[ind])
variant_data$Summary[ind] = variant_data$Lineage[ind]

variant_data$Summary[ind] = plyr::mapvalues(x = variant_data$Summary[ind],
                                            from = c('BA.5.2', 'BA.5.2.1'),
                                            to = c('BA.5','BA.5'))
table(variant_data$Summary, useNA = 'ifany')
####################################################


plot_variant_data=merge(x = variant_data,
                        y = clin_data, by.x = 'ID', by.y = 'Label',all.x = T)

plot_variant_data = dplyr::arrange(plot_variant_data, Rand_date_time)
par(las=1,cex.axis=1.3,cex.lab=1.3)
plot_variant_data$var=as.numeric(factor(plot_variant_data$Summary,
                                        levels=c('Delta','BA.1','BA.2','BA.4','BA.5')))
var_cols = RColorBrewer::brewer.pal(n = length(unique(plot_variant_data$var)), name = 'Set2')
plot(as.POSIXct(plot_variant_data$Rand_date_time), 1:nrow(plot_variant_data), 
     xlab='Randomisation Date', ylab='Genotyped Patient number',
     panel.first=grid(),col = adjustcolor(var_cols[plot_variant_data$var],.5),
     pch = 14+plot_variant_data$var,
     cex=1.5)
legend('bottomright',pch=15:19,col=var_cols,cex=1.5,
       legend = c('Delta','BA.1','BA.2','BA.4','BA.5'),inset=0.03)

abline(v=as.POSIXct('2022-07-01'))
table(variant_data$Summary, useNA = 'ifany')

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




##******************** Log database ***********************
##*********************************************************
##*********************************************************

log_data = haven::read_dta('../Data/InterimSampleLog.dta')
log_data$sl_barc[log_data$sl_barc=='']=NA
log_data = log_data[!is.na(log_data$sl_barc), ]
log_data = log_data[!is.na(log_data$sl_sampdat), ]
log_data = log_data[!is.na(log_data$sl_samptim), ]


##******************** Virus Density PCR database *********
##*********************************************************
##********************************** ***********************
fnames = list.files('../Data/CSV files',full.names = T,recursive = T)
library(readr)

my_specs = cols(
  `SUBJECT ID` = col_character(),
  BARCODE = col_character(),
  INITIALS = col_character(),
  Location = col_character(),
  `TIME-POINT` = col_character(),
  `Time Collected` = col_time(format = ""),
  `COLLECTION DATE` = col_character(),
  VOLUME = col_double(),
  `UNIT (mL)` = col_character(),
  `STORAGE DATE` = col_character(),
  `STORAGE TIME` = col_character(),
  BOX = col_double(),
  POSITION = col_character(),
  NOTE = col_character(),
  `Date Received` = col_character(),
  `Sample ID` = col_character(),
  `N/S Gene` = col_character(),
  RNaseP = col_character(),
  `Target conc. c/mL` = col_number(),
  `Lot no.` = col_character()
)

for(i in 1:length(fnames)){
  writeLines(sprintf('Loading data from file:\n %s \n*********************************************************', fnames[i]))
  if(i==1){
    temp = readr::read_csv(fnames[i])
  } else {
    temp = readr::read_csv(fnames[i],col_types = my_specs)
  }
  # make sure that PCR plates have unique codes across sites
  if(length(grep(pattern = 'Brazil', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste(x,'Brazil',sep='_'))
  }
  
  if(length(grep(pattern = 'Thailand', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste(x,'Thailand',sep='_'))
  }
  if(i==1){
    Res=temp
  } else {
    Res = rbind(Res,temp)
    
  }
  
}

ind_rm = is.na(Res$`SUBJECT ID`) & is.na(Res$`Sample ID`)
Res = Res[!ind_rm, ]
Res$Lab = NA
Res$Lab[grep(pattern = 'Thailand',x = Res$`Lot no.`)]='Thailand'
Res$Lab[grep(pattern = 'Brazil',x = Res$`Lot no.`)]='Brazil'

Res$`Lot no.` = tolower(Res$`Lot no.`)
Res$Plate = as.numeric(as.factor(Res$`Lot no.`))
Res$Lab = as.factor(Res$Lab)


writeLines('\nShowing the number of plates and the number of samples per plate:')
print(table(Res$Plate))
range(table(Res$Plate))
if(max(table(Res$Plate))>98){
  writeLines('**************XXXXXXXXX MORE THAN 98 samples on a single plate!! XXXXXXXX************')
}

writeLines('\nAre there 96 samples per plate?')
print(table(table(Res$Plate)==96))

## Extract standard curve data by plate
ind = grep('std', Res$`Sample ID`)
SC = Res[ind, c('Sample ID','N/S Gene','Target conc. c/mL','Plate','Lab')]

SC$`N/S Gene`[SC$`N/S Gene`=='Undetermined'] = 40
SC$CT_NS = as.numeric(SC$`N/S Gene`)
SC$log10_true_density = log10(as.numeric(SC$`Target conc. c/mL`))
SC$ID = apply(SC[, c('Sample ID','Plate')], 1, function(x) paste(x[1],as.numeric(x[2]), sep='_'))

# Select columns
cols = c('ID','Plate','CT_NS','log10_true_density','Lab')
SC = SC[,cols]

## Extract sample data
Res = Res[!is.na(Res$`SUBJECT ID`), ]

# Res = Res[Res$Location != 'Saliva', ]
table(Res$Location)

writeLines('Clinical data from the following patients not in database:\n')
print(unique(Res$`SUBJECT ID`[!Res$`SUBJECT ID` %in%  clin_data$Label]))
Res = Res[Res$`SUBJECT ID` %in% clin_data$Label, ]

## get missing data from log file
# table(Res$BARCODE %in% log_data$sl_barc)
# xx=Res[!Res$BARCODE %in% log_data$sl_barc, c('BARCODE', 'SUBJECT ID')]
# write.csv(x = xx, file = '~/Downloads/missing_barcodes.csv')

writeLines('Number of samples per patient:')
range(table(Res$`SUBJECT ID`))
hist(table(Res$`SUBJECT ID`))

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
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'SAL',replacement = 'Saliva')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'TLS',replacement = 'TSL')

Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'TSL',replacement = 'Left_tonsil')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'RTS',replacement = 'Right_tonsil')

table(Res$Swab_ID)


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


sampling_time_conflicts = c()
na_sample_times = c()



# manual corrections

log_data$sl_sampdat[log_data$sl_barc=='20SA069'] = '2022-01-23'

log_data$sl_sampdat[log_data$sl_barc=='20QE973'] = '2022-06-22'
log_data$sl_sampdat[log_data$sl_barc=='20QE970'] = '2022-06-22'


log_data$sl_samptim[log_data$sl_barc=='20SA069'] = '09:11:00'

log_data$sl_samptim[log_data$sl_barc=='20RA973'] = '11:03:00'
log_data$sl_samptim[log_data$sl_barc=='20RA979'] = '11:04:00'
log_data$sl_samptim[log_data$sl_barc=='20RA976'] = '11:03:00'
log_data$sl_samptim[log_data$sl_barc=='20RA982'] = '11:04:00'

log_data$sl_barc = toupper(log_data$sl_barc)

Res$`COLLECTION DATE`=gsub(pattern = '/', replacement = '-',x = Res$`COLLECTION DATE`,fixed = T)
Res$`COLLECTION DATE`=gsub(pattern = 'Sept', replacement = 'Sep',x = Res$`COLLECTION DATE`,
                           fixed = T,ignore.case = F)

for(i in 1:nrow(Res)){
  
  id = Res$ID[i]
  ind = which(clin_data$Label==id)
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
  
  
  ## Calculate time since randomisation for the visit
  rand_time = clin_data$Rand_date_time[ind]
  Res$Rand_date[i] = as.character(rand_time)
  
  # PCR data
  if(!any(c(is.na(Res$`COLLECTION DATE`[i]), is.na(Res$`Time Collected`[i])))){
    sample_time = 
      as.character(as.POSIXct(paste(Res$`COLLECTION DATE`[i],
                                    Res$`Time Collected`[i], sep=' '),
                              tryFormats = c('%d-%b-%y %X',
                                             '%d/%b/%Y %X',
                                             '%d-%b-%Y %X',
                                             '%d-%m-%Y %X')))
  } else {
    sample_time=NA
    na_sample_times=c(na_sample_times,i)
  }
  
  # get time of sampling from log file
  barcode_i = Res$BARCODE[i]
  ind_log = which(log_data$sl_barc==barcode_i)
  
  
  if(length(ind_log)>1) {
    writeLines(sprintf('duplicated entries for barcode %s', barcode_i))
    ind_log = ind_log[1]
  }
  # get the sample time from the log file
  if(length(ind_log)>0){
    # check right patient!
    if(!Res$ID[i]==log_data$Label[ind_log]){
      writeLines(sprintf('BARCODE %s is ID %s in sample log and ID %s in PCR data',
                         Res$BARCODE[i], log_data$Label[ind_log], Res$ID[i]))
    }
    sample_time_log = 
      as.character(as.POSIXct(paste(log_data$sl_sampdat[ind_log[1]],
                                    log_data$sl_samptim[ind_log[1]], sep=' ')))
  } else {
    sample_time_log = NA # NA if missing
  }
  
  s_times = c(sample_time, sample_time_log)
  if(all(is.na(s_times))){
    # writeLines(sprintf('No sample time for patient %s at timepoint %s',id, Res$`TIME-POINT`[i]))
    my_sample_time = NA
  } else {
    if(all(!is.na(s_times))){
      if(difftime(sample_time_log, sample_time,units = 'hours')>2){
        writeLines(sprintf('Conflicting sample times: ID %s, Timepoint %s, barcode %s: difference of %s hours',
                           id,
                           Res$Timepoint_ID[i],
                           barcode_i, 
                           round(difftime(sample_time_log, sample_time, units = 'hours'))))
        sampling_time_conflicts=c(sampling_time_conflicts,barcode_i)
        my_sample_time = s_times[2] # trust the log as updated by Padd
      } else {
        my_sample_time = s_times[1]
      }
    } else {
      my_sample_time = s_times[which(!is.na(s_times))]
    }
  }
  
  Res$Time[i] = difftime(my_sample_time,rand_time,units = 'days')
  if(is.na(Res$Time[i])) {
    if(!grepl(pattern = 'TH57', x = Res$`SUBJECT ID`[i])&
       !grepl(pattern = 'TH58', x = Res$`SUBJECT ID`[i])&
       Res$`SUBJECT ID`[i] != 'PLT-TH1-205' & 
       Res$`SUBJECT ID`[i] != 'PLT-TH1-208'){
      writeLines(sprintf('Missing sample time for patient %s at timepoint %s',id,Res$Timepoint_ID[i]))
    }
  } else {
    if(Res$Time[i] < -.1 & 
       !(Res$Site[i] %in% c('th057', 'th058'))){
      writeLines(sprintf('Negative sample time for patient %s at timepoint %s: %s days',
                         id,Res$Timepoint_ID[i],
                         round(Res$Time[i],1)))
    }
  }
  
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

writeLines(sprintf('Missing sample times for %s',unique(Res$ID[na_sample_times])))

ind_na_time = is.na(Res$Time)
writeLines(sprintf('Missing time for following samples: %s',
                   unique(Res$ID[ind_na_time])))
# use protocol time instead of actual time
Res$Time[ind_na_time]=Res$Timepoint_ID[ind_na_time]

ind_neg_time = Res$Time < 0
writeLines(sprintf('Negative time for following samples: %s',
                   unique(Res$ID[ind_neg_time])))
# use protocol time
Res$Time[ind_neg_time]=Res$Timepoint_ID[ind_neg_time]



Res = dplyr::arrange(Res, Rand_date, ID, Time)

writeLines(sprintf('there are a total of %s patients in the PCR database',
                   length(unique(Res$`SUBJECT ID`))))

variant_data$Variant = variant_data$Summary
Res = merge(Res, variant_data[, c('ID','Variant')], by = 'ID', all.x = T)

## Add genotyping data
ind_missing_variant = is.na(Res$Variant)
Res$Variant_Imputed=0 #1: imputed; 0: genotyped
Res$Variant_Imputed[ind_missing_variant]=1

# Imputation based on date if not yet typed
d1 = as.POSIXct('2022-01-01')
d2 = as.POSIXct('2022-02-20')
d3 = as.POSIXct('2022-07-01')
ind_Delta = Res$Rand_date < d1
ind_BA1 = Res$Rand_date >= d1 & Res$Rand_date < d2
ind_BA2 = Res$Rand_date >= d2 & Res$Rand_date < d3
ind_BA5 = Res$Rand_date >= d3
Res$Variant[ind_missing_variant&ind_Delta] = 'Delta'
Res$Variant[ind_missing_variant&ind_BA1] = 'BA.1'
Res$Variant[ind_missing_variant&ind_BA2] = 'BA.2'
Res$Variant[ind_missing_variant&ind_BA5] = 'BA.5'


Res$Epoch = 0
Res$Epoch[Res$Rand_date > as.POSIXct('2022-04-01')] = 1 # stopped ivermectin
Res$Epoch[Res$Rand_date > as.POSIXct('2022-04-18')] = 2 # added fluoxetine
Res$Epoch[Res$Rand_date > as.POSIXct('2022-06-10')] = 3 # stopped remdesivir
table(Res$Epoch[!duplicated(Res$`SUBJECT ID`)], useNA = 'ifany')

##***********************************************
cols = c('ID','Time','Trt','Site','Timepoint_ID',
         'Swab_ID','Rand_date','Any_dose','N_dose',
         'Antibody_test','Weight','BMI','Plate','BARCODE',
         'Age', 'Sex', 'Symptom_onset','Variant','Variant_Imputed',
         'CT_NS','CT_RNaseP','Epoch', 'Per_protocol_sample','Lab')

writeLines('\n column names:')
print(cols)
Res = Res[, cols]
Res$censor = ifelse(Res$CT_NS==40, 'left', 'none')

SC = dplyr::arrange(SC, Lab, Plate, ID)


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
plot(control_dat$CT, jitter(control_dat$log10_true_density), xlim=c(20,40))
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
writeLines(sprintf('out of a total of %s samples, %s are below LLOQ (%s%%)',
                   nrow(Res),
                   sum(Res$CT_NS==40),
                   round(100*mean(Res$CT_NS==40),1)))

table(Res$log10_viral_load == Res$log10_cens_vl)

writeLines('The follwoing IDs have duplicated barcodes:')

# bs_dup = Res$BARCODE[duplicated(Res$BARCODE)]
# View(Res[Res$BARCODE %in% bs_dup, ])
print(unique(Res$ID[duplicated(Res$BARCODE)]))
Res = Res[!duplicated(Res$BARCODE), ]

###### Write csv files
# Overall data files
write.table(x = SC, file = 'interim_control_dat.csv', row.names = F, sep=',')

Res = dplyr::arrange(Res, Rand_date, ID, Time)
write.table(x = Res, file = 'interim_dat.csv', row.names = F, sep=',')

write.table(x = screen_failure, file = 'interim_screening_dat.csv', row.names = F, sep=',')

# remove the saliva samples
Res = Res[Res$Swab_ID != 'Saliva', ]

ind_Thailand = Res$Site %in% c('th001','th057','th058')
ind_Brazil= Res$Site %in% c('br003')

#***********************************************************************#
#************************ Specific analysis data files *****************#

#************************* Ivermectin Analysis *************************#
#* Only Thailand
IDs_Ivermectin = union(unique(Res$ID[Res$Trt %in% 'Ivermectin' & ind_Thailand]),
                       unique(Res$ID[Res$Trt %in%  "No study drug" &
                                       Res$Rand_date < "2022-04-22 00:00:00" & ind_Thailand]))
IDs_pos_control = unique(Res$ID[Res$Rand_date < "2021-12-17 00:00:00" & 
                                  Res$Trt == 'Regeneron' & 
                                  Res$Site == 'th001']) # first 10 randomised to REGN in FTM
Res_ivermectin = Res[Res$ID %in% c(IDs_Ivermectin, IDs_pos_control), ]
pk_summaries = read.csv('../PLATCOV-Ivermectin/PK_summaries.csv')
Res_ivermectin = merge(Res_ivermectin, pk_summaries, by = 'ID', all.x = T)
Res_ivermectin = dplyr::arrange(Res_ivermectin, Rand_date, ID, Time)

write.table(x = Res_ivermectin, file = 'Ivermectin_analysis.csv', sep = ',',row.names = F,quote = F)
# anonymise IDs for publication
Res_ivermectin$ID = as.numeric(as.factor(Res_ivermectin$ID))

write.table(x = Res_ivermectin, file = '../PLATCOV-Ivermectin/Ivermectin_analysis.csv', sep = ',',row.names = F,quote = F)


#************************* Remdesivir Analysis *************************#
#* Thailand and Brazil
IDs_Remdesivir = unique(Res$ID[Res$Trt %in% c('Remdesivir',"No study drug") &
                                 Res$Rand_date < "2022-06-11 00:00:00"])
Res_Remdesivir = Res[Res$ID %in% IDs_Remdesivir, ]
Res_Remdesivir = dplyr::arrange(Res_Remdesivir, Rand_date, ID, Time)

write.table(x = Res_Remdesivir, file = 'Remdesivir_analysis.csv', row.names = F, sep=',',quote = F)


#************************* Favipiravir Analysis *************************#
#* Thailand and Brazil
IDs_Favipiravir = unique(Res$ID[Res$Trt %in% c('Favipiravir',"No study drug")])
Res_Favipiravir = Res[Res$ID %in% IDs_Favipiravir, ]
Res_Favipiravir = dplyr::arrange(Res_Favipiravir, Rand_date, ID, Time)
write.table(x = Res_Favipiravir, file = 'Favipiravir_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Regeneron Analysis *************************#
#* Thailand only
IDs_REGN = unique(Res$ID[Res$Trt %in% c('Regeneron',"No study drug") &
                           ind_Thailand &
                           Res$Rand_date < "2022-08-25 00:00:00"])
Res_REGN = Res[Res$ID %in% IDs_REGN, ]
Res_REGN = dplyr::arrange(Res_REGN, Rand_date, ID, Time)
write.table(x = Res_REGN, file = 'REGN_analysis.csv', row.names = F, sep=',', quote = F)

#************************* Regeneron + Remdesvir Analysis *************************#
Res_Paper2 = Res[Res$ID %in% union(IDs_REGN,IDs_Remdesivir), ]
Res_Paper2 = dplyr::arrange(Res_Paper2, Rand_date, ID, Time)
write.table(x = Res_Paper2, file = 'Paper2_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Fluoxetine Analysis *************************#
#* Thailand - first and then Brazil in June - need to check
IDs_Fluoxetine = unique(Res$ID[Res$Trt %in% c('Fluoxetine',"No study drug") &
                                 Res$Rand_date > "2022-04-01 00:00:00"])
Res_Fluoxetine = Res[Res$ID %in% IDs_Fluoxetine, ]
Res_Fluoxetine = dplyr::arrange(Res_Fluoxetine, Rand_date, ID, Time)
write.table(x = Res_Fluoxetine, file = 'Fluoxetine_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Molnupiravir Analysis *************************#
#* Thailand only
IDs_Molnupiravir = unique(Res$ID[Res$Trt %in% c('Molnupiravir',"No study drug") &
                                   Res$Rand_date > "2022-06-03 00:00:00" & 
                                   ind_Thailand])
Res_Molnupiravir = Res[Res$ID %in% IDs_Molnupiravir, ]
Res_Molnupiravir = dplyr::arrange(Res_Molnupiravir, Rand_date, ID, Time)
write.table(x = Res_Molnupiravir, file = 'Molnupiravir_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Paxlovid Analysis *************************#
#* Thailand only
IDs_Paxlovid = unique(Res$ID[Res$Trt %in% c('Nirmatrelvir + Ritonavir',"No study drug") &
                               Res$Rand_date > "2022-06-03 00:00:00" &
                               ind_Thailand])
Res_Paxlovid = Res[Res$ID %in% IDs_Paxlovid, ]
Res_Paxlovid = dplyr::arrange(Res_Paxlovid, Rand_date, ID, Time)
write.table(x = Res_Paxlovid, file = 'Paxlovid_analysis.csv', row.names = F, sep=',', quote = F)

#************************* Paxlovid v Molnupiravir Analysis *************************#
#* Thailand only
Res_Paxlovid_Molnupiravir = Res[Res$ID %in% union(IDs_Paxlovid,IDs_Molnupiravir), ]
Res_Paxlovid_Molnupiravir = dplyr::arrange(Res_Paxlovid_Molnupiravir, Rand_date, ID, Time)
write.table(x = Res_Paxlovid_Molnupiravir, file = 'Paxlovid_Molnupiravir_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Nitazoxanide Analysis *************************#
#* Brazil only
IDs_Nitazoxanide = unique(Res$ID[Res$Trt %in% c('Nitazoxanide',"No study drug") &
                                   ind_Brazil])
Res_Nitazoxanide = Res[Res$ID %in% IDs_Nitazoxanide, ]
Res_Nitazoxanide = dplyr::arrange(Res_Nitazoxanide, Rand_date, ID, Time)
write.table(x = Res_Nitazoxanide, file = 'Nitazoxanide_analysis.csv', row.names = F, sep=',', quote = F)


save(IDs_Remdesivir,
     IDs_Paxlovid, 
     IDs_Molnupiravir,
     IDs_Fluoxetine, 
     IDs_REGN, 
     IDs_Favipiravir,
     IDs_Nitazoxanide,
     IDs_Ivermectin, file = 'ID_list.RData')


writeLines('Barcodes of conflicting sample times:')
print(sampling_time_conflicts)
rm(list=ls())
