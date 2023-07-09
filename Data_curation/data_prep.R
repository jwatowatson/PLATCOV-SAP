##############################################################################################
# PLATCOV project
# This script prepares data for further analyses
##############################################################################################
library(dplyr)
library(readxl)
library(readr)
library(ggplot2)
library(lme4)
library(lubridate)
##Define user folder path####################################################################
user <- 'james'#"Chang" # Change here

#1 Analysis_data folder
if(user == "Chang"){
  prefix_analysis_data <- "PLATCOV_SAP"
}else{
    prefix_analysis_data <- ".."}

#2 Dropbox folder (PLATCOV_Analysis)
if(user == "Chang"){
  prefix_dropbox <- "C:/Users/Phrutsamon/Dropbox/PLATCOV_Analysis"
}else{
  prefix_dropbox <- "~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/"
}

#3 Downloads folder
if(user == "Chang"){
  prefix_downloads <- "PLATCOV_SAP"
}else{
  prefix_downloads <- "~"
}

#4 Dropbox folder for randomisation
if(user == "Chang"){
  prefix_drop_rand <- "C:/Users/Phrutsamon/Dropbox/PLATCOV" 
}else{
  prefix_drop_rand <- "~/Dropbox/PLATCOV"
}

#5 Data_curation folder
if(user == "Chang"){
  prefix_dat_cur <- "PLATCOV-SAP/Data_curation/" 
}else{
  prefix_dat_cur <- ""
}
#############################################################################################
##******************** Clinical database *******************
##*********************************************************
##*********************************************************
##*
# --- Fever ---
fever_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/temperature_data.csv")) #A dropbox folder shared by James
fever_data = fever_data %>% distinct(Label, .keep_all = T)
# --- Clinical data (symptom onsets) ---
clin_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimEnrolment.dta"))
clin_data = merge(clin_data, fever_data[, c('Label','Fever_Baseline')], by='Label', all = T)
# --- Final status (completed the trials?) ---
final_status = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimFinalStatus.dta"))
final_status = final_status[!is.na(final_status$fs_compyn), ]
# --- Adverse effects ---
AE_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimAE.dta"))
# --- Extract screening failure data ---
table(clin_data$scrpassed, useNA = 'ifany')
ind = !is.na(clin_data$scrpassed) & clin_data$scrpassed==0
screen_failure =
  clin_data[ind,
            c("Trial","Site","scrid","scrdat",
              "scrpassed","reason_failure","scrnote")]

screen_failure$reason_failure = sjlabelled::as_character(screen_failure$reason_failure)
write.csv(x = screen_failure, file = paste0(prefix_downloads, "/Downloads/screening_failures.csv"))
# --- Excluding screening failure data ---
clin_data = clin_data[!ind, ]
sort(unique(clin_data$Label))

clin_data$rangrp = sjlabelled::as_character(clin_data$rangrp) #Randomisation group labels
AE_data = merge(AE_data, clin_data[, c('Label','rangrp')])
# --- check data for missing values ---
#1. Day
ind = !is.na(clin_data$cov_symphr) & is.na(clin_data$cov_sympday)
if(sum(ind)>0) clin_data$cov_sympday[ind] = clin_data$cov_symphr[ind]/24 #If day is missing >>> divide hours by 24
#2. Symptomatic day
ind = is.na(clin_data$cov_sympday)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing time since symptom onset', clin_data$Label[ind]))
  clin_data$cov_sympday[is.na(clin_data$cov_sympday)]=2 # Assign Symptomatic day = 2 for missing data
}
#3. Age in year
clin_data$age_yr[clin_data$Label=='PLT-TH57-002'] = 21 # A special case
ind = is.na(clin_data$age_yr) & !is.na(clin_data$dob_my)
# calculate age at randomisation 
for(i in which(ind)){
  clin_data$age_yr[i] = trunc((as.POSIXct(clin_data$dob_my[i], format='%d/%m/%Y') %--% as.POSIXct(clin_data$randat[i])) / years(1))
}

ind = is.na(clin_data$age_yr)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing age', clin_data$Label[ind]))
}
#4. BMI
clin_data$BMI = clin_data$weight/(clin_data$height/100)^2
clin_data$Weight = clin_data$weight

ind = is.na(clin_data$BMI)
if(sum(ind)>0) {
  writeLines(sprintf('Patient %s has missing BMI', clin_data$Label[ind]))
}

# --- Cross check with online randomisation app data ---
# NOTE: Sites TH57 and TH58 did not use app so cannot cross check
rand_app_data = rbind(read.csv(paste0(prefix_drop_rand, "/data-TH1.csv")),
                      read.csv(paste0(prefix_drop_rand, "/data-BR3.csv")),
                      read.csv(paste0(prefix_drop_rand, "/data-LA08.csv"))) %>%
  filter(!is.na(Treatment))

rand_app_data$ID = paste0('PLT-', rand_app_data$site,'-',
                          stringr::str_pad(rand_app_data$randomizationID, 3, pad = "0"))
rand_app_data$Site = plyr::mapvalues(x = rand_app_data$site, from=c('TH1','BR3'),to=c('th001','br003'))

writeLines('The following randomisation database IDs are not in clinical database:\n')
print(rand_app_data$ID[!rand_app_data$ID %in% clin_data$Label])

# rand_app_data = rand_app_data[rand_app_data$ID%in% clin_data$Label, ]
rand_app_data$Rand_Time = as.POSIXct(rand_app_data$Date, format = '%a %b %d %H:%M:%S %Y',tz = 'GMT')
rand_app_data$tzone = plyr::mapvalues(x = rand_app_data$site,
                                      from = c('TH1','LA08','BR3'),
                                      to = c('Asia/Bangkok','Asia/Bangkok','America/Sao_Paulo'))
rand_app_data$Rand_Time_TZ=NA
for(i in 1:nrow(rand_app_data)){
  rand_app_data$Rand_Time_TZ[i] = as.character(with_tz(rand_app_data$Rand_Time[i], tzone = rand_app_data$tzone[i]))
}

clin_data = merge(clin_data, rand_app_data,
                  by.x = c('Label','Site'),
                  by.y = c("ID",'Site'), all = T)

# --- discrepancies between randomisation database and CRFs? ---
clin_data$sex_char = sjlabelled::as_character(clin_data$sex.x)
clin_data$Rand_date_time = NA
for(i in 1:nrow(clin_data)){
  if(!is.na(clin_data$randat[i]) & !is.na(clin_data$rantim[i])){
    clin_data$Rand_date_time[i] = as.character(as.POSIXct(paste(clin_data$randat[i],
                                                                clin_data$rantim[i], sep=' ')))
  }
}

table(clin_data$rangrp)
clin_data$rangrp[clin_data$rangrp=='Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir'
clin_data$rangrp[clin_data$rangrp=='Molnupiravir and Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir + Molnupiravir'

ind = is.na(clin_data$Treatment) & !is.na(clin_data$rangrp)
clin_data$Treatment[ind] = clin_data$rangrp[ind]

ind = !is.na(clin_data$Treatment) & is.na(clin_data$rangrp)
clin_data$rangrp[ind] = clin_data$Treatment[ind]

# Reporting inconsistency between randomisation database (shiny app; Dropbox) and clinical database (CRF)
if(any(! clin_data$Treatment == clin_data$rangrp)) {
  writeLines(sprintf('Randomisation inconsistent for %s', 
                     clin_data$Label[clin_data$Treatment != clin_data$rangrp]))
  ind_diff = clin_data$Treatment != clin_data$rangrp
  View(clin_data[ind_diff,c('Label','rangrp','Treatment')])
  clin_data$rangrp[ind_diff] = clin_data$Treatment[ind_diff]
}

#######################
# variables we need are:
# * age
# * sex
# * time since symptom onset
# * number of vaccine doses

####### Age #########
ind = is.na(clin_data$age_yr) & !is.na(clin_data$age)
clin_data$age_yr[ind] = clin_data$age[ind]
ind = !is.na(clin_data$age_yr) & is.na(clin_data$age)
clin_data$age[ind] = clin_data$age_yr[ind]
if(any(is.na(clin_data$age))){
  writeLines(sprintf('Age missing for %s', clin_data$Label[is.na(clin_data$age)]))
}

ind = !is.na(clin_data$age) & (!clin_data$age == clin_data$age_yr)
if(any(ind)) {
  writeLines(sprintf('Age inconsistent for %s', clin_data$Label[ind]))
}
print(clin_data[clin_data$age != clin_data$age_yr, c('Label','age','age_yr')])
colnames(clin_data)[which(names(clin_data) == 'age_yr')] <- 'Age'

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
clin_data$Sex = as.numeric(clin_data$sex_char=='Male')

####### Rand Time #########
ind = is.na(clin_data$Rand_date_time) & !is.na(clin_data$Rand_Time_TZ)
clin_data$Rand_date_time[ind] = clin_data$Rand_Time_TZ[ind]
ind = !is.na(clin_data$Rand_date_time) & is.na(clin_data$Rand_Time_TZ)
clin_data$Rand_Time_TZ[ind] = clin_data$Rand_date_time[ind]
Rand_diffs = apply(clin_data[,c('Rand_date_time','Rand_Time_TZ')],1,
                   function(x) difftime(x[1], x[2], units='mins'))
if(any(abs(Rand_diffs)>5)){
  writeLines(sprintf('More than 5 min difference in rand time for %s',
                     clin_data$Label[which(abs(Rand_diffs)>5)]))
}
print(clin_data[which(abs(Rand_diffs)>5), c('Label','Rand_Time_TZ','Rand_date_time') ])
clin_data$Rand_date_time[which(abs(Rand_diffs)>5)] = 
  clin_data$Rand_Time_TZ[which(abs(Rand_diffs)>5)]

clin_data$Symptom_onset = clin_data$cov_sympday
clin_data$Trt = clin_data$rangrp

clin_data = clin_data[, c('Label','Trt','Sex','Age','randat',
                          "Rand_date_time",'BMI','Weight',
                          'Symptom_onset','Site','Fever_Baseline')]

## Per protocol for treatment data
# trt_distcont_data = haven::read_dta('../Data/InterimChangeTreatment.dta')
trt_distcont_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimDrugRescue.dta")) %>% filter(!is.na(dardat))


# ### Variant data
# ## PCR variant data (up until July 2022)
# fnames_var = list.files('~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/Data/PCR genotyping/',pattern = 'variant genotyping',full.names = T)
# dat = lapply(fnames_var,read_excel)
# for(i in 1:length(dat)) dat[[i]] = dat[[i]][, c("SUBJECT ID","Summary")]
# dat = do.call(what = rbind, dat)
# dat = dat[!is.na(dat$Summary),]
# var_data = dat[!duplicated(dat$`SUBJECT ID`), ]
# colnames(var_data)=c('ID', 'Summary_PCR')
# 
# var_data$Summary_PCR = plyr::mapvalues(x = var_data$Summary_PCR,
#                                        from = c('Delta_B1.617.2','Omicron BA.2_B.1.1.529',
#                                                 'Omicron_B.1.1.529','undetermined'),
#                                        to = c('Delta','BA.2','BA.1',NA))
# writeLines(sprintf('We have PCR variant genotyping for %s patients',nrow(var_data)))

## Nanopore data

source(paste0(prefix_dat_cur, "get_nanopore_data.R"))
variant_data = get_nanopore_data(prefix_dropbox = prefix_dropbox)
variant_data = merge(variant_data, clin_data, by.x='ID', by.y = 'Label')
ggplot(variant_data, aes(Rand_date_time, after_stat(count), group=Variant, fill = Variant)) +
  geom_density(position = "fill")

##******************** Vaccine database *******************
##*********************************************************
##*********************************************************
vacc_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimVaccine.dta"))

writeLines(sprintf('No vaccine data for %s',
                   clin_data$Label[!clin_data$Label %in% vacc_data$Label]))

## Some cleaning
# vacc_data$vc_name = tolower(vacc_data$vc_name)
writeLines('\nVaccine names before cleaning:')
print(table(vacc_data$vc_name))

clin_data$Any_dose = NA
clin_data$N_dose = NA
clin_data$Time_since_last_dose = NA
clin_data$N_dose_mRNA = NA
clin_data$Any_dose_mRNA = NA

vacc_date_cols = grep('dos', colnames(vacc_data))
## find any bad dates
all_dates = unlist(unique(vacc_data[, vacc_date_cols]))
all_dates = all_dates[all_dates!='']
bad_dates = all_dates[which(is.na(parse_date_time(all_dates, orders = 'dmy')))]
vacc_data$Label[which(apply(vacc_data[, vacc_date_cols], 1, function(x) length(intersect(x ,bad_dates))>0))]

for(i in 1:nrow(clin_data)){
  id = clin_data$Label[i]
  ind = which(vacc_data$Label==id)
  ind_mRNA = which(vacc_data$Label==id & vacc_data$vc_name %in% c('Moderna','Pfizer'))
  
  vac_dates = unlist(unique(vacc_data[ind, vacc_date_cols]))
  vac_dates = vac_dates[vac_dates!='']
  
  vac_dates_mRNA = unlist(unique(vacc_data[ind_mRNA, vacc_date_cols]))
  vac_dates_mRNA = vac_dates_mRNA[vac_dates_mRNA != '']
  clin_data$N_dose[i] = length(vac_dates)
  clin_data$Any_dose[i] = c('No','Yes')[1+as.numeric(clin_data$N_dose[i]>0)]
  
  clin_data$N_dose_mRNA[i] = length(vac_dates_mRNA)
  clin_data$Any_dose_mRNA[i] = c('No','Yes')[1+as.numeric(clin_data$N_dose_mRNA[i]>0)]
  
  if(clin_data$Any_dose[i]=='Yes'){
    most_recent_vac = max(parse_date_time(vac_dates, orders = 'dmy'))
    clin_data$Time_since_last_dose[i] =  
      difftime(clin_data$Rand_date_time[i],
               most_recent_vac,units = 'days')
  }
  
}




##******************** Log database ***********************
##*********************************************************
##*********************************************************
log_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimSampleLog.dta"))
log_data$sl_barc[log_data$sl_barc=='']=NA
log_data = log_data[!is.na(log_data$sl_barc), ]
log_data = log_data[!is.na(log_data$sl_sampdat), ]
log_data = log_data[!is.na(log_data$sl_samptim), ]
log_data$sl_barc = tolower(log_data$sl_barc)


##******************** Virus Density PCR database *********
##*********************************************************
##********************************** ***********************
fnames = list.files(paste0(prefix_dropbox, "/Data/CSV files"),full.names = T,recursive = T)

# column type specification for each csv file
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
  # writeLines(sprintf('Loading data from file:\n %s \n*********************************************************', fnames[i]))
  
  # read in file
  temp = readr::read_csv(fnames[i],col_types = my_specs)
  
  # check for duplicates
  if(any(duplicated(temp$BARCODE[!is.na(temp$BARCODE)]))){
    writeLines(sprintf('in file %s there are duplicate barcodes',fnames[i]))
  }
  # make sure that PCR plates have unique codes across sites
  if(length(grep(pattern = 'Brazil', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Brazil',x,sep='_'))
  }
  
  if(length(grep(pattern = 'Thailand', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Thailand',x,sep='_'))
  }
  if(i==1){
    Res=temp
  } else {
    Res = rbind(Res,temp)
  }
}

# take out the summary PCR columns
ind_rm = Res$`Sample ID` %in% c('PC','R-squared','Efficiency (%)','Slope (M)')|
  (is.na(Res$`SUBJECT ID`) & is.na(Res$`Sample ID`))
sum(ind_rm)
Res = Res[!ind_rm, ]

# take out duplicate rows
ind_rm = !is.na(Res$BARCODE) & duplicated(Res$BARCODE)
sum(ind_rm)
Res = Res[!ind_rm, ]

sort(unique(Res$`Lot no.`))
Res$`Lot no.`[Res$`Lot no.`=='Thailand_D10 lot1']="Thailand_D10 Lot 1"
Res$`Lot no.`[Res$`Lot no.`=='Thailand_D10 lot2']="Thailand_D10 Lot 2"

# make the plate/lab variables
Res$Lab = NA
Res$Lab[grep(pattern = 'Thailand',x = Res$`Lot no.`)]='Thailand'
Res$Lab[grep(pattern = 'Brazil',x = Res$`Lot no.`)]='Brazil'

Res$`Lot no.` = tolower(Res$`Lot no.`)
sort(unique(Res$`Lot no.`))
Res$Plate = as.numeric(as.factor(Res$`Lot no.`))
Res$Lab = as.factor(Res$Lab)

## Missing samples
ind_missing = grep('not',Res$`Sample ID`,ignore.case = T)
writeLines('Missing data for the following samples:')
Res$`TIME-POINT`[ind_missing]
Res = Res[-ind_missing, ]

writeLines('\nShowing the number of plates and the number of samples per plate:')
xx = table(Res$`Lot no.`)
print(sort(xx))
range(table(Res$Plate))

range(table(Res$`Lot no.`[grep(pattern = 'std', x = Res$`Sample ID`)]))
table(Res$`Lot no.`[grep(pattern = 'std', x = Res$`Sample ID`)])
if(max(table(Res$Plate))>96){
  writeLines('**************XXXXXXXXX MORE THAN 96 samples on a single plate!! XXXXXXXX************')
}

## Extract standard curve data by plate
ind = grep('std', Res$`Sample ID`)
SC = Res[ind, c('Sample ID','N/S Gene','Target conc. c/mL','Plate','Lab','Lot no.')]

SC$`N/S Gene`[SC$`N/S Gene`=='Undetermined'] = 40
SC$CT_NS = as.numeric(SC$`N/S Gene`)
SC$log10_true_density = log10(as.numeric(SC$`Target conc. c/mL`))
SC$ID = apply(SC[, c('Sample ID','Plate')], 1, function(x) paste(x[1],as.numeric(x[2]), sep='_'))

# Select columns
cols = c('ID','Plate','CT_NS','log10_true_density','Lab','Lot no.')
SC = SC[,cols]

## Extract sample data - this was for Liz Batty
# Res = Res[!is.na(Res$`SUBJECT ID`), ]
# write_csv(x = Res[, c("SUBJECT ID","BARCODE","Location","TIME-POINT","Time Collected","COLLECTION DATE")])
# Res = Res[Res$Location != 'Saliva', ]
table(Res$Location, useNA = 'ifany')

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

table(Res$Location)
# Res$Swab_ID = gsub(Res$Swab_ID, pattern = '1',replacement = '')
# Res$Swab_ID = gsub(Res$Swab_ID, pattern = '2',replacement = '')
Res$Swab_ID[grep(Res$Swab_ID, pattern = 'SAL')] = 'Saliva'
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'TLS',replacement = 'TSL')

Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'TSL1',replacement = 'Left_tonsil_1')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'TSL2',replacement = 'Left_tonsil_2')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'TSL',replacement = 'Left_tonsil')

Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'RTS1',replacement = 'Right_tonsil_1')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'RTS2',replacement = 'Right_tonsil_2')
Res$Swab_ID = gsub(Res$Swab_ID, pattern = 'RTS',replacement = 'Right_tonsil')

table(Res$Swab_ID, useNA = 'ifany')


Res$Timepoint_ID = Res$`TIME-POINT`
table(Res$Timepoint_ID, useNA = 'ifany')
Res$Timepoint_ID[Res$Timepoint_ID=='D0H0']='0'
Res$Timepoint_ID[Res$Timepoint_ID=='D0PRE']='0'
Res$Timepoint_ID = gsub(Res$Timepoint_ID,pattern = 'D',replacement = '',fixed = T)
Res$Timepoint_ID = gsub(Res$Timepoint_ID,pattern = 'Tx',replacement = '',fixed = T)
Res$Timepoint_ID = as.numeric(Res$Timepoint_ID)
table(Res$Timepoint_ID, useNA = 'ifany')

include_cols = c('Label','Site','Rand_date_time',
                 'Trt','Age','BMI',
                 'Weight','Sex','Symptom_onset','Fever_Baseline',
                 'Any_dose','N_dose','Time_since_last_dose',
                 'Any_dose_mRNA','N_dose_mRNA')
Res = merge(Res, clin_data[,include_cols], all.x = T, by.x = 'ID', by.y = 'Label')


Res$Time = NA
Res$Per_protocol_all = NA
Res$Per_protocol_sample = NA

Res$ID_log=Res$Time_log=NA

sampling_time_conflicts = c()
na_sample_times = c()



# manual corrections

log_data$sl_sampdat[log_data$sl_barc=='20sa069'] = '2022-01-23'
log_data$sl_samptim[log_data$sl_barc=='20sa069'] = '09:11:00'
log_data$sl_samptim[log_data$sl_barc=='20ra973'] = '11:03:00'
log_data$sl_samptim[log_data$sl_barc=='20ra979'] = '11:04:00'
log_data$sl_samptim[log_data$sl_barc=='20ra976'] = '11:03:00'
log_data$sl_samptim[log_data$sl_barc=='20ra982'] = '11:04:00'

log_data$sl_barc = toupper(log_data$sl_barc)

Res$`COLLECTION DATE`=gsub(pattern = '/', replacement = '-',x = Res$`COLLECTION DATE`,fixed = T)
Res$`COLLECTION DATE`=gsub(pattern = 'Sept', replacement = 'Sep',x = Res$`COLLECTION DATE`,
                           fixed = T,ignore.case = F)

for(i in 1:nrow(Res)){
  
  id = Res$ID[i]
  ## time since randomisation for the visit
  rand_time = Res$Rand_date_time[i]
  
  # PCR data
  if(!any(c(is.na(Res$`COLLECTION DATE`[i]), is.na(Res$`Time Collected`[i])))){
    sample_time =
      as.character(as.POSIXct(paste(Res$`COLLECTION DATE`[i],
                                    Res$`Time Collected`[i], sep=' '),
                              tryFormats = c('%d-%b-%y %X',
                                             '%d/%b/%Y %X',
                                             '%d-%b-%Y %X',
                                             '%d-%m-%y %X',
                                             '%d-%m-%Y %X'
                              )))
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
    Res$ID_log[i]=log_data$Label[ind_log]
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
  Res$Time_log[i] = difftime(s_times[2],Res$Rand_date[Res$ID==Res$ID_log[i]][1],units = 'days')
  
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
    if(sample_time >= as.POSIXct(trt_distcont_data$dardat[ind_discont],
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



Res = dplyr::arrange(Res, Rand_date_time, ID, Time)

writeLines(sprintf('there are a total of %s patients in the PCR database',
                   length(unique(Res$`SUBJECT ID`))))

Res = merge(Res, variant_data[,c('ID','Variant')], by = 'ID', all.x = T)

## Add genotyping data
ind_missing_variant = is.na(Res$Variant) | Res$Variant=='none'
Res$Variant_Imputed=0 #1: imputed; 0: genotyped
Res$Variant_Imputed[ind_missing_variant]=1

# Imputation based on date if not yet typed
d1 = as.POSIXct('2022-01-01')
d2 = as.POSIXct('2022-02-20')
d3 = as.POSIXct('2022-07-01')
d4 = as.POSIXct('2022-11-01')
d5 = as.POSIXct('2023-01-01')

ind_Delta = Res$Rand_date < d1
ind_BA1 = Res$Rand_date >= d1 & Res$Rand_date < d2
ind_BA2 = (Res$Rand_date >= d2 & Res$Rand_date < d3) 
ind_BA5 = Res$Rand_date >= d3 & Res$Rand_date < d4
ind_BA2.75 = Res$Rand_date >= d4 & Res$Rand_date <d5
ind_xbb = Res$Rand_date >= d5 

Res$Variant[ind_missing_variant&ind_Delta] = 'Delta'
Res$Variant[ind_missing_variant&ind_BA1] = 'BA.1'
Res$Variant[ind_missing_variant&ind_BA2] = 'BA.2'
Res$Variant[ind_missing_variant&ind_BA5] = 'BA.5'
Res$Variant[ind_missing_variant&ind_BA2.75] = 'BA.2.75'
Res$Variant[ind_missing_variant&ind_xbb] = 'XBB'


Res$Epoch = 0
Res$Epoch[Res$Rand_date > as.POSIXct('2022-04-01')] = 1 # stopped ivermectin
Res$Epoch[Res$Rand_date > as.POSIXct('2022-04-18')] = 2 # added fluoxetine
Res$Epoch[Res$Rand_date > as.POSIXct('2022-06-10')] = 3 # stopped remdesivir
Res$Epoch[Res$Rand_date > as.POSIXct('2022-10-31')] = 4 # stopped favipiravir
Res$Epoch[Res$Rand_date > as.POSIXct('2023-02-13')] = 5 # stopped molnupiravir
table(Res$Epoch[!duplicated(Res$`SUBJECT ID`)], useNA = 'ifany')

# manually correct error
Res$Swab_ID[Res$BARCODE=='20LH895'] = 'Right_tonsil_1'
Res$Swab_ID[Res$BARCODE=='20RR176'] = 'Left_tonsil_1'


Res$Rand_date = format.Date(x = Res$Rand_date_time, format='%Y-%m-%d')
##***********************************************
cols = c('ID','Time','Trt','Site','Timepoint_ID',
         'Swab_ID','Rand_date','Any_dose','N_dose','Time_since_last_dose',
         'Any_dose_mRNA','N_dose_mRNA',
         'Weight','BMI','Plate','Fever_Baseline','BARCODE',
         'Age', 'Sex', 'Symptom_onset','Variant','Variant_Imputed',
         'CT_NS','CT_RNaseP','Epoch', 'Per_protocol_sample','Lab', 'Lot no.')

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


# random slope and intercept
conv_mod = lmer(log10_true_density ~ 1 + CT + (1+CT|batch),
                data = control_dat,
                control = lmerControl(optimizer ="Nelder_Mead"))

preds = predict(conv_mod)
plot(control_dat$CT, jitter(control_dat$log10_true_density), xlim=c(20,40),
     col = control_dat$Lab)
for(bb in levels(control_dat$batch)){
  ind = control_dat$batch==bb
  lines(control_dat$CT[ind], preds[ind], col = as.numeric(control_dat$Lab[ind]=='Thailand')+1)
}

preds_all =
  predict(conv_mod,
          newdata = data.frame(CT=Res$CT_NS,
                               batch=as.factor(Res$Plate)),
          allow.new.levels = F)
preds_cens =
  predict(conv_mod,
          newdata = data.frame(CT=rep(40,nrow(Res)),
                               batch=as.factor(Res$Plate)),
          allow.new.levels = F)

Res$log10_viral_load = preds_all
Res$log10_cens_vl = preds_cens

Res$log10_viral_load[Res$CT_NS==40]=
  Res$log10_cens_vl[Res$CT_NS==40]
writeLines(sprintf('out of a total of %s samples, %s are below LLOQ (%s%%)',
                   nrow(Res),
                   sum(Res$CT_NS==40),
                   round(100*mean(Res$CT_NS==40),1)))

table(Res$log10_viral_load == Res$log10_cens_vl)

writeLines('The following IDs have duplicated barcodes:')

print(unique(Res$ID[duplicated(Res$BARCODE)]))

Res = Res[!duplicated(Res$BARCODE), ]

###### Write csv files
# Overall data files
write.table(x = SC, file = paste0(prefix_analysis_data, "/Analysis_Data/interim_control_dat.csv"), row.names = F, sep=',')

Res = dplyr::arrange(Res, Rand_date, ID, Time)
write.table(x = Res, file = paste0(prefix_analysis_data, "/Analysis_Data/interim_all_analysis.csv"), row.names = F, sep=',')

write.table(x = screen_failure, file = paste0(prefix_analysis_data, "/Analysis_Data/interim_screening_dat.csv"), row.names = F, sep=',')


Res = 
  Res %>% filter(Swab_ID != 'Saliva') %>% # remove the saliva samples
  mutate(Country= case_when(Site %in% c('th001','th057','th058') ~ 'Thailand',
                            Site == 'br003' ~ 'Brazil',
                            Site == 'la008' ~ 'Laos'))

fever_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/temperature_data.csv"))
fever_data = fever_data %>% 
  mutate(ID = Label,
         ax_temperature = fut_temp)
write.table(x = fever_data[, c('ID','Time','ax_temperature','Fever_Baseline')], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/fever_interim.csv"), 
            row.names = F, sep=',', quote = F)


symptom_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/symptom_data.csv"))
symptom_data = symptom_data %>% 
  mutate(ID = Label,
         Any_symptom = sq_yn)
write.table(x = symptom_data[, c('ID','Timepoint_ID','Any_symptom','heart_rate')], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/symptoms_interim.csv"), 
            row.names = F, sep=',', quote = F)

#***********************************************************************#
#************************ Specific analysis data files *****************#

#************************* Ivermectin Analysis *************************#
#* Only Thailand
# Res_ivermectin = 
#   Res %>% filter((Trt %in% c('Ivermectin',"No study drug")) |
#                    (Trt == 'Regeneron' & Rand_date < "2021-12-17 00:00:00" & Site == 'th001'),
#                  Rand_date < '2022-04-22 00:00:00',
#                  Country == 'Thailand') %>%
#   arrange(Rand_date, ID, Time)
# pk_summaries = read.csv('../PLATCOV-Ivermectin/PK_summaries.csv')
# Res_ivermectin = merge(Res_ivermectin, pk_summaries, by = 'ID', all.x = T)
# 
# write.table(x = Res_ivermectin, file = 'Ivermectin_analysis.csv', sep = ',',row.names = F,quote = F)
# # anonymise IDs for publication
# Res_ivermectin$ID = as.numeric(as.factor(Res_ivermectin$ID))
# 
# write.table(x = Res_ivermectin, file = '../Analysis_Data/Ivermectin_analysis.csv', sep = ',',row.names = F,quote = F)


#************************* Remdesivir Analysis *************************#
#* Thailand and Brazil
Res_Remdesivir =
  Res %>% filter(Trt %in% c('Remdesivir',"No study drug"),
                 Rand_date < '2022-06-11',
                 Country %in% c('Thailand','Brazil')) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_Remdesivir, file = paste0(prefix_analysis_data, "/Analysis_Data/Remdesivir_analysis.csv"), row.names = F, sep=',',quote = F)


fever_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/temperature_data.csv"))
fever_data = fever_data %>% 
  filter(Label %in% unique(Res_Remdesivir$ID)) %>%
  mutate(ID = Label,
         ax_temperature = fut_temp)
write.table(x = fever_data[, c('ID','Time','ax_temperature','Fever_Baseline')], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/Remdesivir_fever.csv"), 
            row.names = F, sep=',', quote = F)


symptom_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/symptom_data.csv"))
symptom_data = symptom_data %>% 
  filter(Label %in% unique(Res_Remdesivir$ID)) %>%
  mutate(ID = Label,
         Any_symptom = sq_yn)
write.table(x = symptom_data[, c('ID','Timepoint_ID','Any_symptom','heart_rate')], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/Remdesivir_symptoms.csv"), 
            row.names = F, sep=',', quote = F)


#************************* Favipiravir Analysis *************************#
#* Thailand and Brazil
# Res_Favipiravir = 
#   Res %>% filter(Trt %in% c('Favipiravir',"No study drug"),
#                  Rand_date < '2022-10-31',
#                  Country %in% c('Thailand','Brazil')) %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_Favipiravir, file = '../Analysis_Data/Favipiravir_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Regeneron Analysis *************************#
#* Thailand only
# Res_REGN = 
#   Res %>% filter(Trt %in% c('Regeneron',"No study drug"),
#                  Country == 'Thailand',
#                  Rand_date < '2022-10-20') %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_REGN, file = '../Analysis_Data/REGN_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Fluoxetine Analysis *************************#
#* Thailand added 1st April 2022; Brazil added 21st June 2022
#*  Stopped in all sites on 8th May 2023
Res_Fluoxetine = 
  Res %>% filter(Trt %in% c('Fluoxetine',"No study drug"),
                 (Country=='Thailand' & Rand_date > "2022-04-01 00:00:00") |
                   (Country=='Brazil' & Rand_date > "2022-06-21 00:00:00") |
                   (Country=='Laos' & Rand_date > "2022-06-21 00:00:00") |
                   (Country=='Pakistan' & Rand_date > "2022-06-21 00:00:00"),
                 Rand_date < "2023-05-09") %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_Fluoxetine, file = paste0(prefix_analysis_data, "/Analysis_Data/Fluoxetine_analysis.csv"), row.names = F, sep=',', quote = F)



Res_Fluoxetine_meta = 
  Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir',
                            'Molnupiravir',
                            "No study drug",
                            'Ivermectin',
                            'Remdesivir',
                            'Favipiravir',
                            'Fluoxetine'),
                 Country %in% c('Thailand','Brazil','Laos','Pakistan'),
                 Rand_date <= "2023-05-09 00:00:00") %>%
  arrange(Rand_date, ID, Time) 

write.table(x = Res_Fluoxetine_meta, 
            file = paste0(prefix_analysis_data, "/Analysis_Data/Fluoxetine_meta_analysis.csv"), 
            row.names = F, sep=',', quote = F)



#************************* Paxlovid v Molnupiravir Analysis *************************#
#* Thailand only
Res_Paxlovid_Molnupiravir = 
  Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir','Molnupiravir',"No study drug"),
                 Rand_date > "2022-06-03 00:00:00",
                 Rand_date < "2023-02-24 00:00:00",
                 Country %in% c('Thailand')) %>%
  arrange(Rand_date, ID, Time) %>% ungroup() 

write.table(x = Res_Paxlovid_Molnupiravir, file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_analysis.csv"), row.names = F, sep=',', quote = F)

## Meta-analysis which we report in the publication using unblinded arms from the same site
Res_Paxlovid_Molnupiravir_meta = 
  Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir',
                            'Molnupiravir',
                            "No study drug",
                            'Ivermectin',
                            'Remdesivir',
                            'Favipiravir'),
                 Rand_date < "2023-02-24 00:00:00",
                 Site =='th001') %>%
  arrange(Rand_date, ID, Time) 

write.table(x = Res_Paxlovid_Molnupiravir_meta, 
            file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_meta_analysis.csv"), 
            row.names = F, sep=',', quote = F)



fever_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/temperature_data.csv"))
fever_data = fever_data %>% 
  filter(Label %in% unique(Res_Paxlovid_Molnupiravir_meta$ID)) %>%
  mutate(ID = Label,
         ax_temperature = fut_temp)
write.table(x = fever_data[, c('ID','Time','ax_temperature','Fever_Baseline')], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_meta_fever.csv"), 
            row.names = F, sep=',', quote = F)



symptom_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/symptom_data.csv"))
symptom_data = symptom_data %>% 
  filter(Label %in% unique(Res_Paxlovid_Molnupiravir_meta$ID)) %>%
  mutate(ID = Label,
         Any_symptom = sq_yn)
write.table(x = symptom_data[, c('ID','Timepoint_ID','Any_symptom','heart_rate')], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_meta_symptoms.csv"), 
            row.names = F, sep=',', quote = F)


#************************* Nitazoxanide Analysis *************************#
#* Brazil and Laos
Res_Nitazoxanide = 
  Res %>% filter(Trt %in% c('Nitazoxanide',"No study drug"),
                 Rand_date > "2022-06-03 00:00:00",
                 Country %in% c('Brazil','Laos'),
                 !ID %in% c("PLT-BR3-006",
                            "PLT-BR3-018",
                            "PLT-BR3-033",
                            "PLT-BR3-043")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_Nitazoxanide, file = paste0(prefix_analysis_data, "/Analysis_Data/Nitazoxanide_analysis.csv"), row.names = F, sep=',', quote = F)


#************************* Evusheld Analysis *************************#
#* Thailand added 2022-09-01; Brazil added 2022-10-31
Res_Evusheld = 
  Res %>% filter(Trt %in% c('Evusheld',"No study drug"),
                 (Country=='Thailand' & Rand_date > "2022-09-01 00:00:00") |
                   (Country=='Brazil' & Rand_date > "2022-10-31 00:00:00")) %>%
  arrange(Rand_date, ID, Time)

evusheld_mutations = read.csv(paste0(prefix_analysis_data, "/Analysis_Data/evusheld_test.csv"))
table(evusheld_mutations$evusheld_resistant, useNA = 'ifany')

Res_Evusheld = merge(Res_Evusheld, evusheld_mutations, by = 'ID', all.x = T)
Res_Evusheld$evusheld_resistant =
  ifelse(is.na(Res_Evusheld$evusheld_resistant),T,
         Res_Evusheld$evusheld_resistant)

Res_Evusheld$evusheld_resistant[Res_Evusheld$Country=='Brazil' & Res_Evusheld$Rand_date<'2022-12-01']=F

write.table(x = Res_Evusheld, file = paste0(prefix_analysis_data, "/Analysis_Data/Evusheld_analysis.csv"), row.names = F, sep=',', quote = F)


#************************* Ensitrelvir Analysis *************************#
#* Thailand and Laos added 2023-03-17

Res_Ensitrelvir = 
  Res %>% filter(Trt %in% c('Ensitrelvir',"No study drug",'Nirmatrelvir + Ritonavir'),
                 (Country=='Thailand' & Rand_date >= "2023-03-17 00:00:00") |
                   (Country=='Laos' & Rand_date >= "2023-03-17 00:00:00")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_Ensitrelvir, file = paste0(prefix_analysis_data, "/Analysis_Data/Ensitrelvir_analysis.csv"), row.names = F, sep=',', quote = F)


#************************* Ineffective Interventions *************************#

Res_ineffective = 
  Res %>% filter(Trt %in% c('Ivermectin',
                            "Favipiravir",
                            "Fluoxetine",
                            "No study drug")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_ineffective, file = paste0(prefix_analysis_data, "/Analysis_Data/Ineffective_analysis.csv"), row.names = F, sep=',', quote = F)

#************************* No Study Drugs *************************#
Res_noStudyDrugs = 
  Res %>% filter(Trt %in% c("No study drug")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_noStudyDrugs, file = paste0(prefix_analysis_data, "/Analysis_Data/Res_noStudyDrugs.csv"), row.names = F, sep=',', quote = F)

#************************* Site TH01 only *************************#
Res_TH1 <-  Res %>% filter(Site %in% c("th001")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_noStudyDrugs, file = paste0(prefix_analysis_data, "/Analysis_Data/Res_TH1.csv"), row.names = F, sep=',', quote = F)

