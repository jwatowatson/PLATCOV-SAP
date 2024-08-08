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
source('user_settings.R')
source('functions.R')
#############################################################################################
##******************** Randomisation data *******************
rand_app_data = rbind(read.csv(paste0(prefix_drop_rand, "/data-TH1.csv")),
                      read.csv(paste0(prefix_drop_rand, "/data-BR3.csv")),
                      read.csv(paste0(prefix_drop_rand, "/data-LA08.csv")),
                      read.csv(paste0(prefix_drop_rand, "/data-PK01.csv"))) %>%
  filter(!is.na(Treatment))
# Defining patient ID
rand_app_data$ID = paste0('PLT-', gsub(x = rand_app_data$site,pattern = '0',replacement = ''),
                          '-',
                          stringr::str_pad(rand_app_data$randomizationID, 3, pad = "0"))
# Converting time zones
rand_app_data$Site = plyr::mapvalues(x = rand_app_data$site, from=c('TH1','BR3','LA08','PK01'),to=c('th001','br003',"la008","pk001"))
rand_app_data$Rand_Time = as.POSIXct(rand_app_data$Date, format = '%a %b %d %H:%M:%S %Y',tz = 'GMT')
rand_app_data$tzone = plyr::mapvalues(x = rand_app_data$site,
                                      from = c('TH1','LA08','BR3','PK01'),
                                      to = c('Asia/Bangkok','Asia/Bangkok','America/Sao_Paulo','Asia/Karachi'))
rand_app_data$Rand_Time_TZ=NA
for(i in 1:nrow(rand_app_data)){
  rand_app_data$Rand_Time_TZ[i] = as.character(with_tz(rand_app_data$Rand_Time[i], tzone = rand_app_data$tzone[i]))
}

# Creating a checklist dataset
# NOTE: Sites TH57 and TH58 did not use app so cannot cross check
check_data <- data.frame(rand_app_data[,c("ID", "sex", "age", "Rand_Time_TZ", "Treatment")])

##*********************************************************
##******************* Clinical database *******************
##*********************************************************
##########  --- Clinical data --- ########## 
clin_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimEnrolment.dta"))
clin_data$scrpassed[clin_data$Label=='PLT-TH1-557']=1

##########  Extract screening failure data
table(clin_data$scrpassed, useNA = 'ifany')
ind = !is.na(clin_data$scrpassed) & clin_data$scrpassed==0
screen_failure =
  clin_data[ind,
            c("Trial","Site","scrid","scrdat",
              "scrpassed","reason_failure","scrnote")]

screen_failure$reason_failure = sjlabelled::as_character(screen_failure$reason_failure)
write.csv(x = screen_failure, file = '../Analysis_Data/screening_failure.csv')

clin_data <- clin_data[!ind,]

##########  Preparing data for checkings 
##### removing patients that not randomised
ind <- clin_data$scrpassed == 1 & is.na(clin_data$rangrp) & clin_data$Label == ""
writeLines(sprintf('This patient %s passed the screening but never randomised', 
                   clin_data$scrid[ind]))
clin_data <- clin_data[!ind,]

clin_data$Sex <- plyr::mapvalues(x = as.numeric(clin_data$sex),
                from = c(1,2),
                to = c('Male','Female'))
### Check which patients has the data entered to the clinical database
check_data$on_macro_yn <- (check_data$ID %in% clin_data$Label)
writeLines(sprintf('Patient %s has no data on MACRO', 
                   check_data$ID[!check_data$on_macro_yn]))
### Check if sex information matched between randomisation database and clinical database
check_data <- merge(check_data, clin_data[,c("Label", "Sex")], by.x = "ID",  by.y = "Label", all.x = T)
check_data$sex_agree_yn <- check_data$sex == check_data$Sex
writeLines(sprintf('Patient %s has mismatched sex data: MACRO = %s and Randomisation = %s', 
                   check_data$ID[!check_data$sex_agree_yn & !is.na(check_data$sex_agree_yn)],
                   check_data$Sex[!check_data$sex_agree_yn & !is.na(check_data$sex_agree_yn)],
                   check_data$sex[!check_data$sex_agree_yn & !is.na(check_data$sex_agree_yn)]))
writeLines(sprintf('Patient %s has missing sex data on MACRO', 
                   check_data$ID[is.na(check_data$sex_agree_yn)]))
sex_problem_ID <- check_data$ID[check_data$on_macro_yn & (!check_data$sex_agree_yn | is.na(check_data$sex_agree_yn))]

# Using sex information from randomisation database in further analyses
for(i in 1:length(sex_problem_ID)){
  clin_data$Sex[which(clin_data$Label == sex_problem_ID[i])] <- rand_app_data$sex[rand_app_data$ID == sex_problem_ID[i]]
}
writeLines("Frequency of sex in all clinical data after corrections:")
print(table(clin_data$Sex, useNA = "always"))
check_data <- check_data[,-which(names(check_data) %in% c("sex", "Sex"))]
### Check if age information matched between randomisation database and clinical database
clin_data$age_yr[clin_data$Label=='PLT-TH57-002'] = 21 # A special case

ind = is.na(clin_data$age_yr) & !is.na(clin_data$dob_my)
# calculate age at randomisation 
for(i in which(ind)){
  clin_data$age_yr[i] = trunc((parse_date_time(clin_data$dob_my[i], c('%d/%m/%Y', '%m/%Y')) %--% parse_date_time(clin_data$randat[i],  c('%Y-%m-%d')))/years(1))
}

check_data <- merge(check_data, clin_data[,c("Label", "dob_my", "age_yr")], by.x = "ID",  by.y = "Label", all.x = T)
check_data$age <- floor(as.numeric(check_data$age))
check_data$age_dob_missing <- is.na(check_data$age_yr)
check_data$age_agree_yn <- check_data$age == check_data$age_yr

writeLines(sprintf('Patient %s has no age or birth date information on MACRO', 
                   check_data$ID[check_data$age_dob_missing]))
writeLines(sprintf('Patient %s has mismatched age data: MACRO = %s and Randomisation = %s', 
                   check_data$ID[!check_data$age_agree_yn & !is.na(check_data$age_agree_yn)],
                   check_data$age_yr[!check_data$age_agree_yn & !is.na(check_data$age_agree_yn)],
                   check_data$age[!check_data$age_agree_yn & !is.na(check_data$age_agree_yn)]
                   ))

check_data <- check_data[,-which(names(check_data) %in% c("dob_my", "age_yr", "age"))]

### Check if symptom onset data is missing
ind = !is.na(clin_data$cov_symphr) & is.na(clin_data$cov_sympday)
if(sum(ind)>0) clin_data$cov_sympday[ind] = clin_data$cov_symphr[ind]/24 #If day is missing >>> divide hours by 24

ind = is.na(clin_data$cov_sympday)

check_data <- merge(check_data, clin_data[,c("Label", "cov_sympday")], by.x = "ID",  by.y = "Label", all.x = T)
check_data$sympday_missing <- is.na(check_data$cov_sympday)
writeLines(sprintf('Patient %s has no information on symptom onset (in hours or days) on MACRO', 
                   check_data$ID[check_data$sympday_missing]))

# Assign Symptomatic day = 2 for missing data 
clin_data$cov_sympday[is.na(clin_data$cov_sympday)]=2

check_data$sympday_exceed_yn <- check_data$cov_sympday > 4
writeLines(sprintf('Patient %s has symptom onset greater than 4 days, which is %s days', 
                   check_data$ID[check_data$sympday_exceed_yn & !is.na(check_data$sympday_exceed_yn)],
                   check_data$cov_sympday[check_data$sympday_exceed_yn & !is.na(check_data$sympday_exceed_yn)]
                   ))

check_data <- check_data[,-which(names(check_data) %in% c("cov_sympday"))]

### Check if weight and height data are missing
clin_data$BMI = clin_data$weight/(clin_data$height/100)^2
clin_data$Weight = clin_data$weight

check_data <- merge(check_data, clin_data[,c("Label", "BMI", "weight", "height")], by.x = "ID",  by.y = "Label", all.x = T)
check_data$weight_height_missing <- is.na(check_data$BMI)
writeLines(sprintf('Patient %s has no information on weights and heights on MACRO', 
                   check_data$ID[check_data$weight_height_missing]))

check_data$weight_height_outlier <- abs(check_data$BMI - mean(check_data$BMI, na.rm = T)) > 3*sd(check_data$BMI, na.rm = T)
writeLines(sprintf('Weights/heights of patient %s is outlier: Weight = %s kg; Height = %s cm', 
                   check_data$ID[check_data$weight_height_outlier & !(check_data$weight_height_missing)],
                   check_data$weight[check_data$weight_height_outlier & !(check_data$weight_height_missing)],
                   check_data$height[check_data$weight_height_outlier & !(check_data$weight_height_missing)]
                   ))

check_data <- check_data[,-which(names(check_data) %in% c("BMI", "weight", "height"))]

### Check if randomisation date and time is correct
clin_data$Rand_date_time = NA
for(i in 1:nrow(clin_data)){
  if(!is.na(clin_data$randat[i]) & !is.na(clin_data$rantim[i])){
    clin_data$Rand_date_time[i] = as.character(as.POSIXct(paste(clin_data$randat[i],
                                                                clin_data$rantim[i], sep=' ')))
  }
}

check_data <- merge(check_data, clin_data[,c("Label", "Rand_date_time")], by.x = "ID",  by.y = "Label", all.x = T)
check_data$Rand_diffs = apply(check_data[,c('Rand_date_time','Rand_Time_TZ')],1, function(x) difftime(x[1], x[2], units='mins'))

# Check missing data
check_data$rand_date_missing <- is.na(check_data$Rand_date_time)
writeLines(sprintf('Patient %s has no information on randomisation date and time on MACRO', 
                   check_data$ID[check_data$rand_date_missing]))

# Check if the differences is more than 5 mins between both databases
check_data$randtime_diff_exceed <- check_data$Rand_diffs > 5
writeLines(sprintf('More than 5 min difference in rand time for %s, MACRO = %s and Randomisation = %s', 
                   check_data$ID[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)],
                   check_data$Rand_date_time[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)],
                   check_data$Rand_Time_TZ[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)]
))
ID_exceed_5mins <- check_data$ID[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)]

# Using randomisation date and time from randomisation database in further analyses
for(i in 1:length(ID_exceed_5mins)){
  clin_data$Rand_date_time[which(clin_data$Label == ID_exceed_5mins[i])] <- rand_app_data$Rand_Time_TZ[rand_app_data$ID == ID_exceed_5mins[i]]
}

check_data <- check_data[,-which(names(check_data) %in% c("Rand_date_time", "Rand_Time_TZ", "Rand_diffs"))]

### Check if randomisation arms matched 
clin_data$rangrp = sjlabelled::as_character(clin_data$rangrp)
clin_data$rangrp[clin_data$rangrp=='Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir'
clin_data$rangrp[clin_data$rangrp=='Molnupiravir and Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir + Molnupiravir'

check_data <- merge(check_data, clin_data[,c("Label", "rangrp")], by.x = "ID",  by.y = "Label", all.x = T)
check_data$rangrp_missing <- is.na(check_data$rangrp)
writeLines(sprintf('Patient %s has no information on treatment arms on MACRO', 
                   check_data$ID[check_data$rangrp_missing]))

check_data$rangrp_agree <- check_data$Treatment == check_data$rangrp
writeLines(sprintf('Randomisation data mismatched for %s, MACRO = %s and Randomisation = %s', 
                   check_data$ID[!check_data$rangrp_agree & !is.na(check_data$rangrp_missing)],
                   check_data$Rand_date_time[!check_data$rangrp_agree & !is.na(check_data$rangrp_missing)],
                   check_data$Rand_Time_TZ[!check_data$rangrp_agree & !is.na(check_data$rangrp_missing)]
))

colnames(clin_data)
check_data <- check_data[,-which(names(check_data) %in% c("rangrp"))]

##########  --- Temperature data --- ########## 
temp_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimFUTemp.dta"))
# Check if there is any patients without follow-up temperature data
check_data$fu_temp_yn <- check_data$ID %in% temp_data$Label
writeLines(sprintf('Patient %s has no information on follow-up tempearature dataset', 
                   check_data$ID[!check_data$fu_temp_yn]))

x1 = temp_data[, c('Site', 'Label', 'visit', 'fut_amdat', "fut_amtim", "fut_amtemp")] #Morning
x2 = temp_data[, c('Site', 'Label', 'visit', 'fut_pmdat', "fut_pmtim", "fut_pmtemp")] #Evening
colnames(x1) = colnames(x2) = 
  c('Site', 'Label', 'visit', 'fut_dat', "fut_tim", "fut_temp")
# Check which patients has missing data in the morning and evening
fu_temp_am_missing <- check_temp_missing(x1, "am")
fu_temp_pm_missing <- check_temp_missing(x2, "pm")

check_data <-  merge(check_data, fu_temp_am_missing, by.x = "ID",  by.y = "Label", all.x = T)
check_data <-  merge(check_data, fu_temp_pm_missing, by.x = "ID",  by.y = "Label", all.x = T)

# Define fever at baseline
fever_data <- prep_tempdata(temp_data, clin_data)
fever_data
# Check which patients has mismatch temperature time and timepoint ID
FUtemp_checktime <- fever_data[abs(fever_data$Timepoint_ID - fever_data$Time) > 1, c("Label", "visit", "Rand_date_time", "temp_time", "visit",
                                                                 "Timepoint_ID", "Time")]
writeLines(sprintf('%s of patient %s: Rand_date = %s; Temp_date = %s; which is %s days post-randomisation', 
                   FUtemp_checktime$visit,
                   FUtemp_checktime$Label,
                   FUtemp_checktime$Rand_date_time,
                   FUtemp_checktime$temp_time,
                   round(FUtemp_checktime$Time,2)
                   ))
FUtemp_checktime_wide <- FUtemp_checktime[,c("Label", "visit", "Time")] %>%
  distinct(Label, visit, .keep_all = T) %>%
  pivot_wider(names_from = visit, values_from = Time, values_fill = NA)
colnames(FUtemp_checktime_wide)[-1] <- paste0("futemp_time_mismatch_", colnames(FUtemp_checktime_wide)[-1])
FUtemp_checktime_wide[,-1] <- !is.na(FUtemp_checktime_wide[,-1])

check_data <-  merge(check_data, FUtemp_checktime_wide, by.x = "ID",  by.y = "Label", all.x = T)
check_data[,colnames(FUtemp_checktime_wide)[-1] ][is.na(check_data[,colnames(FUtemp_checktime_wide)[-1] ])] <- F

# Extract fever at baseline
fever_data = fever_data %>% distinct(Label, .keep_all = T)
clin_data <- merge(clin_data, fever_data[, c('Label','Fever_Baseline')], by='Label', all = T)

# Check if baseline temperature data is available?
check_data <-  merge(check_data, clin_data[, c('Label','Fever_Baseline')], by.x = "ID",  by.y = "Label", all.x = T)
check_data$fever_baseline_missing <- is.na(check_data$Fever_Baseline)
writeLines(sprintf('Patient %s has no information on baseline fever', 
                   check_data$ID[check_data$fever_baseline_missing]))

##########  --- Symptom data --- ########## 
vita_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimVitalSigns.dta"))
check_data$vita_missing <- !check_data$ID %in% vita_data$Label

symp=haven::read_dta(paste0(prefix_dropbox, "/Data/InterimSymptoms.dta"))
check_data$symp_missing <- !check_data$ID %in% symp$Label

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

symp_data = merge(symp, HR_data[, c('Label', 'Timepoint_ID','heart_rate')],
                  all=T, by = c('Label', "Timepoint_ID"))
symp_data$Timepoint_ID = as.numeric(symp_data$Timepoint_ID)
symp_data = symp_data %>% arrange(Label, Timepoint_ID)

write.csv(x = symp_data, file = '../Analysis_Data/symptom_data.csv', row.names = F, quote = F)
#----------------------------------------------------------------------------------------
symptom_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/symptom_data.csv"))
# cleaning
symptom_data = symptom_data %>% 
  mutate(ID = Label,
         Any_symptom = sq_yn) %>%
  rename("sq_soreyn" = "sq_sore")
#----------------------------------------------------------------------------------------
# Listing 'other symptoms'
other_symptoms <- symptom_data %>%
  select(matches("^sq.*des$")) %>%
  unlist() %>%
  unname() %>%
  table()

write.csv(other_symptoms, "../Analysis_Data/other_symptoms.csv", row.names = F)
#----------------------------------------------------------------------------------------
# Check if 'other symptoms' have been filled
col_others <- colnames(symptom_data)[grepl('^sq.*des$', colnames(symptom_data))]

symptom_data <- symptom_data %>%
  mutate(sq_otheryn = NA) %>%
  mutate(across(all_of(col_others), ~ifelse(. == "", NA, .)))

for(i in 1:nrow(symptom_data)){
  if(!all(is.na(symptom_data[i,col_others]))){symptom_data$sq_otheryn[i] <- 1} else {symptom_data$sq_otheryn[i] <- 0}
}

# Check if 1 in sq_yn always provided details
columns_to_check <- colnames(symptom_data)[grepl('^sq.*yn$', colnames(symptom_data))]
columns_to_check <- columns_to_check[-1]
check_sqyn <- symptom_data %>% 
  filter(sq_yn == 1) %>%
  mutate(row_sum = (rowSums(select(., all_of(columns_to_check)), na.rm = TRUE))) %>%
  filter(row_sum == 0)

write.csv(check_sqyn, "symptoms_no_details.csv", row.names = F) # These entries report the presence of symptoms without providing more details
#----------------------------------------------------------------------------------------
col_symptom_gr <- colnames(symptom_data)[grepl('^sq.*gr$', colnames(symptom_data))]
col_symptom_gr <- col_symptom_gr[!grepl('oth', col_symptom_gr)]

summary(symptom_data[,col_symptom_gr]) #Check if there is any dodgy grade data

symptom_data <- symptom_data %>%
  mutate(across(.cols = all_of(col_symptom_gr), 
                .fns = ~ifelse(!is.na(sq_yn) & is.na(.), 0, .)))

write.table(x = symptom_data[, c('ID','Timepoint_ID','Any_symptom','heart_rate', col_symptom_gr)], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/symptoms_interim.csv"), 
            row.names = F, sep=',', quote = F)

#----------------------------------------------------------------------------------------
##########  --- Final status data --- ########## 
final_status = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimFinalStatus.dta"))
final_status = final_status[!is.na(final_status$fs_compyn), ]

# Check if there is any patients without final status data
check_data$fs_missing <- !check_data$ID %in% final_status$Label
writeLines(sprintf('Patient %s has no information on final status dataset', 
                   check_data$ID[check_data$fs_missing]))

check_data <-  merge(check_data, final_status[,c('Label', 'fs_compyn', 'fs_rsn', 'fs_rsnothsp', 'fs_ae',
                                                 'fs_sae', 'fs_diecov')], 
                     by.x = "ID",  by.y = "Label", all.x = T)
# Check if all patients that did not complete the treatment has reasons indicated
check_data$fail_reason_yn <- T
ID_no_fail_rsn <- check_data$ID[!is.na(check_data$fs_compyn) & check_data$fs_compyn == 0 & is.na(check_data$fs_rsn)]
check_data$fail_reason_yn[check_data$ID %in% ID_no_fail_rsn] <- F
writeLines(sprintf('Patient %s does not have a reason indicated for non-completion', 
                   ID_no_fail_rsn))

check_data$fail_reason_other <- T
check_data$fs_rsnothsp[check_data$fs_rsnothsp == ""] <- NA
ID_no_fail_rsn_other <- check_data$ID[!is.na(check_data$fs_rsn) & check_data$fs_rsn == 5 & is.na(check_data$fs_rsnothsp)]
check_data$fail_reason_yn[check_data$ID %in% ID_no_fail_rsn_other] <- F
writeLines(sprintf('Patient %s does not have a reason indicated for non-completion for "other reasons"', 
                   ID_no_fail_rsn_other))

# Check if any patient died
writeLines(sprintf('Patient %s was flagged as dead!!!! Please check!!!', 
                   check_data$ID[check_data$fs_diecov == 1 & !is.na(check_data$fs_diecov)]))
check_data$dead_flag <- F
check_data$dead_flag[check_data$fs_diecov == 1 & !is.na(check_data$fs_diecov)] <- T

# Check if SAEs are flagged as AEs
writeLines(sprintf('Patient %s has SAEs but not flagged as AEs',
                   check_data$ID[check_data$fs_sae == 1 &  check_data$fs_ae == 0 & !(is.na(check_data$fs_sae) & is.na(check_data$fs_ae))]))
check_data$sae_ae_match <- T
check_data$sae_ae_match[check_data$fs_sae == 1 &  check_data$fs_ae == 0 & !(is.na(check_data$fs_sae) & is.na(check_data$fs_ae))] <- F

check_data <- check_data[,-which(names(check_data) %in% c('fs_compyn', 'fs_rsn', 'fs_rsnothsp', 'fs_ae', 'fs_sae', 'fs_diecov'))]
#######################################################################################################################
clin_data$Symptom_onset = clin_data$cov_sympday
clin_data$Trt = clin_data$rangrp
clin_data$Sex = as.numeric(clin_data$sex)
clin_data$Sex[clin_data$Sex == 2] <- 0
clin_data$Age <- clin_data$age_yr

# Select variables of interest
clin_data = clin_data[, c('Label','Trt','Sex','Age','randat',
                          "Rand_date_time",'BMI','Weight',
                          'Symptom_onset','Site','Fever_Baseline')]
#######################################################################################################################
##########  --- Variant data --- ########## 
Sample_ID_map <- extract_FASTA() # This function compiled all FASTA files and saved it.
# Run Nextclade for extracting mutation
# Can't run on MORU internet
# Input: Combined Fasta file "all_fasta.fasta"; SARS-CoV-2 data (downloaded from Nextclade: https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/usage.html)
# Output: Interested only .tsv file
# Need Nextclade installed
##------------------------------------------------------------------------------------
re_download = F
system_used = "mac"

if(re_download){
  arg_download <- "nextclade dataset get --name nextstrain/sars-cov-2/wuhan-hu-1/orfs --output-dir ../Analysis_Data/Nextclade/sars-cov-2"
  arg_download
  if(system_used == "windows"){
    shell(arg_download)
  } else {
    system(arg_download) 
  }
}
##------------------------------------------------------------------------------------
sars_cov_2_data <- "../Analysis_Data/Nextclade/sars-cov-2" # path to the data downloaded from Nextclade
output_folder <- "../Analysis_Data/Nextclade/output" # folder name for outputs
input_fasta <- "../Analysis_Data/all_fasta.fasta" # path to input fasta files

arguments <- paste0("nextclade run --input-dataset ", sars_cov_2_data, " --output-all ", output_folder, " ", input_fasta)
arguments

if(system_used == "windows"){
  shell(arguments)
} else {
  system(arguments) 
}
##------------------------------------------------------------------------------------
# Run python to classify varaints
variant_data = get_nanopore_data(prefix_analysis_data = prefix_analysis_data, run_python = F, system_used = "mac")
variant_data = merge(variant_data, clin_data, by.x='ID', by.y = 'Label', all.y = T)
#To reduce the number of variant groups, as suggested by Liz.
variant_data$Variant2 <- as.character(variant_data$Variant)
variant_data$Variant2[variant_data$Variant2 %in% c("BA.5.2", "BA.5.5", "BQ.1")] <- "BA.5"
variant_data$Variant2[variant_data$Variant2 %in% c("BN.1.2", "BN.1.3", "CH.1.1")] <- "BA.2.75"
variant_data$Variant2[variant_data$Variant2 %in% c("XBB1.5-like with F456L")] <- "XBB.1.5-like"
variant_data$Variant2 <- as.factor(variant_data$Variant2)

check_data <- merge(check_data, variant_data[,c("ID", "Variant")], by = "ID", all.x = T)
check_data$Variant_missing <- is.na(check_data$Variant)
check_data <- check_data[,-which(names(check_data) %in% c('Variant'))]

table(check_data$Variant_missing)
#######################################################################################################################
##########  --- Vaccine data --- ########## 
vacc_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimVaccine.dta"))

writeLines(sprintf('No vaccine data for %s',
                   clin_data$Label[!clin_data$Label %in% vacc_data$Label]))
check_data$vacc_data_missing <- !check_data$ID %in% vacc_data$Label

## Some cleaning
# vacc_data$vc_name = tolower(vacc_data$vc_name)
writeLines('\nVaccine names before cleaning:')
print(table(vacc_data$vc_name))

clin_data$Any_dose = NA
clin_data$N_dose = NA
clin_data$Time_since_last_dose = NA
clin_data$N_dose_mRNA = NA
clin_data$Any_dose_mRNA = NA
clin_data$Any_dose_mRNA = NA

vacc_date_cols = grep('dos', colnames(vacc_data))
## find any bad dates
all_dates = unlist(unique(vacc_data[, vacc_date_cols]))
all_dates = all_dates[all_dates!='']
bad_dates = all_dates[which(is.na(parse_date_time(all_dates, orders = 'dmy')))]
vacc_data$Label[which(apply(vacc_data[, vacc_date_cols], 1, function(x) length(intersect(x ,bad_dates))>0))]

vac_error_ID <- NULL
time_vac_negative_ID <- NULL
clin_data$Label[!clin_data$Label %in% vacc_data$Label]
for(i in 1:nrow(clin_data)){
  id = clin_data$Label[i]
  ind = which(vacc_data$Label==id)
  vacc_data[vacc_data$Label == id, ]
  
  ind_mRNA = which(vacc_data$Label==id & vacc_data$vc_name %in% c('Moderna','Pfizer','Chula-Cov19'))
  
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
    if(clin_data$Time_since_last_dose[i] < 0 & !is.na(clin_data$Time_since_last_dose[i])){
      time_vac_negative_ID <- c(time_vac_negative_ID, id)}
  }
  if(any(vacc_data$vc_statyn[vacc_data$Label == id] == 1) & (length(vac_dates) == 0)){vac_error_ID <- c(vac_error_ID, id)}
}
# No vaccine date details
writeLines(sprintf('Patient %s flagged to have vaccinated but no details',
                   vac_error_ID))
check_data$vacc_detail_missing <- check_data$ID %in% vac_error_ID
# Negative time since last vaccination
writeLines(sprintf('Patient %s has a negative time since last vaccination',
                   time_vac_negative_ID))
check_data$vacc_time_negative <- check_data$ID %in% time_vac_negative_ID
#######################################################################################################################
##########  --- Sample log data --- ########## 
log_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimSampleLog.dta"))
log_data$sl_barc[log_data$sl_barc=='']=NA
log_data = log_data[!is.na(log_data$sl_barc), ]
log_data = log_data[!is.na(log_data$sl_sampdat), ]
log_data = log_data[!is.na(log_data$sl_samptim), ]
log_data$sl_barc = tolower(log_data$sl_barc) # avoids case dependency

## Per protocol for treatment data
trt_distcont_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimDrugRescue.dta")) %>% filter(!is.na(dardat))

##########  --- Viral density data --- ########## 
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
  if(length(grep(pattern = 'Laos', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Laos',x,sep='_'))
  }
  if(length(grep(pattern = 'Pakistan', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Pakistan',x,sep='_'))
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

## 
Res$`SUBJECT ID` = gsub(pattern = '_', replacement = '-', x = Res$`SUBJECT ID`,fixed = T)
Res$`SUBJECT ID` = gsub(pattern = 'PK01', replacement = 'PK1', x = Res$`SUBJECT ID`,fixed = T)

# make the plate/lab variables
Res$Lab = NA
Res$Lab[grep(pattern = 'Thailand',x = Res$`Lot no.`, ignore.case = T)]='Thailand'
Res$Lab[grep(pattern = 'Brazil',x = Res$`Lot no.`, ignore.case = T)]='Brazil'
Res$Lab[grep(pattern = 'Lao',x = Res$`Lot no.`, ignore.case = T)]='Laos'
Res$Lab[grep(pattern = 'Paki',x = Res$`Lot no.`, ignore.case = T)]='Pakistan'

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

for(i in 1:length(unique(Res$Plate))){
  tmp <- Res[Res$Plate == unique(Res$Plate)[i],]
  tmp <- tmp[!duplicated(tmp),]
  if(nrow(tmp) > 96){
    writeLines(sprintf('Sample lot %s has more than 96 samples (%s samples) on a single plate!!!',
                       tmp$`Lot no.`[1],
                       nrow(tmp)),
               )
}
}

### --- Extract standard curve data by plate --- ###
ind = grep('std', Res$`Sample ID`)
SC = Res[ind, c('Sample ID','N/S Gene','Target conc. c/mL','Plate','Lab','Lot no.')]

SC = SC[!duplicated(SC),]

SC$`N/S Gene`[SC$`N/S Gene`=='Undetermined'] = 40
SC$CT_NS = as.numeric(SC$`N/S Gene`)
SC = SC[!is.na(SC$CT_NS), ]
SC$log10_true_density = log10(as.numeric(SC$`Target conc. c/mL`))
SC$ID = apply(SC[, c('Sample ID','Plate')], 1, function(x) paste(x[1],as.numeric(x[2]), sep='_'))

# Select columns
cols = c('ID','Plate','CT_NS','log10_true_density','Lab','Lot no.')
SC = SC[,cols]

## Extract sample data - this was for Liz Batty
Res = Res[!is.na(Res$`SUBJECT ID`), ]
Res$ID_sample = apply(Res, 1, function(x) paste(x[c("SUBJECT ID","Location","TIME-POINT")], collapse = '_'))
#write_csv(x = Res[, c('ID_sample',"SUBJECT ID","BARCODE","Location","TIME-POINT","Time Collected","COLLECTION DATE")],file = '~/Downloads/Liz.csv')

Res = Res[Res$Location != 'Saliva', ]
Res <- Res[rowSums(is.na(Res)) != ncol(Res),]

writeLines('Clinical data from the following patients not in database:\n')
print(unique(Res$`SUBJECT ID`[!Res$`SUBJECT ID` %in%  clin_data$Label]))
########## NEED TO CHECK WHICH PATIENTS HAVE PCR DATA BUT NOT  CLINICAL DATA ###########
# "PLT-LA8-011"
# "PLT-BR3-101"
Res = Res[Res$`SUBJECT ID` %in% clin_data$Label, ]

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
########### NEED TO FINALIZE WHICH LOCATION BELONGS TO WHICH SUB STUDY GROUP #############
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
#Manual correction
Res$Timepoint_ID[Res$ID == "PLT-LA8-015"] <- "D0"

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
na_time_since_rand = NULL
negative_time_since_rand = NULL
timepoint_id_not_matched = NULL
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
        sampling_time_conflicts=rbind(sampling_time_conflicts,
                                      data.frame("ID" = Res$`SUBJECT ID`[i], "BARCODE" = barcode_i,
                                                 "Timepoint_ID" = Res$Timepoint_ID[i],
                                                 "Sample_time_log" =  sample_time_log,
                                                 "Sample_time" = sample_time))
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
      
      na_time_since_rand <- rbind(na_time_since_rand,
                                  data.frame("ID" = Res$`SUBJECT ID`[i], "BARCODE" = barcode_i,
                                             "Timepoint_ID" = Res$Timepoint_ID[i],
                                  "Rand_date_time" =  rand_time,
                                  "Sample_time" = my_sample_time))
    }
  } else {
    if(Res$Time[i] < -.1 &
       !(Res$Site[i] %in% c('th057', 'th058'))){
      writeLines(sprintf('Negative sample time for patient %s at timepoint %s: %s days, BARCODE is %s',
                         id,Res$Timepoint_ID[i],
                         round(Res$Time[i],1),
                         barcode_i))
      
      negative_time_since_rand <- rbind(negative_time_since_rand,
                                  data.frame("ID" = Res$`SUBJECT ID`[i], "BARCODE" = barcode_i,
                                             "Timepoint_ID" = Res$Timepoint_ID[i],
                                             "Time" = Res$Time[i],
                                  "Rand_date_time" =  rand_time,
                                  "Sample_time" = my_sample_time
                                  ))
    }
  }
  
  if(!is.na(Res$Time[i]) & abs(Res$Time[i] - Res$Timepoint_ID[i])>1){
    timepoint_id_not_matched <- rbind(timepoint_id_not_matched,
                                      data.frame("ID" = Res$`SUBJECT ID`[i], "BARCODE" = barcode_i,
                                                 "Timepoint_ID" = Res$Timepoint_ID[i],
                                                 "Time" = Res$Time[i],
                                                  "Rand_date_time" =  rand_time,
                                                 "Sample_time" = my_sample_time
                                                 ))
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


sampling_time_conflicts #between sample log dataset and PRC dataset > 2 hours
na_sample_times #Missing data on sampling time
na_time_since_rand #Missing data on time since randomisation
negative_time_since_rand #Negative time since randomisation
timepoint_id_not_matched


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

Res = merge(Res, variant_data[,c('ID','Variant', 'Variant2')], by = 'ID', all.x = T)

# manually correct error
Res$Swab_ID[Res$BARCODE=='20LH895'] = 'Right_tonsil_1'
Res$Swab_ID[Res$BARCODE=='20RR176'] = 'Left_tonsil_1'


Res$Rand_date = format.Date(x = Res$Rand_date_time, format='%Y-%m-%d')
##***********************************************
cols = c('ID','Time','Trt','Site','Timepoint_ID',
         'Swab_ID','Rand_date','Any_dose','N_dose','Time_since_last_dose',
         'Any_dose_mRNA','N_dose_mRNA',
         'Weight','BMI','Plate','Fever_Baseline','BARCODE',
         'Age', 'Sex', 'Symptom_onset','Variant', 'Variant2',
         'CT_NS','CT_RNaseP', 'Per_protocol_sample','Lab', 'Lot no.')

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
conv_mod = lmer(log10_true_density ~ 1 + CT*Lab + (1+CT|batch),
                data = control_dat,
                control = lmerControl(optimizer ="Nelder_Mead"))

preds = predict(conv_mod)

#Plot standard curve############################################################################
plot(control_dat$CT, jitter(control_dat$log10_true_density), xlim=c(18,40),
     col = control_dat$Lab)
for(bb in levels(control_dat$batch)){
  ind = control_dat$batch==bb
  lines(control_dat$CT[ind], preds[ind], col = as.numeric(control_dat$Lab[ind]=='Thailand')+1)
}
#Plot standard curve (ggplot) ############################################################################
control_dat2 <- control_dat[1:length(preds),]
control_dat2$preds <- preds

G1 <- ggplot(control_dat2, aes(x = CT, y = log10_true_density, col = Lab)) +
  theme_bw() +
  geom_line(aes(x = CT, y = preds, group = batch), alpha = 0.5, linewidth = 0.5) +
  scale_color_manual(values = (c("red", "black", "blue", "#CF4DCE")), name = "") +
  geom_point(width  = 0, height = 0.05, shape = 1, alpha = 1, size = 3.5) +
  xlab("CT values") +
  ylab("Viral densities (log10 genomes/mL)") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")
  ) +
  facet_wrap(Lab~., ncol = 4) +
  scale_y_continuous(breaks = seq(2,10,1)) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"))  +
  ggtitle("A) Standard curves")
G1

# png("../Plots/standard_curves.png", width = 10, height = 4, units = "in", res = 350)
# G1
# dev.off()

control_dat2$resid <- resid(conv_mod)

var <- aggregate(list("resid" = control_dat2$resid), list("Lab" = control_dat2$Lab), var)
var$Labels <- paste0("var = ", sprintf("%.3f", round(var$resid,3)))


G2 <- ggplot(control_dat2, aes(x =resid)) + 
  geom_histogram(aes(y=after_stat(density)), bins = 20) +
  facet_wrap(Lab~., ncol = 4) +
  theme_bw() +
  xlab("Residuals of standard curve") +
  ylab("Density") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")
  )  +
  geom_vline(xintercept = 0, col = "red", linetype = "dashed", linewidth = 0.75) +
  geom_text(data = var, aes(x = 0.25, y = 6, label = Labels)) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold")) +
  ggtitle("B) Residuals of fitted model")
G2

# library(cowplot)
# 
# png("../Plots/standard_curves.png", width = 12, height = 8, units = "in", res = 350)
# plot_grid(
#   plot_grid(
#     G1 + theme(legend.position = "none")
#     , G2
#     , ncol = 1
#     , align = "hv")
#   , plot_grid(
#     get_legend(G1)
#     , ggplot() + theme_minimal()
#     , ncol =1)
#   , rel_widths = c(9,1)
# )
# dev.off()
############################################################################
preds_all =
  predict(conv_mod,
          newdata = data.frame(CT=Res$CT_NS, Lab=Res$Lab,
                               batch=as.factor(Res$Plate)),
          allow.new.levels = F)
preds_cens =
  predict(conv_mod,
          newdata = data.frame(CT=rep(40,nrow(Res)), Lab=Res$Lab,
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
write.table(x = control_dat, file = paste0(prefix_analysis_data, "/Analysis_Data/interim_control_dat.csv"), row.names = F, sep=',')

Res = dplyr::arrange(Res, Rand_date, ID, Time)
write.table(x = Res, file = paste0(prefix_analysis_data, "/Analysis_Data/interim_all_analysis.csv"), row.names = F, sep=',')

write.table(x = screen_failure, file = paste0(prefix_analysis_data, "/Analysis_Data/interim_screening_dat.csv"), row.names = F, sep=',')


Res = 
  Res %>% filter(Swab_ID != 'Saliva') %>% # remove the saliva samples
  mutate(Country= case_when(Site %in% c('th001','th057','th058') ~ 'Thailand',
                            Site == 'br003' ~ 'Brazil',
                            Site == 'la008' ~ 'Laos',
                            Site == 'pk001' ~ 'Pakistan'))
###### Fever data------------------------------------------------------------------
fever_data = read_csv(file = paste0(prefix_analysis_data, "/Analysis_Data/temperature_data.csv"))
fever_data = fever_data %>% 
  mutate(ID = Label,
         ax_temperature = fut_temp)
write.table(x = fever_data[, c('ID','Time','ax_temperature','Fever_Baseline')], 
            file = paste0(prefix_analysis_data, "/Analysis_Data/fever_interim.csv"), 
            row.names = F, sep=',', quote = F)
####################################################################################### 
################################################################################################
#***********************************************************************#
#*************************      Serology       ************************#
serology_data <- read.csv("../Analysis_data/Serology_estimated.csv") # Estimated from "parse_serology_data.R"
baseline_serology_data <- serology_data %>% 
  filter(Day == 0) %>%
  mutate(Baseline_log_IgG = Mean_log_IgG) %>%
  select("ID", "Baseline_log_IgG")

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
# Res_Remdesivir =
#   Res %>% filter(Trt %in% c('Remdesivir',"No study drug"),
#                  Rand_date < '2022-06-11',
#                  Country %in% c('Thailand','Brazil')) %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_Remdesivir, file = paste0(prefix_analysis_data, "/Analysis_Data/Remdesivir_analysis.csv"), row.names = F, sep=',',quote = F)
# 
# 
# fever_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/temperature_data.csv"))
# fever_data = fever_data %>% 
#   filter(Label %in% unique(Res_Remdesivir$ID)) %>%
#   mutate(ID = Label,
#          ax_temperature = fut_temp)
# write.table(x = fever_data[, c('ID','Time','ax_temperature','Fever_Baseline')], 
#             file = paste0(prefix_analysis_data, "/Analysis_Data/Remdesivir_fever.csv"), 
#             row.names = F, sep=',', quote = F)
# 
# 
# symptom_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/symptom_data.csv"))
# symptom_data = symptom_data %>% 
#   filter(Label %in% unique(Res_Remdesivir$ID)) %>%
#   mutate(ID = Label,
#          Any_symptom = sq_yn)
# write.table(x = symptom_data[, c('ID','Timepoint_ID','Any_symptom','heart_rate')], 
#             file = paste0(prefix_analysis_data, "/Analysis_Data/Remdesivir_symptoms.csv"), 
#             row.names = F, sep=',', quote = F)


#************************* Favipiravir Analysis *************************#
#* Thailand and Brazil
# Res_Favipiravir =
#   Res %>% filter(Trt %in% c('Favipiravir',"No study drug"),
#                  Rand_date < '2022-10-31',
#                  Country %in% c('Thailand','Brazil')) %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_Favipiravir, 
#             file =  paste0(prefix_analysis_data, "/Analysis_Data/Favipiravir_analysis.csv"), 
#             row.names = F, sep=',', quote = F)


#************************* Regeneron Analysis *************************#
#* Thailand only
# Res_REGN =
#   Res %>% filter(Trt %in% c('Regeneron',"No study drug"),
#                  Country == 'Thailand',
#                  Rand_date < '2022-10-21 00:00:00') %>%
#   arrange(Rand_date, ID, Time)
# 
# regeneron_mutations = read.csv(paste0(prefix_analysis_data, "/Analysis_Data/regeneron_mutations.csv"))
# regeneron_mutations <- regeneron_mutations %>% select(Patient_ID, test_regeneron)
# colnames(regeneron_mutations)[1] <- "ID"
# table(regeneron_mutations$test_regeneron, useNA = 'ifany')
# 
# Res_REGN = merge(Res_REGN, regeneron_mutations, by = 'ID', all.x = T)
# Res_REGN <- merge(Res_REGN, baseline_serology_data, by = "ID", all.x = T)
# 
# #impute variants
# Res_REGN$test_regeneron[Res_REGN$Rand_date > as.Date("2022-01-01") & Res_REGN$Rand_date < as.Date("2022-04-01") & is.na(Res_REGN$Variant2)] <- "G446S"
# Res_REGN$Variant2[Res_REGN$Rand_date > as.Date("2022-01-01") & Res_REGN$Rand_date < as.Date("2022-04-01") & is.na(Res_REGN$Variant2)] <- "BA.1"
# 
# Res_REGN$test_regeneron[Res_REGN$Rand_date > as.Date("2022-10-01") & is.na(Res_REGN$Variant2)] <- "G446S"
# Res_REGN$Variant2[Res_REGN$Rand_date > as.Date("2022-10-01") & is.na(Res_REGN$Variant2)] <- "BA.2.75"
# 
# Res_REGN$Resistant_test <- "Sensitive"
# Res_REGN$Resistant_test[Res_REGN$test_regeneron != "Wildtype"] <- "Resistant"
# 
# #impute IgG
# mean_baseline_IgG <- Res_REGN %>% ungroup() %>%
#   filter(Timepoint_ID==0, Rand_date > as.Date('2022-04-01')) %>%
#   distinct(ID, .keep_all = T) %>%
#   summarise(mean_baseline_IgG = mean(Baseline_log_IgG, na.rm = T)) %>% as.numeric
# 
# Res_REGN$Baseline_log_IgG[which(is.na(Res_REGN$Baseline_log_IgG))] <- mean_baseline_IgG
# 
# write.table(x = Res_REGN, file = '../Analysis_Data/REGN_analysis.csv', row.names = F, sep=',', quote = F)
# 
# 
# IDs <- Res_REGN %>% select(ID) %>% as.vector() %>% unlist() %>% unique()
# #Temperature data
# fever_REGN <-  fever_data %>% filter(Label %in% IDs)
# write.table(x = fever_REGN, file = '../Analysis_Data/REGN_fever_analysis.csv', row.names = F, sep=',', quote = F)
# #Symptom data
# symptom_REGN <-  symptom_data %>% filter(Label %in% IDs)
# write.table(x = symptom_REGN, file = '../Analysis_Data/REGN_symptom_analysis.csv', row.names = F, sep=',', quote = F)

#************************* Fluoxetine Analysis *************************#
#* Thailand added 1st April 2022; Brazil added 21st June 2022
#*  Stopped in all sites on 8th May 2023
# Res_Fluoxetine =
#   Res %>% filter(Trt %in% c('Fluoxetine',"No study drug"),
#                  (Country=='Thailand' & Rand_date > "2022-04-01 00:00:00") |
#                    (Country=='Brazil' & Rand_date > "2022-06-21 00:00:00") |
#                    (Country=='Laos' & Rand_date > "2022-06-21 00:00:00") |
#                    (Country=='Pakistan' & Rand_date > "2022-06-21 00:00:00"),
#                  Rand_date < "2023-05-09") %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_Fluoxetine, file = paste0(prefix_analysis_data, "/Analysis_Data/Fluoxetine_analysis.csv"), row.names = F, sep=',', quote = F)
#
#



# Res_Fluoxetine_meta =
#   Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir',
#                             'Molnupiravir',
#                             "No study drug",
#                             'Ivermectin',
#                             'Remdesivir',
#                             'Favipiravir',
#                             'Fluoxetine',
#                             'Regeneron'),
#                  Country %in% c('Thailand','Brazil','Laos','Pakistan'),
#                  Rand_date <= "2023-05-09 00:00:00") %>%
#   arrange(Rand_date, ID, Time)
# 
# write.table(x = Res_Fluoxetine_meta,
#             file = paste0(prefix_analysis_data, "/Analysis_Data/Fluoxetine_meta_analysis.csv"),
#             row.names = F, sep=',', quote = F)
# 


#************************* Paxlovid v Molnupiravir Analysis *************************#
#* Thailand only - this is used in the Lancet Infectious Diseases paper 2023
# Res_Paxlovid_Molnupiravir = 
#   Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir','Molnupiravir',"No study drug"),
#                  Rand_date > "2022-06-03 00:00:00",
#                  Rand_date < "2023-02-24 00:00:00",
#                  Country %in% c('Thailand')) %>%
#   arrange(Rand_date, ID, Time) %>% ungroup() 
# 
# write.table(x = Res_Paxlovid_Molnupiravir, file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_analysis.csv"), row.names = F, sep=',', quote = F)
# 
# ## Meta-analysis which we report in the publication using unblinded arms from the same site
# Res_Paxlovid_Molnupiravir_meta = 
#   Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir',
#                             'Molnupiravir',
#                             "No study drug",
#                             'Ivermectin',
#                             'Remdesivir',
#                             'Favipiravir'),
#                  Rand_date < "2023-02-24 00:00:00",
#                  Site =='th001') %>%
#   arrange(Rand_date, ID, Time) 
# 
# write.table(x = Res_Paxlovid_Molnupiravir_meta, 
#             file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_meta_analysis.csv"), 
#             row.names = F, sep=',', quote = F)
# 
# 
# fever_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/temperature_data.csv"))
# fever_data = fever_data %>% 
#   filter(Label %in% unique(Res_Paxlovid_Molnupiravir_meta$ID)) %>%
#   mutate(ID = Label,
#          ax_temperature = fut_temp)
# write.table(x = fever_data[, c('ID','Time','ax_temperature','Fever_Baseline')], 
#             file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_meta_fever.csv"), 
#             row.names = F, sep=',', quote = F)
# 
# 
# 
# symptom_data = read_csv(paste0(prefix_analysis_data, "/Analysis_Data/symptom_data.csv"))
# symptom_data = symptom_data %>% 
#   filter(Label %in% unique(Res_Paxlovid_Molnupiravir_meta$ID)) %>%
#   mutate(ID = Label,
#          Any_symptom = sq_yn)
# write.table(x = symptom_data[, c('ID','Timepoint_ID','Any_symptom','heart_rate')], 
#             file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_Molnupiravir_meta_symptoms.csv"), 
#             row.names = F, sep=',', quote = F)
# 

#************************* Paxlovid v No study drug - recent data only *************************#
#* Thailand only - this is used for internal analyses
# Res_Paxlovid_recent = 
#   Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir',"No study drug"),
#                  Rand_date > "2023-02-24 00:00:00",
#                  Country %in% c('Thailand')) %>%
#   arrange(Rand_date, ID, Time) %>% ungroup() 
# 
# write.table(x = Res_Paxlovid_recent, file = paste0(prefix_analysis_data, "/Analysis_Data/Paxlovid_recent_analysis.csv"), row.names = F, sep=',', quote = F)
# 

## vaccine data
# vacc_data_molnupiravir = vacc_data %>% filter(Label %in% Res_Paxlovid_Molnupiravir$ID) %>%
#   group_by(Label) %>%
#   mutate(vac_combos = paste(sort(unique(vc_name)), collapse = '/')) %>%
#   distinct(Label, .keep_all = T)
# table(vacc_data_molnupiravir$vac_combos)
#************************* Nitazoxanide Analysis *************************#
#* Brazil and Laos
Res_Nitazoxanide = 
  Res %>% filter(Trt %in% c('Nitazoxanide',"No study drug"),
                 Rand_date > "2022-06-03 00:00:00",
                 Country %in% c('Brazil','Laos', 'Pakistan'),
                 !ID %in% c("PLT-BR3-006", ## why are we excluding these patients at this stage - should be at the analysis stage?
                            "PLT-BR3-018",
                            "PLT-BR3-033",
                            "PLT-BR3-043")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_Nitazoxanide, file = paste0(prefix_analysis_data, "/Analysis_Data/Nitazoxanide_analysis.csv"), row.names = F, sep=',', quote = F)


#************************* Evusheld Analysis *************************#
#* Thailand added 2022-09-01; Brazil added 2022-10-31
# Res_Evusheld =
#   Res %>% filter(Trt %in% c('Evusheld',"No study drug"),
#                  (Country=='Thailand' & Rand_date > "2022-09-01 00:00:00" & Rand_date < "2023-07-05 00:00:00") |
#                    (Country=='Brazil' & Rand_date > "2022-10-31 00:00:00" & Rand_date < "2023-01-04 00:00:00")) %>%
#   arrange(Rand_date, ID, Time)
# 
#  evusheld_mutations = read.csv(paste0(prefix_analysis_data, "/Analysis_Data/evusheld_mutations.csv"))
#  evusheld_mutations <- evusheld_mutations %>% select(Patient_ID, test_evusheld)
#  colnames(evusheld_mutations)[1] <- "ID"
#  table(evusheld_mutations$test_evusheld, useNA = 'ifany')
# 
#  Res_Evusheld = merge(Res_Evusheld, evusheld_mutations, by = 'ID', all.x = T)
#  Res_Evusheld <- merge(Res_Evusheld, baseline_serology_data, by = "ID", all.x = T)
# 
#  Res_Evusheld$Rand_date <- as.Date(Res_Evusheld$Rand_date)
# 
#  #impute variants
#  Res_Evusheld$test_evusheld[Res_Evusheld$Site == "br003" & is.na(Res_Evusheld$Variant2)] <- "F486* and (R346* or K444*)"
#  Res_Evusheld$Variant2[Res_Evusheld$Site == "br003"  & is.na(Res_Evusheld$Variant2)] <- "BA.5"
# 
#  Res_Evusheld$test_evusheld[Res_Evusheld$Site == "th001" & is.na(Res_Evusheld$Variant2) & Res_Evusheld$Rand_date <= as.Date("2023-01-01")] <- "R346* or K444*"
#  Res_Evusheld$Variant2[Res_Evusheld$Site == "th001"  & is.na(Res_Evusheld$Variant2) & Res_Evusheld$Rand_date <= as.Date("2023-01-01")] <- "BA.2.75"
# 
#  Res_Evusheld$test_evusheld[Res_Evusheld$Site == "th001" & is.na(Res_Evusheld$Variant2) & Res_Evusheld$Rand_date > as.Date("2023-04-01")] <- "F486* and (R346* or K444*)"
#  Res_Evusheld$Variant2[Res_Evusheld$Site == "th001"  & is.na(Res_Evusheld$Variant2) & Res_Evusheld$Rand_date > as.Date("2023-04-01")] <- "XBB.1.5-like"
# 
#  Res_Evusheld$test_evusheld[is.na(Res_Evusheld$test_evusheld) & Res_Evusheld$Variant2 == "BA.2.75"]  <- "R346* or K444*"
#  Res_Evusheld$test_evusheld[is.na(Res_Evusheld$test_evusheld) & Res_Evusheld$Variant2 == "BA.5"]  <- "F486* and (R346* or K444*)"
# 
#  Res_Evusheld$Resistant_test <- "Partially resistant"
#  Res_Evusheld$Resistant_test[Res_Evusheld$test_evusheld == "F486* and (R346* or K444*)"] <- "Resistant"
# 
#  #impute IgG
#  mean_baseline_IgG <- Res_Evusheld %>% ungroup() %>%
#    filter(Timepoint_ID==0) %>%
#    distinct(ID, .keep_all = T) %>%
#    summarise(mean_baseline_IgG = mean(Baseline_log_IgG, na.rm = T)) %>% as.numeric
# 
#  Res_Evusheld$Baseline_log_IgG[which(is.na(Res_Evusheld$Baseline_log_IgG))] <- mean_baseline_IgG
# 
#  table(Res_Evusheld$Variant2, useNA = 'ifany')
#  table(Res_Evusheld$test_evusheld, useNA = 'ifany')
# 
#  write.table(x = Res_Evusheld, file = paste0(prefix_analysis_data, "/Analysis_Data/Evusheld_analysis.csv"), row.names = F, sep=',', quote = F)
# 
# IDs <- Res_Evusheld %>% select(ID) %>% as.vector() %>% unlist() %>% unique()
# #Temperature data
# fever_Evusheld <-  fever_data %>% filter(Label %in% IDs)
# write.table(x = fever_Evusheld, file = '../Analysis_Data/Evusheld_fever_analysis.csv', row.names = F, sep=',', quote = F)
# #Symptom data
# symptom_Evusheld <-  symptom_data %>% filter(Label %in% IDs)
# write.table(x = symptom_Evusheld, file = '../Analysis_Data/Evusheld_symptom_analysis.csv', row.names = F, sep=',', quote = F)


#************************* Ensitrelvir Analysis *************************#
#* Thailand and Laos added 2023-03-17

Res_Ensitrelvir = 
  Res %>% filter(Trt %in% c('Ensitrelvir',"No study drug",'Nirmatrelvir + Ritonavir'),
                 (Country=='Thailand' &
                    Rand_date >= "2023-03-17 00:00:00" &
                    Rand_date < "2024-04-22") |
                   (Country=='Laos' & 
                      Rand_date >= "2023-03-17 00:00:00"&
                      Rand_date < "2024-04-22")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_Ensitrelvir, file = paste0(prefix_analysis_data, "/Analysis_Data/Ensitrelvir_analysis.csv"), row.names = F, sep=',', quote = F)


# Res_Ensitrelvir_allpax = 
#   Res %>% filter(Trt %in% c('Ensitrelvir',"No study drug",'Nirmatrelvir + Ritonavir'),
#                  (Country=='Thailand') |
#                    (Country=='Laos')) %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_Ensitrelvir_allpax, file = paste0(prefix_analysis_data, "/Analysis_Data/Ensitrelvir_allpax_analysis.csv"), row.names = F, sep=',', quote = F)

#************************* Nirmatrelvir+Molnupiravir Analysis *************************#
#* Thailand, Brazil and Laos added 2023-05-29 
Res_MolPax = 
  Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir + Molnupiravir',"No study drug",'Nirmatrelvir + Ritonavir'),
                 Country %in% c('Thailand','Laos','Brazil'), 
                 Rand_date >= "2023-05-29 00:00:00" ) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_MolPax, file = paste0(prefix_analysis_data, "/Analysis_Data/MolPax_combination_analysis.csv"), row.names = F, sep=',', quote = F)



#************************* Hydroxychloroquine Analysis *************************#
#* Thailand added 2024-01-02
Res_HCQ = 
  Res %>% filter(Trt %in% c('Hydroxychloroquine',"No study drug"),
                 (Country=='Thailand' & Rand_date > "2024-01-01 00:00:00")) %>%
  arrange(Rand_date, ID, Time)
write.table(x = Res_HCQ, file = paste0(prefix_analysis_data, "/Analysis_Data/Hydroxychloroquine_analysis.csv"), row.names = F, sep=',', quote = F)


#************************* Ineffective Interventions *************************#
# Res_ineffective = 
#   Res %>% filter(Trt %in% c('Ivermectin',
#                             "Favipiravir",
#                             "Fluoxetine",
#                             "No study drug")) %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_ineffective, file = paste0(prefix_analysis_data, "/Analysis_Data/Ineffective_analysis.csv"), row.names = F, sep=',', quote = F)

#************************* No Study Drugs *************************#
# Res_noStudyDrugs = 
#   Res %>% filter(Trt %in% c("No study drug")) %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_noStudyDrugs, file = paste0(prefix_analysis_data, "/Analysis_Data/Res_noStudyDrugs.csv"), row.names = F, sep=',', quote = F)

#************************* Site TH01 only *************************#
# Res_TH1 <-  Res %>% filter(Site %in% c("th001")) %>%
#   arrange(Rand_date, ID, Time)
# write.table(x = Res_noStudyDrugs, file = paste0(prefix_analysis_data, "/Analysis_Data/Res_TH1.csv"), row.names = F, sep=',', quote = F)



#************************* Unblinded arm meta-analysis *************************#
Res_Unblinded_meta =
  Res %>% filter(Trt %in% c('Nirmatrelvir + Ritonavir',
                            'Molnupiravir',
                            "No study drug",
                            'Ivermectin',
                            'Remdesivir',
                            'Favipiravir',
                            'Regeneron'),
                 Country %in% c('Thailand','Brazil','Laos','Pakistan'),
                 Rand_date <= "2023-10-20 00:00:00"
  ) %>%
  arrange(Rand_date, ID, Time)

write.table(x = Res_Unblinded_meta,
            file = paste0(prefix_analysis_data, "/Analysis_Data/Unblinded_meta_analysis.csv"),
            row.names = F, sep=',', quote = F)



#************************* Unblinded all *************************#
# Res_Unblinded_all = 
#   Res %>% filter(!Trt %in% c('Nirmatrelvir + Ritonavir + Molnupiravir',
#                             'Nitazoxanide',
#                             'Ensitrelvir')) %>%
#   arrange(Rand_date, ID, Time) 
# 
# 
# write.table(x = Res_Unblinded_all, 
#             file = paste0(prefix_analysis_data, "/Analysis_Data/Unblinded_all_analysis.csv"), 
#             row.names = F, sep=',', quote = F)

