##*********************************************************
##******************* Clinical database *******************
##*********************************************************
library(dplyr)
library(tidyr)


######################################################################
# 1. Loading clinical data
load_clinical_data <- function(){
  # load the clinical data
  writeLines('Reading clinical database from MACRO...\n')
  clin_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimEnrolment.dta"))
  clin_data$scrpassed[clin_data$Label=='PLT-TH1-557']=1
  
  # check duplications
  writeLines('Clinical database: Checking duplicated patient IDs:\n')
  duplicated_clin_data <- names(which(table(clin_data$Label) > 1))
  duplicated_clin_data <- duplicated_clin_data[duplicated_clin_data != ""]
  duplicated_clin_data
  for(i in duplicated_clin_data){
    writeLines(sprintf('- Patient %s has multiple screening IDs: %s \n', 
                       duplicated_clin_data,
                       paste(clin_data$scrid[clin_data$Label == i], collapse = ", ")
                       ))
  }

  # manual correction for duplications
  ind <- which(clin_data$Label %in% duplicated_clin_data & is.na(clin_data$scrdat))
  if(length(ind>0)){clin_data <- clin_data[-ind,]}
  
  return(clin_data)
}

# 2. Check screening failure (Clinical data)
check_screen_failure <- function(clin_data){
  # check screening failure
  writeLines('Clinical database: Checking missing screening failure information:\n')
  
  missing_scr_failure <- clin_data %>% filter(is.na(scrpassed)) %>% select(scrid, Label, rangrp, scrpassed)
  
  writeLines(sprintf('%s patients has missing information on screening results:\n',  # they are currently removed from further analyses
                     nrow(missing_scr_failure)
                     ))
  print(missing_scr_failure)
  
  
  # exporting screening failure data
  ind = !is.na(clin_data$scrpassed) & clin_data$scrpassed==0
  screen_failure = clin_data[ind, c("Trial","Site","scrid","scrdat", "scrpassed","reason_failure","scrnote")]
  screen_failure$reason_failure = sjlabelled::as_character(screen_failure$reason_failure)
  write.csv(x = screen_failure, file = '../Analysis_Data/screening_failure.csv')
  
  # excluding failed patients from the clinical dataset
  if(length(ind > 0)){ clin_data <- clin_data[!ind,]}
  
  return(clin_data)
}



sink("Queries/console_output.txt", split = T)
clin_data <- load_clinical_data()
clin_data <- check_screen_failure(clin_data)
sink()













######################################################################
clin_data <- load_clinical_data()






##########  Preparing data for checking
### Check which screening status are missing
ind <- is.na(clin_data$scrpassed) & clin_data$Label %in% check_data$ID
writeLines(sprintf('This patient %s does not have the information about screening status', 
                   clin_data$scrid[ind]))
check_data$screen_status_missing <- (check_data$ID %in% clin_data$Label[ind])


##### removing patients that not randomised
ind <- clin_data$scrpassed == 1 & is.na(clin_data$rangrp) & clin_data$Label == ""
ind[is.na(ind)] <- T
writeLines(sprintf('This patient %s passed the screening but never randomised', 
                   clin_data$scrid[ind]))

ind <- (is.na(clin_data$scrpassed) | clin_data$scrpassed == 1) & is.na(clin_data$rangrp) & clin_data$Label == ""
writeLines(sprintf('This patient %s does not have the information on screening status and not randomised', 
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
# Using age information from randomisation database in further analyses
age_yr_missing <- clin_data$Label[is.na(clin_data$age_yr)]
for(i in 1:length(age_yr_missing)){
  clin_data$age_yr[clin_data$Label == age_yr_missing[i]] <- check_data$age[check_data$ID == age_yr_missing[i]]
}
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
clin_data$cov_sympday[clin_data$cov_sympday > 4] <- 4

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

# Manual correction 
clin_data$Weight[clin_data$Label == "PLT-TH1-1341"] <- clin_data$height[clin_data$Label == "PLT-TH1-1341"]
clin_data$height[clin_data$Label == "PLT-TH1-1341"] <- clin_data$weight[clin_data$Label == "PLT-TH1-1341"]
clin_data$BMI[clin_data$Label == "PLT-TH1-1341"] <- clin_data$Weight[clin_data$Label == "PLT-TH1-1341"]/(clin_data$height[clin_data$Label == "PLT-TH1-1341"]/100)^2

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
writeLines(sprintf('More than 5 min difference in rand time for %s, MACRO = %s and Randomisation = %s: (%s mins)', 
                   check_data$ID[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)],
                   check_data$Rand_date_time[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)],
                   check_data$Rand_Time_TZ[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)],
                   round(check_data$Rand_diffs[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)],1)
))
ID_exceed_5mins <- check_data$ID[check_data$randtime_diff_exceed & !is.na(check_data$randtime_diff_exceed)]

# Using randomisation date and time from randomisation database in further analyses
for(i in 1:length(ID_exceed_5mins)){
  clin_data$Rand_date_time[which(clin_data$Label == ID_exceed_5mins[i])] <- rand_app_data$Rand_Time_TZ[rand_app_data$ID == ID_exceed_5mins[i]]
}

ID_missing_date <- check_data$ID[is.na(check_data$Rand_date_time)]
for(i in 1:length(ID_missing_date)){
  clin_data$Rand_date_time[which(clin_data$Label == ID_missing_date[i])] <- rand_app_data$Rand_Time_TZ[rand_app_data$ID == ID_missing_date[i]]
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

ID_missing_trt <- check_data$ID[is.na(check_data$rangrp)]
for(i in 1:length(ID_missing_date)){
  clin_data$rangrp[which(clin_data$Label == ID_missing_trt[i])] <- rand_app_data$Treatment[rand_app_data$ID == ID_missing_trt[i]]
}

check_data <- check_data[,-which(names(check_data) %in% c("rangrp"))]

##########  --- Temperature data --- ########## 