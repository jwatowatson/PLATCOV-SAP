##*********************************************************
##******************* Clinical database *******************
##*********************************************************
######################################################################
# 1. Loading clinical data
load_clinical_data <- function(){
  # load the clinical data
  writeLines('##########################################################################')
  writeLines('##########################################################################')
  
  writeLines('Reading clinical database from MACRO...')
  clin_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimEnrolment.dta"))
  clin_data$scrpassed[clin_data$Label=='PLT-TH1-557']=1
  
  writeLines(sprintf('Clinical dataset contains %s rows with %s unique Screening IDs and %s unique Patient IDs',
                     nrow(clin_data),
                     length(unique(clin_data$scrid)),
                     length(unique(clin_data$Label))
  ))
  writeLines('##########################################################################')
  
  # check duplications
  writeLines('### Clinical database: Checking multiple screening IDs:')
  duplicated_clin_data <- names(which(table(clin_data$Label) > 1))
  duplicated_clin_data <- duplicated_clin_data[duplicated_clin_data != ""]
  duplicated_clin_data
  writeLines(sprintf('%s patients have multiple screening ID',
                     length(duplicated_clin_data)))
  for(i in duplicated_clin_data){
    writeLines(sprintf('- Patient %s has multiple screening IDs: %s', 
                       duplicated_clin_data,
                       paste(clin_data$scrid[clin_data$Label == i], collapse = ", ")
                       ))
  writeLines('##########################################################################')
    
  }

  # manual correction for duplications
  ind <- which(clin_data$Label %in% duplicated_clin_data & is.na(clin_data$scrdat))
  if(length(ind>0)){clin_data <- clin_data[-ind,]}
  
  return(clin_data)
}
# -----------------------------------------------------------------------------------------------
# 2. Check data entries
check_MACRO_clinical_database <- function(clin_data, rand_app_data){
  writeLines('### Clinical database: Checking MACRO data entry progress for baseline information:')
  on_macro_yn <- rand_app_data$ID %in% clin_data$Label
  
  writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient baseline data on MACRO: %s', 
                     nrow(rand_app_data) + 19, # 19 = number of patients from site 057 and 058
                     clin_data %>% filter(!is.na(Label)) %>% distinct(Label) %>% nrow()
                     )
             )
  
  IDs_pending <- rand_app_data %>% filter(!on_macro_yn) %>% pull(ID)
  writeLines(sprintf('Baseline data entries for the following %s patients are pending:',
                     length(IDs_pending)))
  print(IDs_pending)
  writeLines('##########################################################################')
  
  return(IDs_pending)
}

# -----------------------------------------------------------------------------------------------
# 3. Check screening failure (Clinical data) + Exporting the patients who failed the screening
check_screen_failure <- function(clin_data){
  # check screening failure
  writeLines('### Clinical database: Checking missing screening failure information:')
  
  missing_scr_failure <- clin_data %>% filter(is.na(scrpassed)) %>% select(scrid, Label, rangrp, scrpassed)
  
  writeLines(sprintf('%s patients have missing information on screening results:\n',  # they are currently removed from further analyses
                     nrow(missing_scr_failure)
                     ))
  print(missing_scr_failure)
  writeLines('##########################################################################')
  
  # exporting screening failure data
  ind = !is.na(clin_data$scrpassed) & clin_data$scrpassed==0
  screen_failure = clin_data[ind, c("Trial","Site","scrid","scrdat", "scrpassed","reason_failure","scrnote")]
  screen_failure$reason_failure = sjlabelled::as_character(screen_failure$reason_failure)
  write.csv(x = screen_failure, file = '../Analysis_Data/screening_failure.csv')
  
  # excluding failed patients from the clinical dataset
  if(length(ind > 0)){ clin_data <- clin_data[!ind,]}
  
  return(clin_data)
}
# -----------------------------------------------------------------------------------------------
# 4. Check randomisation information
check_randomisation_info <- function(clin_data){
  # Passed the screening but not randomised?
  writeLines('### Clinical database: Checking randomisation information:')
  ind <- clin_data$scrpassed == 1 & is.na(clin_data$rangrp) & clin_data$Label == ""
  ind[is.na(ind)] <- T
  writeLines(sprintf('%s patients marked as passed the screening but have not been randomised:\n', 
                    sum(ind)))
  clin_data[ind,] %>% select(scrid, scrpassed, Label, rangrp) %>% print()
  clin_data <- clin_data[!ind,]
  writeLines('##########################################################################')
  
  # Not randomised with missing screening status?
  writeLines('### Clinical database: Checking randomisation information:')
  ind <- (is.na(clin_data$scrpassed) | clin_data$scrpassed == 1) & is.na(clin_data$rangrp) & clin_data$Label == ""
  ind[is.na(ind)] <- T
  writeLines(sprintf('Patient %s does not have the information on screening status and not randomised:\n', 
                     sum(ind)))
  clin_data[ind,] %>% select(scrid, scrpassed, Label, rangrp) %>% print()
  clin_data <- clin_data[!ind,]
  writeLines('##########################################################################')
  
  return(clin_data)
}
# -----------------------------------------------------------------------------------------------
# 5. Check sex data
check_sex <- function(clin_data, IDs_pending, rand_app_data){
  clin_data$Sex <- plyr::mapvalues(x = as.numeric(clin_data$sex),
                                   from = c(1,2),
                                   to = c('Male','Female'))
  # Missing Sex data?
  ind <- is.na(clin_data$sex) #& (clin_data$Label %in% IDs_pending)
  writeLines('### Clinical database: Checking missing sex information:')
  writeLines(sprintf('%s patients with data on MACRO have missing sex information:', 
                     sum(ind)))
  clin_data[ind,] %>% pull(Label) %>% print()# select(scrid, Label, rangrp, sex) %>% print(n = Inf)
  
  writeLines('##########################################################################')
  # Mismatched Sex data?
  writeLines('### Clinical database: Checking mismatched sex information from Randomisation database:')
  check_data <- merge(clin_data[,c("Label", "Sex")], rand_app_data, by.y = "ID",  by.x = "Label", all.x = T)
  check_data <- check_data %>% filter(!check_data$Label %in% IDs_pending)
  
  ind <- (check_data$Sex != check_data$sex) & !is.na(check_data$sex)  & !is.na(check_data$Sex)
  writeLines(sprintf('%s patients have mismatched sex data:', sum(ind))) 
  id_sex_mismatched <- check_data[ind, ] %>% pull(Label) 
  check_data %>% filter(Label %in% id_sex_mismatched) %>% 
    mutate(MACRO = Sex, Randomisation_DB = sex) %>% 
    select(Label, Treatment, MACRO, Randomisation_DB) %>% 
    print()
  writeLines('##########################################################################')
  #manual correction
  writeLines('### [MANUAL CORRECTIONS]: Subsequent analyses use sex information from the randomisation database')
  for(i in id_sex_mismatched){clin_data$Sex[clin_data$Label == i] <- rand_app_data$sex[rand_app_data$ID == i]}
  writeLines('##########################################################################')
  
  return(clin_data)
}
######################################################################
# 6. Check age data
check_age <- function(clin_data, IDs_pending, rand_app_data){
  # A special case
  clin_data$age_yr[clin_data$Label=='PLT-TH57-002'] = 21 # A special case
  
  writeLines('### Clinical database: Checking missing age information:')
  # Missing Age data?
  ind = is.na(clin_data$age_yr) & !is.na(clin_data$dob_my)
  
  writeLines('##########################################################################')
  # calculate age at randomisation 
  for(i in which(ind)){
    clin_data$age_yr[i] = trunc((parse_date_time(clin_data$dob_my[i], c('%d/%m/%Y', '%m/%Y')) %--% parse_date_time(clin_data$randat[i],  c('%Y-%m-%d')))/years(1))
  }
  
  check_data <- merge(clin_data[,c("Label", "randat", "dob_my", "age_yr")], rand_app_data, by.y = "ID",  by.x = "Label", all.x = T)
  check_data <- check_data %>% filter(!check_data$Label %in% IDs_pending)
  
  check_data$age <- floor(as.numeric(check_data$age))
  check_data$age_dob_missing <- is.na(check_data$age_yr)
  check_data$age_agree_yn <- check_data$age == check_data$age_yr
  id_age_mismatched <- check_data %>% filter(!age_agree_yn) %>% pull(Label)
  
  writeLines(sprintf('%s patients have no information on age/date of birth on MACRO:', 
                     sum(check_data$age_dob_missing)))
  check_data %>% filter(age_dob_missing ) %>% pull(Label) %>%  as.character() %>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s patients have mismatched information on age/date with the randomisation database:', 
                     sum(!check_data$age_agree_yn, na.rm = T)               
  ))
  check_data %>% filter(Label %in% id_age_mismatched) %>% 
    mutate(dob_my = dob_my, MACRO = age_yr, Randomisation_DB = age) %>% 
    select(Label, Treatment, randat, dob_my, MACRO, Randomisation_DB) %>% 
    print()
  writeLines('##########################################################################')
  #manual correction
  writeLines('### [MANUAL CORRECTIONS]: Subsequent analyses use age information from the randomisation database')
  for(i in id_age_mismatched){clin_data$age_yr[clin_data$Label == i] <- rand_app_data$age[rand_app_data$ID == i]}
  writeLines('##########################################################################')
  clin_data$age_yr <- as.numeric(clin_data$age_yr)
  
  
  writeLines('### Clinical database: Checking distributions of age:')
  summary(clin_data$age_yr) %>% print()
  writeLines('##########################################################################')
  
  return(clin_data)
}
######################################################################
# 6. Check symptom onset data
check_symptom_onset <- function(clin_data, IDs_pending, rand_app_data){
  writeLines('### Clinical database: Checking missing symptom onset information:')
  ### Check if symptom onset data is missing
  ind = !is.na(clin_data$cov_symphr) & is.na(clin_data$cov_sympday)
  if(sum(ind)>0) clin_data$cov_sympday[ind] = clin_data$cov_symphr[ind]/24 #If day is missing >>> divide hours by 24
  
  ind = is.na(clin_data$cov_sympday)
  writeLines(sprintf('%s patients have no information on symptom onset (in hours or days) on MACRO:', 
                     sum(ind)))
  clin_data[ind,] %>% pull(Label) %>% as.character() %>% print()
  check_data <- merge(clin_data[,c("Label", "cov_sympday")], rand_app_data,  by.y = "ID",  by.x = "Label", all.x = T)
  writeLines('##########################################################################')
  
  
  ind_4days <- check_data$cov_sympday > 4 & !is.na(check_data$cov_sympday)
  writeLines(sprintf('%s patients have symptom onset on MACRO greater than 4 days:', 
                     sum(ind_4days, na.rm = T)))
  check_data[ind_4days,] %>% select(Label, Treatment, Rand_Time, cov_sympday) %>% print()
  writeLines('##########################################################################')
  
  # Assign Symptomatic day = 2 for missing data 
  writeLines('### [MANUAL CORRECTIONS]: Symptom onset of 2 days are assigned for missing data')
  clin_data$cov_sympday[clin_data$Label %in% (clin_data[ind,] %>% pull(Label))] <- 2
  
  
  # Assign Symptomatic day = 4 for data greater than 4 days 
  writeLines('### [MANUAL CORRECTIONS]: Symptom onset of 4 days are assigned for onset > 4 days')
  clin_data$cov_sympday[clin_data$Label %in% (check_data[ind_4days,] %>% pull(Label))] <- 4
  
  writeLines('##########################################################################')
  
  writeLines('### Clinical database: Checking distributions of symptom onset:')
  summary(clin_data$cov_sympday) %>% print()
  writeLines('##########################################################################')
  
  return(clin_data)
}
######################################################################
# 6. Check weight and height data
check_weight_height <- function(clin_data, IDs_pending){
  writeLines('### Clinical database: Checking missing data for weight/height:')
  ### Check if weight and height data are missing
  clin_data$BMI = clin_data$weight/(clin_data$height/100)^2
  clin_data$Weight = clin_data$weight
  
  ind_weight_height_missing <- is.na(clin_data$BMI)
  writeLines(sprintf('%s patients have no information on weights or heights on MACRO:', 
                     sum(ind_weight_height_missing)))
  clin_data[ind_weight_height_missing, ] %>% pull(Label) %>% as.character() %>% print()
  writeLines('##########################################################################')
  
  
  writeLines('### Clinical database: Checking outliers for weight/height:')
  weight_height_outlier <- abs(clin_data$BMI - mean(clin_data$BMI, na.rm = T)) > 6*sd(clin_data$BMI, na.rm = T) & !is.na(clin_data$BMI)
  writeLines(sprintf('Weights/heights of %s patients are outliers:', 
                     sum(weight_height_outlier)))
  clin_data[weight_height_outlier, ] %>% select(Label, Sex, weight, height,BMI) %>% print()
  writeLines('##########################################################################')
  writeLines('### [MANUAL CORRECTIONS]: Outliers are likely from switching weights and heights')
  for(i in clin_data[weight_height_outlier, ] %>% pull(Label) %>% as.character()){
    weight_correct <- clin_data$height[clin_data$Label == i]
    height_correct  <- clin_data$weight[clin_data$Label == i]
    clin_data$height[clin_data$Label == i] <- height_correct
    clin_data$weight[clin_data$Label == i] <- weight_correct
  }
  
  clin_data$BMI = clin_data$weight/(clin_data$height/100)^2
  writeLines('##########################################################################')
  writeLines('### Clinical database: Checking distributions of BMI after manual corrections:')
  summary(clin_data$BMI) %>% print()
  hist(clin_data$BMI)
  writeLines('##########################################################################')
  
 return(clin_data)
  
}

######################################################################
#7. Check randomisation date and time data
check_rand_date_time <- function(clin_data, IDs_pending, rand_app_data){
  writeLines('### Clinical database: Checking randomisation date and time:')
  ### Check if randomisation date and time is correct
  clin_data$Rand_date_time = NA
  for(i in 1:nrow(clin_data)){
    if(!is.na(clin_data$randat[i]) & !is.na(clin_data$rantim[i])){
      clin_data$Rand_date_time[i] = as.character(as.POSIXct(paste(clin_data$randat[i],
                                                                  clin_data$rantim[i], sep=' ')))
    }
  }
  
  check_data <- merge(rand_app_data, clin_data[,c("Label", "Rand_date_time")], by.x = "ID",  by.y = "Label", all.x = T)
  check_data$Rand_diffs = apply(check_data[,c('Rand_date_time','Rand_Time_TZ')],1, function(x) difftime(x[1], x[2], units='mins'))
  
  # Check missing data
  writeLines('### Clinical database: Checking if randomisation date and time are missing:')
  rand_date_missing <- is.na(check_data$Rand_date_time)
  
  writeLines(sprintf('%s patients have no information on randomisation date and time on MACRO:', 
                     sum(rand_date_missing)))
  check_data[rand_date_missing, ] %>% pull(ID) %>% as.character() %>% print()
  writeLines('##########################################################################')
  
  # Check if the differences is more than 5 mins between both databases
  writeLines('### Clinical database: Checking if more than 5 min difference in randomisation date/time between MACRO and Randomisation database')
  randtime_diff_exceed <- check_data$Rand_diffs > 5 & !is.na(check_data$Rand_diffs)
  writeLines(sprintf('%s patients have more than 5 minutes difference in randomisation date/time:', 
                     sum(randtime_diff_exceed)))
  check_data[randtime_diff_exceed, ] %>% 
    mutate(MACRO = Rand_date_time, Randomisation_DB = Rand_Time_TZ) %>% 
    select(ID, Treatment, MACRO, Randomisation_DB, Rand_diffs) %>% print()
  writeLines('##########################################################################')
  writeLines('### [MANUAL CORRECTIONS]: Using randomisation date and time from randomisation database in further analyses')
  ID_time_errors <- check_data[rand_date_missing|randtime_diff_exceed,] %>% pull(ID) %>% as.character()
  for(i in ID_time_errors){clin_data$Rand_date_time[clin_data$Label == i] <- rand_app_data$Rand_Time_TZ[rand_app_data$ID == i]}
  writeLines('##########################################################################')
  
  return(clin_data)
}
######################################################################
#8. Check randomisation arms
check_rand_arms <- function(clin_data, IDs_pending, rand_app_data){
  
  ### Check if randomisation arms matched 
  clin_data$rangrp = sjlabelled::as_character(clin_data$rangrp)
  clin_data$rangrp[clin_data$rangrp=='Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir'
  clin_data$rangrp[clin_data$rangrp=='Molnupiravir and Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir + Molnupiravir'
  
  check_data <- merge(rand_app_data, clin_data[,c("Label", "rangrp")], by.x = "ID",  by.y = "Label", all.x = T)
  
  writeLines('### Clinical database: Checking if randomisation arms are missing:')
  rangrp_missing <- is.na(check_data$rangrp)
  writeLines(sprintf('%s patients have no information on treatment arms on MACRO:', 
                     sum(rangrp_missing)))
  check_data[rangrp_missing, ] %>% pull(ID) %>% as.character() %>% print()
  writeLines('##########################################################################')
  
  
  writeLines('### Clinical database: Checking if randomisation arms matched with the randomisation database:')
  rangrp_disagree <- (check_data$Treatment != check_data$rangrp) & !is.na(check_data$rangrp)
  writeLines(sprintf('%s patients have mismatched treatment arms:', 
                     sum(rangrp_disagree)))
  check_data[rangrp_disagree, ] %>% mutate(MACRO = rangrp, randomisation_DB = Treatment) %>% 
    select(ID, sex, MACRO, randomisation_DB)  %>% print()
  writeLines('##########################################################################')
  writeLines('### [MANUAL CORRECTIONS]: Using randomisation arms from randomisation database in further analyses')
  ID_rangrp_errors <- check_data[rangrp_missing|rangrp_disagree,] %>% pull(ID) %>% as.character()
  for(i in ID_rangrp_errors){clin_data$rangrp[clin_data$Label == i] <- rand_app_data$Treatment[rand_app_data$ID == i]}
  writeLines('##########################################################################')
  
  writeLines('### Clinical database: Checking frequencies of randomisation arms after manual corrections:')
  table(clin_data$rangrp, useNA = "always") %>% as.data.frame() %>% print()
  writeLines('##########################################################################')
  return(clin_data)
}
