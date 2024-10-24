##*********************************************************
##******************* Temperature database *******************
##*********************************************************
###########################################################################
# 1. Load temperature data
load_temp_data <- function(rand_app_data){
  ##########  --- Temperature data --- ########## 
  writeLines('##########################################################################')
  writeLines('##########################################################################')
  writeLines('Reading temperature database from MACRO...')
  writeLines('##########################################################################')
  temp_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimFUTemp.dta"))
  id_data = temp_data %>% distinct(Label) %>% pull(Label) %>% as.character()
  
  writeLines('### Temperature database: Checking MACRO data entry progress:')
  writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient temperature data on MACRO: %s', 
                     nrow(rand_app_data) + 19, # 19 = number of patients from site 057 and 058
                     length(id_data)
  )
  )
  writeLines('##########################################################################')
  id_missing <- rand_app_data$ID[which(!rand_app_data$ID %in% id_data)]
  
  writeLines(sprintf('%s patients have no temperature data on any visits:', 
                     length(id_missing))
  )
  print(id_missing)
  writeLines('##########################################################################')
  writeLines('### Temperature database: Missing data by visits:')
  temp_data %>%
    group_by(visit) %>%
    summarise(missing_AM = sum(is.na(fut_amdat)),
              missing_PM = sum(is.na(fut_pmdat))) %>% print()
  writeLines('##########################################################################')
  
  return(temp_data)
  
}
# -----------------------------------------------------------------------------------------------
# 2. Preparing temperature data
prep_tempdata <- function(temp_data, clin_data){
  x1 = temp_data[, c('Site', 'Label', 'visit', 'fut_amdat', "fut_amtim", "fut_amtemp")]
  x2 = temp_data[, c('Site', 'Label', 'visit', 'fut_pmdat', "fut_pmtim", "fut_pmtemp")]
  
  colnames(x1) = colnames(x2) = 
    c('Site', 'Label', 'visit', 'fut_dat', "fut_tim", "fut_temp")
  temp_data = rbind(x1, x2)
  
  temp_data = temp_data[(!is.na(temp_data$fut_dat) & 
                           !is.na(temp_data$fut_tim) &
                           !is.na(temp_data$fut_temp)), ]
  
  temp_data$temp_time = 
    as.character(as.POSIXct(apply(temp_data[, c('fut_dat','fut_tim')],
                                  1, paste, collapse =' '),
                            format = '%F %H:%M:%S'))
  temp_data = merge(temp_data, 
                    clin_data, 
                    by = c('Label','Site'),
                    all.x=T)
  
  temp_data = temp_data %>% filter(!is.na(scrpassed), 
                                   scrpassed==1, 
                                   !is.na(randat),
                                   rantim != '') %>% 
    ungroup() %>%
    mutate(rangrp = sjlabelled::as_character(rangrp),
           Rand_date_time = as.character(as.POSIXct(paste(randat,
                                                          rantim, sep=' '))),
           Timepoint_ID = gsub(x = visit,replacement = '',pattern='D'),
           Timepoint_ID = gsub(x = Timepoint_ID, replacement = '', pattern = 'H1'),
           Timepoint_ID = as.numeric(Timepoint_ID)) %>%  filter(!is.na(temp_time)) %>%
    arrange(Label, temp_time, Rand_date_time)
  
  temp_data$Label[temp_data$Rand_date_time > temp_data$temp_time]
  for(i in 1:nrow(temp_data)){
    temp_data$Time[i]=as.numeric(difftime(temp_data$temp_time[i],
                                          temp_data$Rand_date_time[i], 
                                          units = 'days'))
  }
  
  temp_data = temp_data %>% group_by(Label) %>%
    mutate(Fever_Baseline = ifelse(any(fut_temp>37 & Time<1), 1, 0))
  
  temp_data$temp_time_date <- sub(" .*", "", temp_data$temp_time)
  temp_data$temp_time_time <- sub(".+? ", "", temp_data$temp_time)
  
  writeLines('##########################################################################')
  writeLines('### Temperature database: Summarize temperature on each visits:')
  
  temp_data %>%
    group_by(visit) %>%
    summarise(min = min(fut_temp),
              Q1 = quantile(fut_temp, 0.25),
              Q2 = quantile(fut_temp, 0.5),
              Q3 = quantile(fut_temp, 0.75),
              max = max(fut_temp)) %>%
    print()
  
  
  write.csv(x = temp_data, file = '../Analysis_Data/temperature_data.csv', row.names = F, quote = F)
  return(temp_data)
}

# -----------------------------------------------------------------------------------------------
#3. Check timepoint ID
check_time_temp <- function(fever_data){
  writeLines('##########################################################################')
  writeLines('### Temperature database: Checking temperature data and time')
  
  FUtemp_checktime <- fever_data[abs(fever_data$Timepoint_ID - fever_data$Time) > 1, c("Label", "visit", "Rand_date_time", "temp_time", 
                                                                                       "Timepoint_ID", "Time")]
  
  FUtemp_checktime <- FUtemp_checktime %>% mutate(expected_date = as.Date(Rand_date_time) + Timepoint_ID)
  
  writeLines('##########################################################################')
  
  negative_time <-   FUtemp_checktime %>% filter(Time < 0)
  writeLines(sprintf('%s patients have negative time of temperature measurement:',
                     nrow(negative_time)
  ))
  negative_time %>% arrange(Label) %>% print()
  
  writeLines('##########################################################################')
  D0 <-  FUtemp_checktime %>% filter(Timepoint_ID  == 0 & Time > 0)
  writeLines(sprintf('On visit D0, %s patients have mismatched time of temperature measurement:',
                     nrow(D0)
  ))
  D0 %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf)
  
  writeLines('##########################################################################')
  D1_D7 <-  FUtemp_checktime %>% filter(Timepoint_ID  %in% 1:7 & Time > 0)
  writeLines(sprintf('On visit D1 to D7, %s patients have mismatched time of temperature measurement:',
                     nrow(D1_D7)
  ))
  D1_D7 %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf)
  
  writeLines('##########################################################################')
  D10_D14 <-  FUtemp_checktime %>% filter(Timepoint_ID   >=10)
  writeLines(sprintf('On visit D10 and D14, %s patients have mismatched time of temperature measurement:',
                     nrow(D10_D14)
  ))
  D10_D14 %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf)
  
  writeLines('##########################################################################')
  missing_temp_time <-  FUtemp_checktime %>% filter(is.na(Time))
  writeLines(sprintf(' %s patients have missing time of temperature measurement:',
                     nrow(missing_temp_time)
  ))
  missing_temp_time %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf)
  
}
# -----------------------------------------------------------------------------------------------
#4. Merging baseline fever
add_baseline_fever <- function(clin_data, fever_data) {
  # Extract fever at baseline
  fever_data_baseline = fever_data %>% distinct(Label, .keep_all = T)
  clin_data <- merge(clin_data, fever_data_baseline[, c('Label','Fever_Baseline')], by='Label', all = T)
  
  writeLines('##########################################################################')
  writeLines('### Temperature database: Checking missing baseline fever information:')
  missing_baseline_temp <- clin_data %>% filter(is.na(Fever_Baseline))
  writeLines(sprintf('%s patients has no information on baseline fever:', 
                     nrow(missing_baseline_temp)))
  missing_baseline_temp %>% pull(Label) %>% print()
  
  return(clin_data)
}