##*********************************************************
##******************* Vital database *********************
##*********************************************************
###########################################################################
# -----------------------------------------------------------------------------------------------
# 1. Load vital sign data
load_vital_data <- function(rand_app_data){
  writeLines('##########################################################################')
  writeLines('##########################################################################')
  writeLines('Reading vital sign database from MACRO...')
  writeLines('##########################################################################')
  vita_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimVitalSigns.dta"))
  
  id_data = vita_data %>% distinct(Label) %>% pull(Label) %>% as.character()
  
  writeLines('### Vital sign database: Checking MACRO data entry progress:')
  writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient vital sign data on MACRO: %s', 
                     nrow(rand_app_data) + 19, # 19 = number of patients from site 057 and 058
                     length(id_data)
  )
  )
  
  writeLines('##########################################################################')
  id_missing <- rand_app_data$ID[which(!rand_app_data$ID %in% id_data)]
  
  writeLines(sprintf('%s patients have no vital sign data on any visits:', 
                     length(id_missing))
  )
  print(id_missing)

  writeLines('##########################################################################')
  writeLines('### Vital sign database: Missing data by visits:')
  vita_data %>%
    group_by(visit) %>%
    summarise(missing = sum(is.na(vs_temp))) %>% print()
  writeLines('##########################################################################')
  
  return(vita_data)
}
# -----------------------------------------------------------------------------------------------
# 2. Preparing vital sign data
prep_vitadata <- function(vita_data, clin_data){
  
  vita_data = vita_data[(!is.na(vita_data$vd_dat) & 
                           !is.na(vita_data$vs_tim) &
                           !is.na(vita_data$vs_temp)), ]
  
  vita_data$vita_time = 
    as.character(as.POSIXct(apply(vita_data[, c('vd_dat','vs_tim')],
                                  1, paste, collapse =' '),
                            format = '%F %H:%M:%S'))
  
  vita_data = merge(vita_data, 
                    clin_data[,c('Label', "Site", 'rangrp', 'scrpassed', 'randat', 'rantim', 
                                 'Rand_date_time')], 
                    by = c('Label','Site'),
                    all.x=T)
  
  
  vita_data = vita_data %>% filter(!is.na(scrpassed), 
                                   scrpassed==1, 
                                   !is.na(randat),
                                   rantim != '') %>% 
    ungroup() %>%
    mutate(Timepoint_ID = gsub(x = visit,replacement = '',pattern='D'),
           Timepoint_ID = gsub(x = Timepoint_ID, replacement = '', pattern = 'H1'),
           Timepoint_ID = as.numeric(Timepoint_ID)) %>%  filter(!is.na(vita_time)) %>%
    arrange(Label, vita_time, Rand_date_time)
  
  vita_data$Label[vita_data$Rand_date_time > vita_data$vita_time]
  
  for(i in 1:nrow(vita_data)){
    vita_data$Time[i]=as.numeric(difftime(vita_data$vita_time[i],
                                          vita_data$Rand_date_time[i], 
                                          units = 'days'))
  }
  
  return(vita_data)
}
# -----------------------------------------------------------------------------------------------
# 3. Check vital sign data
check_vita_time <- function(vita_data){
  writeLines('##########################################################################')
  writeLines('### Vital sign database: Checking temperature data and time')
  
  VS_checktime <- vita_data[abs(vita_data$Timepoint_ID - vita_data$Time) > 1, c("Label", "visit", "Rand_date_time", "vita_time", 
                                                                                       "Timepoint_ID", "Time")]
  
  VS_checktime <- VS_checktime %>% mutate(expected_date = as.Date(Rand_date_time) + Timepoint_ID)
  
  writeLines('##########################################################################')
  
  negative_time <-   VS_checktime %>% filter(Time < 0)
  writeLines(sprintf('%s patients have negative time of vital sign measurement:',
                     nrow(negative_time)
  ))
  negative_time %>% arrange(Label) %>% print(n = Inf, na.print = "NA")
  
  writeLines('##########################################################################')
  D0 <-  VS_checktime %>% filter(Timepoint_ID  == 0 & Time > 0)
  writeLines(sprintf('On visit D0, %s patients have mismatched time of vital sign measurement:',
                     nrow(D0)
  ))
  D0 %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf, na.print = "NA")
  
  writeLines('##########################################################################')
  D1_D7 <-  VS_checktime %>% filter(Timepoint_ID  %in% 1:7 & Time > 0)
  writeLines(sprintf('On visit D1 to D7, %s patients have mismatched time of vital sign measurement:',
                     nrow(D1_D7)
  ))
  D1_D7 %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf, na.print = "NA")
  
  writeLines('##########################################################################')
  D10_D14 <-  VS_checktime %>% filter(Timepoint_ID   >=10)
  writeLines(sprintf('On visit D10 and D14, %s patients have mismatched time of vital sign measurement:',
                     nrow(D10_D14)
  ))
  D10_D14 %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf, na.print = "NA")
  
  writeLines('##########################################################################')
  missing_temp_time <-  VS_checktime %>% filter(is.na(Time))
  writeLines(sprintf(' %s patients have missing time of temperature measurement:',
                     nrow(missing_temp_time)
  ))
  missing_temp_time %>% arrange(Label, Timepoint_ID)  %>%  print(n = Inf, na.print = "NA")
  
  G <- vita_data %>%
    ggplot() +
    geom_jitter(aes(x = Timepoint_ID, y = Timepoint_ID-Time), size = 3, alpha = 0.5, width = 0.05) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    theme_bw(base_size = 13) +
    scale_x_continuous(breaks = 0:14)
  
return(G)
  
}
# -----------------------------------------------------------------------------------------------
# 4. Check vital sign data

check_vital_signs <- function(vita_data){
  writeLines('##########################################################################')
  writeLines('### Vital sign database: Checking outliers at 6 sd')
  RR_outliers <- abs(vita_data$vs_rr - mean(vita_data$vs_rr, na.rm = T)) > 6*sd(vita_data$vs_rr, na.rm = T) & !is.na(vita_data$vs_rr)
  HR_outliers <- abs(vita_data$vs_hr - mean(vita_data$vs_hr, na.rm = T)) > 6*sd(vita_data$vs_hr, na.rm = T) & !is.na(vita_data$vs_hr)
  SBP_outliers <- abs(vita_data$vs_sbp - mean(vita_data$vs_sbp, na.rm = T)) > 6*sd(vita_data$vs_sbp, na.rm = T) & !is.na(vita_data$vs_sbp)
  DBP_outliers <- abs(vita_data$vs_dbp - mean(vita_data$vs_dbp, na.rm = T)) > 6*sd(vita_data$vs_dbp, na.rm = T) & !is.na(vita_data$vs_dbp)
  VStemp_outliers <- abs(vita_data$vs_temp - mean(vita_data$vs_temp, na.rm = T)) > 6*sd(vita_data$vs_temp, na.rm = T) & !is.na(vita_data$vs_temp)
  
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows contain an outlier in respiratory rate:', 
                     nrow(vita_data[RR_outliers,]))
  )
  vita_data[RR_outliers,] %>% select(Label, visit, rangrp, vs_rr) %>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows contain an outlier in heart rate:', 
                     nrow(vita_data[HR_outliers,]))
  )
  vita_data[HR_outliers,] %>% select(Label, visit, rangrp, vs_hr ) %>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows contain an outlier in systolic blood pressure:', 
                     nrow(vita_data[SBP_outliers,]))
  )
  vita_data[SBP_outliers,] %>% select(Label, visit, rangrp, vs_sbp ) %>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows contain an outlier in diastolic blood pressure:', 
                     nrow(vita_data[DBP_outliers,]))
  )
  vita_data[DBP_outliers,] %>% select(Label, visit, rangrp,  vs_dbp ) %>% print()

  DBP_outrace_SBP <- vita_data$vs_dbp > vita_data$vs_sbp & !is.na(vita_data$vs_sbp) & !is.na(vita_data$vs_dbp)
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows have a higher diastolic than systolic blood pressure:', 
                     nrow(vita_data[DBP_outrace_SBP,]))
  )
  vita_data[DBP_outrace_SBP,] %>% select(Label, visit, rangrp,  vs_dbp, vs_sbp ) %>% print()
  
  
  G1 <- ggplot(vita_data, aes(x = vs_sbp, y = vs_dbp)) +
    geom_point(alpha = 0.5, size = 3) +
    geom_abline(slope=1, intercept = 0, col = "red", linetype = "dashed") +
    xlim(20,200) +
    ylim(20,200) +
    theme_bw(base_size = 13) +
    xlab("Systolic blood pressure") +
    ylab("Diastolic blood pressure")
  
  
  G2 <- ggplot(vita_data, aes(x = vs_hr, y = vs_rr)) +
    geom_point(alpha = 0.5, size = 3) +
    theme_bw(base_size = 13) +
    xlab("Heart rate") +
    ylab("Respiratory rate")
  
return(ggarrange(G1, G2, ncol = 2) %>% suppressWarnings())  

}

