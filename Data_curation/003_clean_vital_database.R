##*********************************************************
##******************* Vital database *********************
##*********************************************************
###########################################################################
# -----------------------------------------------------------------------------------------------
# 1. Load vital sign data
load_vital_data <- function(rand_app_data, query_file_name){
  data_name <- 'InterimVitalSigns.dta'
  file_name <- paste0(prefix_dropbox, "/Data/", data_name)
  version <- file.info(file_name)$mtime %>% as.Date()
  today <- Sys.Date() 
  
  vita_data = haven::read_dta(file_name)
  
  sink(query_file_name, split = T)
  writeLines(sprintf('PLATCOV data queries\nData: %s\nReceived date: %s\nQuery date: %s\n',
                     data_name,
                     version,
                     today)
  )
  sink()
  writeLines('##########################################################################')
  write.table(data.frame('Dataset' = "",
                         'CRF form/Topic' = "",
                         'Question/Variable' = "",
                         'Query message' = "",
                         'Example data' = ""), 
              query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
  
  id_data = vita_data %>% distinct(Label) %>% pull(Label) %>% as.character()
  
  writeLines('### Vital sign database: Checking MACRO data entry progress:')
  writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient vital sign data on MACRO: %s', 
                     nrow(rand_app_data) + 19, # 19 = number of patients from site 057 and 058
                     length(id_data)
  )
  )
  writeLines('##########################################################################')
  
  visits_ID <- unique(vita_data$visit)
  visits_ID <- visits_ID[visits_ID != "D0H1"]
  
  cat("\n", file=query_file_name, append=TRUE)
  for(i in 1:length(visits_ID)){
    missing_info <- vita_data %>% filter(visit == visits_ID[i], 
                                         is.na(vd_dat) & vs_tim == "" & is.na(vs_temp)) %>% select(Label, visit)
    if(nrow(missing_info) > 0){
      report <- missing_info
      report <- data.frame("Dataset" = data_name, #Dataset
                           "CRF form/Topic" = visits_ID[i], #'CRF form/Topic'
                           "Question/Variable" = "Vital signs", #'Question/Variable'
                           "Query message" = paste0(nrow(report), " patients have no data on MACRO about vital signs on ", visits_ID[i]), #'Query message'
                           report #'Example data'
      ) %>%
        mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
      
      write.table(report, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
      cat("\n", file=query_file_name, append=TRUE)
      
      writeLines(sprintf('Query: %s\n',
                         report$Query.message[1])
      )
      print(report[,-(1:4)])
      writeLines('##########################################################################')
    }
  }
  
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
  writeLines('### Vital sign database: Checking vital sign data and time')
  
  vita_data <- vita_data %>%
    mutate(time_diff = abs(Timepoint_ID - Time)) %>%
    mutate(expected_date = as.Date(Rand_date_time) + Timepoint_ID) 
  
  VS_checktime <- vita_data %>%
    filter(time_diff > 1) %>%
    select(Label, visit, Timepoint_ID, Time, time_diff, randat, rantim, vd_dat, vs_tim, expected_date) %>%
    rename("vs_time" = "Time") %>%
    mutate(vs_time = round(vs_time, 1),
           time_diff = round(time_diff, 1))
  
  # Negative time of temperature measurement
  negative_time <-   VS_checktime %>% filter(vs_time < 0)
  
  if(nrow(negative_time) > 0){
    negative_time <- data.frame("Dataset" = data_name, #Dataset
                                "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                "Question/Variable" = "Vital signs", #'Question/Variable'
                                "Query message" = paste0(nrow(negative_time), " patients have negative time of vital sign measurement. [randomised time (randat and rantime) after vital sign measurement (vd_dat and vs_tim)]. Please see expected date in the 'expected_date' column"), #'Query message'
                                negative_time #'Example data'
    ) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(negative_time, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
    cat("\n", file=query_file_name, append=TRUE)
    
    writeLines(sprintf('Query: %s\n',
                       negative_time$Query.message[1])
    )
    print(negative_time[,-(1:4)])
    writeLines('##########################################################################')
    
  }
  
  # Check time mismatched
  for(i in c(0:7, 10, 14)){
    vs_time_mismatched <- VS_checktime %>% filter(Timepoint_ID  == i & vs_time > 0)
    if(nrow(vs_time_mismatched) > 0){
      vs_time_mismatched <- data.frame("Dataset" = data_name, #Dataset
                                         "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                         "Question/Variable" = "Vital sign", #'Question/Variable'
                                         "Query message" = paste0(nrow(vs_time_mismatched), " patients have mismatched time of vital sign measurement (vd_dat and vs_tim) with randomisation time (randat and rantime), which is greater than 1 day. Please see expected date in the 'expected_date' column"), #'Query message'
                                       vs_time_mismatched #'Example data'
      ) %>%
        mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
      
      write.table(vs_time_mismatched, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
      cat("\n", file=query_file_name, append=TRUE)
      
      writeLines(sprintf('Query: %s\n',
                         vs_time_mismatched$Query.message[1])
      )
      print(vs_time_mismatched[,-(1:4)])
      writeLines('##########################################################################')
      
    }
  }
  
  G <- vita_data %>%
    ggplot() +
    geom_jitter(aes(x = Timepoint_ID, y = time_diff), size = 3, alpha = 0.5, width = 0.05) +
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
  
  outliers <- list("Temperature" = VStemp_outliers, 
                   "Respiratory rate" = RR_outliers, 
                   "Heart rate" = HR_outliers, 
                    "Systolic blood pressure" = SBP_outliers, 
                    "Diastolic blood pressure" = DBP_outliers)
  means = list("Temperature" =  mean(vita_data$vs_temp, na.rm = T), 
              "Respiratory rate" = mean(vita_data$vs_rr, na.rm = T), 
              "Heart rate" =  mean(vita_data$vs_hr, na.rm = T), 
              "Systolic blood pressure" = mean(vita_data$vs_sbp, na.rm = T), 
              "Diastolic blood pressure" = mean(vita_data$vs_dbp, na.rm = T))
  colnamesss = c("vs_temp", "vs_rr", "vs_hr", "vs_sbp", "vs_dbp", "vs_temp")
  
  
  for(i in 1:length(outliers)){
    outlier_i <- outliers[[i]]
    outlier_i <- vita_data[outlier_i,]
    mean_i <- round(means[[i]], 1)
    
    
    
    if(nrow(outlier_i) > 0){
         outlier_i <- outlier_i[, c("Label", "visit", colnamesss[i])]
      
         outlier_i <- data.frame("Dataset" = data_name, #Dataset
                                       "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                       "Question/Variable" = "Vital sign", #'Question/Variable'
                                       "Query message" = paste0(names(outliers)[i]," (", colnamesss[i], ") of ", nrow(outlier_i), " patients are outliers (greater than 6 sd). The mean value = ", mean_i, ". Please check the source documents"), #'Query message'
                                 outlier_i #'Example data'
      ) %>%
        mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
      
      write.table(outlier_i, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
      cat("\n", file=query_file_name, append=TRUE)
      
      writeLines(sprintf('Query: %s\n',
                         outlier_i$Query.message[1])
      )
      print(outlier_i[,-(1:4)])
      writeLines('##########################################################################')
    }
  }
  
  
  DBP_outrace_SBP <- vita_data$vs_dbp > vita_data$vs_sbp & !is.na(vita_data$vs_sbp) & !is.na(vita_data$vs_dbp)
  DBP_outrace_SBP <-   vita_data[DBP_outrace_SBP, c("Label", "visit", "vs_sbp", "vs_dbp")]
  
  if(nrow(DBP_outrace_SBP) > 0){
    DBP_outrace_SBP <- data.frame("Dataset" = data_name, #Dataset
                                "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                "Question/Variable" = "Vital signs", #'Question/Variable'
                                "Query message" = paste0(nrow(DBP_outrace_SBP), " patients have higher diastolic blood pressure (vs_dbp) than systolic blood pressure (vs_sbp). Please check the source documents."), #'Query message'
                                DBP_outrace_SBP #'Example data'
    ) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(DBP_outrace_SBP, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
    cat("\n", file=query_file_name, append=TRUE)
    
    writeLines(sprintf('Query: %s\n',
                       DBP_outrace_SBP$Query.message[1])
    )
    print(DBP_outrace_SBP[,-(1:4)])
    writeLines('##########################################################################')
    
  }
  
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

