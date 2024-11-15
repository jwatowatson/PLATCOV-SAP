##*********************************************************
##******************* Temperature database *******************
##*********************************************************
###########################################################################
# 1. Load temperature data
  load_temp_data <- function(rand_app_data, query_file_name){
    # load the clinical data
    data_name <- 'InterimFUTemp.dta'
    file_name <- paste0(prefix_dropbox, "/Data/", data_name)
    version <- file.info(file_name)$ctime %>% as.Date()
    today <- Sys.Date()
    temp_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimFUTemp.dta"))
    
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
    
    id_data = temp_data %>% distinct(Label) %>% pull(Label) %>% as.character()
    
    writeLines('### Temperature database: Checking MACRO data entry progress:')
    writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient temperature data on MACRO: %s', 
                       nrow(rand_app_data) + 19, # 19 = number of patients from site 057 and 058
                       length(id_data)
    )
    )
    writeLines('##########################################################################')
    
    # Missing temperature follow-up data?
    id_missing <- rand_app_data$ID[which(!rand_app_data$ID %in% id_data)]
    
    if(length(id_missing) > 0){
      id_missing_matrix <- matrix(id_missing, ncol = 5)
      
      id_missing_matrix <- data.frame("Dataset" = data_name, #Dataset
                                      "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                      "Question/Variable" = "Follow-up temperature", #'Question/Variable'
                                      "Query message" = paste0(length(id_missing), " patients have no data on MACRO about follow-up temperature (D1 to D7, D10, and D14). IDs are listed here:"), #'Query message'
                                      id_missing_matrix #'Example data'
      ) %>%
        mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
      
      write.table(id_missing_matrix, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
      cat("\n", file=query_file_name, append=TRUE)
      
      writeLines(sprintf('Query: %s\n',
                         id_missing_matrix$Query.message[1])
      )
      print(id_missing_matrix[,-(1:4)])
      writeLines('##########################################################################')
    }
    
    
    # Missing at any visit??
    for(i in paste0("D", c(0:7, 10, 14))){
      
      fu_temp_missing <- temp_data %>% filter(visit == i) %>% mutate(AM = is.na(fut_amdat), PM = is.na(fut_pmdat)) %>% 
        filter(AM & PM) %>%
        select(Label, visit, fut_amdat, fut_amtim, fut_amtemp,
               fut_pmdat, fut_pmtim, fut_pmtemp)
      
      fu_temp_missing <- data.frame("Dataset" = data_name, #Dataset
                                    "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                    "Question/Variable" = "Follow-up temperature", #'Question/Variable'
                                    "Query message" = paste0(nrow(fu_temp_missing), " patients have missing follow-up temperature on ", i), #'Query message'
                                    fu_temp_missing #'Example data'
      ) %>%
        mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
      
      write.table(fu_temp_missing, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
      cat("\n", file=query_file_name, append=TRUE)
      
      writeLines(sprintf('Query: %s\n',
                         fu_temp_missing$Query.message[1])
      )
      print(fu_temp_missing[,-(1:4)])
      writeLines('##########################################################################')
      
    }
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
  
  fever_data <- fever_data %>%
    mutate(time_diff = abs(Timepoint_ID - Time)) %>%
    mutate(expected_date = as.Date(Rand_date_time) + Timepoint_ID) 
  
  FUtemp_checktime <- fever_data %>%
    filter(time_diff > 1) %>%
    select(Label, visit, Timepoint_ID, Time, time_diff, randat, rantim, fut_dat, fut_tim, expected_date) %>%
    rename("temp_time" = "Time") %>%
    mutate(temp_time = round(temp_time, 1),
           time_diff = round(time_diff, 1))
  
  writeLines('##########################################################################')
  
  # Negative time of temperature measurement
  negative_time <-   FUtemp_checktime %>% filter(temp_time < 0)
  
  if(nrow(negative_time) > 0){
    negative_time <- data.frame("Dataset" = data_name, #Dataset
                                "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                "Question/Variable" = "Follow-up temperature", #'Question/Variable'
                                "Query message" = paste0(nrow(negative_time), " patients have negative time of temperature measurement. [randomised time (randat and rantime) after temperature measurement (fut_dat and fut_tim; combined am and pm)]. Please see expected date in the 'expected_date' column"), #'Query message'
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
    temp_time_mismatched <- FUtemp_checktime %>% filter(Timepoint_ID  == i & temp_time > 0)
    if(nrow(temp_time_mismatched) > 0){
      temp_time_mismatched <- data.frame("Dataset" = data_name, #Dataset
                                         "CRF form/Topic" = "Day 0 to Day 14", #'CRF form/Topic'
                                         "Question/Variable" = "Follow-up temperature", #'Question/Variable'
                                         "Query message" = paste0(nrow(temp_time_mismatched), " patients have mismatched time of temperature measurement (fut_dat and fut_tim) with randomisation time (randat and rantime), which is greater than 1 day. Please see expected date in the 'expected_date' column"), #'Query message'
                                         temp_time_mismatched #'Example data'
      ) %>%
        mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
      
      write.table(temp_time_mismatched, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
      cat("\n", file=query_file_name, append=TRUE)
      
      writeLines(sprintf('Query: %s\n',
                         temp_time_mismatched$Query.message[1])
      )
      print(temp_time_mismatched[,-(1:4)])
      writeLines('##########################################################################')
      
    }
  }
  
  # Check the mismatches between randomisation date and temperature date
  G <- fever_data %>%
    ggplot() +
    geom_jitter(aes(x = Timepoint_ID, y = time_diff), size = 3, alpha = 0.5, width = 0.05) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    theme_bw(base_size = 13) +
    scale_x_continuous(breaks = 0:14)
  
  return(G)
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