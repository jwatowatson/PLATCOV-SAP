# 1. Load vital sign data
load_vital_data <- function(rand_app_data){
  # load the clinical data
  data_name <- 'InterimVitalSigns.dta'
  file_name <- paste0(prefix_dropbox, "/Data/", data_name)
  version <- file.info(file_name)$ctime %>% as.Date()
  today <- Sys.Date()
  vita_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimVitalSigns.dta"))

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
 
   # Missing vital sign data?
  id_missing <- rand_app_data$ID[which(!rand_app_data$ID %in% id_data)]
  
  if(length(id_missing) > 0){
    n_rows <- ceiling(length(id_missing) / 5)
    id_missing <- c(id_missing, rep("", n_rows * 5 - length(id_missing)))
    id_missing_matrix <- matrix(id_missing, ncol = 5)
    
    id_missing_matrix <- data.frame("Dataset" = data_name, #Dataset
                                    "CRF form/Topic" = "Baseline and follow-up", #'CRF form/Topic'
                                    "Question/Variable" = "Vital signs", #'Question/Variable'
                                    "Query message" = paste0(length(id_missing), " patients have no data on MACRO about vital signs. IDs are listed here:"), #'Query message'
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
  
  writeLines('##########################################################################')
  writeLines('### Vital sign database: Missing data by visits:')
  
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
  
  
  
  vita_data %>%
    group_by(visit) %>%
    summarise(missing = sum(is.na(vs_temp))) %>% print()
  writeLines('##########################################################################')
  
  return(vita_data)
}
# ----------------------------------------------------------
# ---------------------------------------------