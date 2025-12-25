##*********************************************************
##******************** Vaccine database *******************
##*********************************************************
###########################################################################
# 1. Load Vaccine Data

## ------------------------------------------------------------------------------
## Function: load_vaccine_data(rand_app_data, query_file_name)
##
## Description:
## Loads the MACRO vaccination dataset and checks completeness against the
## SHINY randomisation database. Identifies patients with missing vaccination entries.
##
## Key Steps:
## 1. Reads `.dta` file from predefined Dropbox path
## 2. Extracts unique patient IDs (`Label`) from vaccine data
## 3. Compares with SHINY randomised IDs to identify pending data
## 4. Outputs progress summary and missing IDs to console
## 5. Converts `vc_statyn` to character for consistency
##
## Parameters:
## - rand_app_data: SHINY randomisation dataset (must contain `ID`)
## - query_file_name: File path for general temperature query CSV
##
## Output:
## - Returns `vacc_data`: the MACRO vaccine dataset with harmonised formats
##
## Notes:
## - File path constructed using global variable `prefix_dropbox`
## - Assumes `Label` uniquely identifies patients in vaccine data
## ------------------------------------------------------------------------------

load_vaccine_data <- function(rand_app_data, query_file_name) {
  writeLines('##########################################################################')
  writeLines('Reading vaccine database from MACRO...')
  writeLines('##########################################################################')
  
  # Define dataset name and file path
  data_name <- 'InterimVaccine.dta'
  file_name <- paste0(prefix_dropbox, "/Data/PLATCOV_22Sep2025/", data_name)
  version <- file.info(file_name)$ctime %>% as.Date()
  today <- Sys.Date()
  
  # Load vaccine data using haven package
  vacc_data <- haven::read_dta(file_name)
  
  # Start new query file with header and blank table
  sink(query_file_name, split = TRUE)
  writeLines(sprintf(
    "PLATCOV data queries\nData: %s\nReceived date: %s\nQuery date: %s\n",
    data_name, version, today))
  sink()
  
  writeLines('##########################################################################')
  
  write.table(data.frame('Dataset' = "", 
                         'CRF form/Topic' = "", 
                         'Question/Variable' = "",
                         'Query message' = "", 
                         'Example data' = ""),
              query_file_name, col.names = TRUE, sep = ",", append = TRUE, row.names = FALSE) %>% suppressWarnings()
  
  id_data <- vacc_data %>% distinct(Label) %>% pull(Label) %>% as.character()
  writeLines('### Vaccination database: Checking MACRO data entry progress:')
  writeLines(sprintf('Number of randomised patients: %s\nNumber of randomised patients with vaccine data on MACRO: %s', 
                     nrow(rand_app_data),
                     length(id_data)))
  
  writeLines('##########################################################################')
  
  id_missing <- rand_app_data$ID[which(!rand_app_data$ID %in% id_data)]
  writeLines(sprintf('%s patients\' vaccination information is still pending:', 
                     length(id_missing)))
  if (length(id_missing) > 0) {
    print(id_missing)
    
    missing_ids_df <- data.frame("Dataset" = data_name,
                                 "CRF form/Topic" = "Vaccination",
                                 "Question/Variable" = "SARS-CoV-2 Vaccine History",
                                 "Query message" = sprintf("%s patients are missing from MACRO vaccine database", length(id_missing)),
                                 "Label" = id_missing) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(missing_ids_df, query_file_name, col.names = TRUE, sep = ",", append = TRUE, row.names = FALSE) %>% suppressWarnings()
    cat("\n", file = query_file_name, append = TRUE)
  }
  
  vacc_data$vc_statyn <- sjlabelled::as_character(vacc_data$vc_statyn)
  
  return(vacc_data)
}



# -----------------------------------------------------------------------------------------------
# 2. Check Vaccine Data

## ------------------------------------------------------------------------------
## Function: check_vaccine_data(clin_data, vacc_data, rand_app_data, query_file_name)
##
## Description:
## Merges clinical and vaccination datasets and verifies completeness and consistency 
## of vaccination status across sources (MACRO vs. SHINY).
##
## Key Steps:
## 1. Merges `clin_data` with `vacc_data` on patient ID, site, and trial
## 2. Merges `clin_data` with `rand_app_data` to obtain SHINY vaccination status
## 3. Checks for missing vaccination status in MACRO
## 4. Checks for mismatch between MACRO and SHINY vaccination fields
## 5. Logs inconsistencies and applies manual corrections
##
## Parameters:
## - clin_data: Clinical dataset with patient-level data (must include 'Label', 'Site', 'Trial')
## - vacc_data: Vaccine dataset from MACRO (must include 'Label', 'vc_statyn')
## - rand_app_data: Randomisation dataset (must include 'ID', 'vaccine')
## - query_file_name: File path for general temperature query CSV
##
## Output:
## - Updated `clin_data` with harmonised vaccination status
##
## Notes:
## - Manual correction: missing or mismatched `vc_statyn` values are overwritten with SHINY `vaccine`
## ------------------------------------------------------------------------------

check_vaccine_data <- function(clin_data, vacc_data, rand_app_data, query_file_name) {
  clin_data <- merge(clin_data, vacc_data,
                     by.x = c('Label', 'Site', 'Trial'),
                     by.y = c('Label', 'Site', 'Trial'),
                     all.x = TRUE)
  
  # clin_data <- merge(clin_data, rand_app_data[, c("ID", "vaccine")],
  #                    by.x = 'Label',
  #                    by.y = 'ID',
  #                    all.x = TRUE)
  
  writeLines('##########################################################################')
  
  writeLines('### Vaccination database: Checking if vaccination status are missing:')
  ind <- is.na(clin_data$vc_statyn)
  writeLines(sprintf('%s patients have no information on vaccination status on MACRO:', 
                     sum(ind)))
  if (sum(ind) > 0) {
    clin_data[ind, ] %>% pull(Label) %>% as.character() %>% print()
    
    no_vaccine_df <- data.frame("Dataset" = "ASTInterimVaccine.dta",
                                "CRF form/Topic" = "Vaccination",
                                "Question/Variable" = "Influenza Vaccine History",
                                "Query message" = sprintf("%s patients are missing vaccination status in MACRO", sum(ind)),
                                clin_data[ind, c("Label", "vc_statyn")]) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(no_vaccine_df, query_file_name, col.names = TRUE, sep = ",", append = TRUE, row.names = FALSE) %>% suppressWarnings()
    cat("\n", file = query_file_name, append = TRUE)
  }
  
  writeLines('##########################################################################')
  
  writeLines('### Vaccination database: Checking mismatched vaccination status:')
  ind_mismatched <- (clin_data$vc_statyn != clin_data$vaccine) &
    !is.na(clin_data$vc_statyn) &
    !is.na(clin_data$vaccine)
  
  if (sum(ind_mismatched) > 0) {
    writeLines(sprintf('%s patients have mismatched vaccination status:', 
                       sum(ind_mismatched)))
    clin_data[ind_mismatched, ] %>%
      select(Label, rangrp, vaccine, vc_statyn) %>%
      print()
    
    mismatch_df <- data.frame("Dataset" = "ASTInterimVaccine.dta",
                              "CRF form/Topic" = "Vaccination",
                              "Question/Variable" = "Influenza Vaccine History",
                              "Query message" = sprintf("%s patients are mismatched in vaccination status between SHINY and MACRO", sum(ind_mismatched)),
                              clin_data[ind_mismatched, c("Label", "rangrp", "vaccine", "vc_statyn")]) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(mismatch_df, query_file_name, col.names = TRUE, sep = ",", append = TRUE, row.names = FALSE) %>% suppressWarnings()
    cat("\n", file = query_file_name, append = TRUE)
    
    writeLines('##########################################################################')
    
    writeLines('### [MANUAL CORRECTIONS]: Using vaccination status from randomisation database in further analyses')
    ID_errors <- clin_data[ind | ind_mismatched, ] %>%
      pull(Label) %>% as.character()
    
    for (i in ID_errors) {
      clin_data$vc_statyn[clin_data$Label == i] <- clin_data$vaccine[clin_data$Label == i]
    }
    
    writeLines('##########################################################################')
  }
  
  return(clin_data)
}