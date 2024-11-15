check_age <- function(clin_data, IDs_pending, rand_app_data, query_file_name){
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
  
  age_missing <- check_data %>% filter(age_dob_missing) %>% select(Label, dob_my, age_yr)

  if(sum(check_data$age_dob_missing) > 0){
    age_missing <- data.frame("Dataset" = data_name, #Dataset
                              "CRF form/Topic" = "Baseline", #'CRF form/Topic'
                              "Question/Variable" = "Demographics", #'Question/Variable'
                              "Query message" = paste0(nrow(age_missing), " patients have missing age (age_yr) AND date of birth (dob_my) information."), #'Query message'
                              age_missing #'Example data'
    ) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(age_missing, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
    cat("\n", file=query_file_name, append=TRUE)
    
    writeLines(sprintf('Query: %s\n',
                       age_missing$Query.message[1])
    )
    print(age_missing[,-(1:4)])
    writeLines('##########################################################################')
    #manual correction
    writeLines('### [MANUAL CORRECTIONS]: Subsequent analyses use age information from the randomisation database')
    for(i in age_missing$Label){clin_data$age_yr[clin_data$Label == i] <- rand_app_data$age[rand_app_data$ID == i]}
    writeLines('##########################################################################')
    
  }
  
  
  # mismatched age?
  age_mismatched <- check_data %>% filter(!age_agree_yn) %>% select(Label, age_yr, age)
  colnames(age_mismatched) <- c("Label", "age_yr_MACRO", "age_yr_SHINY")
  if(nrow(age_mismatched) > 0){
    age_mismatched <- data.frame("Dataset" = data_name, #Dataset
                              "CRF form/Topic" = "Baseline", #'CRF form/Topic'
                              "Question/Variable" = "Demographics", #'Question/Variable'
                              "Query message" = paste0(nrow(age_mismatched), " patients have mismatched age information with randomisation database (SHINY). Please check with the source documents"), #'Query message'
                              age_mismatched #'Example data'
    ) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(age_mismatched, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
    cat("\n", file=query_file_name, append=TRUE)
    
    writeLines(sprintf('Query: %s\n',
                       age_mismatched$Query.message[1])
    )
    print(age_mismatched[,-(1:4)])
    writeLines('##########################################################################')
    #manual correction
    writeLines('### [MANUAL CORRECTIONS]: Subsequent analyses use age information from the randomisation database')
    for(i in age_mismatched$Label){clin_data$age_yr[clin_data$Label == i] <- rand_app_data$age[rand_app_data$ID == i]}
    writeLines('##########################################################################')
  }
  
  clin_data$ v <- as.numeric(clin_data$age_yr)
  
  # writeLines('### Clinical database: Checking distributions of age:')
  # ggplot(clin_data, aes(x = randat, y =  age_yr)) +
  #   geom_point(size = 3, alpha = 0.25) +
  #   theme_bw(base_size = 13) +
  #   xlab("Randomisation date") +
  #   ylab("Age (years)")
  writeLines('##########################################################################')
  
  return(clin_data)
}
