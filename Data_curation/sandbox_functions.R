check_rand_arms <- function(clin_data, IDs_pending, rand_app_data, query_file_name){
  clin_data$rangrp = sjlabelled::as_character(clin_data$rangrp)
  clin_data$rangrp[clin_data$rangrp=='Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir'
  clin_data$rangrp[clin_data$rangrp=='Molnupiravir and Nirmatrelvir/ritonavir']='Nirmatrelvir + Ritonavir + Molnupiravir'
  
  check_data <- merge(rand_app_data, clin_data, by.x = "ID",  by.y = "Label", all.x = T)
  
  ### Check missing data
  writeLines('### Clinical database: Checking if randomisation arms are missing:')
  rangrp_missing <- is.na(check_data$rangrp)
  rangrp_missing <- check_data %>% filter(rangrp_missing) %>% select(ID, rangrp) %>% rename("Label" = "ID")
  
  if(nrow(rangrp_missing) > 0){
    rangrp_missing <- data.frame("Dataset" = data_name, #Dataset
                                       "CRF form/Topic" = "Baseline", #'CRF form/Topic'
                                       "Question/Variable" = "Randomisation arm", #'Question/Variable'
                                       "Query message" = paste0(nrow(rangrp_missing), " patients have missing randomisation arm information"), #'Query message'
                                 rangrp_missing #'Example data'
    ) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(rangrp_missing, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
    cat("\n", file=query_file_name, append=TRUE)
    
    writeLines(sprintf('Query: %s\n',
                       rangrp_missing$Query.message[1])
    )
    print(rangrp_missing[,-(1:4)])
    writeLines('##########################################################################')
    writeLines('### [MANUAL CORRECTIONS]: Using randomisation arms from randomisation database in further analyses')
    for(i in rangrp_missing$Label){clin_data$rangrp[clin_data$Label == i] <- rand_app_data$Treatment[rand_app_data$ID == i]}
    writeLines('##########################################################################')
  }
  
  ### Check mismatched data
  writeLines('### Clinical database: Checking if randomisation arms matched with the randomisation database:')
  rangrp_disagree <- (check_data$Treatment != check_data$rangrp) & !is.na(check_data$rangrp)
  rangrp_disagree <- check_data %>% filter(rangrp_disagree) %>% select(ID, rangrp, Treatment) %>% rename("Trt_SHINY" = "Treatment")
  
  if(nrow(rangrp_disagree) > 0){
    rangrp_disagree <- data.frame("Dataset" = data_name, #Dataset
                                 "CRF form/Topic" = "Baseline", #'CRF form/Topic'
                                 "Question/Variable" = "Randomisation arm", #'Question/Variable'
                                 "Query message" = paste0(nrow(rangrp_disagree), " patients mismatched randomisation arm information with the randomisation database (SHINY). Please check."), #'Query message'
                                 rangrp_disagree #'Example data'
    ) %>%
      mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
    
    write.table(rangrp_disagree, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
    cat("\n", file=query_file_name, append=TRUE)
    
    writeLines(sprintf('Query: %s\n',
                       rangrp_disagree$Query.message[1])
    )
    print(rangrp_disagree[,-(1:4)])
    writeLines('##########################################################################')
    writeLines('### [MANUAL CORRECTIONS]: Using randomisation arms from randomisation database in further analyses')
    for(i in rangrp_disagree$Label){clin_data$rangrp[clin_data$Label == i] <- rand_app_data$Treatment[rand_app_data$ID == i]}
    writeLines('##########################################################################')
  }
  
  return(clin_data)
}
