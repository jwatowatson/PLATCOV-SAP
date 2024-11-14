load_clinical_data <- function(){
  # load the clinical data
  writeLines('##########################################################################')
  data_name <- 'InterimEnrolment.dta'
  file_name <- paste0(prefix_dropbox, "/Data/", data_name)
  version <- file.info(file_name)$ctime %>% as.Date()
  today <- Sys.Date()
  clin_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimEnrolment.dta"))
  clin_data$scrpassed[clin_data$Label=='PLT-TH1-557']=1
  
  
  sink('Queries/test.csv')
  writeLines(sprintf('PLATCOV data queries\nData: %s\nReceived date: %s\nQuery date: %s\n',
                     data_name,
                     version,
                     today)
  )
  sink()
  
  write.table(B, 'Queries/test.csv', col.names=T, sep=",", append=TRUE, row.names = F)
  
  
  
  writeLines('##########################################################################')
  
  # writeLines(sprintf('Clinical dataset contains %s rows with %s unique Screening IDs and %s unique Patient IDs',
  #                    nrow(clin_data),
  #                    length(unique(clin_data$scrid)),
  #                    length(unique(clin_data$Label))
  # ))
  
  # Check if a patient has multiple screening ID
  duplicated_clin_data <- names(which(table(clin_data$Label) > 1))
  duplicated_clin_data <- duplicated_clin_data[duplicated_clin_data != ""]
  duplicated_clin_data
  duplicated_clin_queries <- clin_data %>%
    filter(Label %in% duplicated_clin_data) %>%
    select(Site, scrid, Label)
  
  A <- data.frame(Dataset = data_name,
             CRF = "Screening",
             Variable = "Screening Number",
             Query = "Different Screening Number ('scrid') were assigned to the same Subject Number ('Label')",
             duplicated_clin_queries
  )
  
  B <- A %>%
    mutate(i = 1:nrow(A),
           Dataset = if_else(i > 1, "", Dataset),
           CRF = if_else(i > 1, "", CRF),
           Variable = if_else(i > 1, "", Variable),
           Query = if_else(i > 1, "", Query)) %>%
    select(-i)
           
  
  C <- list(B,B)
  lapply(C, function(x) write.table( data.frame(x), 'Queries/test.csv'  , append= T, sep=',', row.names = F ))
  
  duplicated_clin_queries

  
  
  
  
  
  
  
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
