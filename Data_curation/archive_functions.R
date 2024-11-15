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