##*********************************************************
##******************* Symptom database *********************
##*********************************************************
###########################################################################
# 1. Load symptom  data
load_symptom_data <- function(rand_app_data){
  writeLines('##########################################################################')
  writeLines('Reading symptom database from MACRO...')
  writeLines('##########################################################################')
  symp=haven::read_dta(paste0(prefix_dropbox, "/Data/InterimSymptoms.dta"))
  
  id_data = symp %>% distinct(Label) %>% pull(Label) %>% as.character()
  writeLines('### Symptom database: Checking MACRO data entry progress:')
  writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient symptom data on MACRO: %s', 
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
  writeLines('### Symptom database: Missing data by visits:')
  symp %>%
    group_by(visit) %>%
    summarise(missing = sum(is.na(vd_dat))) %>% print()
  writeLines('##########################################################################')
  
  
  return(symp)
  
}
###########################################################################
# 2. Load symptom  data

prep_symptom_data <- function(symp, vita_data){
  
  symp = symp[!is.na(symp$sq_yn),]
  sort(table(symp$Label),decreasing = T)
  ids_symp_data = names(which(table(symp$Label)>4)) ## need at least 4 records to be included
  
  symp = symp %>% filter(Label %in% ids_symp_data, visit != 'D0H1') %>%
    mutate(Timepoint_ID = gsub(x=visit, pattern='D',replacement=''))
  
  HR_data = vita_data[, c('Label','visit','vs_hr')]
  HR_data$Timepoint_ID = gsub(pattern = 'D',replacement = '',x = HR_data$visit)
  HR_data$Timepoint_ID = as.numeric(gsub(pattern = 'H1',replacement = '',x = HR_data$Timepoint_ID))
  HR_data = HR_data[!is.na(HR_data$vs_hr), ]
  HR_data$heart_rate = HR_data$vs_hr
  
  HR_data = aggregate(heart_rate ~ Label + Timepoint_ID, data = HR_data, mean)
  
  symp_data = merge(symp, HR_data[, c('Label', 'Timepoint_ID','heart_rate')],
                    all=T, by = c('Label', "Timepoint_ID"))
  symp_data$Timepoint_ID = as.numeric(symp_data$Timepoint_ID)
  symp_data = symp_data %>% arrange(Label, Timepoint_ID)
  
  write.csv(x = symp_data, file = '../Analysis_Data/symptom_data.csv', row.names = F, quote = F)
  
  return(symp_data)
  
}




