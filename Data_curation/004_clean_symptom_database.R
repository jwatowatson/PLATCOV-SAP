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
###########################################################################
# 3. Check symptom yes/no
check_symptom_data <- function(symp_data){
  symptom_data = symp_data %>% 
    mutate(ID = Label, Any_symptom = sq_yn) %>%
    rename("sq_soreyn" = "sq_sore")
  
  col_others <- colnames(symptom_data)[grepl('^sq.*des$', colnames(symptom_data))]
  # Checking if other symptom is presented ?
  symptom_data <- symptom_data %>%
    mutate(sq_otheryn = NA) %>%
    mutate(across(all_of(col_others), ~ifelse(. == "", NA, .)))
  for(i in 1:nrow(symptom_data)){
    if(!all(is.na(symptom_data[i,col_others]))){symptom_data$sq_otheryn[i] <- 1} else {symptom_data$sq_otheryn[i] <- 0}
  }
  
  # Check if 1 in sq_yn always provided details
  columns_to_check <- colnames(symptom_data)[grepl('^sq.*yn$', colnames(symptom_data))]
  columns_to_check <- columns_to_check[-1]
  check_sqyn <- symptom_data %>% 
    mutate(symptom_counts = (rowSums(select(., all_of(columns_to_check)), na.rm = TRUE))) %>%
    select(all_of(c("Label", "Timepoint_ID", "sq_yn", columns_to_check, "symptom_counts")))
  
  writeLines('### Symptom database: Checking symptom data mismatched:')
  writeLines(sprintf('%s rows have sq_yn column = 1, but no details explained:', 
                     check_sqyn %>% filter(sq_yn == 1 & symptom_counts == 0) %>% nrow)
  )
  check_sqyn %>% filter(sq_yn == 1 & symptom_counts == 0) %>% select(Label, Timepoint_ID, sq_yn, symptom_counts) %>% print()
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows have sq_yn column = 0, but indicated of having at least 1 symptom:', 
                     check_sqyn %>% filter(sq_yn == 0 & symptom_counts != 0) %>% nrow)
  )
  check_sqyn %>% filter(sq_yn == 0 & symptom_counts != 0) %>% select(Label, Timepoint_ID, sq_yn, symptom_counts) %>% print()
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows have missing data on the sq_yn column:', 
                     check_sqyn %>% filter(is.na(sq_yn)) %>% nrow)
  )
  check_sqyn %>% filter(is.na(sq_yn)) %>% select(Label, Timepoint_ID, sq_yn, symptom_counts) %>% print()
  
  writeLines('##########################################################################')
  writeLines('### [Export] A list of free-text other symptoms...')
  other_symptoms <- symptom_data %>%
    select(matches("^sq.*des$")) %>%
    unlist() %>%
    unname() %>%
    table()
  
  write.csv(other_symptoms, "../Analysis_Data/other_symptoms.csv", row.names = F)
  writeLines('##########################################################################')
  
  return(symptom_data)
}
###########################################################################
# 4. Check symptom grades
check_symptom_grades <- function(symptom_data){
  col_symptom_gr <- colnames(symptom_data)[grepl('^sq.*gr$', colnames(symptom_data))]
  col_symptom_gr <- col_symptom_gr[!grepl('oth', col_symptom_gr)]
  
  col_symptom_yn <- colnames(symptom_data)[grepl('^sq.*yn$', colnames(symptom_data))]
  col_symptom_yn <-  col_symptom_yn[-c(1, length(col_symptom_yn))]
  writeLines('##########################################################################')
  for(i in 1:length(col_symptom_gr)){
    yn <- symptom_data[,col_symptom_yn[i]] %>% as.numeric()
    gr <- symptom_data[,col_symptom_gr[i]]
    
    writeLines(sprintf('Checking symptom %s: %s rows reported this symptom without grades:', 
                       col_symptom_yn[i],
                       symptom_data[(yn == 1 & !is.na(yn)) & is.na(gr),] %>% nrow()
                       )
    )
    symptom_data[(yn == 1 & !is.na(yn)) & is.na(gr),] %>% select(all_of(c("Label", "Timepoint_ID", col_symptom_yn[i], col_symptom_gr[i]))) %>% print()
    writeLines('#--------------------------------------------------------------------------')
    writeLines(sprintf('Checking symptom %s: %s rows reported this symptom grades without checking y/n:', 
                       col_symptom_yn[i],
                       symptom_data[(yn == 0 | is.na(yn)) & !is.na(gr),] %>% nrow()
    )
    )
    symptom_data[(yn == 0 | is.na(yn)) & !is.na(gr),] %>% select(all_of(c("Label", "Timepoint_ID", col_symptom_yn[i], col_symptom_gr[i]))) %>% print()
    writeLines('#--------------------------------------------------------------------------')
    writeLines(sprintf('Checking symptom %s: %s rows reported this symptom with grade not in the range:', 
                       col_symptom_yn[i],
                       symptom_data[!gr %in% 0:3 & !is.na(gr),] %>% nrow()
    )
    )
    symptom_data[!gr %in% 0:3 & !is.na(gr),] %>% select(all_of(c("Label", "Timepoint_ID", col_symptom_yn[i], col_symptom_gr[i]))) %>% print()
    writeLines('##########################################################################')
    
  }
  
  write.table(x = symptom_data[, c('ID','Timepoint_ID','Any_symptom','heart_rate', col_symptom_gr)], 
              file = paste0(prefix_analysis_data, "/Analysis_Data/symptoms_interim.csv"), 
              row.names = F, sep=',', quote = F)
}
###########################################################################
# 5. Check ohter symptom grades
check_other_symptom_grades <- function(symptom_data){
  col_symptom_gr <- colnames(symptom_data)[grepl('^sq.*gr$', colnames(symptom_data))]
  col_symptom_gr <- col_symptom_gr[grepl('oth', col_symptom_gr)]
  
  col_symptom_des <- colnames(symptom_data)[grepl('^sq.*des$', colnames(symptom_data))]
  writeLines('##########################################################################')
  
  for(i in 1:length(col_symptom_gr)){
    des <- symptom_data[,col_symptom_des[i]]
    gr <- symptom_data[,col_symptom_gr[i]]
    
    writeLines(sprintf('Checking other symptom %s: %s rows reported a symptom without grades:', 
                       col_symptom_des[i],
                       symptom_data[(des != "" & !is.na(des)) & is.na(gr),] %>% nrow()
    )
    )
    symptom_data[(des != "" & !is.na(des)) & is.na(gr),] %>% select(all_of(c("Label", "Timepoint_ID", col_symptom_des[i], col_symptom_gr[i]))) %>% print()
    writeLines('#--------------------------------------------------------------------------')
    
    writeLines(sprintf('Checking other symptom %s: %s rows reported grades with no symptom details:', 
                       col_symptom_des[i],
                       symptom_data[(des %in% c("", ".") | is.na(des)) & !is.na(gr),] %>% nrow()
    )
    )
    symptom_data[(des %in% c("", ".") | is.na(des)) & !is.na(gr),] %>% select(all_of(c("Label", "Timepoint_ID", col_symptom_des[i], col_symptom_gr[i]))) %>% print()
    writeLines('##########################################################################')
    
  }

}
  