##*********************************************************
##******************* Final status database ***************
##*********************************************************
###########################################################################
# 1. Load final status  data
load_final_status_data <- function(rand_app_data){
  writeLines('##########################################################################')
  writeLines('Reading final status database from MACRO...')
  writeLines('##########################################################################')
  final_status = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimFinalStatus.dta"))
  
  id_data = final_status %>% distinct(Label) %>% pull(Label) %>% as.character()
  writeLines('### Final status database: Checking MACRO data entry progress:')
  writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient final status data on MACRO: %s', 
                     nrow(rand_app_data) + 19, # 19 = number of patients from site 057 and 058
                     length(id_data)
  )
  )

  writeLines('##########################################################################')
  id_missing <- rand_app_data$ID[which(!rand_app_data$ID %in% id_data)]
  writeLines(sprintf('%s patients have no final status data on any visits:', 
                     length(id_missing))
  )
  print(id_missing)
  
  return(final_status)
    
}
#######################################################################################################################
check_final_status <- function(final_status){
  writeLines('##########################################################################')
  writeLines('### Final status database: Checking missing data on final status')
  writeLines('##########################################################################')
  
  missing_fs <-   final_status %>% filter(is.na(fs_compyn))
  writeLines(sprintf('%s patients do not have data on if they have completed the trial:',
                     nrow(missing_fs)
  ))
  missing_fs %>% arrange(Label) %>% select(Label, rangrp, fs_compyn) %>% print(n = Inf, na.print = "NA")
  
  writeLines('##########################################################################')
  
  not_complete <- final_status %>% filter(fs_compyn == 0)
  writeLines(sprintf('%s patients have been flagged of not completing the trial:',
                     nrow(not_complete)
  ))
  table(not_complete$fs_rsn, useNA = 'always')  %>% print()
  haven::print_labels(not_complete$fs_rsn)
  writeLines('##########################################################################')
  no_details <-   not_complete %>% filter(fs_rsn == 5, 
                                          fs_rsnothsp == "" & fs_comm == "") %>% select(Label, rangrp, fs_compyn, fs_rsn, fs_rsnothsp, fs_comm)
  writeLines(sprintf('%s patients who stopped the trial for other reason do not have further details:',
                     nrow(no_details)
  ))
  no_details %>% print(n = Inf, na.print = "NA")
  writeLines('##########################################################################')
  writeLines('### Final status database: Checking if any patients indicated dead:')
  dead_flag <- final_status %>% filter(fs_diecov == 1)
  writeLines(sprintf('%s patients were flagged as dead!!!!', 
                     nrow(dead_flag)))
  dead_flag %>% select(Label, rangrp, fs_diecov) %>% print(n = Inf, na.print = "NA")
  
  writeLines('##########################################################################')
  writeLines('### Final status database: Checking adversed events:')
  writeLines(sprintf('%s patients have missing AE information:', 
                     sum(is.na(final_status$fs_ae))))
  writeLines('##########################################################################')
  AE <- final_status %>% filter(fs_ae == 1)
  writeLines(sprintf('%s patients reported AEs, %s of which were SAEs and %s do not have information:', 
                    nrow(AE),
                    AE %>% filter(fs_sae == 1) %>% nrow(),
                    AE %>% filter(is.na(fs_sae)) %>% nrow()
                    ))
  AE %>% filter(is.na(fs_sae)) %>% select(Label, rangrp, fs_ae, fs_sae) %>% print(n = Inf, na.print = "NA")
  writeLines('##########################################################################')
  mismatched_AE <- final_status %>% filter(fs_sae == 1 & (fs_ae == 0|is.na(fs_ae)))
  writeLines(sprintf('%s patients reported SAEs but not indicated of having an AE', 
                     nrow(mismatched_AE)
  ))
  mismatched_AE %>% select(Label, rangrp, fs_ae, fs_sae) %>% print(n = Inf, na.print = "NA")
  writeLines('##########################################################################')
  writeLines('### Final status database: Checking treatment change:')
  writeLines(sprintf('%s patients do not have information on treatment change', 
                     nrow(final_status %>% filter(is.na(fs_tretchange)))
  ))
  final_status %>% filter(is.na(fs_tretchange)) %>% pull(Label) %>% as.character() %>% print()
  writeLines('##########################################################################')
  trt_changed <-   final_status %>% filter(fs_tretchange  == 1)
  writeLines(sprintf('%s patients do not have information on treatment change, %s of which do not have information on the date of treatment change', 
                     nrow(trt_changed),
                     nrow(trt_changed %>% filter(is.na(fs_tretchandat)))))
  
  trt_changed %>% filter(is.na(fs_tretchandat)) %>% select(Label, rangrp, fs_tretchange, fs_tretchandat) %>% print(n = Inf, na.print = "NA")
  
  final_status%>%filter(fs_conmed == 1)
  writeLines('##########################################################################')
  writeLines('### Final status database: checking unscheduled visits:')
  writeLines(sprintf('%s patients do not have information on unscheduled visits', 
                     nrow(final_status %>% filter(is.na(fs_unsche)))))
  final_status %>% filter(is.na(fs_unsche)) %>% pull(Label) %>% as.character() %>% print()
}



