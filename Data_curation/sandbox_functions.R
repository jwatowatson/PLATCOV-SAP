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
# ---------------------------------------------