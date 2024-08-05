library(dplyr)
library(tidyr)
temp_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimFUTemp.dta"))

prep_tempdata <- function(temp_data, clin_data){
  x1 = temp_data[, c('Site', 'Label', 'visit', 'fut_amdat', "fut_amtim", "fut_amtemp")]
  x2 = temp_data[, c('Site', 'Label', 'visit', 'fut_pmdat', "fut_pmtim", "fut_pmtemp")]
  
  colnames(x1) = colnames(x2) = 
    c('Site', 'Label', 'visit', 'fut_dat', "fut_tim", "fut_temp")
  temp_data = rbind(x1, x2)
  
  temp_data = temp_data[(!is.na(temp_data$fut_dat) & 
                           !is.na(temp_data$fut_tim) &
                           !is.na(temp_data$fut_temp)), ]
  
  temp_data$temp_time = 
    as.character(as.POSIXct(apply(temp_data[, c('fut_dat','fut_tim')],
                                  1, paste, collapse =' '),
                            format = '%F %H:%M:%S'))
  temp_data = merge(temp_data, 
                    clin_data, 
                    by = c('Label','Site'),
                    all.x=T)
  
  temp_data = temp_data %>% filter(!is.na(scrpassed), 
                                   scrpassed==1, 
                                   !is.na(randat),
                                   rantim != '') %>% 
    ungroup() %>%
    mutate(rangrp = sjlabelled::as_character(rangrp),
           Rand_date_time = as.character(as.POSIXct(paste(randat,
                                                          rantim, sep=' '))),
           Timepoint_ID = gsub(x = visit,replacement = '',pattern='D'),
           Timepoint_ID = gsub(x = Timepoint_ID, replacement = '', pattern = 'H1'),
           Timepoint_ID = as.numeric(Timepoint_ID)) %>%  filter(!is.na(temp_time)) %>%
    arrange(Label, temp_time, Rand_date_time)
  
  temp_data$Label[temp_data$Rand_date_time > temp_data$temp_time]
  for(i in 1:nrow(temp_data)){
    temp_data$Time[i]=as.numeric(difftime(temp_data$temp_time[i],
                                          temp_data$Rand_date_time[i], 
                                          units = 'days'))
  }
  

  
  temp_data = temp_data %>% group_by(Label) %>%
    mutate(Fever_Baseline = ifelse(any(fut_temp>37 & Time<1), 1, 0))
  
  write.csv(x = temp_data, file = '../Analysis_Data/temperature_data.csv', row.names = F, quote = F)
  return(temp_data)
}

check_temp_missing <- function(dat, ampm){
  long_temp <- dat[,c("Label", "visit", "fut_temp")]
  wide_temp <- long_temp %>%
    pivot_wider(names_from = visit, values_from = fut_temp, values_fill = NA)
  colnames(wide_temp)[-1] <- paste0("missing_", ampm, "_temp_", colnames(wide_temp)[-1])
  wide_temp[,-1] <- is.na(wide_temp[,-1])
  wide_temp
}
