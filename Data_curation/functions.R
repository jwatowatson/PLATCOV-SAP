library(dplyr)
library(tidyr)
library(seqinr)
library(readxl)
library(stringr)
####################################################################################
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
####################################################################################
check_temp_missing <- function(dat, ampm){
  long_temp <- dat[,c("Label", "visit", "fut_temp")]
  wide_temp <- long_temp %>%
    pivot_wider(names_from = visit, values_from = fut_temp, values_fill = NA)
  colnames(wide_temp)[-1] <- paste0("missing_", ampm, "_temp_", colnames(wide_temp)[-1])
  wide_temp[,-1] <- is.na(wide_temp[,-1])
  wide_temp
}
####################################################################################
# Combine all FASTA files for Nextclade analysis
extract_FASTA <- function(){
  listfile <- list.files(paste0(prefix_dropbox, "/DATA/PLATCOV_FASTA/"), full.names = T,recursive = T)
  listfile_fasta <- listfile[-grep(".xlsx", listfile)]
  
  all_fasta <- do.call(c, lapply(listfile_fasta, read.fasta))
  
  write.fasta(all_fasta, names = names(all_fasta), file = paste0(prefix_analysis_data, "/Analysis_Data/", "all_fasta.fasta"))
  
  names <-  listfile[grep(".xlsx", listfile)]
  # Extract Sequence ID map
  Sample_ID_map <- NULL
  for(i in 1:length(names)){
    if(grepl("Brazil", names[i])){
      temp <- read_xlsx(names[i])  
      cols <- c("Subject ID", "sample")
      temp <- temp[, cols]
      colnames(temp) <- c("Patient_ID", "Sequence_ID")
      
      Sample_ID_map <- rbind(Sample_ID_map, temp)
    } else if (grepl("Laos", names[i])) {
      temp <- read_xlsx(names[i])  
      cols <- c("Patient ID", "GISAID name")
      temp <- temp[, cols]
      colnames(temp) <- c("Patient_ID", "Sequence_ID")
      Sample_ID_map <- rbind(Sample_ID_map, temp)
    }  else if (grepl("Thailand", names[i])) {
      temp <- read_xlsx(names[i])
      if("Patient ID" %in% temp[1,]){colnames(temp) <- temp[1,]; temp <- temp[-1,]}
      cols <- c("Patient_ID", "Sequence_ID")
      temp <- temp[, cols]
      colnames(temp) <- c("Patient_ID", "Sequence_ID")
      temp$Patient_ID <- sapply(str_split(str_trim(temp$Patient_ID, "both"), "-|_"), function(x) paste(x, collapse = "-"))
      
      Sample_ID_map <- rbind(Sample_ID_map, temp)
    }
  }
  
  #Sample_ID_map$Patient_ID<- str_extract(Sample_ID_map$Patient_ID, '.+?\\d{3,4}')
  write.csv(Sample_ID_map, file = paste0(prefix_analysis_data, "/Analysis_Data/", "sequencing_ID_map.csv"), row.names = F)
  Sample_ID_map
}
####################################################################################
##### run the python script to convert lineage names into a usable set

get_nanopore_data = function(prefix_analysis_data, run_python=T, system_used){
  if(run_python){  
    if(system_used == "windows"){
      shell(run_lineage_classifier) 
    } else {
      system("python3 lineage_classifier.py --input ../Analysis_Data/lineages.csv --output ../Analysis_Data/newlineagelist.csv")
    }
  }
  
  res_update = read.csv(paste0(prefix_analysis_data, "/Analysis_Data/newlineagelist.csv"))
  res_update = res_update[!duplicated(res_update$Original), ]
  write.csv(x = res_update, file = paste0(prefix_analysis_data, "/Analysis_Data/newlineagelist.csv"),row.names = F)
  print(unique(res_update$Original))
  
  res <- read.csv("../Analysis_Data/lineages.csv")
  colnames(res) <- c("ID", "Lineages")
  
  res$Variant = plyr::mapvalues(res$Lineage, from = res_update$Original, to = res_update$VariantClass)
  res = res[,c('ID','Variant')]
  write.csv(x = res, file = paste0(prefix_analysis_data, "/Analysis_Data/variant_data.csv"),row.names = F)
  res
  
}
####################################################




