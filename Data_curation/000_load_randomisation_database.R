##*********************************************************
##******************** Randomisation database *******************
##*********************************************************
load_randomisation_data <- function(){
  rand_app_data = rbind(read.csv(paste0(prefix_drop_rand, "/data-TH1.csv")),
                        read.csv(paste0(prefix_drop_rand, "/data-BR3.csv")),
                        read.csv(paste0(prefix_drop_rand, "/data-LA08.csv")),
                        read.csv(paste0(prefix_drop_rand, "/data-PK01.csv"))) %>%
    filter(!is.na(Treatment))
  # Defining patient ID
  rand_app_data$ID = paste0('PLT-', gsub(x = rand_app_data$site,pattern = '0',replacement = ''),
                            '-',
                            stringr::str_pad(rand_app_data$randomizationID, 3, pad = "0"))
  # Converting time zones
  rand_app_data$Site = plyr::mapvalues(x = rand_app_data$site, from=c('TH1','BR3','LA08','PK01'),to=c('th001','br003',"la008","pk001"))
  rand_app_data$Rand_Time = as.POSIXct(rand_app_data$Date, format = '%a %b %d %H:%M:%S %Y',tz = 'GMT')
  rand_app_data$tzone = plyr::mapvalues(x = rand_app_data$site,
                                        from = c('TH1','LA08','BR3','PK01'),
                                        to = c('Asia/Bangkok','Asia/Bangkok','America/Sao_Paulo','Asia/Karachi'))
  rand_app_data$Rand_Time_TZ=NA
  for(i in 1:nrow(rand_app_data)){
    rand_app_data$Rand_Time_TZ[i] = as.character(with_tz(rand_app_data$Rand_Time[i], tzone = rand_app_data$tzone[i]))
  } 
  
  return(rand_app_data)
}

