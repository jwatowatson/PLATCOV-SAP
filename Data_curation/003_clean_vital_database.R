##*********************************************************
##******************* Vital database *********************
##*********************************************************
###########################################################################

load_vital_data <- function(rand_app_data){
  writeLines('##########################################################################')
  writeLines('##########################################################################')
  writeLines('Reading vital sign database from MACRO...')
  writeLines('##########################################################################')
  vita_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimVitalSigns.dta"))
  
  
  
  
}

glimpse(vita_data)

vita_data %>%
  group_by(visit) %>%
  filter(!is.na(vs_temp)) %>%
  summarise(min = min(vs_temp),
            Q1 = quantile(vs_temp, 0.25),
            Q2 = quantile(vs_temp, 0.5),
            Q3 = quantile(vs_temp, 0.75),
            max = max(vs_temp))

hist(vita_data$vs_temp)
hist(vita_data$vs_hr)
hist(vita_data$vs_rr)
hist(vita_data$vs_sbp)
hist(vita_data$vs_dbp)

plot(vita_data$vs_hr, vita_data$vs_rr)
plot(vita_data$vs_sbp, vita_data$vs_dbp)

glimpse(vita_data)



##########  --- Symptom data --- ########## 
vita_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimVitalSigns.dta"))
check_data$vita_missing <- !check_data$ID %in% vita_data$Label