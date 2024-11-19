##############################################################################################
# PLATCOV project
# This script prepares data for further analyses
##############################################################################################
library(dplyr)
library(readxl)
library(readr)
library(ggplot2)
library(lme4)
library(lubridate)
library(anytime)
library(ggpubr)
#library(plyr)
##Define user folder path####################################################################
source('user_settings.R')
source('functions.R')

source("000_load_randomisation_database.R")
source("001_clean_clinical_database.R")
source("002_clean_temperature_database.R")
source("003_clean_vital_database.R")
source("004_clean_symptom_database.R")
source("005_clean_final_status_database.R")
source("006_clean_CBC_database.R")


options(max.print = 5000)
#############################################################################################
# 000 Randomization database
rand_app_data <- load_randomisation_data()
#############################################################################################
## 001 Clinical database
query_file_name <-  'Queries/01_PLATCOV_queries_clinical_database.csv'

sink("Queries/01_queries_clinical_database.txt", split = T)
clin_data <- load_clinical_data(query_file_name)
IDs_pending <- check_MACRO_clinical_database(clin_data, rand_app_data)
clin_data <- check_screen_failure(clin_data, query_file_name)
clin_data <- check_randomisation_info(clin_data, query_file_name)
clin_data <- check_sex(clin_data, IDs_pending, rand_app_data, query_file_name)
clin_data <- check_age(clin_data, IDs_pending, rand_app_data, query_file_name)
clin_data <- check_symptom_onset(clin_data, IDs_pending, rand_app_data, query_file_name)
clin_data <- check_weight_height(clin_data, IDs_pending, query_file_name)
clin_data <- check_rand_date_time(clin_data, IDs_pending, rand_app_data, query_file_name)
clin_data <- check_rand_arms(clin_data, IDs_pending, rand_app_data, query_file_name)
sink()

#############################################################################################

#############################################################################################
## 002 Vital sign database

query_file_name <-  'Queries/02_PLATCOV_queries_vital_sign_database.csv'







sink("Queries/03_queries_vital_sign_database.txt", split = T)
vita_data <- load_vital_data(rand_app_data)
vita_data <- prep_vitadata(vita_data, clin_data)
check_vita_time(vita_data)
check_vital_signs(vita_data)
sink()





## 002 Temperature database
query_file_name <-  'Queries/02_PLATCOV_queries_temperature_database.csv'


sink("Queries/02_queries_temperature_database.txt", split = T)
temp_data <- load_temp_data(rand_app_data, query_file_name)
fever_data <- prep_tempdata(temp_data, clin_data)
check_time_temp(fever_data)
sink()




clin_data <- add_baseline_fever(clin_data, fever_data)


#############################################################################################
## 004 Symptom database
sink("Queries/04_queries_symptom_database.txt", split = T)
symp <- load_symptom_data(rand_app_data)
symp_data <- prep_symptom_data(symp, vita_data)
symptom_data <- check_symptom_data(symp_data)
check_symptom_grades(symptom_data)
check_other_symptom_grades(symptom_data)
sink()
## ------------------------------------------------------------------------------------------
write.table(x = merge(clin_data[,c("Label", "rangrp", "cov_sympday", "Sex", "Weight", "height", "age_yr")], symptom_data, by = "Label", all.y = T), 
            file = paste0(prefix_analysis_data, "/Analysis_Data/symptoms_interim_details.csv"), 
            row.names = F, sep=',', quote = F)
#############################################################################################
## 005 Final status database
sink("Queries/05_queries_final_status_database.txt", split = T)
final_status <- load_final_status_data(rand_app_data)
check_final_status(final_status)
sink()

#############################################################################################
## 006 Final status database
sink("Queries/06_queries_cbc_database.txt", split = T)
CBC_data <- load_CBC_data(rand_app_data)
CBC_data <-  check_date_time_cbc(CBC_data, clin_data)
CBC_data <- check_RBC_HB_PLT(CBC_data)
CBC_data <- check_WBCs(CBC_data)
sink()














# #############################################################################################
# clin_data$Symptom_onset = clin_data$cov_sympday
# clin_data$Trt = clin_data$rangrp
# clin_data$Sex = as.numeric(as.factor(clin_data$Sex))
# clin_data$Sex <- clin_data$Sex-1
# clin_data$Age <- clin_data$age_yr
# 
# clin_data <- clin_data[!is.na(clin_data$rangrp) &!(is.na(clin_data$randat) & is.na(clin_data$Rand_date_time)),]
# 
# # Select variables of interest
# clin_data = clin_data[, c('Label','Trt','Sex','Age','randat',
#                           "Rand_date_time",'BMI','Weight',
#                           'Symptom_onset','Site','Fever_Baseline')]
# #############################################################################################


