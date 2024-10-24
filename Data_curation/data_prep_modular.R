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
##Define user folder path####################################################################
source('user_settings.R')
source('functions.R')

source("000_load_randomisation_database.R")
source("001_clean_clinical_database.R")
source("002_clean_temperature_database.R")
source("003_clean_vital_database.R")

options(max.print = 5000)
#############################################################################################
# 000 Randomization database
rand_app_data <- load_randomisation_data()
#############################################################################################
# 001 Clinical database
sink("Queries/01_queries_clinical_database.txt", split = T)
clin_data <- load_clinical_data()
IDs_pending <- check_MACRO_clinical_database(clin_data, rand_app_data)
clin_data <- check_screen_failure(clin_data)
clin_data <- check_randomisation_info(clin_data)
clin_data <- check_sex(clin_data, IDs_pending, rand_app_data)
clin_data <- check_age(clin_data, IDs_pending, rand_app_data)
clin_data <- check_symptom_onset(clin_data, IDs_pending, rand_app_data)
clin_data <- check_symptom_onset(clin_data, IDs_pending, rand_app_data)
clin_data <- check_weight_height(clin_data, IDs_pending)
clin_data <- check_rand_date_time(clin_data, IDs_pending, rand_app_data)
clin_data <- check_rand_arms(clin_data, IDs_pending, rand_app_data)
sink()
#############################################################################################
## 002 Temperature database
sink("Queries/02_queries_temperature_database.txt", split = T)
temp_data <- load_temp_data(rand_app_data)
fever_data <- prep_tempdata(temp_data, clin_data)
check_time_temp(fever_data)
clin_data <- add_baseline_fever(clin_data, fever_data)
sink()
#############################################################################################
## 003 Vital sign database
sink("Queries/03_queries_vital_sign_database.txt", split = T)

vita_data <- load_vital_data(rand_app_data)
vita_data <- prep_vitadata(vita_data, clin_data)
check_vita_time(vita_data)


sink()



