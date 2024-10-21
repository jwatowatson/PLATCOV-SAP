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
#############################################################################################
rand_app_data <- load_randomisation_data()




sink("Queries/console_output.txt", split = T)
clin_data <- load_clinical_data()

IDs_pending <- check_MACRO_clinical_database(clin_data, rand_app_data)

clin_data <- check_screen_failure(clin_data)
clin_data <-check_randomisation_info(clin_data)






sink()





#############################################################################################









