################################################################################
# Data prep for Nirmatrelvir PK analysis
################################################################################
library(tidyverse)
################################################################################
source('user_settings.R')
################################################################################
# Dosing database ###### ###### ###### ###### ###### ###### ######
data_name <- 'PLATCOV_Nirmatrelvir.dta'
file_name <- paste0(prefix_dropbox, "/Data/", data_name)
version <- file.info(file_name)$mtime %>% as.Date()
today <- Sys.Date()
query_file_name <- 'Queries/PLATCOV_queries_dosing_nirmatrelvir.csv'

# Loading data
dosing_data = haven::read_dta(file_name)
dosing_data = dosing_data[!is.na(dosing_data$da_dat),]
dosing_data <- dosing_data %>%
  mutate(Date_Time = (as.POSIXct(paste(da_dat, da_sttim, sep=' '), format = "%Y-%m-%d %H:%M:%S"))) %>%
  group_by(Label) %>%
  mutate(schedule_time = str_extract(da_9tp,"(?<=H)\\d+") %>% as.numeric(),
         Time = difftime(Date_Time, Date_Time[which(schedule_time==0)], units = 'hours') %>% as.numeric(),
         error = abs(Time-schedule_time) > 6,
         missing = (da_sttim == "")|is.na(da_sttim)) %>%
  ungroup()

ggplot(dosing_data, aes(x = schedule_time, y = Time)) +
  geom_point(aes(col = error)) +
  geom_abline(slope = 1) +
  theme_bw()

# Creating query file
sink(query_file_name, split = T)
writeLines(sprintf('PLATCOV data queries\nData: %s\nReceived date: %s\nQuery date: %s\n',
                   data_name,
                   version,
                   today)
)
sink()
write.table(data.frame('Dataset' = "",
                       'CRF form/Topic' = "",
                       'Question/Variable' = "",
                       'Query message' = "",
                       'Example data' = ""), 
            query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()

# Missing time points
flag_missing <-dosing_data %>%
  filter(missing)

if(nrow(flag_missing)>0){
  queries <- flag_missing %>% select(Label, da_9tp, da_dat, da_sttim)
  queries <- data.frame("Dataset" = data_name, #Dataset
                                        "CRF form/Topic" = "Drug administration/Nirmetrelvir", #'CRF form/Topic'
                                        "Question/Variable" = "Actual time", #'Question/Variable'
                        "Query message" = "Acutual time of these time points are missing. If the missing values are at D0H0, please check if the later time points also make sense.", #'Query message'
                        queries #'Example data'
  ) %>%
    mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
  
  write.table(queries, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
  cat("\n", file=query_file_name, append=TRUE)
  print(queries)
}  

# Mismatched time points
flag_error <-dosing_data %>%
  group_by(Label) %>%
  mutate(flag = any(error)) %>%
  ungroup() %>%
  filter(flag)

if(nrow(flag_error)>0){
  queries <- flag_error %>% filter(error) %>% select(Label, da_dat, da_sttim, da_9tp, Time)
  queries <- data.frame("Dataset" = data_name, #Dataset
                        "CRF form/Topic" = "Drug administration/Nirmetrelvir", #'CRF form/Topic'
                        "Question/Variable" = "Actual time", #'Question/Variable'
                        "Query message" = "Mismatched expected time (da_9tp) and the recorded time (Time, calculated from da_dat and da_sttim compared to D0H0). Please check recorded date and time.", #'Query message'
                        queries #'Example data'
  ) %>%
    mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
  
  write.table(queries, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
  cat("\n", file=query_file_name, append=TRUE)
  print(queries)
}  

# Dosing Check
flag_dose <- dosing_data %>%
  filter(da_tabs != 2 | da_totabs != 300)

if(nrow(flag_dose)>0){
  queries <- flag_dose %>% select(Label,  da_9tp, da_tabs, da_totabs)
  queries <- data.frame("Dataset" = data_name, #Dataset
                        "CRF form/Topic" = "Drug administration/Nirmetrelvir", #'CRF form/Topic'
                        "Question/Variable" = "No. of tabs/Total dose", #'Question/Variable'
                        "Query message" = "Either the number of tablets is not 2 or the total dose is not 300. Please check.", #'Query message'
                        queries #'Example data'
  ) %>%
    mutate(across(1:4, ~ if_else(row_number() > 1, "", .)))
  
  write.table(queries, query_file_name, col.names=T, sep=",", append=TRUE, row.names = F) %>% suppressWarnings()
  cat("\n", file=query_file_name, append=TRUE)
  print(queries)
}
################################################################################
################################################################################
# Pk data ###### ###### ###### ###### ###### ###### ######
pk_data = read.csv(paste0(prefix_dropbox, '/Data/PK/AR_23011_NIV.csv'))
pk_data = pk_data %>%
  mutate(Date_time = paste(Collection.date, Collection.time, sep = ' '),
         Date_time = as.POSIXct(Date_time, format = "%d-%b-%y %I:%M %p"))

# manual correction for sample barcodes
pk_data = pk_data %>%
  mutate(Sample.barcode = if_else(Patient.ID == "PLT-TH1-1350" & Protocol.time == "D3", "20MO137", Sample.barcode))

# Sample log data ###### ###### ###### ###### ###### ###### ######
sample_log = haven::read_dta(paste0(prefix_dropbox, '/Data/InterimSampleLog.dta'))
sample_log_pk = haven::read_dta(paste0(prefix_dropbox, '/Data/InterimSampleLog_Unscheduled_PK.dta'))

# Check if there is a sample with missing sample log
sample_log_check <- pk_data$Sample.barcode %in% sample_log_pk$sl3_barc | pk_data$Sample.barcode %in% sample_log$sl_barc 
pk_data %>% filter(!sample_log_check)

sample_log = sample_log %>%
  filter(sl_barc %in% pk_data$Sample.barcode) %>%
  mutate(sampling_date_time = paste(sl_sampdat, sl_samptim, sep = " ")) #%>%
#   select(Trial, Site, Label, RepeatNo, sl_barc, sl_sampdat, sl_samptim)
# 
sample_log_pk = sample_log_pk%>%
  filter(sl3_barc %in% pk_data$Sample.barcode) #%>%
#   select(Trial, Site, Label, RepeatNo, sl3_barc, sl3_sampdat, s3_samptim)

#colnames(sample_log_pk) <- colnames(sample_log)
#sample_log_PK <- rbind(sample_log, sample_log_pk)
#left_join(pk_data, sample_log_PK, by = c("Sample.barcode" = "sl_barc"))
################################################################################
clin_data <- read.csv("../Analysis_Data/cleaned_clinical_data.csv")
clin_data = clin_data %>% 
  filter(Label %in% pk_data$Patient.ID) %>% 
  select(Label, rangrp, Rand_date_time, age_yr, height, weight, cov_sympday, Sex)
################################################################################
pk_data_merged = left_join(clin_data, pk_data, by = c("Label" = "Patient.ID"))
pk_data_merged <- pk_data_merged %>%
  mutate(sampling_time_pr = difftime(Date_time, Rand_date_time, units = 'hours') %>% as.numeric(),
         D = str_extract(Protocol.time,"(?<=D)\\d+") %>% as.numeric(),
         H = str_extract(Protocol.time,"(?<=H)\\d+") %>% as.numeric(),
         H = if_else(is.na(H), 0, H),
         schedule_time = (24*D)+H,
         drug_level = if_else(`Nirmatrelvir..NIV...ng.mL.` == "< LLOQ", "0", `Nirmatrelvir..NIV...ng.mL.`),
         drug_level = as.numeric(drug_level))

ggplot(pk_data_merged, aes(x = schedule_time, y = sampling_time_pr)) +
  geom_point()

pk_data_merged %>% select(Label, Sample.barcode, Protocol.time, schedule_time, sampling_time_pr, Rand_date_time, Date_time) %>%
  filter(abs(schedule_time-sampling_time_pr) > 6 | sampling_time_pr < 0)
################################################################################
dosing_data_merged = left_join(clin_data, dosing_data, by = c("Label"))
dosing_data_merged <- dosing_data_merged %>%
  mutate(dosing_time_pr = difftime(Date_Time, Rand_date_time, units = 'hours') %>% as.numeric())

library(ggforce)

IDs <- pk_data_merged$Label %>% unique()
ind_plots <- list()
for(i in 1:length(IDs)){
  ID = IDs[i]
  ind_plots[[i]] <- ggplot() +
    geom_point(data = dosing_data_merged %>% filter(Label == ID), mapping = aes(x = dosing_time_pr, y = 0)) +
   # geom_vline(data = pk_data_merged  %>% filter(Label == ID), mapping = aes(xintercept = sampling_time_pr), col = 'red', linetype = "dashed") +
    geom_point(data = pk_data_merged  %>% filter(Label == ID), mapping = aes(x = sampling_time_pr, y = drug_level), col = 'red') +
    geom_line(data = pk_data_merged  %>% filter(Label == ID), mapping = aes(x = sampling_time_pr, y = drug_level), col = 'red') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = seq(0,100,24), limits = c(0,100)) +
    ylim(0, 20000) +
    ggtitle(ID) +
    #facet_wrap_paginate(Label~., ncol = 4, nrow = 5, page = 2) +
    theme_bw(base_size = 11) +
    xlab("") +
    ylab("")
}



ind_plot_all <- ggarrange(plotlist =  ind_plots, nrow = 7, ncol = 4, common.legend = T, legend = "right")
library(grid)

for(i in 1:length(ind_plot_all)){
  png(paste0("~/Project/PLATCOV-Nirmatrelvir_levels/Plots/Individual_plots_nir_level/Ind_plot_", i, ".png"), width = 12, height = 12, units = "in", res = 350)
  print(annotate_figure( ind_plot_all[[i]], bottom = textGrob("Time since randomisation (days)", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
                         left = textGrob("Drug plasma concentration (mg/mL)", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))))
  dev.off()
} 

