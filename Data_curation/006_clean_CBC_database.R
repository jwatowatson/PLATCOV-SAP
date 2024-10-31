##*********************************************************
##******************* Complete blood count database *********
##*********************************************************
###########################################################################
# 1. Load CBC  data
load_CBC_data <- function(rand_app_data){
  writeLines('##########################################################################')
  writeLines('Reading CBC database from MACRO...')
  writeLines('##########################################################################')
  CBC_data = haven::read_dta(paste0(prefix_dropbox, "/Data/InterimCBC.dta"))
  
  id_data = CBC_data %>% distinct(Label) %>% pull(Label) %>% as.character()
  writeLines('### CBC database: Checking MACRO data entry progress:')
  writeLines(sprintf('Number of randomised patient: %s\nNumber of randomised patient CBC data on MACRO: %s', 
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
  return(CBC_data)
}

###########################################################################
# 2. Check date and time
check_date_time_cbc <- function(CBC_data, clin_data){
  writeLines('##########################################################################')
  writeLines('CBC database: Checking missing data by time points')
  missing_tp <- CBC_data %>%
    filter(is.na(fb_dat)) %>%
    select(Label, visit, fb_tp, fb_dat)
  writeLines(sprintf('%s rows have no information on CBC data:', 
                     nrow(missing_tp))
  )
  missing_tp %>% print(n = Inf, na.print = "NA")
  writeLines('##########################################################################')
  missing_fb_tim <- CBC_data %>%
    filter(!is.na(fb_dat) & fb_tim == "") %>%
    select(Label, visit, fb_tp, fb_dat, fb_tim)
  writeLines(sprintf('%s rows have no information on blood collection time:', 
                     nrow(missing_fb_tim))
  )
  missing_fb_tim %>% print(n = Inf, na.print = "NA")
 
  
  # Covert to time point ID
  CBC_data$Timepoint_ID <-  plyr::mapvalues(CBC_data$visit,  from = 0:3,
                                            to = c(0,3,7,14)) %>% as.factor()
  
  writeLines('##########################################################################')
  writeLines('CBC database: Check date and time')
  
  CBC_data = CBC_data[(!is.na(CBC_data$fb_dat) & 
                           !is.na(CBC_data$fb_tim)), ]
  CBC_data$cbc_time = 
    as.character(as.POSIXct(apply(CBC_data[, c('fb_dat','fb_tim')],
                                  1, paste, collapse =' '),
                            format = '%F %H:%M:%S'))
  CBC_data = merge(CBC_data, 
                    clin_data, 
                    by = c('Label','Site', 'Trial'),
                    all.x=T)
  CBC_data$time <- (difftime(as.Date(CBC_data$fb_dat), as.Date(CBC_data$Rand_date_time), units = "days")) %>% as.numeric()
  G <- ggplot(CBC_data, aes(x = Timepoint_ID, y = time-as.numeric(as.character(Timepoint_ID)))) +
    geom_jitter(size = 3, alpha = 0.5, width = 0.1) +
    theme_bw(base_size = 13) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") 
  print(G)

  CBC_time_mismatched <- CBC_data %>% 
    mutate(time_diff = abs(time - as.numeric(as.character(Timepoint_ID)))) %>%
    filter(time_diff > 1) %>%
    select(Label, fb_tp, fb_dat, Rand_date_time, Timepoint_ID, time, time_diff)

  writeLines(sprintf('%s rows blood collection time mismatched:', 
                     nrow(CBC_time_mismatched)))
  CBC_time_mismatched %>% print(n = Inf, na.print = "NA")
  
  return(CBC_data)
  }
###########################################################################
# 3. Check RBC//HB//PLT
check_RBC_HB_PLT <- function(CBC_data){
  writeLines('##########################################################################')
  writeLines('### CBC database: Missing data (RBC/HB/PLT):')
  ind_missing_RBC <- is.na(CBC_data$fb_rbc)
  ind_missing_HB <- is.na(CBC_data$fb_hb)
  ind_missing_PLT <- is.na(CBC_data$fb_plt)
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Red blood cell counts:', 
                     sum(ind_missing_RBC)))
  CBC_data$Label[ind_missing_RBC] %>% print()
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Haemoglobin:', 
                     sum(ind_missing_HB)))
  CBC_data$Label[ind_missing_HB] %>% print()
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Platelet counts:', 
                     sum(ind_missing_PLT)))
  CBC_data$Label[ind_missing_PLT] %>% print()
  
  writeLines('##########################################################################')
  writeLines('### CBC database: Checking outliers at 6 sd (RBC/HB/PLT):')
  writeLines('##########################################################################')
  CBC_data = CBC_data %>%
    mutate(
      outlier_RBC = abs(fb_rbc - mean(fb_rbc, na.rm = T)) > 6 * sd(fb_rbc, na.rm = T),
      outlier_HB = abs(fb_hb - mean(fb_hb, na.rm = T)) > 6 * sd(fb_hb, na.rm = T),
      outlier_PLT = abs(fb_plt - mean(fb_plt, na.rm = T)) > 6 * sd(fb_plt, na.rm = T)
           )
  
  writeLines(sprintf('Red blood cell counts of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_RBC, na.rm = T),
                     mean(CBC_data$fb_rbc, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_rbc, na.rm = T) %>% round(digits = 1)
                     ))
  CBC_data %>% filter(outlier_RBC) %>% select(Label, fb_tp, fb_rbc) %>% print()
  writeLines('## ---------------------------------------------------------------------')
  writeLines(sprintf('Haemoglobin of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_HB, na.rm = T),
                     mean(CBC_data$fb_hb, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_hb, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_HB) %>% select(Label, fb_tp, fb_hb) %>% print()
  writeLines('## ---------------------------------------------------------------------')
  writeLines(sprintf('Platelet counts of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_PLT, na.rm = T),
                     mean(CBC_data$fb_plt, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_plt, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_PLT) %>% select(Label, fb_tp, fb_plt) %>% print()
  
  A <- ggplot(CBC_data) +
    geom_point(aes(x = fb_rbc, y = fb_hb, col = Site), size = 2.5, alpha = 0.5) +
    theme_bw(base_size = 13) +
    xlab("Red blood cell counts") +
    ylab("Haemoglobin levels")
  
  B <- ggplot(CBC_data) +
    geom_point(aes(x = fb_rbc, y = fb_plt, col = Site), size = 2.5, alpha = 0.5) +
    theme_bw(base_size = 13) +
    xlab("Red blood cell counts") +
    ylab("Platelet counts")
  
  C <- ggarrange(A, B, ncol = 2, align = "h", common.legend = T, legend = "right")
  print(C)
  
  return(CBC_data)
}

###########################################################################
# 4. Check WBCs
check_WBCs <- function(CBC_data){
  writeLines('##########################################################################')
  writeLines('##########################################################################')
  writeLines('### CBC database: Missing data (WBCs):')
  
  ind_missing_WBC <- is.na(CBC_data$fb_wbc)
  # ind_missing_NEU <- is.na(CBC_data$fb_neu)
  # ind_missing_LYMP <- is.na(CBC_data$fb_lymp)
  # ind_missing_EOS <- is.na(CBC_data$fb_eos)
  # ind_missing_MONO <- is.na(CBC_data$fb_mono)
  # ind_missing_BASO <- is.na(CBC_data$fb_baso)
  
  ind_missing_NEUPER <- is.na(CBC_data$fb_neuper)
  ind_missing_LYMPPER <- is.na(CBC_data$fb_lympper)
  ind_missing_EOSPER <- is.na(CBC_data$fb_eosper)
  ind_missing_MONOPER <- is.na(CBC_data$fb_monoper)
  ind_missing_BASOPER <- is.na(CBC_data$fb_basoper)
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Total white blood cell counts:', 
                     sum(ind_missing_WBC)))
  CBC_data$Label[ind_missing_WBC] %>% unique() %>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Neutrophil percentage:', 
                     sum(ind_missing_NEUPER)))
  CBC_data$Label[ind_missing_NEUPER] %>% unique() %>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Lymphocyte percentage:', 
                     sum(ind_missing_LYMPPER)))
  CBC_data$Label[ind_missing_LYMPPER] %>% unique()%>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Eosinophil percentage:', 
                     sum(ind_missing_EOSPER)))
  CBC_data$Label[ind_missing_EOSPER] %>% unique()%>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Monocyte percentage:', 
                     sum(ind_missing_MONOPER)))
  CBC_data$Label[ind_missing_MONOPER] %>% unique()%>% print()
  
  writeLines('##########################################################################')
  writeLines(sprintf('%s rows are missing for Basophil percentage:', 
                     sum(ind_missing_BASOPER)))
  CBC_data$Label[ind_missing_BASOPER] %>% unique()%>% print()
  
  
  writeLines('##########################################################################')
  writeLines('### CBC database: Checking outliers at 6 sd (WBCs):')
  writeLines('##########################################################################')
  writeLines('## [Manual correction]: 80 patients from Brazil from August 2023 to March 2024 used the wrong unit')
  writeLines('##########################################################################')
  ind_WBC_wrong <- CBC_data$Site == "br003" & CBC_data$fb_wbc > 1000 & !is.na(CBC_data$fb_wbc)
  print(CBC_data[ind_WBC_wrong,c('Label', 'fb_tp', 'fb_dat', 'fb_wbc')], n = Inf, na.print = "NA")
  
  CBC_data$fb_wbc[ind_WBC_wrong] <- CBC_data$fb_wbc[ind_WBC_wrong]/1000
  
  CBC_data = CBC_data %>%
    mutate(
      outlier_WBC = abs(fb_wbc - mean(fb_wbc, na.rm = T)) > 6 * sd(fb_wbc, na.rm = T),
      outlier_NEUPER = abs(fb_neuper - mean(fb_neuper, na.rm = T)) > 6 * sd(fb_neuper, na.rm = T),
      outlier_LYMPPER = abs(fb_lympper - mean(fb_lympper, na.rm = T)) > 6 * sd(fb_lympper, na.rm = T),
      outlier_EOSPER = abs(fb_eosper - mean(fb_eosper, na.rm = T)) > 6 * sd(fb_eosper, na.rm = T),
      outlier_MONOPER = abs(fb_monoper - mean(fb_monoper, na.rm = T)) > 6 * sd(fb_monoper, na.rm = T),
      outlier_BASOPER  = abs(fb_basoper - mean(fb_basoper, na.rm = T)) > 6 * sd(fb_basoper, na.rm = T)
    )
  #---------------------------------------------------------------------------------------------------------
  writeLines(sprintf('Total white blood cell counts of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_WBC, na.rm = T),
                     mean(CBC_data$fb_wbc, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_wbc, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_WBC) %>% select(Label, fb_tp, fb_wbc) %>% print()
  
  G <- ggplot(CBC_data, aes(x = fb_dat, y = fb_wbc, col = Site)) +
    geom_point(size = 3, alpha = 0.5) +
    theme_bw() +
    facet_wrap(.~ Site)
  print(G)
  writeLines('## ---------------------------------------------------------------------')
  writeLines(sprintf('Lymphocyte percentage of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_LYMPPER, na.rm = T),
                     mean(CBC_data$fb_lympper, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_lympper, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_LYMPPER) %>% select(Label, fb_tp, fb_lympper) %>% print()
  
  G_lymp <- ggplot(CBC_data, aes(x = fb_dat, y = fb_lympper, col = Site)) +
    geom_point(size = 3, alpha = 0.5) +
    theme_bw() +
    facet_wrap(.~ Site)
  print(G_lymp)
  
  writeLines('## ---------------------------------------------------------------------')
  writeLines(sprintf('Neutrophil percentage of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_NEUPER, na.rm = T),
                     mean(CBC_data$fb_neuper, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_neuper, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_NEUPER) %>% select(Label, fb_tp, fb_neuper) %>% print()
  
  G_neu <- ggplot(CBC_data, aes(x = fb_dat, y = fb_neuper, col = Site)) +
    geom_point(size = 3, alpha = 0.5) +
    theme_bw() +
    facet_wrap(.~ Site)
  print(G_neu)
  
  writeLines('## ---------------------------------------------------------------------')
  writeLines(sprintf('Eosinophil percentage of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_EOSPER, na.rm = T),
                     mean(CBC_data$fb_eosper, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_eosper, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_EOSPER) %>% select(Label, fb_tp, fb_eosper) %>% print()
  
  G_eos <- ggplot(CBC_data, aes(x = fb_dat, y = fb_eosper, col = Site)) +
    geom_point(size = 3, alpha = 0.5) +
    theme_bw() +
    facet_wrap(.~ Site)
  print(G_eos)
  
  writeLines('## ---------------------------------------------------------------------')
  writeLines(sprintf('Monocyte percentage of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_MONOPER, na.rm = T),
                     mean(CBC_data$fb_monoper, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_monoper, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_MONOPER) %>% select(Label, fb_tp, fb_monoper) %>% print()
  
  G_mono <- ggplot(CBC_data, aes(x = fb_dat, y = fb_monoper, col = Site)) +
    geom_point(size = 3, alpha = 0.5) +
    theme_bw() +
    facet_wrap(.~ Site)
  print(G_mono)
  
  writeLines('## ---------------------------------------------------------------------')
  writeLines(sprintf('Basophil percentage of %s patients are outliers from mean = %s, sd = %s:', 
                     sum(CBC_data$outlier_BASOPER, na.rm = T),
                     mean(CBC_data$fb_basoper, na.rm = T) %>% round(digits = 1),
                     sd(CBC_data$fb_basoper, na.rm = T) %>% round(digits = 1)
  ))
  CBC_data %>% filter(outlier_BASOPER) %>% select(Label, fb_tp, fb_basoper) %>% print()
  
  G_baso <- ggplot(CBC_data, aes(x = fb_dat, y = fb_basoper, col = Site)) +
    geom_point(size = 3, alpha = 0.5) +
    theme_bw() +
    facet_wrap(.~ Site)
  print(G_baso)
  
  return(CBC_data)
}



# 
# ind_missing_NEUPER <- is.na(CBC_data$fb_neuper)
# ind_missing_LYMPPER <- is.na(CBC_data$fb_lympper)
# ind_missing_EOSPER <- is.na(CBC_data$fb_eosper)
# ind_missing_MONOPER <- is.na(CBC_data$fb_monoper)
# ind_missing_BASOPER <- is.na(CBC_data$fb_basoper)
# 
CBC_data <- CBC_data %>%
  mutate(fb_allWBCper = CBC_data[,c('fb_neuper', 'fb_lympper', 'fb_eosper', 'fb_monoper', 'fb_basoper')] %>% rowSums())
  
ggplot(CBC_data) +
  geom_point(aes(x = fb_dat, y = fb_allWBCper)) +
  facet_grid(.~Site)

colnames(CBC_data)

which(round(CBC_data$fb_allWBCper) < 100)

table(CBC_data$fb_bandper, useNA = 'always')

CBC_data$wbc_percentage_sum <- CBC_data[,c('fb_neuper', 'fb_lympper', 'fb_eosper', 'fb_monoper', 'fb_basoper')] %>% rowSums()

CBC_data %>% 
  filter(Label %in% c('PLT-TH1-495', 'PLT-TH1-496', 'PLT-TH1-638', 'PLT-TH1-639')) %>%
  select(Label, visit, fb_neuper, fb_lympper, fb_eosper, fb_monoper, , fb_basoper, wbc_percentage_sum)


CBC_data[6627,]


# 
# ggplot(CBC_data, aes(x = fb_basoper, y = fb_monoper, col = Site)) +
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw()
