library(tidyverse)
library(gridExtra)
library(pheatmap)
library(scales)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggtext)
library(forcats)
library(randomForest)

################################################################################
# 1. Data Loading
################################################################################

### Path Specification
# Mac
# abs_path <- "~/Dropbox/PLATCOV_Analysis/Data/Cytokines"
# Windows
abs_path <- path.expand("~/../Dropbox/PLATCOV_Analysis/Data/Cytokines")

create_df <- F # True if want to recreate dataframe

# Choices of csv data

# Don't use this because some batch has deviated standard curve
# ff_combo0 = list.files(path=abs_path,pattern = '*.csv',full.names = T, recursive = T)

# Use this
ff_combo1 = list.files(path=file.path(abs_path, 'Viral combo 1'),pattern = '*.csv',full.names = T, recursive = T)
ff_combo2 = list.files(path=file.path(abs_path, 'Viral combo 2'),pattern = '*.csv',full.names = T, recursive = T)
ff = c(ff_combo1, ff_combo2)

if (create_df) {
    #ff = ff_combo0
    dat = list()
    for(i in 1:length(ff)){
      plate_name = gsub(pattern = ',',replacement = '',x = readLines(ff[i], n = 1))
      dat[[i]]=readr::read_csv(file = ff[i], skip = 1)
      dat[[i]]$plate=plate_name
    }
    dat_all = dplyr::bind_rows(dat)
    
    # Fill in day of samples
    dat_all$day = NA
    dat_all$day[grep('D0',x = dat_all$Sample)]='D0'
    dat_all$day[grep('D3',x = dat_all$Sample)]='D3'
    dat_all$day[grep('D7',x = dat_all$Sample)]='D7'
    dat_all$day[grep('D14',x = dat_all$Sample)]='D14'
    
    cytokines = unique(dat_all$Assay)
    #dat_all = dat_all %>% filter(Sample!='Control sample')
    control_data = dat_all %>% filter(Sample=='Control')
    standards_data = dat_all %>% filter(!is.na(Concentration))
    
    dat_all = dat_all %>% filter(Sample!='Control', is.na(Concentration))
    write_csv(x = control_data, file = 'Analysis_Data/cytokine_control_merged.csv')
    write_csv(x = dat_all, file = 'Analysis_Data/cytokine_data_merged.csv')
    write_csv(x = standards_data, file = 'Analysis_Data/cytokine_standards_merged.csv')
}

# ff_combo1 = list.files(path='~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/Data/Cytokines/Viral combo 1/',pattern = '*.csv',full.names = T, recursive = T)
# ff_combo2 = list.files(path='~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/Data/Cytokines/Viral combo 2/',pattern = '*.csv',full.names = T, recursive = T)
# ff = c(ff_combo1, ff_combo2)
dat = list()
for(i in 1:length(ff)){
  plate_name = gsub(pattern = ',',replacement = '',x = readLines(ff[i], n = 1))
  dat[[i]]=readr::read_csv(file = ff[i], skip = 1)
  dat[[i]]$plate=plate_name
}
dat_all = dplyr::bind_rows(dat)

dat_all$day = NA
dat_all$day[grep('D0',x = dat_all$Sample)]='D0'
dat_all$day[grep('D3',x = dat_all$Sample)]='D3'
dat_all$day[grep('D7',x = dat_all$Sample)]='D7'
dat_all$day[grep('D14',x = dat_all$Sample)]='D14'

cytokines = unique(dat_all$Assay)
control_data = dat_all %>% filter(Sample=='Control')
standards_data = dat_all %>% filter(!is.na(Concentration))

dat_all = dat_all %>% filter(Sample!='Control', is.na(Concentration))
write_csv(x = dat_all, file = 'Analysis_Data/cytokine_data_merged.csv')
write_csv(x = standards_data, file = 'Analysis_Data/cytokine_standards_merged.csv')


pdf('cytokine_data.pdf')

for(cytk in cytokines){
  cytk_dat = dat_all %>% filter(Assay==cytk, Concentration>0 | is.na(Concentration))
  cytk_std_dat = standards_data %>% filter(Assay==cytk)
  my_lims_x = range(c(cytk_dat$Signal, cytk_std_dat$Signal))
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

  h1 = cytk_dat %>% filter(is.na(Concentration)) %>% ggplot(aes(x=Signal, colour = day))+
    geom_density() + scale_x_log10(limits = my_lims_x)+theme_minimal()

  p1=cytk_std_dat %>% ggplot(aes(y=Concentration, x = Signal, colour = plate))+
    scale_x_log10()+
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,)+
    geom_point()+theme_minimal()+ylab('pg/ml')+
    geom_smooth(se = F)+
    annotation_logticks(sides='l')+ggtitle(cytk)
  print(grid.arrange(h1,p1))
  cytk_dat %>% filter(is.na(Concentration)) %>% ggplot(aes(x=Signal))+
    geom_histogram() + scale_x_log10(limits=c(1,10^6))
}

# CBC Data
file_path_plt <- paste0("C:/Users/Kantapong/Dropbox/PLATCOV_Analysis", 
                        "/Data/")
cbc_data_plt <- haven::read_dta(paste0(file_path_plt, "InterimCBC.dta"))

file_path_ast <- paste0("C:/Users/Kantapong/Dropbox/ADASTRA_Analysis", 
                        "/Data/")
cbc_data_ast <- haven::read_dta(paste0(file_path_ast, "ASTInterimCBC.dta"))

# Process wbc count
cbc_data_plt$fb_neu <- cbc_data_plt$fb_wbc * cbc_data_plt$fb_neuper / 100
cbc_data_plt$fb_lymp <- cbc_data_plt$fb_wbc * cbc_data_plt$fb_lympper / 100
cbc_data_plt$fb_eos <- cbc_data_plt$fb_wbc * cbc_data_plt$fb_eosper / 100
cbc_data_plt$fb_mono <- cbc_data_plt$fb_wbc * cbc_data_plt$fb_monoper / 100
cbc_data_plt$fb_baso <- cbc_data_plt$fb_wbc * cbc_data_plt$fb_basoper / 100

cbc_data_ast$fb_neu  <- cbc_data_ast$fb_wbc * cbc_data_ast$fb_neuper / 100
cbc_data_ast$fb_lymp <- cbc_data_ast$fb_wbc * cbc_data_ast$fb_lympper / 100
cbc_data_ast$fb_eos  <- cbc_data_ast$fb_wbc * cbc_data_ast$fb_eosper / 100
cbc_data_ast$fb_mono <- cbc_data_ast$fb_wbc * cbc_data_ast$fb_monoper / 100
cbc_data_ast$fb_baso <- cbc_data_ast$fb_wbc * cbc_data_ast$fb_basoper / 100
################################################################################
# 2. Write PDF of Cytokine's Calculated Concentration Distribution
################################################################################
# As it turns out, using all csv files might not be viable due to batch effect
# Some standard curve deviates largely from
# control_data <- read_csv(file = 'Analysis_Data/cytokine_control_merged2.csv')
# dat_all <- read_csv(file = 'Analysis_Data/cytokine_data_merged2.csv')

control_data <- read_csv(file = 'Analysis_Data/cytokine_control_merged.csv')
dat_all <- read_csv(file = 'Analysis_Data/cytokine_data_merged.csv')
standards_data <- read_csv(file = 'Analysis_Data/cytokine_standards_merged.csv')

control_data2 <- read_csv(file = 'Analysis_Data/cytokine_control_merged2.csv')
dat_all2 <- read_csv(file = 'Analysis_Data/cytokine_data_merged2.csv')
standards_data2 <- read_csv(file = 'Analysis_Data/cytokine_standards_merged2.csv')

write_pdf <- F # True if want to write PDF

if (write_pdf) {
    jan2025_plates <- setdiff(
      unique(standards_data2$plate),
      unique(standards_data$plate)
    )
    
    standards_both <- bind_rows(
        standards_data  %>% mutate(batch = "Others"),
        standards_data2 %>% mutate(batch = if_else(plate %in% jan2025_plates, "Jan2025", "Others"))
    )
    
    # Optional: set nice palette for the two groups
    pal_batch <- c("Jan2025" = "#d62728",  # red
                   "Others"  = "#1f77b4")  # blue
    
    pdf("cytokine_data_batch.pdf")
    
    for (cytk in unique(dat_all$Assay)) {
      # All samples for this cytokine (as before)
      cytk_dat <- dat_all %>%
        filter(Assay == cytk, Concentration > 0 | is.na(Concentration))
      
      # Standards for this cytokine from BOTH tables, labeled by batch
      cytk_std_dat <- standards_both %>% filter(Assay == cytk)
      
      my_lims_x <- range(c(cytk_dat$Signal, cytk_std_dat$Signal), na.rm = TRUE)
      breaks        <- 10^(-10:10)
      minor_breaks  <- rep(1:9, 21) * (10^rep(-10:10, each = 9))
      
      # Density of unknowns (per day) – unchanged
      h1 <- cytk_dat %>%
        filter(is.na(Concentration)) %>%
        ggplot(aes(x = Signal, colour = day)) +
        geom_density() +
        scale_x_log10(limits = my_lims_x) +
        theme_minimal()
      
      # Standard curve: color by batch (Jan2025 vs Others)
      p1 <- cytk_std_dat %>%
        ggplot(aes(y = Concentration, x = Signal, colour = batch)) +
        scale_x_log10() +
        scale_y_log10(breaks = breaks, minor_breaks = minor_breaks) +
        geom_point(alpha = 0.8) +
        theme_minimal() +
        ylab("pg/ml") +
        geom_smooth(se = FALSE) +
        annotation_logticks(sides = "l") +
        ggtitle(cytk) +
        scale_colour_manual(values = pal_batch, name = "Standards batch")
      
      print(grid.arrange(h1, p1))
      
      # (Optional) Histogram of unknowns
      cytk_dat %>%
        filter(is.na(Concentration)) %>%
        ggplot(aes(x = Signal)) +
        geom_histogram() +
        scale_x_log10(limits = c(1, 10^6))
    }
    dev.off()
}

if (write_pdf) {
    pdf('cytokine_data.pdf')
    
    for(cytk in cytokines){
      cytk_dat = dat_all %>% filter(Assay==cytk, Concentration>0 | is.na(Concentration))
      cytk_std_dat = standards_data %>% filter(Assay==cytk)
      my_lims_x = range(c(cytk_dat$Signal, cytk_std_dat$Signal))
      breaks <- 10^(-10:10)
      minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
      
      # All samples that belong to such cytokine, for each day
      h1 = cytk_dat %>% filter(is.na(Concentration)) %>% ggplot(aes(x=Signal, colour = day))+
        geom_density() + scale_x_log10(limits = my_lims_x)+theme_minimal()
      
      # Standard curve for such cytokine, for each plate
      p1=cytk_std_dat %>% ggplot(aes(y=Concentration, x = Signal, colour = plate))+
        scale_x_log10()+
        scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,)+
        geom_point()+theme_minimal()+ylab('pg/ml')+
        # theme(legend.position = 'none')+
        geom_smooth(se = F)+
        annotation_logticks(sides='l')+ggtitle(cytk)
      print(grid.arrange(h1,p1))
      cytk_dat %>% filter(is.na(Concentration)) %>% ggplot(aes(x=Signal))+
        geom_histogram() + scale_x_log10(limits=c(1,10^6))
    }
    dev.off()
}

################################################################################
# 3. Data Preprocessing
################################################################################
dat_all = dat_all %>% filter(is.na(Concentration))
dat_all$ID=NA
for (i in 1:nrow(dat_all)) {
  xx = strsplit(x = dat_all$Sample[i],split = '-',fixed = T)[[1]]
  dat_all$ID[i] = paste(xx[1],xx[2],xx[3],sep = '-')
  #dat_all$ID[i] = paste(xx[1],xx[2],xx[4],sep = '-')
}

all_cytokines <- unique(dat_all$Assay)

# Pick cytokines with high variation among days
# cytokines = c("IFN-g","IL-1RA","IP-10") -- previous analysis: use to plot 
cytokines = c("IFN-g","IL-10", "IL-6", "IL-8", "G-CSF", "IFN-a2a",
              "IL-1RA","IP-10", "MCP-1")

# Duplicated rows (some patients belong in more than one plate)
# The issue is that for this patient AST-TH1-205, D14 appear twice when it is likely D7 and D14
dup_IDs <- dat_all %>%
  group_by(ID, day, plate, Assay) %>%
  filter(n() > 1) #%>%
#  distinct(ID, day)

# Remove those rows
# 508 patients in total
dat_all <- dat_all %>%
  anti_join(dup_IDs, by = c("ID", "day"))

# Retrieve PCR data (PLATCOV) and join with our cytokine data
# Nirmatrelvir+Ritonavir // NSD // Regeneron
# 324 patients
pcr_dat_plt = read_csv('Analysis_Data/interim_all_analysis.csv')
pcr_dat_plt = pcr_dat_plt %>% filter(ID %in% dat_all$ID)#, Time<5)

# Retrieve PCR data (AD ASTRA) and join with our cytokine data
# Favipiravir // NSD
# 184 patients
pcr_dat_ast = read_csv('../ADASTRA-SAP/Analysis_Data/interim_all_analysis.csv')
pcr_dat_ast = pcr_dat_ast %>% filter(ID %in% dat_all$ID)#, Time<5)

# Retrieve leukocytes data for PLATCOV and AD ASTRA
# Also 324 for platcov and 184 for ad astra --> no missing
cbc_data_plt = cbc_data_plt %>% filter(Label %in% dat_all$ID)
cbc_data_ast = cbc_data_ast %>% filter(Label %in% dat_all$ID)

cbc_data_plt <- cbc_data_plt %>%
  mutate(fb_tp = factor(fb_tp, levels = c("Baseline", "D3", "D7", "D14"), ordered = TRUE))
cbc_data_ast <- cbc_data_ast %>%
  mutate(fb_tp = factor(fb_tp, levels = c("Baseline", "D3", "D7", "D14"), ordered = TRUE))

# Timeline of recruitment
plot_timeline <- function(df, study_name = "Study") {
  df2 <- df %>%
    mutate(Rand_date = as.Date(Rand_date))
  
  # Extract only one row for each id, will contain ID, Trt, Rand_date columns
  id_rand <- df2 %>%
    filter(!is.na(Rand_date), !is.na(ID), !is.na(Trt)) %>%
    group_by(ID, Trt) %>%
    summarise(rand_date = min(Rand_date), .groups = "drop") %>%
    group_by(ID) %>%   # if an ID appears under multiple Trt, keep earliest overall
    slice_min(rand_date, with_ties = FALSE) %>%
    ungroup()
  
  # Plot
  p <- ggplot(id_rand, aes(x = rand_date, y = Trt, color = Trt)) +
    geom_point(alpha = 0.8, size = 2) +
    labs(
      title = paste0("Recruitment timeline — ", study_name),
      x = "Randomisation date",
      y = "Treatment"
    ) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
  
  return(p)
}

timeline_ast <- plot_timeline(pcr_dat_ast, "AD ASTRA")
timeline_plt <- plot_timeline(pcr_dat_plt, "PLATCOV")

# Enrollment timeline plot
timeline_ast
timeline_plt

# Enrollment table count
pcr_dat_ast %>%
  group_by(Trt) %>%
  summarise(n_enroll = n_distinct(ID), .groups = "drop")

pcr_dat_plt %>%
  group_by(Trt) %>%
  summarise(n_enroll = n_distinct(ID), .groups = "drop")

##### PLATCOV Cytokine Data #####
# Keep PLATCOV patients only and transform "D0" to numeric days
dat_all_plt <- dat_all %>%
  filter(grepl("^PLT", ID)) %>% 
  mutate(day_num = as.numeric(gsub("D", "", day)))

# Number of patients that have each day_num
patients_per_day <- dat_all_plt %>%
  group_by(day_num) %>%
  summarise(n_patients = n_distinct(ID)) %>%
  arrange(day_num)

print(patients_per_day)

# Number of patients that have # of days of data (1–4)
patients_by_num_days <- dat_all_plt %>%
  group_by(ID) %>%
  summarise(n_days = n_distinct(day_num)) %>%
  count(n_days, name = "n_patients") %>%
  arrange(n_days)

print(patients_by_num_days)

cbc_data_plt <- cbc_data_plt %>%
  mutate(
    day_num = case_when(
      grepl("^baseline$", fb_tp, ignore.case = TRUE) ~ 0,
      grepl("^D\\s*\\d+$", fb_tp, ignore.case = TRUE) ~ as.numeric(gsub("[^0-9]", "", fb_tp)),
      grepl("^\\d+$", fb_tp) ~ as.numeric(fb_tp),   # handles plain numbers like "0","3","7","14"
      TRUE ~ NA_real_
    )
  )

# Convert dat_all_plt into wide format
# capture all assays so pivot_wider creates ALL columns even if absent in some groups
all_assays <- sort(unique(dat_all_plt$Assay))

dat_all_plt_wide <- dat_all_plt %>%
  mutate(Assay = factor(Assay, levels = all_assays)) %>%  # force full set of assay cols
  pivot_wider(
    id_cols     = c(ID, day_num),
    names_from  = Assay,
    values_from = `Calc. Concentration`,
    values_fill = NA_real_,
    # if duplicates exist per (ID, day_num, Assay), pick a rule:
    values_fn   = list(`Calc. Concentration` = ~ dplyr::first(na.omit(.x)))
  )

# Join dat_all_plt with pcr_dat_plt and cbc_data_plt by ID and Label
# left join on pcr_dat_plt which offers fullest history
pcr_cytokine_dat_plt <- pcr_dat_plt %>%
  left_join(
    dat_all_plt_wide,
    by = c("ID" = "ID", "Timepoint_ID" = "day_num"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    cbc_data_plt,
    by = c("ID" = "Label", "Timepoint_ID" = "day_num"),
    relationship = "many-to-many"
  )



# inner join; will only contain D0, D3
# pcr_cytokine_dat_plt <- pcr_dat_plt %>%
#   inner_join(dat_all_plt, by = c("ID" = "ID", "Timepoint_ID" = "day_num"),
#              relationship = "many-to-many") %>% 
#   inner_join(cbc_data_plt, by = c("ID" = "Label", "Timepoint_ID" = "day_num"),
#              relationship = "many-to-many")

# 50 patients are lost due to cytokine samples not obtained on D0-D5
lost_plt_ids <- setdiff(
  unique(pcr_dat_plt$ID),
  unique(pcr_cytokine_dat_plt$ID)
)

View(dat_all_plt %>% filter(ID %in% lost_plt_ids))

# Keep only patients with 120 rows (80 samples from D0, 40 samples from D3) -- 269 Patients
# For context,
# Patient 114, 227, 309, 359, 390 only have value from D0
# Patient 1677 only has value from D3
pcr_cytokine_dat_plt_filtered <- pcr_cytokine_dat_plt %>%
  group_by(ID) %>%
  filter(n() == 120) %>%
  ungroup()

setdiff(unique(pcr_cytokine_dat_plt$ID), unique(pcr_cytokine_dat_plt_filtered$ID))

##### AD ASTRA Cytokine Data #####
# Keep AD ASTRA patients only and transform "D0" to numeric days
dat_all_ast <- dat_all %>%
  filter(grepl("^AST", ID)) %>% 
  mutate(day_num = as.numeric(gsub("D", "", day)))

# Number of patients that have each day_num
patients_per_day <- dat_all_ast %>%
  group_by(day_num) %>%
  summarise(n_patients = n_distinct(ID)) %>%
  arrange(day_num)

print(patients_per_day)

# Number of patients that have # of days of data (1–4)
patients_by_num_days <- dat_all_ast %>%
  group_by(ID) %>%
  summarise(n_days = n_distinct(day_num)) %>%
  count(n_days, name = "n_patients") %>%
  arrange(n_days)

print(patients_by_num_days)

cbc_data_ast <- cbc_data_ast %>%
  mutate(
    day_num = case_when(
      grepl("^baseline$", fb_tp, ignore.case = TRUE) ~ 0,
      grepl("^D\\s*\\d+$", fb_tp, ignore.case = TRUE) ~ as.numeric(gsub("[^0-9]", "", fb_tp)),
      grepl("^\\d+$", fb_tp) ~ as.numeric(fb_tp),   # handles plain numbers like "0","3","7","14"
      TRUE ~ NA_real_
    )
  )

# Pivot cytokine from ad astra into wide format
all_assays <- sort(unique(dat_all_ast$Assay))

dat_all_ast_wide <- dat_all_ast %>%
  mutate(Assay = factor(Assay, levels = all_assays)) %>%  # force full set of assay cols
  pivot_wider(
    id_cols     = c(ID, day_num),
    names_from  = Assay,
    values_from = `Calc. Concentration`,
    values_fill = NA_real_,
    # if duplicates exist per (ID, day_num, Assay), pick a rule:
    values_fn   = list(`Calc. Concentration` = ~ dplyr::first(na.omit(.x)))
  )

# Join dat_all_ast with pcr_dat_ast and cbc_data_ast by ID and Label
# left join on pcr_dat_ast which offers fullest history
pcr_cytokine_dat_ast <- pcr_dat_ast %>%
  left_join(
    dat_all_ast_wide,
    by = c("ID" = "ID", "Timepoint_ID" = "day_num"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    cbc_data_ast,
    by = c("ID" = "Label", "Timepoint_ID" = "day_num"),
    relationship = "many-to-many"
  )

ast_fever <- read_csv("D:/Data-Work/ADASTRA-SAP/Analysis_Data/interim_all_fever_analysis.csv")
ast_symp <- read_csv("D:/Data-Work/ADASTRA-SAP/Analysis_Data/interim_all_symptom_analysis.csv")

ast_fever <- ast_fever %>% select(Label, visit, fut_dat, fut_tim, fut_temp, temp_time)
ast_symp <- ast_symp %>% select(-Trial, -Site)

# Join dat_all_ast with pcr_dat_ast and cbc_data_ast by ID and Label
# will only contain D0, D3
# pcr_cytokine_dat_ast <- pcr_dat_ast %>%
#   inner_join(dat_all_ast, by = c("ID" = "ID", "Timepoint_ID" = "day_num"),
#              relationship = "many-to-many") %>% 
#   inner_join(cbc_data_ast, by = c("ID" = "Label", "Timepoint_ID" = "day_num"),
#              relationship = "many-to-many")

# No patients are lost due to cytokine samples not obtained on D0-D5
lost_ast_ids <- setdiff(
  unique(pcr_dat_ast$ID),
  unique(pcr_cytokine_dat_ast$ID)
)

#View(dat_all_ast %>% filter(ID %in% lost_ast_ids))

# Keep only patients with 120 rows (80 samples from D0, 40 samples from D3) -- 183 Patients
# For context,
# Patient 320 only have value from D0
pcr_cytokine_dat_ast_filtered <- pcr_cytokine_dat_ast %>%
  group_by(ID) %>%
  filter(n() == 120) %>%
  ungroup()

setdiff(unique(pcr_cytokine_dat_ast$ID), unique(pcr_cytokine_dat_ast_filtered$ID))

################################################################################
# 4. EDA for Cytokine and Relevance to PCR Data
################################################################################

# Figure 0: Look at PLT/AST studies, and plot WBC levels across time
cbc_data_long <- cbc_data_ast %>%
  mutate(fb_tp = factor(fb_tp, levels = c("Baseline", "D3", "D7", "D14"), ordered = TRUE)) %>%
  pivot_longer(
    cols = c(fb_neu, fb_lymp, fb_eos, fb_mono, fb_baso),
    names_to = "Cell_type",
    values_to = "Count"
  )

ggplot(
  cbc_data_long %>% filter(fb_tp %in% c("Baseline", "D3", "D7", "D14")),
  aes(x = fb_tp, y = Count)
) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.6, fill = "lightblue") +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8, color = "gray40") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 2) +
  facet_wrap(~Cell_type, scales = "free_y") +
  labs(
    x = "Timepoint",
    y = "Cell count",
    title = "Trajectory of white blood cell subtypes across timepoints",
    subtitle = "Box-and-whisker + mean trajectory (red line) per subtype"
  ) +
  theme_minimal(base_size = 14)

# Figure 1: Look at PLT/AST studies, and plot cytokine levels across time


# Median across patients (for overlay lines) (D0, D3, D7, D14)
median_across_patients_plt <- dat_all_plt %>%
  group_by(Assay, day_num) %>%
  summarise(median_conc = median(`Calc. Concentration`, na.rm = TRUE), .groups = "drop")

# Median across patients (for overlay lines) (D0, D3, D7, D14)
median_across_patients_ast <- dat_all_ast %>%
  group_by(Assay, day_num) %>%
  summarise(median_conc = median(`Calc. Concentration`, na.rm = TRUE), .groups = "drop")

### PLATCOV
# Make sure day_num is a factor so axis shows exactly 0, 3, 7, 14
median_across_patients_plt <- median_across_patients_plt %>%
  mutate(
    day_num = factor(day_num, levels = c(0, 3, 7, 14)),
    Assay   = factor(Assay, levels = all_cytokines)  # control facet order
  )

stats_plt <- median_across_patients_plt %>%
  group_by(Assay) %>%
  summarise(
    n               = sum(!is.na(median_conc)),
    mean_conc       = mean(median_conc, na.rm = TRUE),
    sd_conc         = sd(median_conc,   na.rm = TRUE),
    sd_over_mean    = ifelse(mean_conc == 0, NA_real_, sd_conc / mean_conc),
    .groups = "drop"
  ) %>%
  arrange(desc(sd_over_mean))

# filter to selected cytokines only
df_plot <- median_across_patients_plt %>% filter(Assay %in% cytokines)

# compute dynamic figure height (inches): ~2.4" per row
n_panels <- length(cytokines)
n_rows   <- ceiling(n_panels / 2)
fig_h    <- n_rows * 2.4

# Build facet labels from stats_plt (format to taste)
lab_tbl <- stats_plt %>%
  filter(Assay %in% unique(df_plot$Assay)) %>%
  mutate(
    Assay_label = sprintf(
      "%s  (mean=%.2f, sd=%.2f, SD/Mean=%.2f)",
      Assay, mean_conc, sd_conc, sd_over_mean
    )
  ) %>%
  select(Assay, Assay_label)

# Join labels to your plotting data
df_plot_labeled <- df_plot %>%
  left_join(lab_tbl, by = "Assay") %>%
  # keep panels in your chosen order
  mutate(Assay_label = factor(Assay_label, levels = lab_tbl$Assay_label[match(levels(Assay), lab_tbl$Assay)]))

# Plot with the label in facet strips
p <- ggplot(df_plot_labeled, aes(x = day_num, y = median_conc, group = Assay, color = Assay)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Assay_label, ncol = 2, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +
  labs(title = "PLATCOV Cytokines", x = "Day", y = "Median concentration (pg/ml)") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p)

# compute the complementary set (others not in cytokines)
others <- sort(as.character(setdiff(all_cytokines, cytokines)))

# prepare data
df_other <- median_across_patients_plt %>% filter(Assay %in% others)

if (nrow(df_other) == 0) {
  message("No assays remain outside the selected cytokines list.")
} else {
  # dynamic figure height: ~2.4\" per row, 2 columns
  n_panels <- length(others)
  n_rows   <- ceiling(n_panels / 2)
  fig_h    <- n_rows * 2.4
  
  # Build facet labels from stats_plt (format to taste)
  lab_tbl <- stats_plt %>%
    filter(Assay %in% unique(df_other$Assay)) %>%
    mutate(
      Assay_label = sprintf(
        "%s  (mean=%.2f, sd=%.2f, SD/Mean=%.2f)",
        Assay, mean_conc, sd_conc, sd_over_mean
      )
    ) %>%
    select(Assay, Assay_label)
  
  # Join labels to your plotting data
  df_other_labeled <- df_other %>%
    left_join(lab_tbl, by = "Assay") %>%
    # keep panels in your chosen order
    mutate(Assay_label = factor(Assay_label, levels = lab_tbl$Assay_label[match(levels(Assay), lab_tbl$Assay)]))
  
  p_other <- ggplot(df_other_labeled, aes(x = day_num, y = median_conc, group = Assay, color = Assay)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~ Assay_label, ncol = 2, scales = "free_y") +
    scale_x_discrete(drop = FALSE) +
    labs(title = "PLATCOV Cytokines (Other)", x = "Day", y = "Median concentration (pg/ml)") +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )

  print(p_other)
}


### AD ASTRA
# Make sure day_num is a factor so axis shows exactly 0, 3, 7, 14
median_across_patients_ast$day_num <- factor(
  median_across_patients_ast$day_num,
  levels = c(0, 3, 7, 14)
)

stats_ast <- median_across_patients_ast %>%
  group_by(Assay) %>%
  summarise(
    n               = sum(!is.na(median_conc)),
    mean_conc       = mean(median_conc, na.rm = TRUE),
    sd_conc         = sd(median_conc,   na.rm = TRUE),
    sd_over_mean    = ifelse(mean_conc == 0, NA_real_, sd_conc / mean_conc),
    .groups = "drop"
  ) %>%
  arrange(desc(sd_over_mean))

# filter to selected cytokines only
df_plot2 <- median_across_patients_ast %>% filter(Assay %in% cytokines)

# compute dynamic figure height (inches): ~2.4" per row
n_panels <- length(cytokines)
n_rows   <- ceiling(n_panels / 2)
fig_h    <- n_rows * 2.4

# Build facet labels from stats_plt (format to taste)
lab_tbl2 <- stats_ast %>%
  filter(Assay %in% unique(df_plot2$Assay)) %>%
  mutate(
    Assay_label = sprintf(
      "%s  (mean=%.2f, sd=%.2f, SD/Mean=%.2f)",
      Assay, mean_conc, sd_conc, sd_over_mean
    )
  ) %>%
  select(Assay, Assay_label)

# Join labels to your plotting data
df_plot2_labeled <- df_plot2 %>%
  left_join(lab_tbl2, by = "Assay")

p2 <- ggplot(df_plot2_labeled, aes(x = day_num, y = median_conc, group = Assay, color = Assay)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Assay_label, ncol = 2, scales = "free_y") +  # free y per subfigure
  scale_x_discrete(drop = FALSE) +
  labs(title = "AD ASTRA Cytokines", x = "Day", y = "Median concentration (pg/ml)") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p2)

# compute the complementary set (others not in cytokines)
others2 <- sort(as.character(setdiff(all_cytokines, cytokines)))

# prepare data
df_other2 <- median_across_patients_ast %>% filter(Assay %in% others2)# %>%

if (nrow(df_other2) == 0) {
  message("No assays remain outside the selected cytokines list.")
} else {
  # dynamic figure height: ~2.4\" per row, 2 columns
  n_panels <- length(others2)
  n_rows   <- ceiling(n_panels / 2)
  fig_h    <- n_rows * 2.4
  
  # Build facet labels from stats_ast (format to taste)
  lab_tbl2 <- stats_ast %>%
    filter(Assay %in% unique(df_other2$Assay)) %>%
    mutate(
      Assay_label = sprintf(
        "%s  (mean=%.2f, sd=%.2f, SD/Mean=%.2f)",
        Assay, mean_conc, sd_conc, sd_over_mean
      )
    ) %>%
    select(Assay, Assay_label)
  
  # Join labels to your plotting data
  df_other2_labeled <- df_other2 %>%
    left_join(lab_tbl2, by = "Assay")
  
  p2_other <- ggplot(df_other2_labeled, aes(x = day_num, y = median_conc, group = Assay, color = Assay)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~ Assay_label, ncol = 2, scales = "free_y") +
    scale_x_discrete(drop = FALSE) +
    labs(title = "AD ASTRA Cytokines (Other)", x = "Day", y = "Median concentration (pg/ml)") +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )
  
  print(p2_other)
}


# Figure 2: Plot Cytokine Concentration Box-and-Whisker across D0, D3, D7, D14

# 1) Filter to IP-10 and compute log10 concentration; average duplicates per (ID, Timepoint_ID)
ip10 <- dat_all_plt %>%
  filter(Assay == "IP-10") %>%
  mutate(log10_conc = ifelse(`Calc. Concentration` > 0, log10(`Calc. Concentration`), NA_real_)) %>%
  group_by(ID, day_num) %>%
  summarise(log10_conc = mean(log10_conc, na.rm = TRUE), .groups = "drop")

# 2) Keep only patients with all four timepoints (0, 3, 7, 14)
ip10_ok <- ip10 %>%
  group_by(ID) %>%
  filter(all(c(0, 3, 7, 14) %in% day_num)) %>%
  ungroup()

# 3) Plot box-and-whisker along with each patient's trajectory
ggplot() +
  # Bold box & whisker per day (use all patients)
  geom_boxplot(data = ip10, aes(x = factor(day_num), y = log10_conc, group = factor(day_num)),
               width = 0.65, fill = "yellow", alpha = 0.9, outlier.alpha = 0.4) +
  geom_line(data = ip10_ok,
            aes(x = factor(day_num), y = log10_conc, group = ID),
            alpha = 0.05, linewidth = 0.3, color = "black") +
  # Faint individual points (only complete cases)
  geom_point(data = ip10_ok, aes(x = factor(day_num), y = log10_conc, group = ID),
             alpha = 0.3, size = 1.2) +
  scale_x_discrete(breaks = c("0", "3", "7", "14")) +
  labs(title = "PLATCOV — IP-10 per-day distribution (boxes) with patient data (points)",
       x = "Day", y = "log10(Calc. Concentration) (pg/ml)") +
  theme_minimal()

## PLATCOV
# Full
days <- c(0, 3, 7, 14)

# 1) Keep only assays you care about; compute log10 and average duplicates per (Assay, ID, day)
all_plt_log0 <- dat_all_plt %>%
  filter(Assay %in% cytokines) %>%
  mutate(log10_conc = ifelse(`Calc. Concentration` > 0, log10(`Calc. Concentration`), NA_real_)) %>%
  group_by(Assay, ID, day_num) %>%
  summarise(log10_conc = mean(log10_conc, na.rm = TRUE), .groups = "drop")

# 2) Mark complete-case patients (have all 0,3,7,14) *per assay*
complete_ids <- all_plt_log0 %>%
  group_by(Assay, ID) %>%
  summarise(has_all = all(days %in% day_num), .groups = "drop") %>%
  filter(has_all) %>%
  select(Assay, ID)

all_plt_log_ok <- all_plt_log0 %>% inner_join(complete_ids, by = c("Assay","ID"))

# 3) Order facets as in your `cytokines` vector + fix x-axis labels
all_plt_log  <- all_plt_log0  %>%
  mutate(
    day_num = factor(day_num, levels = days),
    Assay   = factor(Assay, levels = cytokines)
  )

all_plt_log_ok <- all_plt_log_ok %>%
  mutate(
    day_num = factor(day_num, levels = days),
    Assay   = factor(Assay, levels = cytokines)
  )

# 4) Dynamic height (inches)
n_panels <- length(cytokines)
n_rows   <- ceiling(n_panels / 2)
fig_h    <- n_rows * 2.4  # tweak 2.4 if you want tighter/looser panels

# 5) Plot: boxplots per day (all pts) + faint trajectories for complete pts; 2 columns; free y per panel
p3 <- ggplot() +
  # per-day distribution (all patients)
  geom_boxplot(
    data = all_plt_log,
    aes(x = day_num, y = log10_conc, group = day_num),
    width = 0.65, alpha = 0.9, outlier.alpha = 0.4
  ) +
  geom_point(
    data = all_plt_log_ok,
    aes(x = day_num, y = log10_conc, group = ID),
    alpha = 0.4, size = 0.9
  ) +
  facet_wrap(~ Assay, ncol = 2, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "PLATCOV — Per-day distributions (boxes) with complete-patient trajectories (lines/points)",
    x = "Day", y = "log10(Calc. Cytokine Concentration) (pg/ml)"
  ) +
  theme_minimal()

print(p3)


# 1) Keep only assays you care about; compute log10 and average duplicates per (Assay, ID, day)
all_ast_log0 <- dat_all_ast %>%
  filter(Assay %in% cytokines) %>%
  mutate(log10_conc = ifelse(`Calc. Concentration` > 0, log10(`Calc. Concentration`), NA_real_)) %>%
  group_by(Assay, ID, day_num) %>%
  summarise(log10_conc = mean(log10_conc, na.rm = TRUE), .groups = "drop")

# 2) Mark complete-case patients (have all 0,3,7,14) *per assay*
complete_ids <- all_ast_log0 %>%
  group_by(Assay, ID) %>%
  summarise(has_all = all(days %in% day_num), .groups = "drop") %>%
  filter(has_all) %>%
  select(Assay, ID)

all_ast_log_ok <- all_ast_log0 %>% inner_join(complete_ids, by = c("Assay","ID"))

# 3) Order facets as in your `cytokines` vector + fix x-axis labels
all_ast_log  <- all_ast_log0  %>%
  mutate(
    day_num = factor(day_num, levels = days),
    Assay   = factor(Assay, levels = cytokines)
  )

all_ast_log_ok <- all_ast_log_ok %>%
  mutate(
    day_num = factor(day_num, levels = days),
    Assay   = factor(Assay, levels = cytokines)
  )

# 4) Dynamic height (inches)
n_panels <- length(cytokines)
n_rows   <- ceiling(n_panels / 2)
fig_h    <- n_rows * 2.4  # tweak 2.4 if you want tighter/looser panels

# 5) Plot: boxplots per day (all pts) + faint trajectories for complete pts; 2 columns; free y per panel
p4 <- ggplot() +
  # per-day distribution (all patients)
  geom_boxplot(
    data = all_ast_log,
    aes(x = day_num, y = log10_conc, group = day_num),
    width = 0.65, alpha = 0.9, outlier.alpha = 0.4
  ) +
  geom_point(
    data = all_ast_log_ok,
    aes(x = day_num, y = log10_conc, group = ID),
    alpha = 0.4, size = 0.9
  ) +
  facet_wrap(~ Assay, ncol = 2, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "AD ASTRA — Per-day distributions (boxes) with complete-patient trajectories (lines/points)",
    x = "Day", y = "log10(Calc. Cytokine Concentration) (pg/ml)"
  ) +
  theme_minimal()

print(p4)


# Figure 3: Baseline Viral Load vs Cytokine Concentration (log-log plot)

### PLATCOV and AD ASTRA
# --- 1) Filter to baseline (TP0) and assays of interest
dat0_plt <- pcr_cytokine_dat_plt_filtered %>%
  filter(Timepoint_ID == 0, Assay %in% all_cytokines)
dat3_plt <- pcr_cytokine_dat_plt_filtered %>%
  filter(Timepoint_ID == 3, Assay %in% all_cytokines)
dat0_ast <- pcr_cytokine_dat_ast_filtered %>%
  filter(Timepoint_ID == 0, Assay %in% all_cytokines)
dat3_ast <- pcr_cytokine_dat_ast_filtered %>%
  filter(Timepoint_ID == 3, Assay %in% all_cytokines)

# --- 2) Keep one unique read per (ID, Assay), prioritize complete rows; add log10 conc
dat0_plt_unique <- dat0_plt %>%
  mutate(
    has_both   = !is.na(log10_viral_load) & !is.na(`Calc. Concentration`),
    log10_conc = ifelse(`Calc. Concentration` > 0, log10(`Calc. Concentration`), NA_real_)
  ) %>%
  arrange(desc(has_both)) %>%
  distinct(ID, Assay, .keep_all = TRUE) %>%
  select(-has_both)

dat0_ast_unique <- dat0_ast %>%
  mutate(
    has_both   = !is.na(log10_viral_load) & !is.na(`Calc. Concentration`),
    log10_conc = ifelse(`Calc. Concentration` > 0, log10(`Calc. Concentration`), NA_real_)
  ) %>%
  arrange(desc(has_both)) %>%
  distinct(ID, Assay, .keep_all = TRUE) %>%
  select(-has_both)

# --- 3) Split the assays into groups of 4 (will make 2x2 panels per figure)
assay_groups <- split(all_cytokines, ceiling(seq_along(all_cytokines) / 4))

# --- 4) Plot each group as a 2x2 figure (5 figures total if 20 assays)
plot_list <- vector("list", length(assay_groups))
for (i in seq_along(assay_groups)) {
  grp <- assay_groups[[i]]
  df_plot <- dat0_plt_unique %>% filter(Assay %in% grp)
  
  # --- compute Spearman correlation per Assay ---
  corr_tbl <- df_plot %>%
    group_by(Assay) %>%
    summarise(
      corr = cor(log10_conc, log10_viral_load,
                 method = "spearman", use = "complete.obs"),
      x = max(log10_conc, na.rm = TRUE),
      y = max(log10_viral_load, na.rm = TRUE),
      .groups = "drop"
    )
  
  # format labels (rounded to 2 decimals)
  corr_tbl <- corr_tbl %>%
    mutate(label = paste0("corr = ", round(corr, 2)))
  
  # --- plot ---
  p <- ggplot(df_plot, aes(x = log10_conc, y = log10_viral_load)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE,
                level = 0.95,
                fill = "grey60",
                alpha = 0.25,
                linewidth = 0.8) +
    facet_wrap(~ Assay, scales = "free", ncol = 2) +
    geom_text(data = corr_tbl,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
    labs(
      title = paste0("PLATCOV: Cytokine Level vs Baseline Viral Load — Set ", i),
      x = "log10(Calc. Cytokine Concentration) (pg/ml)",
      y = "log10(Viral Load) (pg/ml)"
    ) +
    theme_minimal()
  
  print(p)
  plot_list[[i]] <- p
}
# Very low correlation for baseline viral load VS cytokine concentration in D0
# For PLATCOV, highest in IFN-g (0.21), IL-6 (0.18), IFN-a2a (0.31), IL-1RA (0.31), IP-10 (0.27)

# --- 4) Plot each group as a 2x2 figure (5 figures total if 20 assays)
plot_list <- vector("list", length(assay_groups))
for (i in seq_along(assay_groups)) {
  grp <- assay_groups[[i]]
  df_plot <- dat0_ast_unique %>% filter(Assay %in% grp)
  
  # --- compute Spearman correlation per Assay ---
  corr_tbl <- df_plot %>%
    group_by(Assay) %>%
    summarise(
      corr = cor(log10_conc, log10_viral_load,
                 method = "spearman", use = "complete.obs"),
      x = max(log10_conc, na.rm = TRUE),
      y = max(log10_viral_load, na.rm = TRUE),
      .groups = "drop"
    )
  
  # format labels (rounded to 2 decimals)
  corr_tbl <- corr_tbl %>%
    mutate(label = paste0("corr = ", round(corr, 2)))
  
  # --- plot ---
  p <- ggplot(df_plot, aes(x = log10_conc, y = log10_viral_load)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE,
                level = 0.95,
                fill = "grey60",
                alpha = 0.25,
                linewidth = 0.8) +
    facet_wrap(~ Assay, scales = "free", ncol = 2) +
    geom_text(data = corr_tbl,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
    labs(
      title = paste0("AD ASTRA: Cytokine Level vs Baseline Viral Load — Set ", i),
      x = "log10(Calc. Cytokine Concentration) (pg/ml)",
      y = "log10(Viral Load) (pg/ml)"
    ) +
    theme_minimal()
  
  print(p)
  plot_list[[i]] <- p
}
# For AD ASTRA, correlation is much less pronounced

# Focus on these fewer cytokines which show high correlation
cytokines2 = c("IFN-g", "IL-6", "IFN-a2a", "IL-1RA", "IP-10")

## PLATCOV
df_plot <- dat0_plt_unique %>% filter(Assay %in% cytokines2)

# --- compute Spearman correlation per Assay ---
corr_tbl <- df_plot %>%
  group_by(Assay) %>%
  summarise(
    corr = cor(log10_conc, log10_viral_load,
               method = "spearman", use = "complete.obs"),
    x = max(log10_conc, na.rm = TRUE),
    y = max(log10_viral_load, na.rm = TRUE),
    .groups = "drop"
  )

# format labels (rounded to 2 decimals)
corr_tbl <- corr_tbl %>%
  mutate(label = paste0("corr = ", round(corr, 2)))

# --- plot ---
p <- ggplot(df_plot, aes(x = log10_conc, y = log10_viral_load)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              level = 0.95,
              fill = "grey60",
              alpha = 0.25,
              linewidth = 0.8) +
  facet_wrap(~ Assay, scales = "free", ncol = 2) +
  geom_text(data = corr_tbl,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
  labs(
    title = paste0("PLATCOV: Cytokine Level vs Baseline Viral Load"),
    x = "log10(Calc. Cytokine Concentration) (pg/ml)",
    y = "log10(Viral Load) (pg/ml)"
  ) +
  theme_minimal()

print(p)


df_plot <- dat0_ast_unique %>% filter(Assay %in% cytokines2)

# --- compute Spearman correlation per Assay ---
corr_tbl <- df_plot %>%
  group_by(Assay) %>%
  summarise(
    corr = cor(log10_conc, log10_viral_load,
               method = "spearman", use = "complete.obs"),
    x = max(log10_conc, na.rm = TRUE),
    y = max(log10_viral_load, na.rm = TRUE),
    .groups = "drop"
  )

# format labels (rounded to 2 decimals)
corr_tbl <- corr_tbl %>%
  mutate(label = paste0("corr = ", round(corr, 2)))

# --- plot ---
p <- ggplot(df_plot, aes(x = log10_conc, y = log10_viral_load)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              level = 0.95,
              fill = "grey60",
              alpha = 0.25,
              linewidth = 0.8) +
  facet_wrap(~ Assay, scales = "free", ncol = 2) +
  geom_text(data = corr_tbl,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
  labs(
    title = paste0("AD ASTRA: Cytokine Level vs Baseline Viral Load"),
    x = "log10(Calc. Cytokine Concentration) (pg/ml)",
    y = "log10(Viral Load) (pg/ml)"
  ) +
  theme_minimal()

print(p)


# --- Helper to plot by an arbitrary grouping column (as string) ---
plot_by_group <- function(dat, assay_groups, group_col, title_prefix = "", pal = NULL) {
  stopifnot(group_col %in% names(dat))
  plot_list <- vector("list", length(assay_groups))
  
  for (i in seq_along(assay_groups)) {
    grp <- assay_groups[[i]]
    df_plot <- dat %>%
      filter(Assay %in% grp)
    
    # ensure factor for clean legend/lines
    df_plot[[group_col]] <- factor(df_plot[[group_col]])
    
    # Spearman per Assay × group
    corr_tbl <- df_plot %>%
      group_by(Assay, .data[[group_col]]) %>%
      summarise(
        corr = suppressWarnings(cor(log10_conc, log10_viral_load,
                                    method = "spearman", use = "complete.obs")),
        .groups = "drop"
      )
    
    # per-facet ranges for placing labels
    ranges <- df_plot %>%
      group_by(Assay) %>%
      summarise(
        xmax = max(log10_conc, na.rm = TRUE),
        xmin = min(log10_conc, na.rm = TRUE),
        ymax = max(log10_viral_load, na.rm = TRUE),
        ymin = min(log10_viral_load, na.rm = TRUE),
        .groups = "drop"
      )
    
    # stacked labels (top-right) per group within each Assay facet
    corr_lbl <- corr_tbl %>%
      left_join(ranges, by = "Assay") %>%
      group_by(Assay) %>%
      arrange(.data[[group_col]], .by_group = TRUE) %>%
      mutate(
        label = paste0(as.character(.data[[group_col]]), ": \u03C1 = ", sprintf("%.2f", corr)),
        y = ymax - (row_number() - 1) * 0.08 * (ymax - ymin),
        x = xmax
      ) %>%
      ungroup()
    
    p <- ggplot(df_plot, aes(x = log10_conc, y = log10_viral_load)) +
      geom_point(aes(color = .data[[group_col]]), alpha = 0.7) +
      geom_smooth(
        aes(color = .data[[group_col]]),
        method = "lm", se = TRUE, level = 0.95, linewidth = 0.8, alpha = 0.2
      ) +
      facet_wrap(~ Assay, scales = "free", ncol = 2) +
      geom_text(
        data = corr_lbl,
        aes(x = x, y = y, label = label, color = .data[[group_col]]),
        inherit.aes = FALSE, hjust = 1.02, vjust = 1.1, size = 3.5, fontface = "bold"
      ) +
      labs(
        title = paste0(title_prefix, i),
        x = "log10(Calc. Cytokine Concentration) (pg/ml)",
        y = "log10(Viral Load) (pg/ml)",
        color = group_col
      ) +
      theme_minimal()
    
    if (!is.null(pal)) {
      p <- p + scale_color_manual(values = pal, name = group_col)
    }
    
    print(p)
    plot_list[[i]] <- p
  }
  
  invisible(plot_list)
}

# ---------- USE IT ----------

# 1) Fever_Baseline as-is
plot_list_fever <- plot_by_group(
  dat  = dat0_plt_unique,
  assay_groups = assay_groups,
  group_col = "Fever_Baseline",
  title_prefix = "PLATCOV (TP0): VL vs log10(conc) by Fever — Set "
)

plot_list_fever <- plot_by_group(
  dat  = dat0_ast_unique,
  assay_groups = assay_groups,
  group_col = "Fever_Baseline",
  title_prefix = "AD ASTRA (TP0): VL vs log10(conc) by Fever — Set "
)

# 2) FluType (A vs B). Provide a palette if you like.
pal_flu <- c("A" = "#1f77b4", "B" = "#d62728")
plot_list_flu <- plot_by_group(
  dat  = dat0_ast_unique,
  assay_groups = assay_groups,
  group_col = "fluType",
  title_prefix = "AD ASTRA: Cytokine Level & Baseline VL by FluType — Set ",
  pal = pal_flu
)

plot_by_group(
  dat = dat0_ast_unique,
  assay_groups = list(cytokines2),   # or split list if you want multiple batches
  group_col = "fluType",
  title_prefix = "AD ASTRA: Cytokine Level vs Baseline Viral Load — ",
  pal = c("A" = "#E41A1C", "B" = "#377EB8", "Other" = "#4DAF4A")
)

#######################
# Figure 4: Box-and-whisker plot of all cytokines (D0 only)
ggplot(dat0_plt_unique, aes(x = Assay, y = log10_conc)) +
  geom_boxplot(outlier.shape = 21, fill = "lightgreen", alpha = 0.6) +
  labs(
    title = "PLATCOV: D0 Cytokine Concentrations",
    x = "Cytokine (Assay)",
    y = "log10(Concentration) (pg/ml)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(dat0_ast_unique, aes(x = Assay, y = log10_conc)) +
  geom_boxplot(outlier.shape = 21, fill = "steelblue", alpha = 0.6) +
  labs(
    title = "AD ASTRA: D0 Cytokine Concentrations",
    x = "Cytokine (Assay)",
    y = "log10(Concentration) (pg/ml)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


####################
# Figure 5: Box-and-whisker plot of all cytokines (D0 only)
dat_d0d3_plt <- pcr_cytokine_dat_plt_filtered %>%
  filter(Timepoint_ID %in% c(0, 3), Assay %in% all_cytokines) %>%
  mutate(
    log10_conc = ifelse(`Calc. Concentration` > 0,
                        log10(`Calc. Concentration`),
                        NA_real_)
  )

# 2) Keep only one read per (ID, Assay, Timepoint_ID) if duplicates exist
dat_d0d3_plt_unique <- dat_d0d3_plt %>%
  arrange(!is.na(log10_conc)) %>%
  distinct(ID, Assay, Timepoint_ID, .keep_all = TRUE)

# Reorder Assay by median log10_conc (descending)
dat_d0d3_plt_unique <- dat_d0d3_plt_unique %>%
  mutate(
    Assay = fct_reorder(Assay, log10_conc, .fun = median, na.rm = TRUE, .desc = TRUE)
  )

# 3) Box-and-whisker plot with side-by-side boxes for Day 0 vs Day 3
ggplot(dat_d0d3_plt_unique,
       aes(x = Assay, y = log10_conc, fill = factor(Timepoint_ID))) +
  geom_boxplot(width = 0.5, outlier.shape = 21, alpha = 0.6, position = position_dodge(width = 0.5)) +
  labs(
    title = "PLATCOV: Cytokine Concentrations at Day 0 vs Day 3",
    x = "Cytokine (Assay)",
    y = "log10(Calc. Concentration)",
    fill = "Day"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 1) AD ASTRA
dat_d0d3_ast <- pcr_cytokine_dat_ast_filtered %>%
  filter(Timepoint_ID %in% c(0, 3), Assay %in% all_cytokines) %>%
  mutate(
    log10_conc = ifelse(`Calc. Concentration` > 0,
                        log10(`Calc. Concentration`),
                        NA_real_)
  )

# 2) Keep only one read per (ID, Assay, Timepoint_ID) if duplicates exist
dat_d0d3_ast_unique <- dat_d0d3_ast %>%
  arrange(!is.na(log10_conc)) %>%
  distinct(ID, Assay, Timepoint_ID, .keep_all = TRUE)

# Reorder Assay by median log10_conc (descending)
dat_d0d3_ast_unique <- dat_d0d3_ast_unique %>%
  mutate(
    Assay = fct_reorder(Assay, log10_conc, .fun = median, na.rm = TRUE, .desc = TRUE)
  )

# 3) Box-and-whisker plot with side-by-side boxes for Day 0 vs Day 3
ggplot(dat_d0d3_ast_unique,
       aes(x = Assay, y = log10_conc, fill = factor(Timepoint_ID))) +
  geom_boxplot(width = 0.5, outlier.shape = 21, alpha = 0.6, position = position_dodge(width = 0.5)) +
  labs(
    title = "AD ASTRA: Cytokine Concentrations at Day 0 vs Day 3",
    x = "Cytokine (Assay)",
    y = "log10(Calc. Concentration)",
    fill = "Day"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#################################

# Fever
# 1) Filter to Day 0, cytokines of interest, compute log10 concentration
dat_d0_plt_fever <- pcr_cytokine_dat_plt_filtered %>%
  filter(Timepoint_ID == 0, Assay %in% all_cytokines) %>%
  mutate(
    log10_conc = ifelse(`Calc. Concentration` > 0,
                        log10(`Calc. Concentration`),
                        NA_real_),
    Fever_Baseline = factor(Fever_Baseline, levels = c(0, 1),
                            labels = c("No fever", "Fever"))
  ) %>%
  filter(!is.na(Fever_Baseline))

# 2) Keep one read per (ID, Assay) if duplicates exist
dat_d0_plt_fever_unique <- dat_d0_plt_fever %>%
  arrange(!is.na(log10_conc)) %>%
  distinct(ID, Assay, .keep_all = TRUE)

# 3) Box-and-whisker plot with side-by-side boxes for Fever vs No fever
dat_d0_plt_fever_unique <- dat_d0_plt_fever_unique %>%
  mutate(Assay = factor(Assay))  # keeps current order

assay_levels <- levels(dat_d0_plt_fever_unique$Assay)
assay_labels <- ifelse(
  assay_levels %in% cytokines,
  paste0("<span style='color:red;'>", assay_levels, "</span>"),
  assay_levels
)

ggplot(dat_d0_plt_fever_unique,
       aes(x = Assay, y = log10_conc, fill = Fever_Baseline)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.6,
               position = position_dodge(width = 0.8)) +
  scale_x_discrete(labels = setNames(assay_labels, assay_levels)) +
  labs(
    title = "PLATCOV: Cytokine Concentrations at Day 0 by Baseline Fever",
    x = "Cytokine (Assay)",
    y = "log10(Calc. Concentration)",
    fill = "Baseline fever"
  ) +
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1))


# 1) Filter to Day 0, cytokines of interest, compute log10 concentration
dat_d0_ast_fever <- pcr_cytokine_dat_ast_filtered %>%
  filter(Timepoint_ID == 0, Assay %in% all_cytokines) %>%
  mutate(
    log10_conc = ifelse(`Calc. Concentration` > 0,
                        log10(`Calc. Concentration`),
                        NA_real_),
    Fever_Baseline = factor(Fever_Baseline, levels = c(0, 1),
                            labels = c("No fever", "Fever"))
  ) %>%
  filter(!is.na(Fever_Baseline))

# 2) Keep one read per (ID, Assay) if duplicates exist
dat_d0_ast_fever_unique <- dat_d0_ast_fever %>%
  arrange(!is.na(log10_conc)) %>%
  distinct(ID, Assay, .keep_all = TRUE)

# 3) Box-and-whisker plot with side-by-side boxes for Fever vs No fever
dat_d0_ast_fever_unique <- dat_d0_ast_fever_unique %>%
  mutate(Assay = factor(Assay))  # keeps current order

assay_levels <- levels(dat_d0_ast_fever_unique$Assay)
assay_labels <- ifelse(
  assay_levels %in% cytokines,
  paste0("<span style='color:red;'>", assay_levels, "</span>"),
  assay_levels
)

ggplot(dat_d0_ast_fever_unique,
       aes(x = Assay, y = log10_conc, fill = Fever_Baseline)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.6,
               position = position_dodge(width = 0.8)) +
  scale_x_discrete(labels = setNames(assay_labels, assay_levels)) +
  labs(
    title = "AD ASTRA: Cytokine Concentrations at Day 0 by Baseline Fever",
    x = "Cytokine (Assay)",
    y = "log10(Calc. Concentration)",
    fill = "Baseline fever"
  ) +
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1))



# # 1) Filter to baseline and the assay of interest
# dat0 <- pcr_cytokine_dat_plt_filtered %>%
#   filter(Timepoint_ID == 0, Assay %in% c("IP-10"))
# 
# # 2) Keep one unique read per (ID, Assay), prioritizing complete rows
# dat0_unique <- dat0 %>%
#   mutate(
#     has_both = !is.na(log10_viral_load) & !is.na(`Calc. Concentration`),
#     # safe log10 transform; non-positive concentrations become NA
#     log10_conc = ifelse(`Calc. Concentration` > 0, log10(`Calc. Concentration`), NA_real_)
#   ) %>%
#   arrange(desc(has_both)) %>%           # complete rows first
#   distinct(ID, Assay, .keep_all = TRUE) %>%
#   select(-has_both)
# 
# # (Optional) sanity check
# dat0_unique %>% count(Assay) %>% arrange(desc(n)) %>% print()
# 
# # 3) Plot: y = log10_viral_load, x = log10(Calc. Concentration)
# ggplot(dat0_unique, aes(x = log10_conc, y = log10_viral_load)) +
#   geom_point(alpha = 0.7) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ Assay, scales = "free_x") +
#   labs(
#     title = "Baseline (Timepoint 0): Viral load vs log10 cytokine concentration",
#     x = "log10(Calc. Concentration)",
#     y = "log10(Viral Load)"
#   ) +
#   theme_minimal()


#############################
# Leukocytes
## PLOT LEUKOCYTES LEVEL ACROSS DIFFERENT CYTOKINE
View(pcr_cytokine_dat_ast)



########## Hierarchical Modeling --> Takes very long ##############
library(brms)
mod_naive = brm(log10_viral_load | cens(censor) ~ Time+CT_RNaseP+(1+Time|ID),
                data = pcr_dat,
                prior = c(prior(lkj(2), class = "cor"),
                          prior(normal(5, 5), class = Intercept),
                          prior(normal(0,1), class='b'),
                          prior(constant(7), class = 'nu')),
                family = 'student', cores = 4, chains = 4)

summary(mod_naive)
pcr_dat$preds = colMeans(posterior_epred(mod_naive,newdata = data.frame(Time=pcr_dat$Time, CT_RNaseP=mean(pcr_dat$CT_RNaseP), ID=pcr_dat$ID)))

pcr_dat_summary = pcr_dat %>% group_by(ID) %>%
  mutate(
    baseline_VL = mean(log10_viral_load[Timepoint_ID==0]),
    slope = coef(lm(preds~Time))['Time']) %>%
  distinct(ID, .keep_all = T)

hist(pcr_dat_summary$slope, breaks = 15)

dat_all_2 = merge(dat_all, pcr_dat_summary[, c('ID','slope','baseline_VL')], by = 'ID')
dat_all_2 = dat_all_2 %>% filter(Assay %in% cytokines)
dat_all_2 = dat_all_2 %>% group_by(Assay, day) %>%
  mutate(Conc_standardised = scale(log10(`Calc. Concentration pg/ml`)))

dat_all_2 %>% ggplot(aes(x=`Calc. Concentration pg/ml`, y = slope))+
  geom_point() + geom_smooth(method = lm) + scale_x_log10()+
  facet_wrap(vars(Assay, day), nrow = 4, ncol = 3)

p_all_vc=dat_all_2 %>% ggplot(aes(x=Conc_standardised, y = slope))+
  geom_point() + geom_smooth(method = lm) + theme_minimal() +
  xlab('Standardised concentration')+ylab('Viral clearance (log10 change per day)')+
  facet_wrap(vars(Assay, day), nrow = 4, ncol = 3)

p_all_vl=dat_all_2 %>% filter(day=='D0') %>%
  ggplot(aes(x=Conc_standardised, y = baseline_VL))+
  geom_point() + geom_smooth(method = lm) + theme_minimal() +
  xlab('Standardised concentration')+ylab('Baseline log10 viral load')+
  facet_wrap(vars(Assay), nrow = 4, ncol = 3)

dat_all_2 %>% filter(day=='D3') %>%
  ggplot(aes(x=Conc_standardised, y = baseline_VL))+
  geom_point() + geom_smooth(method = lm) + theme_minimal() +
  xlab('Standardised concentration')+ylab('Baseline log10 viral load')+
  facet_wrap(vars(Assay), nrow = 4, ncol = 3)

ggsave(filename =  '~/Downloads/first_100_baseline_VL.pdf', plot = p_all_vl)
ggsave(filename =  '~/Downloads/first_100_viral_clearance.pdf', plot = p_all_vc)


summary(lm(baseline_VL ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IFN-g')))
summary(lm(baseline_VL ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IL-1RA')))
summary(lm(baseline_VL ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IP-10')))

summary(lm(slope ~ log10(`Calc. Concentration pg/ml`), data = dat_all_2%>%filter(day=='D0',Assay=='IFN-g')))

summary(lm(slope ~ Conc_standardised, data = dat_all_2%>%filter(day=='D3',Assay=='IL-1RA')))
summary(lm(slope ~ log10(`Calc. Concentration pg/ml`), data = dat_all_2%>%filter(day=='D3',Assay=='IL-1RA')))

summary(lm(slope ~ Conc_standardised, data = dat_all_2%>%filter(day=='D0',Assay=='IP-10')))
summary(lm(slope ~ log10(`Calc. Concentration pg/ml`), data = dat_all_2%>%filter(day=='D0',Assay=='IP-10')))


xx=dat_all_2 %>% filter(day=='D0') %>%
  pivot_wider(names_from = Assay,
              values_from = Conc_standardised,
              id_cols = c(ID, slope, baseline_VL))
mod1 = lm(slope ~ . , data = xx[, !colnames(xx) %in% c('ID','baseline_VL')])
summary(mod1)

mod2 = lm(baseline_VL ~ . , data = xx[, !colnames(xx) %in% c('ID','slope')])
summary(mod2)
