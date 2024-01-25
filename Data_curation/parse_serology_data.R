library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

source('user_settings.R')
RERUN=F
if(RERUN){
  library(readxl)
  library(stringr)
  path = paste0(prefix_dropbox, '/Data/Serology/csv_files/450 nm/')
  
  f_meta = paste0(prefix_dropbox, '/Data/Serology/csv_files/PLATCOV sample list up to 18 Oct 23_updated10Nov23.xlsx')
  all_sheets = readxl::excel_sheets(f_meta)
  dat_list = list()
  for(my_sh in all_sheets){
    xx=unlist(readxl::read_excel(path = f_meta, sheet = my_sh, skip = 0,n_max = 1,col_names = F))
    
    plate_id = str_match(xx, "\\(\\s*(.*?)\\s*\\)")[2]
    dat = readxl::read_excel(path = f_meta, sheet = my_sh, skip=1)
    if(ncol(dat) == 13){
      dat = dat[, -13]
    }
    colnames(dat)= c("ID","SUBJECT NO","Time_ID","SAMPLE","COLLECTION","day",
                     "VOLUME","BARCODE","BOX/BAG","POSITION","Remark","Dilution")
    dat$ID = paste0('un', str_pad(dat$ID, width = 3, pad = '0',side = 'left'))
    
    dat$plate_ID = paste('plate', as.numeric(plate_id), sep = '_')
    dat_list[[which(my_sh==all_sheets)]]=dat
  }
  
  for(i in 1:length(dat_list)){
    print(c(i, nrow(dat_list[[1]]), ncol(dat_list[[i]])))
  }
  meta_dat_all = do.call("rbind", dat_list)
  unique(meta_dat_all$plate_ID)
  unique(meta_dat_all$Dilution)
  
  ##### PROBLEM IN EXCEL FORMATTING ###########
  meta_dat_all$Dilution = ifelse(meta_dat_all$Dilution %in% c('0.118055555555556','7.6388888888888895E-2'),
                                 '1:50', meta_dat_all$Dilution)
  table(meta_dat_all$Dilution)
  
  plate_data_list=list()
  f_names = list.files(path = path, pattern=".csv")
  for(i in 1:length(f_names)) {
    plate_ID = (strsplit(f_names[i], "\\_|\\."))[[1]][1]
    plate_data <- read_csv(paste(path, f_names[i], sep=""), col_names = TRUE, skip = 9)
    plate_data$plate_ID = paste('plate', 
                                as.numeric(gsub(pattern = 'PLATCOV Plate ',
                                                replacement = '',x = plate_ID)),
                                sep='_')
    plate_data_list[[i]] = plate_data
  }
  output_df <- do.call("rbind", plate_data_list)
  unique(output_df$plate_ID)
  
  
  
  std_data = output_df %>% filter(Type=='Standard') %>% arrange(plate_ID)
  std_data$plate_ID_int = as.numeric(as.factor(std_data$plate_ID))
  unique(std_data$Sample)
  CAL_VALS = c(200,200,150,150,100,100,80,80,40,40,20,20,10,10)
  names(CAL_VALS) = c("CAL10001","CAL10002","CAL20001","CAL20002","CAL30001","CAL30002",
                      "CAL40001","CAL40002","CAL50001","CAL50002","CAL60001","CAL60002",
                      "CAL70001","CAL70002")
  std_data$true_conc = CAL_VALS[std_data$Sample]
  std_data = std_data %>% 
    filter(! (Sample=='CAL70002' & plate_ID=='plate_34')) %>%
    arrange(plate_ID_int,true_conc)
  
  output_df_unk = output_df %>% filter(Type=='Unknown') %>% arrange(plate_ID)
  output_df_unk$Unique_sample_ID = apply(output_df_unk[, c('Sample','plate_ID')], 1, paste, collapse='_')
  meta_dat_all$Unique_sample_ID = apply(meta_dat_all[, c('ID','plate_ID')], 1, paste, collapse='_')
  output_df_unk$Unique_sample_ID=tolower(output_df_unk$Unique_sample_ID)
  output_df_unk$plate_ID_int = as.numeric(as.factor(output_df_unk$plate_ID))
  
  output_df_unk = merge(output_df_unk, 
                        meta_dat_all[, c('SUBJECT NO','Time_ID',
                                         'BARCODE','Dilution',
                                         'Unique_sample_ID')],
                        by = 'Unique_sample_ID', all = T)
  table(output_df_unk$Dilution)
  output_df_unk$ID_patient_timepoint = apply(output_df_unk[,c('SUBJECT NO','Time_ID')],1,paste,collapse='_')
  output_df_unk = output_df_unk %>%
    mutate(Dilution_factor = as.numeric(gsub(x=Dilution,pattern='1:',replacement=''))/50)
  output_df_unk$ID = as.numeric(as.factor(output_df_unk$ID_patient_timepoint))
  save(output_df_unk, std_data, file = '../Analysis_Data/serology_all.RDS')
  write_csv(output_df_unk, file = '../Analysis_Data/serology_unknowns.csv')
  write_csv(std_data, file = '../Analysis_Data/serology_controls.csv')
} else {
  load('../Analysis_Data/serology_all.RDS')
}

par(las=1, family='serif')
plot(jitter(std_data$true_conc), std_data$Abs,panel.first=grid(),
     xlab='True concentration', ylab='Absorbance at 450 nm',
     ylim = c(0,3))

# output_df_unk = output_df_unk %>% filter(ID %in% 1:2000)
## Run Bayesian model ##
mod = stan_model(file = 'serology_standard_curves.stan')
stan_data = data = list(N_controls=nrow(std_data),
                        K_plates=max(std_data$plate_ID_int),
                        conc_controls = std_data$true_conc/200,
                        absorb_controls=std_data$Abs,
                        ind_plate_controls=std_data$plate_ID_int,
                        N_sample = max(output_df_unk$ID),
                        N_obs = nrow(output_df_unk),
                        dilution_factor = 1/output_df_unk$Dilution_factor,
                        ind_sample = output_df_unk$ID,
                        absorb_obs = output_df_unk$Abs,
                        ind_plate_obs = output_df_unk$plate_ID_int)
sero_all_fit = sampling(mod, stan_data, iter=2000, chain=4, thin=1,
                        pars=c('L_Omega','g_controls','tau_controls',
                               'g_unk','tau_unk','theta_rand'), # we don't save this as it takes up lots of memory!
                        include=FALSE)
save(sero_all_fit, file = '../Rout/sero_fit.RData')

load(file = '../Rout/sero_fit.RData')

traceplot(sero_all_fit, pars = c('beta','sigmasq_u','sigma'))
summary(sero_all_fit, pars='beta')$summary


log_x_init = rstan::extract(sero_all_fit, pars='log_x_init')$log_x_init
log_IgG_ests = t(apply(log_x_init, 2, function(x) c(mean(x), sd(x))))
colnames(log_IgG_ests) = c('Mean_log_IgG', 'SE_IgG')
output_df_unk_unique = output_df_unk %>% arrange(ID) %>%
  distinct(ID, .keep_all = T)
output_df_unk_unique = cbind(output_df_unk_unique, log_IgG_ests)
hist(output_df_unk_unique$SE_IgG, breaks = 100)
hist(output_df_unk_unique$Mean_log_IgG, breaks = 100)

output_df_unk_unique$Day = ifelse(output_df_unk_unique$Time_ID=='D0H0','D0',output_df_unk_unique$Time_ID)
output_df_unk_unique$Day = as.numeric(gsub(x = output_df_unk_unique$Day,pattern = 'D',replacement = ''))
unique(output_df_unk_unique$Day)

sero_old = read_csv('../Analysis_Data/Serology_IgG.csv')
unique(sero_old$Day)

output_df_unk_unique$ID = output_df_unk_unique$`SUBJECT NO`
output_df_unk_unique = merge(output_df_unk_unique,
                             sero_old,
                             by=c('ID','Day'),
                             all.x = T)

plot((output_df_unk_unique$Mean_log_IgG+log(200)+log(50))/log(10),
     output_df_unk_unique$log10_IgG,
     xlab='Bayesian model mean estimate (log10)',
     ylab = 'Previous value - taken from thresholding (log10)',
     panel.first=grid())
lines(-4:20, -4:20, col='grey', lwd=3)

xx_low_low = output_df_unk_unique %>% filter(log10_IgG==0,Mean_log_IgG < -2)
xx_low_high = output_df_unk_unique %>% filter(log10_IgG==0,Mean_log_IgG > 0)

par(mfrow=c(1,2))
hist(output_df_unk_unique$SE_IgG, breaks=100, xlab='Standard error on log IgG estimate',
     main='')

plot((output_df_unk_unique$Mean_log_IgG+log(200)+log(50))/log(10),
     output_df_unk_unique$SE_IgG,
     panel.first=grid(), xlab='Bayesian model mean estimate (log10)',
     ylab = 'Standard error on log IgG estimate')
xx = output_df_unk_unique %>%
       filter(SE_IgG>0.4, Mean_log_IgG>2)
View(output_df_unk %>% filter(output_df_unk$ID_patient_timepoint %in% xx$ID_patient_timepoint))
write.csv(output_df_unk_unique, file = '../Analysis_Data/Serology_estimated.csv')
