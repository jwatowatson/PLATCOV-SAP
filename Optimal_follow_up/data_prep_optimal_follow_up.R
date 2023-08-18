library(dplyr)
library(tidyr)
############################################################################

############################################################################
source('../priors.R')
source('../functions.R')
############################################################################
mITT <- function(platcov_dat){
  platcov_dat = platcov_dat %>% group_by(ID) %>%
    mutate(
      mITT = any(Per_protocol_sample==1 & Timepoint_ID>=3) & !all(CT_NS==40))
  return(as.data.frame(platcov_dat))
}

############################################################################
# Loading data
remdesivir_data = read.csv('../Analysis_Data/Remdesivir_analysis.csv')
paxlovid_data = read.csv( '../Analysis_Data/Paxlovid_Molnupiravir_analysis.csv')
recent_paxlovid_data = read.csv('../Analysis_Data/Paxlovid_recent_analysis.csv')

remdesivir_data <- mITT(remdesivir_data)
paxlovid_data <- mITT(paxlovid_data)
recent_paxlovid_data <- mITT(recent_paxlovid_data)

data_list <- list(remdesivir_data, paxlovid_data, recent_paxlovid_data)
############################################################################
interventions_all <- c("Remdesivir", "Molnupiravir",  "Molnupiravir", "Nirmatrelvir + Ritonavir", "Nirmatrelvir + Ritonavir")
ref_arms_all <- c("No study drug", "No study drug", "Nirmatrelvir + Ritonavir", "No study drug",  "No study drug")
data_ID <- c(1,2,2,2,3)
pairs <- paste(interventions_all,ref_arms_all, data_ID, sep = "_")

Dmax_all <- c(2:7, 14)

bootstrap_rep = 20
model_settings <-  expand.grid(mod = '../Stan_models/Linear_model_RNaseP.stan',
                               prior = 1,
                               cov_matrices = 1,
                               Dmax = Dmax_all,
                               pairs = pairs,
                               boot_rep = 1:bootstrap_rep)
model_settings <- model_settings %>%
  separate(pairs, c("intervention", "ref_arm", "data_ID"), "_")

model_settings$Niter = 6000
model_settings$Nwarmup = 3000
model_settings$Nthin = 10
model_settings$Nchain = 4

model_settings$data_ID <- as.numeric(model_settings$data_ID)

############################################################################
model_setup_f_name = "model_settings.RData"

save(model_settings, 
     data_list,
     all_priors,
     file = model_setup_f_name)

