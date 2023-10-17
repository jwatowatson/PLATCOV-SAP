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
ineffective_data = read.csv('../Analysis_Data/Ineffective_analysis.csv')

remdesivir_data <- mITT(remdesivir_data)
paxlovid_data <- mITT(paxlovid_data)
recent_paxlovid_data <- mITT(recent_paxlovid_data)
ineffective_data <- mITT(ineffective_data)

data_list <- list(remdesivir_data, paxlovid_data, recent_paxlovid_data, ineffective_data)
############################################################################
interventions_all = c("Remdesivir", "Molnupiravir",  "Nirmatrelvir + Ritonavir", "Nirmatrelvir + Ritonavir", "Nirmatrelvir + Ritonavir", "Ivermectin", "Favipiravir", "Fluoxetine")
ref_arms_all = c("No study drug",   "No study drug", "Molnupiravir",              "No study drug",            "No study drug", "No study drug", "No study drug", "No study drug")
data_ID <- c(1,2,2,2,3,4,4,4)
pairs <- paste(interventions_all,ref_arms_all, data_ID, sep = "_")

Dmax_all <- c(2:7, 14)

bootstrap_rep = 20

pair_ID <- 6:8 # For ineffective drugs (should be 1:5 for success drugs)

model_settings <-  expand.grid(mod = '../Stan_models/Linear_model_RNaseP.stan',
                               prior = 1,
                               cov_matrices = 1,
                               Dmax = Dmax_all,
                               pairs = pairs[pair_ID],
                               boot_rep = 1:bootstrap_rep)

model_settings <- model_settings %>%
  separate(pairs, c("intervention", "ref_arm", "data_ID"), "_")

model_settings$Niter = 4000
model_settings$Nwarmup = 2000
model_settings$Nthin = 8
model_settings$Nchain = 4

writeLines(sprintf('Numbers of posterior samples in total is %s',
                   unique((model_settings$Niter-model_settings$Nwarmup)*model_settings$Nchain/model_settings$Nthin)))
model_settings$data_ID <- as.numeric(model_settings$data_ID)

############################################################################
model_setup_f_name = "model_settings_ineffective.RData"# change here for ineffective drugs

save(model_settings, 
     data_list,
     all_priors,
     file = model_setup_f_name)
