library(seqinr)
library(readxl)
library(stringr)
####################################################################################
source("user_settings.R")
####################################################################################
# Combine all FASTA files for Nextclade analysis
listfile <- list.files(paste0(prefix_dropbox, "/DATA/PLATCOV_FASTA/"), full.names = T,recursive = T)
listfile_fasta <- listfile[-grep(".xlsx", listfile)]

all_fasta <- do.call(c, lapply(listfile_fasta, read.fasta))

write.fasta(all_fasta, names = names(all_fasta), file = paste0(prefix_analysis_data, "/Analysis_Data/", "all_fasta.fasta"))
####################################################################################
names <-  listfile[grep(".xlsx", listfile)]
# Extract Sequence ID map
Sample_ID_map <- NULL
for(i in 1:length(names)){
  if(grepl("Brazil", names[i])){
    temp <- read_xlsx(names[i])  
    cols <- c("Subject ID", "sample")
    temp <- temp[, cols]
    colnames(temp) <- c("Patient_ID", "Sequence_ID")
    
    Sample_ID_map <- rbind(Sample_ID_map, temp)
  } else if (grepl("Laos", names[i])) {
    temp <- read_xlsx(names[i])  
    cols <- c("Patient ID", "GISAID name")
    temp <- temp[, cols]
    colnames(temp) <- c("Patient_ID", "Sequence_ID")
    #temp$Patient_ID <- sapply(str_split(str_trim(temp$Patient_ID, "both"), "-"), function(x) paste(x[1:3], collapse = "-"))
    #temp$Sequence_ID <- sapply(str_split(str_trim(temp$Sequence_ID , "both"), "/"), function(x) paste0(x[3]))
    
    Sample_ID_map <- rbind(Sample_ID_map, temp)
  }  else if (grepl("Thailand", names[i])) {
    temp <- read_xlsx(names[i])
    if("Patient ID" %in% temp[1,]){colnames(temp) <- temp[1,]; temp <- temp[-1,]}
    temp <- temp[, cols]
    colnames(temp) <- c("Patient_ID", "Sequence_ID")
    temp$Patient_ID <- sapply(str_split(str_trim(temp$Patient_ID, "both"), "-|_"), function(x) paste(x, collapse = "-"))
    
    Sample_ID_map <- rbind(Sample_ID_map, temp)
  }
}
Sample_ID_map

write.csv(Sample_ID_map, file = paste0(prefix_analysis_data, "/Analysis_Data/", "sequencing_ID_map.csv"), row.names = F)
####################################################################################
# Run Nextclade for extracting mutation
# Can't run on MORU internet
# Input: Combined Fasta file "all_fasta.fasta"; SARS-CoV-2 data (downloaded from Nextclade: https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/usage.html)
# Output: Interested only .tsv file
# Need Nextclade installed
##------------------------------------------------------------------------------------
re_download = T
if(re_download){
  arg_download <- "nextclade dataset get --name nextstrain/sars-cov-2/wuhan-hu-1/orfs --output-dir ../Analysis_Data/Nextclade/sars-cov-2"
  arg_download
  #for Windows
  #shell(arg_download)
  ### For Mac #### TO BE CHANGED #######
   system(arg_download) 
}
##------------------------------------------------------------------------------------
sars_cov_2_data <- "../Analysis_Data/Nextclade/sars-cov-2" # path to the data downloaded from Nextclade
output_folder <- "../Analysis_Data/Nextclade/output" # folder name for outputs
input_fasta <- "../Analysis_Data/all_fasta.fasta" # path to input fasta files

arguments <- paste0("nextclade run --input-dataset ", sars_cov_2_data, " --output-all ", output_folder, " ", input_fasta)
arguments
### For Windows
#shell(arguments) 

### For Mac #### TO BE CHANGED #######
 system(arguments) 
####################################################################################










