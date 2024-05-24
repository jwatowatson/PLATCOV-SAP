prefix_dat_cur <- "D:/PLATCOV-SAP/Data_curation/"
prefix_dropbox <- "C:/Users/Phrutsamon/Dropbox/PLATCOV_Analysis"

## Variant data
source(paste0(prefix_dat_cur, "get_nanopore_data.R"))
variant_data = get_nanopore_data(prefix_dropbox = prefix_dropbox, run_python = T) #MORU firewalls blocks this

write.csv(variant_data, paste0(prefix_dropbox, "/Data/variant_data.csv"), row.names = F)



