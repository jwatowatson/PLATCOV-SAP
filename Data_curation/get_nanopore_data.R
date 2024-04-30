get_nanopore_data = function(prefix_dropbox, run_python=F){
  # fnames_var = list.files(paste0(prefix_dropbox,"/Data/Nanopore sequencing"),full.names = T,pattern = '.csv')
  # res_nano = lapply(fnames_var, read.csv)
  # for(i in 1:length(res_nano)) res_nano[[i]] = res_nano[[i]][, c('Sequence.name','Scorpio.call','Lineage')]
  # res_nano = do.call(rbind, res_nano )
  # res_nano$ID = sapply(res_nano$Sequence.name, FUN = function(x) strsplit(x,split = '_')[[1]][1])
  # if(any(duplicated(res_nano$ID))) print('warning - some duplicates in the nanopore output - will be removed')
  # res = res_nano[!duplicated(res_nano$ID), ]
  # write.csv(x = res[, c('ID','Lineage')], file = paste0(prefix_analysis_data, "/Analysis_Data/lineages.csv"),row.names = F)

  if(run_python){  
  ##### run the python script to convert lineage names into a usable set
  ### For Windows
   shell(run_lineage_classifier) #run_lineage_classifier is defined in user_setting.R
  ### For Mac
  #system("python3 lineage_classifier.py --input ../Analysis_Data/lineages.csv --output ../Analysis_Data/newlineagelist.csv")
}
  res_update = read.csv(paste0(prefix_analysis_data, "/Analysis_Data/newlineagelist.csv"))
  #res_update <- res_update[-1,] # Remove header
  res_update = res_update[!duplicated(res_update$Original), ]
  write.csv(x = res_update, file = paste0(prefix_analysis_data, "/Analysis_Data/newlineagelist.csv"),row.names = F)
  print(unique(res_update$Original))
  
  res <- read.csv("../Analysis_Data/lineages.csv")
  colnames(res) <- c("ID", "Lineages")
  
  res$Variant = plyr::mapvalues(res$Lineage, from = res_update$Original, to = res_update$VariantClass)
  res = res[,c('ID','Variant')]
  res
}
####################################################


