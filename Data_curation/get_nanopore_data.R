get_nanopore_data = function(prefix_dropbox){
  fnames_var = list.files(paste0(prefix_dropbox,"/Data/Nanopore sequencing"),full.names = T,pattern = '.csv')
  res_nano = lapply(fnames_var, read.csv)
  for(i in 1:length(res_nano)) res_nano[[i]] = res_nano[[i]][, c('Sequence.name','Scorpio.call','Lineage')]
  res_nano = do.call(rbind, res_nano )
  res_nano$ID = sapply(res_nano$Sequence.name, FUN = function(x) strsplit(x,split = '_')[[1]][1])
  if(any(duplicated(res_nano$ID))) print('warning - some duplicates in the nanopore output - will be removed')
  res = res_nano[!duplicated(res_nano$ID), ]
  write.csv(x = res[, c('ID','Lineage')], file = paste0(prefix_analysis_data, "/Analysis_Data/lineages.csv"),row.names = F)
  
  ## run the python script to convert lineage names into a usable set
  
  res_update = read.csv(paste0(prefix_analysis_data, "/Analysis_Data/output.csv"))
  res_update = res_update[!duplicated(res_update$Original), ]
  print(unique(res_update$Original))
  res$Variant = plyr::mapvalues(res$Lineage, from = res_update$Original, to = res_update$VariantClass)
  res = res[,c('ID','Variant')]
  res
}
####################################################



