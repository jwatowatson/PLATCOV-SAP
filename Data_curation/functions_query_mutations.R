Reformat_Mutations <- function(nextclade_file, naming_file) {
  #import the nextclade data
  nextclade <- read_tsv(nextclade_file)
  nextclade$seqName2 <- nextclade$seqName
  nextclade$seqName2[grep("/", nextclade$seqName2)] <- sapply(str_split(str_trim(nextclade$seqName2[grep("/", nextclade$seqName2)], "both"), "/"), function(x) paste0(x[3]))
  
  naming <- read.csv(naming_file)
  naming$Sequence_ID2 <- naming$Sequence_ID
  naming$Sequence_ID2[grep("/",  naming$Sequence_ID2)] <- sapply(str_split(str_trim(naming$Sequence_ID2[grep("/",  naming$Sequence_ID2)], "both"), "/"), function(x) paste0(x[3]))
  
  
  missing_IDs <- nextclade$seqName2[!nextclade$seqName2 %in% naming$Sequence_ID2]
  writeLines(sprintf('This sample lacks of Patient ID or Sequence ID %s', missing_IDs))
  
  for(i in 1:length(missing_IDs)){
    if(grepl("MTM",  missing_IDs[i])){naming <- rbind(naming, c(NA, missing_IDs[i], missing_IDs[i]))}
    if(grepl("PLT",  missing_IDs[i])){naming <- rbind(naming, c(missing_IDs[i], NA, missing_IDs[i]))}
  }
  
  joined <- merge(nextclade, naming, by.x = "seqName2", by.y = "Sequence_ID2")
  
  #select amino acid changes and missing nucleotide ranges only and sort names
  aachanges <- joined %>% select(c(Patient_ID, aaSubstitutions, aaDeletions, aaInsertions, missing))
  
  #split up substitutions, insertions, and deletions into separate columns in three separate operations
  #and change from a CSV of subs/dels/insertions into a wide form TRUE/FALSE for every mutation
  aasubs <- aachanges %>%
    mutate(items = str_split(aaSubstitutions, ",")) %>%
    unnest(cols = items, keep_empty = TRUE) %>%
    mutate(value = TRUE) %>%
    spread(key = items, value = value, fill = FALSE)
  
  aadeletions <- aachanges %>%
    mutate(items = str_split(aaDeletions, ",")) %>%
    unnest(cols = items, keep_empty = TRUE) %>%
    mutate(value = TRUE) %>%
    spread(key = items, value = value, fill = FALSE)
  
  aainsertions <- aachanges %>%
    mutate(items = str_split(aaInsertions, ",")) %>%
    unnest(cols = items, keep_empty = TRUE) %>%
    mutate(value = TRUE) %>%
    spread(key = items, value = value, fill = FALSE)
  
  #now merge those three tables back together, removing the redundant columns
  aa_matrix <- full_join(aasubs, aadeletions, keep=FALSE, by=join_by(Patient_ID, aaSubstitutions,    aaDeletions,aaInsertions,missing))
  aa_matrix <- full_join(aa_matrix, aainsertions, keep=FALSE, by=join_by(Patient_ID,aaSubstitutions,    aaDeletions,aaInsertions,missing))
  
  
  return(aa_matrix)
  
}

Query_Mutations(aa_matrix, querystring = "S:G446S", nucl_positions = "22898,22899,22900")



Query_Mutations<- function(aa_matrix, querystring, nucl_positions) {
  
  #first, add a column to determine whether the queried amino acid is missing based on the nucleotide missingness
  positions <- as.numeric(unlist(strsplit(nucl_positions, ',')))
  #probably a more elegant way to do this, but check if each of the 3 bases is in the missing range and return true/false
  
  aa_matrix <- aa_matrix %>%
    mutate( base1_missing = map_lgl(missing, ~ Check_Missing_Range(positions[1],.x)))
  aa_matrix <- aa_matrix %>%
    mutate( base2_missing = map_lgl(missing, ~ Check_Missing_Range(positions[2],.x)))
  aa_matrix <- aa_matrix %>%
    mutate( base3_missing = map_lgl(missing, ~ Check_Missing_Range(positions[3],.x)))
  
  #then merge into a single column which is TRUE if any of the three positions are missing
  aa_matrix <- aa_matrix %>%
    mutate(any_base_missing = rowSums(across(c(base1_missing,base2_missing, base3_missing)) > 0))
  
  #rename the column using the queryname
  colnames(aa_matrix)[colnames(aa_matrix) == "any_base_missing"] =paste0(querystring,"_missing") 
  
  #extract the mutation columns which match the query string
  #if there is a wildcard and we accept any mutation at the final position
  if (endsWith(querystring, "*")) {
    #knock off the final * for querying
    querystring_short <- str_sub(querystring, 1, -2)
    #select all matching columns
    aa_matches <- aa_matrix %>% select(c(Patient_ID,starts_with(querystring_short)))
    #merge column if any of the columns are true
    aa_matches <- aa_matches %>%
      mutate(any_true = rowSums(select(aa_matches, -c(paste0(querystring,"_missing"),Patient_ID))) > 0)
    #drop the columns for each individual mutation
    aa_matches <- aa_matches %>% select(c('Patient_ID', paste0(querystring,"_missing"),'any_true'))
    #rename the column using the queryname
    colnames(aa_matches)[colnames(aa_matches) == "any_true"] =paste0(querystring_short,"_any") 
  } else {
    aa_matches <- aa_matrix %>% select(c(Patient_ID,starts_with(querystring)))  
  }
  
  aa_matches[,grep("missing", colnames(aa_matches))] <- aa_matches[,grep("missing", colnames(aa_matches))] > 0
  
  
  return(aa_matches)
}

Check_Missing_Range <- function(number, range_string) {
  range_segments <- strsplit(range_string, ",")[[1]]
  
  within_range <- any(map_lgl(range_segments, function(segment) {
    if (grepl("-", segment)) {
      range_bounds <- strsplit(segment, "-")[[1]]
      lower_bound <- as.numeric(range_bounds[1])
      upper_bound <- as.numeric(range_bounds[2])
      number >= lower_bound && number <= upper_bound
    } else {
      as.numeric(segment) == number
    }
  }))
  return(within_range)
}

Check_Missing_Range <- function(number, range_string) {
  range_segments <- strsplit(range_string, ",")[[1]]
  
  within_range <- any(map_lgl(range_segments, function(segment) {
    if (grepl("-", segment)) {
      range_bounds <- strsplit(segment, "-")[[1]]
      lower_bound <- as.numeric(range_bounds[1])
      upper_bound <- as.numeric(range_bounds[2])
      number >= lower_bound && number <= upper_bound
    } else {
      as.numeric(segment) == number
    }
  }))
  return(within_range)
}


call_mutations <- function(mutation_data){
  mutation_list <- apply(mutation_data, 1, function(x) Query_Mutations(aa_matrix, querystring = x[2], nucl_positions = x[3]))
  mutation_summary <- mutation_list %>% reduce(left_join, by = "Patient_ID")
  missing_ind <- grep("missing", colnames(mutation_summary))
  mutation_summary[,missing_ind] <- mutation_summary[,missing_ind]>0
  return(mutation_summary)
}















