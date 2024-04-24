library(tidyverse)
library(stringr)
library(purrr)
library(readxl)
#############################################################
source("functions_query_mutations.R")
#############################################################
# Need to run "prep_fasta_files.R" before (with Nextclade installed). 
nextclade_file <- "../Analysis_Data/Nextclade/output/nextclade.tsv"
naming_file <- "../Analysis_Data/sequencing_ID_map.csv"
#############################################################
# Transforming data for calling mutations
aa_matrix <- Reformat_Mutations(nextclade_file, naming_file)
#############################################################
# Mutations of interest
# Evusheld
evusheld_mutations <- data.frame(
  "position" = c(346, 444, 486),
  "querystring" = c("S:R346*", "S:K444*", "S:F486*"),
  "nucl_positions" = c("22598,22599,22600", "22892,22893,22894", "23018,23019,23020")
)
# Regeneron
regeneron_mutations <- data.frame(
  "position" = c(446, 444),
  "querystring" = c("S:G446S", "S:K444*"),
  "nucl_positions" = c("22898,22899,22900", "22892,22893,22894")
)
# Ensitrelvir
ensitrelvir_mutations <- data.frame(
  "position" = c(21, 49, 49, 50, 144, 166, 166, 167, 173, 252, 304, 168),
  "querystring" = c("S:T21I", "S:M49I", "S:M49L", "S:L50F", "S:S144A", "S:E166A", "S:E166V", "S:L167F", "S:A173V", "S:P252L", "S:T304I", "S:P168-"),
  "nucl_positions" = c("22898,22899,22900", "22892,22893,22894")
)
#############################################################
# Calling mutations
evusheld_muts <- call_mutations(mutation_data = evusheld_mutations)
regeneron_muts <- call_mutations(mutation_data = regeneron_mutations)
#ensitrelvir_muts <- call_mutations(mutation_data = ensitrelvir_mutations)
#############################################################
# Infering resistance
# Evusheld
evusheld_muts <- evusheld_muts %>% mutate(
  evusheld_test = case_when(
    (`S:F486_any` == TRUE & `S:F486*_missing` == FALSE) & 
      ((`S:R346_any` == TRUE & `S:R346*_missing` == FALSE)|(`S:K444_any` == TRUE & `S:K444*_missing` == FALSE)) ~ TRUE,
    (`S:F486_any` == FALSE & `S:F486*_missing` == FALSE) ~ FALSE,
    ((`S:R346_any` == FALSE & `S:R346*_missing` == FALSE)|(`S:K444_any` == FALSE & `S:K444*_missing` == FALSE)) ~ FALSE,
    TRUE ~ NA
    )
  )
#############################################################
# Resistance
regeneron_muts <- regeneron_muts %>% mutate(
  regeneron_test = case_when(
    (`S:G446S` == TRUE & `S:G446S_missing` == FALSE) ~ TRUE,
    (`S:K444_any` == TRUE & `S:K444*_missing` == FALSE) ~ TRUE,
    (`S:G446S` == FALSE & `S:G446S_missing` == FALSE) & (`S:K444_any` == FALSE & `S:K444*_missing` == FALSE) ~ FALSE,
    TRUE ~ NA
  )
)
#############################################################
write.csv(evusheld_muts, '../Analysis_Data/evusheld_mutations.csv', row.names = F)
write.csv(regeneron_muts, '../Analysis_Data/regeneron_mutations.csv', row.names = F)
#############################################################