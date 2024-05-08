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
joined_data <- match_nextclade_ID(nextclade_file, naming_file)
aa_matrix <- Reformat_Mutations(joined_data)
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
# evusheld_muts <- evusheld_muts %>% mutate(
#   evusheld_test = case_when(
#     (`S:F486_any` == TRUE & `S:F486*_missing` == FALSE) &
#       ((`S:R346_any` == TRUE & `S:R346*_missing` == FALSE)|(`S:K444_any` == TRUE & `S:K444*_missing` == FALSE)) ~ "Resistant",
#     (`S:F486_any` == FALSE & `S:F486*_missing` == FALSE) ~ "Wild type",
#     ((`S:R346_any` == FALSE & `S:R346*_missing` == FALSE)|(`S:K444_any` == FALSE & `S:K444*_missing` == FALSE)) ~ "Wild type",
#     TRUE ~ NA
#     )
#   )

evusheld_muts <- evusheld_muts %>% mutate(
  test_F486_any = case_when(
    (`S:F486_any` == TRUE & `S:F486*_missing` == FALSE) ~ TRUE,
    (`S:F486_any` == FALSE & `S:F486*_missing` == FALSE) ~ FALSE,
    TRUE ~ NA
  ),
  test_R346_any = case_when(
    (`S:R346_any` == TRUE & `S:R346*_missing` == FALSE) ~ TRUE,
    (`S:R346_any` == FALSE & `S:R346*_missing` == FALSE) ~ FALSE,
    TRUE ~ NA
  ),
  test_K444_any = case_when(
    (`S:K444_any` == TRUE & `S:K444*_missing` == FALSE) ~ TRUE,
    (`S:K444_any` == FALSE & `S:K444*_missing` == FALSE) ~ FALSE,
    TRUE ~ NA
  ),
  test_evusheld = case_when(
    test_F486_any == FALSE & (test_R346_any == FALSE & test_K444_any == FALSE) ~ "Wild type",
    test_F486_any == TRUE & (test_R346_any == FALSE & test_K444_any == FALSE) ~ "F846*",
    test_F486_any == FALSE & (test_R346_any == TRUE | test_K444_any == TRUE) ~ "R346* or K444*",
    test_F486_any == TRUE & (test_R346_any == TRUE | test_K444_any == TRUE) ~ "F846* and (R346* or K444*)",
    TRUE ~ NA
  )
)
evusheld_muts
#############################################################
# Resistance
# regeneron_muts <- regeneron_muts %>% mutate(
#   regeneron_test = case_when(
#     (`S:G446S` == TRUE & `S:G446S_missing` == FALSE) ~ "Resistant",
#     (`S:K444_any` == TRUE & `S:K444*_missing` == FALSE) ~ "Resistant",
#     (`S:G446S` == FALSE & `S:G446S_missing` == FALSE) & (`S:K444_any` == FALSE & `S:K444*_missing` == FALSE) ~ "Wildtype",
#     TRUE ~ NA
#   )
# )

regeneron_muts <- regeneron_muts %>% mutate(
  test_G446S = case_when(
    (`S:G446S` == TRUE & `S:G446S_missing` == FALSE) ~ TRUE,
    (`S:G446S` == FALSE & `S:G446S_missing` == FALSE) ~ FALSE,
    TRUE ~ NA
  ),
  test_K444_any = case_when(
    (`S:K444_any` == TRUE & `S:K444*_missing` == FALSE) ~ TRUE,
    (`S:K444_any` == FALSE & `S:K444*_missing` == FALSE) ~ FALSE,
    TRUE ~ NA
  ),
  test_regeneron = case_when(
    test_G446S == TRUE & test_K444_any == FALSE ~ "G446S",
    test_G446S == FALSE & test_K444_any == TRUE ~ "K444*",
    test_G446S == TRUE & test_K444_any == TRUE ~ "G446S and K444*", 
    test_G446S == FALSE & test_K444_any == FALSE ~ "Wildtype",
    TRUE ~ NA
  )
)
regeneron_muts

# regeneron_muts <- regeneron_muts %>% mutate(
#   test_G446S = case_when(
#     (`S:G446S` == TRUE & `S:G446S_missing` == FALSE) ~ TRUE,
#     (`S:G446S` == FALSE & `S:G446S_missing` == FALSE) ~ FALSE,
#     TRUE ~ NA
#   ),
#   test_K444_any = case_when(
#     (`S:K444_any` == TRUE & `S:K444*_missing` == FALSE) ~ TRUE,
#     (`S:K444_any` == FALSE & `S:K444*_missing` == FALSE) ~ FALSE,
#     TRUE ~ NA
#   ),
#   regeneron_test_number = as.numeric(test_G446S) + as.numeric(test_K444_any),
#   regeneron_test = case_when(
#     regeneron_test_number == 2 ~ "Fully resistant",
#     regeneron_test_number == 1 ~ "Partially resistant",
#     regeneron_test_number == 0 ~ "Wild type",
#     TRUE ~ NA
#   )
# )
# regeneron_muts

#############################################################
#Export
write.csv(evusheld_muts, '../Analysis_Data/evusheld_mutations.csv', row.names = F)
write.csv(regeneron_muts, '../Analysis_Data/regeneron_mutations.csv', row.names = F)
#############################################################
joined_data_export <- joined_data %>% select(Patient_ID, Nextclade_pango)  

# Patient_ID, seqName, seqName2, Nextclade_pango, clade, partiallyAliased,
# clade_nextstrain, clade_who, clade_display, qc.overallScore, qc.overallStatus,
# coverage

joined_data_export$ID <- lapply(str_split(str_trim(joined_data_export$Patient_ID, "both"), "-"), function(x) paste(x[1], x[2], x[3], sep = "-")) %>% unlist()
colnames(joined_data_export)[2] <- "Lineage"
joined_data_export <- joined_data_export %>% select(ID, Lineage)

write.table(joined_data_export, file = "../Analysis_Data/lineages.csv", sep=",", row.names = F, col.names=FALSE)


