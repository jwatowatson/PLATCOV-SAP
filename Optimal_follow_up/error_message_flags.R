# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#######################################################################################
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)

#######################################################################################
load('model_settings.RData')
#######################################################################################
res = array(dim = c(nrow(model_settings), 2+7))

error_1 <- "Tail Effective Samples Size \\(ESS\\) is too low, indicating posterior variances and tail quantiles may be unreliable"
error_2 <- "There were [0-9]+ transitions after warmup that exceeded the maximum treedepth"
error_3 <- "Examine the pairs\\(\\) plot to diagnose sampling problems"
error_4 <- "Bulk Effective Samples Size \\(ESS\\) is too low, indicating posterior means and medians may be unreliable"
error_5 <- "There were [0-9]+ chains where the estimated Bayesian Fraction of Missing Information was low"
error_6 <- "The largest R-hat is [0-9]+\\.[0-9]+, indicating chains have not mixed"
error_7 <- "There were [0-9]+ divergent transitions after warmup"
#######################################################################################
for(i in 1:nrow(model_settings)){
  error_text <- paste(readLines(paste('~/Downloads/o_and_e_files/output.e24445707_',i,'.out',sep='')), collapse="\n")
  output_text <- paste(readLines(paste('~/Downloads/o_and_e_files/output.o24445707_',i,'.out',sep='')), collapse="\n")
  (output_text_v <- (str_split(output_text, "\n")[[1]]))
  (error_text_v <- (str_split(error_text, "\n")[[1]]))
  
  error_1_check <- as.numeric(any(grepl(error_1, error_text_v)))
  error_2_check <- as.numeric(any(grepl(error_2, error_text_v)))
  error_3_check <- as.numeric(any(grepl(error_3, error_text_v)))
  error_4_check <- as.numeric(any(grepl(error_4, error_text_v)))
  error_5_check <- as.numeric(any(grepl(error_5, error_text_v)))
  error_6_check <- as.numeric(any(grepl(error_6, error_text_v)))
  error_7_check <- as.numeric(any(grepl(error_7, error_text_v)))
  
  res[i, ] = c(length(output_text_v), length(error_text_v),
             error_1_check, error_2_check, error_3_check, error_4_check, 
             error_5_check, error_6_check, error_7_check)

}

res <- as.data.frame(res)
colnames(res) <- c("output_len", "error_len", 
                   "error_1", "error_2", "error_3", "error_4",
                   "error_5", "error_6", "error_7")
#######################################################################################
res_conditions <- data.frame(model_settings[,c("Dmax", "intervention", "ref_arm", "boot_rep", "data_ID")], res)

error_summary <- subset(res_conditions, select = -c(output_len, boot_rep, error_len)) %>% 
  group_by(Dmax, intervention, ref_arm, data_ID) %>% 
  summarise_each(list(sum)) %>% 
  as.data.frame() %>%
  arrange(data_ID, intervention, ref_arm, Dmax)
#######################################################################################
error_summary <- gather(error_summary, error, counts, error_1:error_7, factor_key=TRUE)
error_summary$Dmax <- as.factor(error_summary$Dmax)
error_summary$Dmax <- factor(error_summary$Dmax, levels = as.character(1:14))

error_summary$intervention[error_summary$intervention == "Nirmatrelvir + Ritonavir"] <- "Paxlovid"
error_summary$ref_arm[error_summary$ref_arm == "Nirmatrelvir + Ritonavir"] <- "Paxlovid"
error_summary$intervention[error_summary$data_ID == 3 & error_summary$intervention == "Paxlovid"] <- "Recent Paxlovid"
error_summary$pairs <- paste0(error_summary$intervention, " vs ", error_summary$ref_arm)
error_summary$pairs <- as.factor(error_summary$pairs)
error_summary$pairs <- factor(error_summary$pairs, levels = c("Remdesivir vs No study drug", "Molnupiravir vs No study drug",
                                                "Paxlovid vs Molnupiravir", "Paxlovid vs No study drug", "Recent Paxlovid vs No study drug"))
error_summary$color <- "black"
error_summary$color[error_summary$counts < 10] <- "white"

error_summary$error <- as.factor(error_summary$error)
levels(error_summary$error) <- c(paste("Error message", 1:7, sep = " "))
error_summary$error <- factor(error_summary$error, levels = rev(levels(error_summary$error)))
#######################################################################################
ggplot(error_summary, aes(x = Dmax, y = error, fill = counts)) +
  geom_tile(col = "white", linewidth = 0.25) +
  geom_text(aes(label = counts, col = color), size = 3) +
  facet_wrap(.~pairs) +
  theme_bw() + 
  scale_fill_viridis_c(guide = "none") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_color_manual(values = c("black","white"), guide = "none") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 14, face = "bold")) +
  xlab("Follow-up day") +
  ylab("Error counts of 20 bootstraps")
####################################################################################### 