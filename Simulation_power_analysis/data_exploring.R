library(ggplot2)
############################################################################
load("Rout/model_settings.RData")
############################################################################
i <- 1
print(model_settings[i,])
############################################################################
dataplot <- data_list[[model_settings$data_ID[i]]]
dataplot <- dataplot[dataplot$Timepoint_ID %in% 0:7,]
intervention <- model_settings$intervention[i]
ref_arm <- model_settings$ref_arm[i]

writeLines(paste('Number of all patients:', length(unique(dataplot$ID)), sep = " "))
writeLines(paste('Number of patients with', intervention,':', length(unique(dataplot$ID[dataplot$Trt == intervention])), sep = " "))
writeLines(paste('Number of patients with', ref_arm,':', length(unique(dataplot$ID[dataplot$Trt == ref_arm])), sep = " "))
############################################################################
dataplot$Timepoint_ID <- as.factor(dataplot$Timepoint_ID)
levels(dataplot$Timepoint_ID) <- c("Day 0", "Day 1", "Day 2", "Day 3", 
                                   "Day 4", "Day 5", "Day 6", "Day 7")   

dataplot2 <- aggregate(log10_viral_load~ID+Timepoint_ID+Trt+Site+BMI+Plate+Age+Sex+Symptom_onset, 
                       data = dataplot, FUN = median)

dataplot3<- aggregate(log10_viral_load~Timepoint_ID+Trt, data = dataplot, FUN = quantile, c(0.25, 0.5, 0.75))
dataplot3[,3:5] <- as.data.frame(as.matrix(dataplot3[,3]))
colnames(dataplot3)[3:5] <- c("Q1", "Q2", "Q3")

G2 <- ggplot() +
  geom_point(data = dataplot2, aes(x = Timepoint_ID, y = log10_viral_load, col = Trt), alpha = 0.4, size = 1.75, shape = 16) +
  geom_line(data = subset(dataplot2, as.numeric(Timepoint_ID) <= 8), aes(x =  Timepoint_ID, y = log10_viral_load, group = ID, col = Trt), 
            alpha = 0.5, linewidth = 0.5, linetype = 1) +
  scale_color_manual(values = c("#FA9494", "#87CEEB"), name = "") +
  ggnewscale::new_scale_color() +
  geom_point(data = dataplot3, aes(x = Timepoint_ID, y = Q2, col = Trt), size = 3, shape = 17) +
  geom_line(data = subset(dataplot3, as.numeric(Timepoint_ID) <= 8), aes(x =  Timepoint_ID, y = Q2, group = Trt, col = Trt), linewidth = 1, linetype = 1) +
  scale_color_manual(values = c("#B31312", "#332FD0"), name = "") +
  theme_bw() +
  scale_x_discrete(drop=F) +
  xlab("") +
  ylab("Log viral loads") + 
  theme(axis.title  = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  ggtitle(paste0(intervention, " vs ", ref_arm)) +
  annotate("text", x = 4.5, y = 8, label = paste(intervention,':', length(unique(dataplot$ID[dataplot$Trt == intervention])), "patients", sep = " "), hjust = 0) +
  annotate("text", x = 4.5, y = 7.5, label = paste(ref_arm,':', length(unique(dataplot$ID[dataplot$Trt == ref_arm])), "patients", sep = " "), hjust = 0)    

G2    
############################################################################
png(paste0("plots/observed_data_", gsub(" \\+.*", "", intervention),"_", gsub(" \\+.*", "", ref_arm), ".png"), width = 9, height = 6, units = "in", res = 350)
G2
dev.off()
