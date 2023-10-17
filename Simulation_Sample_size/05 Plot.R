load('Rout/sequential_sim_params.RData')
######################################################################################################
library(ggplot2)
library(stringr)
library(ggpubr)
######################################################################################################
list_files <- list.files("Sim_out")

simResults <- do.call(rbind, lapply(paste0("Sim_out/",list_files), read.csv))
simResults$sim_id <-  as.numeric(str_extract(list_files, "[0-9]+"))

simResults <- simResults[order(simResults$sim_id),]
simResults$N_futility_success[is.na(simResults$N_futility_success)] <- max(simResults$N_futility_success, na.rm = T)
######################################################################################################
simResults$decision <- "Inconclusive"
simResults$decision[simResults$success == 1] <- "Success"
simResults$decision[simResults$futility == 1] <- "Futility"
simResults$decision <- as.factor(simResults$decision)

simResults$intervention_effect <- as.factor(simResults$intervention_effect)
levels(simResults$intervention_effect) <- c("trt_effects = 0.8", "trt_effects = 1.0", "trt_effects = 1.2",
                                            "trt_effects = 1.4", "trt_effects = 1.6")

# simResults$N_swabs_per_day <- as.factor(simResults$N_swabs_per_day )
# levels(simResults$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")

simResults$day_plans <- sim_settings$day_plans
######################################################################################################
res_plots <- list()

for(i in 1:length(levels(simResults$intervention_effect))){
  res_plots[[i]] <- ggplot(simResults[simResults$intervention_effect == levels(simResults$intervention_effect)[i],], aes(as.factor(N_swabs_per_day))) +
    geom_bar(aes(fill = decision), col = "black") +
    facet_grid(. ~ day_plans) +
    theme_bw() +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab("Number of swabs per day") +
    ylab("Counts") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          panel.spacing = unit(1, "lines")) +
    scale_fill_manual(values = c("#BB2525", "#F4D160", "#219C90"), drop = F, name = "Decisions") +
    ggtitle(levels(simResults$intervention_effect)[i])
} 
######################################################################################################
png(paste0("plots/decisions.png"), width = 15, height = 8, units = "in", res = 400)
ggarrange(plotlist = res_plots, nrow = 2, ncol = 3, labels = "AUTO", common.legend = T, legend = "bottom")
dev.off()
######################################################################################################
res_hist <- list()

for(i in 1:length(levels(simResults$intervention_effect))){
subdat <- simResults[simResults$intervention_effect == levels(simResults$intervention_effect)[i],]

subdat$N_swabs_per_day <- as.factor(subdat$N_swabs_per_day )
 levels(subdat$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")

res_hist[[i]] <-  ggplot(subdat, aes(x = N_futility_success, fill = decision)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_histogram(binwidth = 20) +
  facet_grid(N_swabs_per_day ~ day_plans) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.3, "lines")) +
  scale_fill_manual(values = c("#BB2525", "#F4D160", "#219C90"), drop = F, name = "Decisions") +
  xlab("Number of patients recruited per arm") +
  ggtitle(levels(simResults$intervention_effect)[i]) +
  #xlim(0,140) +
  ylim(0,50)

} 

png(paste0("plots/decisions_histogram.png"), width = 16, height = 8, units = "in", res = 400)
ggarrange(plotlist = res_hist, nrow = 2, ncol = 3, labels = "AUTO", common.legend = T, legend = "bottom")
dev.off()
######################################################################################################


subdat





