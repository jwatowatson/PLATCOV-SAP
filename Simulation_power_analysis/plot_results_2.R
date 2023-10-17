load('sim_settings_extended.RData')

library(ggplot2)
library(fitdistrplus)
library(tidyr)

stop_threslod <- 1.2

Summary <- NULL

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out/sim_out_extended_',i,'.csv',sep=''))
  quantile <- as.data.frame(t(quantile(res_data$trt_effect, c(0.5, 0.025, 0.975, 0.1, 0.9))))
  success <- ifelse(sum(res_data$trt_effect > log(stop_threslod))/nrow(res_data) > 0.9, 1, 0) 
  futility <- ifelse(sum(res_data$trt_effect < log(stop_threslod))/nrow(res_data) > 0.9, 1, 0) 
  inconclusive <- ifelse((success + futility == 0), 1, 0)
  summary <- data.frame(quantile, success, futility, inconclusive)
  Summary <- rbind(Summary, summary)
}

Summary
colnames(Summary)[1:5] <- c("Median", "Lower95", "Upper95", "Lower80", "Upper80")

res_all <- cbind(sim_settings, Summary)
res_all$decisions <- "inconclusive"
res_all$decisions[res_all$success == 1] <- "success"
res_all$decisions[res_all$futility == 1] <- "futility"
res_all$decisions <- as.factor(res_all$decisions)

###################################################################################################
# res_agg <- aggregate(list(res_all$success, res_all$futility, res_all$inconclusive), 
#                      by=list(res_all$N, res_all$day_plans, res_all$N_swabs_per_day, res_all$trt_effects), FUN=sum)

# colnames(res_agg) <- c("N", "day_plans", "N_swabs_per_day", "trt_effects",
#                        "success", "futility", "inconclusive")
# res_plot <- gather(res_agg, "decisions", "counts", success:inconclusive, factor_key = T)


res_all$trt_effects <- as.factor(res_all$trt_effects)
levels(res_all$trt_effects) <- c("trt_effects = 0.8", "trt_effects = 1.0", "trt_effects = 1.2",
                               "trt_effects = 1.4", "trt_effects = 1.6")

res_all$N_swabs_per_day <- as.factor(res_all$N_swabs_per_day)
levels(res_all$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")

library(ggplot2)  
library(ggpubr)

res_plots <- list()

for(i in 1:length(levels(res_all$trt_effects))){
  res_plots[[i]] <- ggplot(res_all[res_all$trt_effects == levels(res_all$trt_effects)[i],], aes(as.factor(N))) +
  geom_bar(aes(fill = decisions), col = "black") +
  facet_grid(N_swabs_per_day ~ day_plans) +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Number of patients recruited") +
  ylab("Counts") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  scale_fill_manual(values = c("#BB2525", "#F4D160", "#219C90"), drop = F, name = "Decisions") +
  ggtitle(levels(res_all$trt_effects)[i])
} 

png(paste0("plots/decisions.png"), width = 15, height = 8, units = "in", res = 400)
ggarrange(plotlist = res_plots, nrow = 2, ncol = 3, labels = "AUTO", common.legend = T, legend = "bottom")
dev.off()





