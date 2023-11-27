library(rstan)
library(dplyr)
library(ggpubr)
#####################################################################################
load('model_settings.RData')
model_settings_effective<- model_settings
model_settings_effective <- model_settings_effective %>%
  filter(data_ID != 3)

res_1 = array(dim = c(nrow(model_settings_effective), 6))

for(i in 1:nrow(model_settings_effective)){
  load(paste('output/Rout_effective/model_fits_',i,'.RData',sep=''))
  thetas_trt = rstan::extract(out, pars='trt_effect')$trt_effect
  res_1[i, ] = c(quantile(thetas_trt, probs = c(0.025, 0.5, 0.975)),
               quantile(thetas_trt, probs = 0.975)-
                 quantile(thetas_trt, probs = 0.025),
               sd(thetas_trt),
               mean(thetas_trt)/sd(thetas_trt))
}
#####################################################################################
load('model_settings_bootstraps.RData')
model_settings_bootstraps <- model_settings

res_2 = array(dim = c(nrow(model_settings_bootstraps), 6))

for(i in 1:nrow(model_settings_bootstraps)){
  load(paste('Rout/model_fits_bootstraps_',i,'.RData',sep=''))
  thetas_trt = rstan::extract(out, pars='trt_effect')$trt_effect
  res_2[i, ] = c(quantile(thetas_trt, probs = c(0.025, 0.5, 0.975)),
                 quantile(thetas_trt, probs = 0.975)-
                   quantile(thetas_trt, probs = 0.025),
                 sd(thetas_trt),
                 mean(thetas_trt)/sd(thetas_trt))
}
#####################################################################################
model_settings <- rbind(model_settings_effective, model_settings_bootstraps)
res <- rbind(res_1, res_2)
res <- as.data.frame(res)
colnames(res) <- c("low", "med", "up", "IC_width", "sd", "z_score")

results <- cbind(model_settings, res)
results$contrasts <- paste0(results$intervention, results$ref_arm, results$data_ID)

unique_contrast <- unique(results$contrasts)
unique_contrast <- unique_contrast[c(1,2,6,3,4,5)]
#####################################################################################
results_mean <- results %>%
  group_by(contrasts, Dmax) %>%
  summarise(med_z = mean(z_score)) %>%
  as.data.frame()

tags <- c("A", "B", "C", "D", "E", "F")

plot_list <- list()

for(i in 1:length(unique_contrast)){
  data_plot_all <- results[results$contrasts == unique_contrast[i], ]
  data_plot_mean <- results_mean[results_mean$contrasts == unique_contrast[i], ]
  
  intervention <- data_plot_all$intervention[1]
  ref_arm <- data_plot_all$ref_arm[1]
  
  if(intervention == "Nirmatrelvir + Ritonavir"){intervention <- "Nirmatrelvir"}
  

  if(intervention == "Nirmatrelvir" & data_plot_all$data_ID[1] == 2 &
     ref_arm == "No study drug"){ref_arm <- "No study drug (before Feb 2023)"}
  if(intervention == "Nirmatrelvir" & data_plot_all$data_ID[1] == 3){ref_arm <- "No study drug (after Feb 2023)"}
  if(intervention == "Regeneron"){intervention <- "Casirivimab/imdevimab"}

  day_max_z <- data_plot_mean$Dmax[data_plot_mean$med_z == max(data_plot_mean$med_z)]
  
  lab <- paste0(intervention, " vs \n", ref_arm)
  
  tag <- tags[i]
  
  plot_list[[i]] <-
  
  ggplot() +
    geom_jitter(data = data_plot_all, aes(x = Dmax, y = z_score), width = 0.15, size = 1.2, shape = 1,
                stroke =1, col = "#7D7C7C", alpha = 0.5) +
    geom_vline(xintercept = day_max_z, linetype = "dashed", linewidth = 0.55) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.55) +
    geom_line(data = data_plot_mean,  aes(x = Dmax, y = med_z),
              col = "red", linewidth = 1) +
    geom_point(data = data_plot_mean,  aes(x = Dmax, y = med_z), size = 2.5, shape = 21,
               stroke = 1, col = "white", fill = "red") +
    theme_bw() +
    scale_x_continuous(breaks = seq(2,14,2), minor_breaks = 2:14) +
    ylim(-3,13) +
    xlab("") +
    ylab("") +
    ggtitle(lab) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
          axis.text = element_text(size = 8),
          plot.tag.position = c(0.175,0.99),
          plot.tag = element_text(size = 9, face = "bold")) +
    labs(tag = tag)
  
  
  
}


G <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 2) #, labels = "AUTO", hjust = -5.5,
              # font.label = list(size = 10))
G

png("plots/z_scores.png", width = 8, height = 6, units = "in", res = 350)
annotate_figure(G, left = textGrob("Z-score", rot = 90, vjust = 0.5, gp = gpar(cex = 0.9, fontface="bold")),
                bottom = textGrob("Days follow-up included", gp = gpar(cex = 0.9, fontface="bold")))
dev.off()



