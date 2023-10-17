dat_plot <- res_all

dat_plot$N_swabs_per_day <- as.factor(dat_plot$N_swabs_per_day)
levels(dat_plot$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")

dat_plot$trt_effects2 <- as.factor(dat_plot$trt_effects)
levels(dat_plot$trt_effects2) <-c("trt = 0.8","trt = 1.0","trt = 1.2",
                                 "trt = 1.4", "trt = 1.6")

ggplot(dat_plot, aes(x = as.factor(N), y = exp(Median), col = N_swabs_per_day, group = N_swabs_per_day)) +
  geom_point( position = position_dodge(width = 0.8)) +
  facet_grid(trt_effects2~ day_plans) +
  geom_hline(aes(yintercept = trt_effects), col = 'black', linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("#B70404", "#5A96E3"), name = "") +
  xlab("Number of patients recruited") +
  ylab("Estimated mu") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.5, "lines")) 
  


res_all
