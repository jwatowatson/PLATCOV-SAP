load('sim_settings_extended.RData')

library(ggplot2)
library(fitdistrplus)

res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out/sim_out_extended_',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)

res_all$conditions <- as.numeric(as.factor(paste(res_all$N, res_all$day_plans, res_all$N_swabs_per_day, res_all$trt_effects, sep = "_")))
res_all$conditions <- factor(res_all$conditions, levels = unique(res_all$conditions))
levels(res_all$conditions) <- as.character(1:length(unique(res_all$conditions)))
res_all$conditions <- as.numeric(res_all$conditions)

res_all$N_swabs_per_day <- as.factor(res_all$N_swabs_per_day)
###################################################################################################
threshlods <- c(1, 1.125, 1.15, 1.2)

threshlod_i <- threshlods[4]

power <- NULL
for(j in 1:length(unique(res_all$conditions))){
  subdata <- res_all[res_all$conditions == j,]
  non_sig <- subdata$Upper < log(threshlod_i) | 
         (subdata$Lower < log(threshlod_i) & (subdata$Upper > log(threshlod_i)))
  
  percent_sig <- 1-(sum(non_sig)/nrow(subdata))
  
  summary <- data.frame("conditions" = subdata$condition[1], "trt_effects" = subdata$trt_effects[1], "N" = subdata$N[1], day_plans = subdata$day_plans[1],
                        "N_swabs_per_day" = subdata$N_swabs_per_day[1],"percent_sig" = percent_sig, "threshold" = threshlod_i)
  power <- rbind(power, summary)
}

power

power$trt_effects <- as.factor(power$trt_effects)
levels(power$trt_effects) <- c("trt_effects = 0.8", "trt_effects = 1.0", "trt_effects = 1.2",
                               "trt_effects = 1.4", "trt_effects = 1.6")

power$N_swabs_per_day <- as.factor(power$N_swabs_per_day)
levels(power$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")


png("plots/power.png", width = 10, height = 5, units = "in", res = 400)


ggplot(power, aes(x = N, y = percent_sig, col = day_plans , gr = day_plans, linetype = N_swabs_per_day)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.75) +
  facet_grid( N_swabs_per_day ~ trt_effects) +
  theme_bw() +
  xlab("Number of patients recruited") +
  ylab("Proportion of significant results") +
  scale_x_continuous(breaks = seq(0,100, 25)) +
  scale_color_manual(values = c("#B70404", "#5A96E3", "#5FD068", "#F79327", "black"), name = "Sampling schedule")




dev.off()

###################################################################################################
trt_effects <- 1

plot_intervals <- function(res_all, trt_effects){
  dat_plot <- res_all[res_all$trt_effects == trt_effects,]
  
  dat_plot$conditions <- factor(dat_plot$day_plans, levels = unique(dat_plot$day_plans))
  
  dat_plot$N_swabs_per_day <- as.factor(dat_plot$N_swabs_per_day)
  levels(dat_plot$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")
  
  dat_plot$N <- as.factor(dat_plot$N)
  levels(dat_plot$N) <-c("N = 25","N = 50", "N = 75", "N = 100")
  
  lab <- as.character(sprintf("%.1f", trt_effects))
  
  dat_plot <- dat_plot[with(dat_plot, order(N,conditions, N_swabs_per_day, Median)), ]
  dat_plot$ID <- rep(1:100,16)
  
  
  G <- ggplot(dat_plot, aes(x = ID, y = Median, col = N_swabs_per_day)) +
    geom_point(size = 0.5) +
    facet_grid(day_plans ~ N) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    geom_hline(yintercept = log(1.2), linetype = "dashed", col = "blue") +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), alpha = 0.5)  +
    scale_color_manual(values = c("#B70404", "#5A96E3", "#5FD068", "#F79327"), name = "") +
    coord_flip() +
    ylab("Log effect size") +
    xlab("Simulation ID") +
    ggtitle(paste0("trt_effects = ", lab)) +
    ylim(-1.6,1.6)
  
  G
  
  png(paste0("plots/trt_effect_",lab, ".png"), width = 10, height = 5, units = "in", res = 400)
  print(G)
  dev.off()
}
###################################################################################################
plot_intervals(res_all = res_all, trt_effects = 0.8)
plot_intervals(res_all = res_all, trt_effects = 1)
plot_intervals(res_all = res_all, trt_effects = 1.2)
plot_intervals(res_all = res_all, trt_effects = 1.4)
plot_intervals(res_all = res_all, trt_effects = 1.6)
###################################################################################################
fit_sd <- function(res_all, trt_effects){
  dat_plot <- res_all[res_all$trt_effects == trt_effects,]
  
  dat_plot$conditions <- factor(dat_plot$day_plans, levels = unique(dat_plot$day_plans))
  dat_plot <- dat_plot[with(dat_plot, order(N,conditions, N_swabs_per_day, Median)), ]


  dat_plot$combinations <- as.factor(paste0(dat_plot$N, dat_plot$day_plans, dat_plot$N_swabs_per_day))
  dat_plot$combinations <- factor(dat_plot$combinations, levels = unique(dat_plot$combinations))
  levels(dat_plot$combinations) <- as.character(1:length(levels(dat_plot$combinations)))
  dat_plot$combinations <- as.numeric(dat_plot$combinations)

  summary <- NULL
  for (k in unique(dat_plot$combinations)) { # <- here
    sub_d <- dat_plot[dat_plot$combinations == k,]
    fit1 = fitdist(sub_d$Median, "norm", fix.arg = list(mean = log(trt_effects))) 
  
    summary_sd <- data.frame("N" = sub_d$N[1], "day_plans" = sub_d$day_plans[1], "N_swabs_per_day" = sub_d$N_swabs_per_day[1],
                             "trt_effects" = trt_effects, "sd" = fit1$estimate)
  
    summary <- rbind(summary, summary_sd)
  }
  
  
  return(summary)
}
###################################################################################################
sd_all  <-  rbind(fit_sd(res_all = res_all, trt_effects = 0.8),
                  fit_sd(res_all = res_all, trt_effects = 1),
                  fit_sd(res_all = res_all, trt_effects = 1.2),
                  fit_sd(res_all = res_all, trt_effects = 1.4),
                  fit_sd(res_all = res_all, trt_effects = 1.6))

sd_all$N_swabs_per_day <- as.factor(sd_all$N_swabs_per_day)
levels(sd_all$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")

 sd_all$N <- as.factor(sd_all$N)

 sd_all$trt_effects <- as.factor(sd_all$trt_effects)
 levels(sd_all$trt_effects) <- c("trt_effects = 0.8", "trt_effects = 1.0", "trt_effects = 1.2",
                                "trt_effects = 1.4", "trt_effects = 1.6")
 
png(paste0("plots/std_err.png"), width = 10, height = 8, units = "in", res = 400)
ggplot(sd_all, aes(x = N, y = sd, shape = N_swabs_per_day, col = day_plans)) +
  geom_point(size = 3, alpha = 0.80)+
  facet_grid(.~day_plans) +
  scale_color_manual(values = c("#B70404", "#5A96E3", "#5FD068", "#F79327"),  name = "Sampling schedule") +
  scale_shape_manual(values = c(16, 17), name = "") +
  theme_bw() +
  ylim(0,0.2) +
  xlab("Number of patients recruited") +
  ylab("Standard errors for treatment effects\n across simulations") +
  theme(axis.title = element_text(size = 14, face = "bold"))
dev.off()
###################################################################################################








