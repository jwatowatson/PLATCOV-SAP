library(rstan)
library(ggpubr)
library(ggplot2)
############################################################################################
load('sim_settings_extended.RData')
############################################################################################
Summary <- NULL
stop_threslod <- 1.2
############################################################################################
for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out/sim_out_extended_',i,'.csv',sep=''))
  quantile <- as.data.frame(t(quantile(res_data$trt_effect, c(0.5, 0.025, 0.975, 0.1, 0.9))))
  success <- ifelse(sum(res_data$trt_effect > log(stop_threslod))/nrow(res_data) > 0.9, 1, 0) 
  futility <- ifelse(sum(res_data$trt_effect < log(stop_threslod))/nrow(res_data) > 0.9, 1, 0) 
  inconclusive <- ifelse((success + futility == 0), 1, 0)
  
  mu_hat <- mean(res_data$trt_effect)
  sigma_hat <- sd(res_data$trt_effect)
  
  summary <- data.frame(quantile, success, futility, inconclusive, mu_hat, sigma_hat)
  Summary <- rbind(Summary, summary)
}

Summary
colnames(Summary)[1:5] <- c("Median", "Lower95", "Upper95", "Lower80", "Upper80")
############################################################################################
res_all <- cbind(sim_settings, Summary)
res_all$decisions <- "inconclusive"
res_all$decisions[res_all$success == 1] <- "success"
res_all$decisions[res_all$futility == 1] <- "futility"
res_all$decisions <- as.factor(res_all$decisions)

res_all$conditions <- as.numeric(as.factor(paste(res_all$N, res_all$day_plans, res_all$N_swabs_per_day, sep = "_")))
res_all$conditions <- factor(res_all$conditions, levels = unique(res_all$conditions))
levels(res_all$conditions) <- as.character(1:length(unique(res_all$conditions)))
res_all$conditions <- as.numeric(res_all$conditions)

N_conditions <- length(unique(res_all$conditions))

res_all$sim_ID <- row.names(res_all)
############################################################################################
data_for_stan <- list(
  N_sim = nrow(res_all),
  mu_hat = res_all$mu_hat,
  sigma_hat = res_all$sigma_hat,
  mu = log(res_all$trt_effects),
  condition_i = res_all$conditions,
  N_conditions = N_conditions
)
############################################################################################
mod = stan_model(file = "model_precisions.stan") # compile 

out <- sampling(mod, 
         data=data_for_stan,
         iter=4000,
         chain=4,
         thin=8,
         warmup=2000,
         save_warmup = FALSE,
         cores = 4,
         seed=1,
         include=T)
out
############################################################################################
post <- extract(out)

condition_data <- res_all[!duplicated(res_all$conditions),c("N", "day_plans", "N_swabs_per_day", "conditions")]

condition_data$sigma_mu <- apply(post$sigma_mu,2,median)
condition_data$sigma_mu_low <- apply(post$sigma_mu,2,quantile, 0.025)
condition_data$sigma_mu_up <- apply(post$sigma_mu,2,quantile, 0.975)

condition_data$mu_sigma <- apply(post$mu_sigma,2,median)
condition_data$mu_sigma_low <- apply(post$mu_sigma,2,quantile, 0.025)
condition_data$mu_sigma_up <- apply(post$mu_sigma,2,quantile, 0.975)

condition_data$sigma_sigma <- apply(post$sigma_sigma,2,median)
condition_data$sigma_sigma_low <- apply(post$sigma_sigma,2, quantile, 0.025)
condition_data$sigma_sigma_up <- apply(post$sigma_sigma,2, quantile, 0.975)
############################################################################################
condition_data$N_swabs_per_day <- as.factor(condition_data$N_swabs_per_day)
levels(condition_data$N_swabs_per_day) <-c("2 swabs per day","4 swabs per day")

condition_data$N <- as.factor(condition_data$N)

#Accuracy
G1 <- ggplot(condition_data, aes(x = N, y = sigma_mu, col = N_swabs_per_day, group = N_swabs_per_day)) +
  geom_errorbar(aes(ymin = sigma_mu_low, ymax = sigma_mu_up), col = "black", width = 0.5, linewidth = 0.2, 
                position = position_dodge(width = 0.8)) +
  geom_point(size = 2.5, alpha = 0.75,position = position_dodge(width = 0.8)) +
  facet_grid(.~ day_plans) +
  theme_bw() +
  scale_color_manual(values = c("#B70404", "#5A96E3"), name = "") +
  xlab("Number of patients recruited") +
  ylab("Variance of mu") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.5, "lines")) +
  ggtitle("Variance of mu")

#Precision
G2 <- ggplot(condition_data, aes(x = N, y = mu_sigma, col = N_swabs_per_day, group = N_swabs_per_day)) +
  geom_errorbar(aes(ymin = mu_sigma_low, ymax = mu_sigma_up), col = "black", width = 0.5, linewidth = 0.2, 
                position = position_dodge(width = 0.8)) +  
  facet_grid(.~ day_plans) +
  geom_point(size = 2.5, alpha = 0.75,position = position_dodge(width = 0.8)) +
  theme_bw() +
  scale_color_manual(values = c("#B70404", "#5A96E3"), name = "") +
  xlab("Number of patients recruited") +
  ylab("Mean of sigma") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.5, "lines")) +
  ggtitle("Mean of sigma")
G2

G3 <- ggplot(condition_data, aes(x = N, y = sigma_sigma, col = N_swabs_per_day, group = N_swabs_per_day)) +
  geom_errorbar(aes(ymin = sigma_sigma_low, ymax = sigma_sigma_up), col = "black", width = 0.5, linewidth = 0.2, 
                position = position_dodge(width = 0.8)) +  
  facet_grid(.~ day_plans) +
  geom_point(size = 2.5, alpha = 0.75,position = position_dodge(width = 0.8)) +
  theme_bw() +
  scale_color_manual(values = c("#B70404", "#5A96E3"), name = "") +
  xlab("Number of patients recruited") +
  ylab("Variance of sigma") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.5, "lines"))  +
  ggtitle("Variance of sigma")
G3
############################################################################################
png(paste0("plots/precision_accuracy.png"), width = 8, height = 8, units = "in", res = 400)
ggarrange(G1, G2, G3, nrow = 3, labels = "AUTO", common.legend = T, legend = "right")
dev.off()
