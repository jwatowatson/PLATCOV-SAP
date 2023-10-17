load('Rout/model_fits1.RData')
########################################################################################
post <- extract(out)
########################################################################################
png(paste0("plots/posterior_", 
           gsub(" \\+.*", "", intervention),"_", gsub(" \\+.*", "", ref_arm), ".png"), 
    width = 9, height = 6, units = "in", res = 350)

par(mfrow = c(2,2))
#1
hist(exp(post$trt_effect), main = "Treatment effects of \nPaxlovid vs No study drug (After March 2023)", 
     xlab = "Treatment effects", col = "#3085C3")
#2
hist((post$beta_0), main = "Baseline viral clearance rate (After March 2023)", 
     xlab = "Baseline viral clearance rate", col = "#3085C3")
#3
hist(post$sigma_logvl, main = "Sigma parameter of t distribution", 
     xlab = "Sigma parameter of t distribution", col = "#3085C3")
#4
hist(post$t_dof, main = "Degree of freedom of t distribution", 
     xlab = "Degree of freedom of t distribution", col = "#3085C3")

par(mfrow = c(1,1))
dev.off()
########################################################################################
