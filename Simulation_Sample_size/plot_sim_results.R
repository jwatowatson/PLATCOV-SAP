#
# f_sims = list.files('sims_out1',full.names = T, pattern = '*.csv')
# names(f_sims) = as.numeric(gsub(gsub(f_sims, pattern = 'sims_out1/sim_out',replacement = ''),pattern = '.csv',replacement = ''))
# f_sims = f_sims[as.character(sort(as.numeric(names(f_sims))))]
#
# sim_res = lapply(f_sims, read.csv)
# sim_all1 = dplyr::bind_cols(sim_res)
# colnames(sim_all1) = names(f_sims)
# load('sims_out1/sim_settings.RData')
# sim_settings1 = sim_settings
#
# f_sims2 = list.files('sims_out2',full.names = T, pattern = '*.csv')
# names(f_sims2) = as.numeric(gsub(gsub(f_sims2, pattern = 'sims_out2/sim_out',replacement = ''),pattern = '.csv',replacement = ''))
# f_sims2 = f_sims2[as.character(sort(as.numeric(names(f_sims2))))]
#
# sim_res2 = lapply(f_sims2, read.csv)
# sim_all2 = dplyr::bind_cols(sim_res2)
# colnames(sim_all2) = names(f_sims2)
# load('sims_out2/sim_settings2.RData')
# sim_settings2 = sim_settings
#
# sim_settings = rbind(sim_settings1, sim_settings2)
# sim_all = cbind(sim_all1, sim_all2)
# save(sim_settings,  file = 'sim_settings.RData')
# write.csv(x = sim_all, file = 'all_sim_results.csv')

sim_all = readr::read_csv(file = 'all_sim_results.csv')[,-1]
load('sim_settings.RData')
cols = RColorBrewer::brewer.pal(n = length(unique(sim_settings$trt_effect_comp)),name = 'Set2')
names(cols)=unique(sim_settings$trt_effect_comp)
## Superiority simulations - check matches what we have observed thus far
ind_superiority = which(sim_settings$trt_control==1)
sims_neg_control = sim_settings[sim_settings$trt_control==1, ]
sims_out_superiority = sim_all[, ind_superiority]

lambda_stop = 1.2
sims_neg_control$test_outcome = apply(sims_out_superiority, 2,
                                      function(x) ifelse(mean(x > log(lambda_stop)) > .9, 1, 0))

sims_neg_control$test_outcome_reject = apply(sims_out_superiority, 2,
                                      function(x) ifelse(mean(x > log(lambda_stop)) < .1, 1, 0))


par(las=1)
par(mfrow=c(1,2))


xx=aggregate(test_outcome ~ N + trt_effect_comp, data = sims_neg_control, mean)
plot(sqrt(xx$N), xx$test_outcome, pch=16, ylim=c(0,1), xaxt='n',
     panel.first=grid(),xlab='N per arm', ylab='Probability of stopping',
     col=cols[as.character(xx$trt_effect_comp)])
title('Success if P[Effect>12.5%]>0.9')
axis(1, at = sqrt(c(25,50,75,100)), labels = c(25,50,75,100))
for(kk in xx$trt_effect_comp){
  ind = xx$trt_effect_comp==kk
  lines(sqrt(xx$N[ind]), xx$test_outcome[ind],lwd=2,col=cols[as.character(xx$trt_effect_comp[ind])])
}

xx=aggregate(test_outcome_reject ~ N + trt_effect_comp, data = sims_neg_control, mean)
plot(sqrt(xx$N), xx$test_outcome_reject, pch=16, ylim=c(0,1), xaxt='n',
     panel.first=grid(),xlab='N per arm', ylab='Probability of stopping',
     col=cols[as.character(xx$trt_effect_comp)])
title('Futility if P[Effect>12.5%]<0.1')
axis(1, at = sqrt(c(25,50,75,100)), labels = c(25,50,75,100))
for(kk in xx$trt_effect_comp){
  ind = xx$trt_effect_comp==kk
  lines(sqrt(xx$N[ind]), xx$test_outcome[ind],lwd=2,col=cols[as.character(xx$trt_effect_comp[ind])])
}

legend('right',title = 'Increase relative to no drug', pch=16,
       inset=0.03,cex=1.3,
       legend = paste0((as.numeric(names(cols))-1)*100,'%'), col=cols, lwd=3)


## Superiority simulations - check matches what we have observed thus far
ind_inferiority = which(sim_settings$trt_control==1.6)
sims_pos_control = sim_settings[sim_settings$trt_control==1.6, ]
sims_out_NI = sim_all[, ind_inferiority]

lambda1= -log(1.6)+log(1.3)
lambda1 = log(0.8)
exp(lambda1)
# sims_pos_control$mean_x = colMeans(sims_out_NI)
# boxplot(mean_x ~ trt_effect_comp, data = sims_pos_control)
sims_pos_control$test_outcome_lambda1 =
  apply(sims_out_NI, 2,
        function(x) ifelse(mean(x > lambda1) > .9, 1, 0))

sims_pos_control$test_outcome_lambda1_reject =
  apply(sims_out_NI, 2,
        function(x) ifelse(mean(x > lambda1) < .1, 1, 0))

sims_pos_control$test_outcome_lower =
  apply(sims_out_NI, 2,
        function(x) ifelse(mean(x < 0) > .99, 1, 0))

par(mfrow=c(1,3),cex.lab=1.5, cex.axis=1.5)
xx=aggregate(test_outcome_lambda1 ~ N + trt_effect_comp, data = sims_pos_control, mean)
plot(sqrt(xx$N), xx$test_outcome_lambda1, pch=16, ylim=c(0,1), xaxt='n',
     panel.first=grid(),xlab='N per arm', ylab='Probability of stopping',
     col=cols[as.character(xx$trt_effect_comp)])
title('Non-inferiority if P[Effect> -20%]>0.9')
axis(1, at = sqrt(c(25,50,75,100)), labels = c(25,50,75,100))
for(kk in xx$trt_effect_comp){
  ind = xx$trt_effect_comp==kk
  lines(sqrt(xx$N[ind]), xx$test_outcome_lambda1[ind], lwd=2,
        col=cols[as.character(xx$trt_effect_comp[ind])])
}


xx=aggregate(test_outcome_lambda1_reject ~ N + trt_effect_comp, data = sims_pos_control, mean)
plot(sqrt(xx$N), xx$test_outcome_lambda1_reject, pch=16, ylim=c(0,1), xaxt='n',
     panel.first=grid(),xlab='N per arm', ylab='Probability of stopping',
     col=cols[as.character(xx$trt_effect_comp)])
title('Not non-inferiority if P[Effect> -20%]<0.1')
axis(1, at = sqrt(c(25,50,75,100)), labels = c(25,50,75,100))
for(kk in xx$trt_effect_comp){
  ind = xx$trt_effect_comp==kk
  lines(sqrt(xx$N[ind]), xx$test_outcome_lambda1_reject[ind], lwd=2,
        col=cols[as.character(xx$trt_effect_comp[ind])])
}
legend('topleft',title = 'Increase relative to no drug', pch=16,
       inset=0.03,cex=1.4,
       legend = paste0((as.numeric(names(cols))-1)*100,'%'), col=cols, lwd=3)


xx=aggregate(test_outcome_lower ~ N + trt_effect_comp,
             data = sims_pos_control,
             mean)
plot(sqrt(xx$N), xx$test_outcome_lower, pch=16, ylim=c(0,1), xaxt='n',
     panel.first=grid(),xlab='N per arm', ylab='Probability of stopping',
     col=cols[as.character(xx$trt_effect_comp)])
title('Stop if P[Effect<0]>0.99')
axis(1, at = sqrt(c(25,50,75,100)), labels = c(25,50,75,100))
for(kk in xx$trt_effect_comp){
  ind = xx$trt_effect_comp==kk
  lines(sqrt(xx$N[ind]), xx$test_outcome_lower[ind], lwd=2,
        col=cols[as.character(xx$trt_effect_comp[ind])])
}
