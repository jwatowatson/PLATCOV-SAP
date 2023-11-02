library(rstan)

load('model_settings.RData')
model_settings_effective<- model_settings

load('model_settings_ineffective.RData')
model_settings_ineffective <- model_settings

res = array(dim = c(nrow(model_settings_effective) + nrow(model_settings_ineffective), 6))

for(i in 1:nrow(model_settings_effective)){
  load(paste('output/Rout_effective/model_fits_',i,'.RData',sep=''))
  thetas_trt = rstan::extract(out, pars='trt_effect')$trt_effect
  res[i, ] = c(quantile(thetas_trt, probs = c(0.025, 0.5, 0.975)),
               quantile(thetas_trt, probs = 0.975)-
                 quantile(thetas_trt, probs = 0.025),
               sd(thetas_trt),
               mean(thetas_trt)/sd(thetas_trt))
}

for(i in (nrow(model_settings_effective)+1):nrow(res)){
  load(paste('output/Rout_ineffective/model_fits_ineffective_',(i-nrow(model_settings_effective)),'.RData',sep=''))
  thetas_trt = rstan::extract(out, pars='trt_effect')$trt_effect
  #model_settings_ineffective[(i-nrow(model_settings_effective)),]
  #hist(exp(thetas_trt))
  
  res[i, ] = c(quantile(thetas_trt, probs = c(0.025, 0.5, 0.975)),
               quantile(thetas_trt, probs = 0.975)-
                 quantile(thetas_trt, probs = 0.025),
               sd(thetas_trt),
               mean(thetas_trt)/sd(thetas_trt))
}

model_settings <- rbind(model_settings_effective, model_settings_ineffective)

unique_contrast = unique(model_settings[, c('intervention',
                                            'ref_arm', 
                                            'data_ID')])

par(mfrow=c(3,3), las=1)
for(j in 1:nrow(unique_contrast)){
  ind = model_settings$intervention==unique_contrast$intervention[j] &
    model_settings$ref_arm==unique_contrast$ref_arm[j] &
    model_settings$data_ID==unique_contrast$data_ID[j]
  dat_out = data.frame(trt_effect = res[ind, 2], day_FUP=model_settings$Dmax[ind])
  plot(jitter(model_settings$Dmax[ind]), res[ind, 2],
       panel.first=grid(), ylim = range(res[,2]),
       xlab='Days follow-up included', ylab = 'Treatment effect (log)')
  xx=aggregate(trt_effect ~ day_FUP, FUN = mean, data = dat_out)
  lines(xx$day_FUP, xx$trt_effect,lwd=2,col='red')
  abline(v=7, h=0)
  title(paste(unique_contrast$intervention[j], unique_contrast$ref_arm[j], sep='\n'))
}



par(mfrow=c(3,3), las=1)
for(j in 1:nrow(unique_contrast)){
  ind = model_settings$intervention==unique_contrast$intervention[j] &
    model_settings$ref_arm==unique_contrast$ref_arm[j] &
    model_settings$data_ID==unique_contrast$data_ID[j]
  dat_out = data.frame(sd_size = res[ind, 4], 
                       day_FUP=model_settings$Dmax[ind])
  plot(jitter(model_settings$Dmax[ind]), res[ind, 4],
       panel.first=grid(), ylim = c(0,max(res[,4])),
       xlab='Days follow-up included', 
       ylab = 'Standard error width (log)')
  xx=aggregate(sd_size ~ day_FUP, FUN = mean, data = dat_out)
  lines(xx$day_FUP, xx$sd_size,lwd=2,col='red')
  abline(v=7, h=0)
  title(paste(unique_contrast$intervention[j], unique_contrast$ref_arm[j], sep='\n'))
}


png("plots/z_scores.png", width = 8, height = 6, units = "in", res = 350)
labs <- c("A", "B", "C", "D", "E")
notes <- c("", "", "", "Before Feb 2023", "After Feb 2023")
par(mfrow=c(2,3), las=1)
for(j in 1:5){#nrow(unique_contrast)){
  
  ind = model_settings$intervention==unique_contrast$intervention[j] &
    model_settings$ref_arm==unique_contrast$ref_arm[j] &
    model_settings$data_ID==unique_contrast$data_ID[j]
  dat_out = data.frame(zscore = res[ind, 6],
                       day_FUP=model_settings$Dmax[ind])
  
  xx=aggregate(zscore ~ day_FUP, FUN = mean,data = dat_out)
  
  plot(xx$day_FUP, xx$zscore,
       panel.first=grid(),ylim=c(min(res[,6]), max(res[,6])),
       xlab='Days follow-up included',
       ylab = 'Z-score estimate',lwd=3,type='b')
  points(jitter(dat_out$day_FUP), dat_out$zscore,col='grey')
  lines(xx$day_FUP, xx$zscore,type='b',lwd=3,col='red')
  abline(v = xx$day_FUP[which.max(xx$zscore)],lty=2)
  abline(h = 0,lty=2)
  title(paste(unique_contrast$intervention[j], unique_contrast$ref_arm[j], sep='\n'))
  mtext(labs[j], 2, adj=3, las=1, padj = -10)
  mtext(notes[j], 2, adj=-0.25, las=1, padj = 6)
  
}
dev.off()
