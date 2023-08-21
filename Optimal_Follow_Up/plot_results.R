library(rstan)

load('model_settings.RData')

res = array(dim = c(nrow(model_settings), 6))


for(i in 1:nrow(model_settings)){
  load(paste('~/Downloads/Rout/model_fits_',i,'.RData',sep=''))
  thetas_trt = rstan::extract(out, pars='trt_effect')$trt_effect
  res[i, ] = c(quantile(thetas_trt, probs = c(0.025, 0.5, 0.975)),
               quantile(thetas_trt, probs = 0.975)-
                 quantile(thetas_trt, probs = 0.025),
               sd(thetas_trt),
               mean(thetas_trt)/sd(thetas_trt))
}

unique_contrast = unique(model_settings[, c('intervention',
                                            'ref_arm', 
                                            'data_ID')])

par(mfrow=c(2,3), las=1)
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



par(mfrow=c(2,3), las=1)
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




par(mfrow=c(2,3), las=1)
for(j in 1:nrow(unique_contrast)){
  ind = model_settings$intervention==unique_contrast$intervention[j] &
    model_settings$ref_arm==unique_contrast$ref_arm[j] &
    model_settings$data_ID==unique_contrast$data_ID[j]
  dat_out = data.frame(zscore = res[ind, 6],
                       day_FUP=model_settings$Dmax[ind])
  
  xx=aggregate(zscore ~ day_FUP, FUN = mean,data = dat_out)
  
  plot(xx$day_FUP, xx$zscore,
       panel.first=grid(),ylim=c(0, max(res[,6])),
       xlab='Days follow-up included',
       ylab = 'Z-score estimate',lwd=3,type='b')
  points(jitter(dat_out$day_FUP), dat_out$zscore,col='grey')
  lines(xx$day_FUP, xx$zscore,type='b',lwd=3,col='red')
  abline(v = xx$day_FUP[which.max(xx$zscore)],lty=2)
  title(paste(unique_contrast$intervention[j], unique_contrast$ref_arm[j], sep='\n'))
}
