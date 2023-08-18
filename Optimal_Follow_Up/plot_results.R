library(rstan)

load('model_settings.RData')

res = array(dim = c(nrow(model_settings), 4))


for(i in 1:nrow(model_settings)){
  load(paste('Rout/model_fits_',i,'.RData',sep=''))
  thetas_trt = extract(out, pars='trt_effect')$trt_effect
  res[i, ] = c(quantile(thetas_trt, probs = c(0.025, 0.5, 0.975)), quantile(thetas_trt, probs = 0.975)-
                 quantile(thetas_trt, probs = 0.025))
}

unique_contrast = unique(model_settings[, c('intervention','ref_arm', 'data_ID')])

par(mfrow=c(2,3), las=1)
for(j in 1:nrow(unique_contrast)){
  ind = model_settings$intervention==unique_contrast$intervention[j] &
    model_settings$ref_arm==unique_contrast$ref_arm[j] &
    model_settings$data_ID==unique_contrast$data_ID[j]
  plot(model_settings$Dmax[ind], res[ind, 2],
       panel.first=grid(), ylim = c(-1,1))
  abline(v=7, h=0)
  title(paste(unique_contrast$intervention[j], unique_contrast$ref_arm[j], sep='\n'))
}
