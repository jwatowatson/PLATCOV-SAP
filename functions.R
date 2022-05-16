plot_effect_estimates = function(effect_ests, plot_models, names_plot_models, place='topright'){
  
  par(bty='n', cex.lab=1.3, cex.axis=1.3,family='serif',
      mar=c(5,11,2,2),las=1, mfrow=c(1,1))
  mod_cols = brewer.pal(n = length(plot_models), name = 'Set1')
  
  K_treatments = nrow(effect_ests[[plot_models[1]]])
  
  xlims = range(sapply(effect_ests[plot_models], rbind))
  plot(NA, NA, xlim = xlims, 
       ylim = c(0.5,K_treatments+.5), 
       panel.first=grid(), ylab='', yaxt='n', type='n', 
       xlab = 'Multiplicative change in slope relative to control',
       xaxt = 'n')
  axis(1, at = log(c(0.7,1,1.5,2)), labels = c(0.7,1,1.5,2))
  axis(2, at = 1:K_treatments,
       labels = rownames(effect_ests[[plot_models[1]]]),
       tick = F)
  index_p = rev(seq(-.2,.2, length.out = length(plot_models)))
  abline(v=0,lwd=2)
  polygon(c(-10, log(study_threshold), log(study_threshold), -10), 
          c(-100, -100, 100, 100), border = NA, 
          col = adjustcolor('grey',.4))
  for(i in 1:length(plot_models)){
    points(effect_ests[[plot_models[i]]][,'mean'], 
           1:K_treatments+index_p[i],
           pch=16, col=mod_cols[i])
    for(j in 1:K_treatments){
      lines(c(effect_ests[[plot_models[i]]][j,'2.5%'],
              effect_ests[[plot_models[i]]][j,'97.5%']),
            rep(j+index_p[i],2),col=mod_cols[i],lwd=1)
      lines(c(effect_ests[[plot_models[i]]][j,'10%'], 
              effect_ests[[plot_models[i]]][j,'90%']),
            rep(j+index_p[i],2),col=mod_cols[i],lwd=3)
    }
  }
  
  legend(place, col=mod_cols,lwd = 2, 
         title = 'Model',cex=1,bty='n',
         legend = names_plot_models)
  
}