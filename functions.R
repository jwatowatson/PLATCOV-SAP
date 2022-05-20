plot_effect_estimates = function(effect_ests, #list of stan outputs
                                 plot_models, # indices of models to plot in list
                                 names_plot_models, # names for plotting
                                 place='topright', 
                                 mylty=1,
                                 mod_cols = NULL,
                                 left_space=8,
                                 study_threshold){
  
  if (length(mylty)==1) mylty = rep(1,length(plot_models))
  if (length(mylty)!=length(plot_models)) stop('length of mylty needs to be the same as the number of input models')
  if (is.null(mod_cols)){
    mod_cols = brewer.pal(n = length(plot_models), name = 'Set1')
    if(length(plot_models)>9) writeLines('too many models - supply user defined colors')
  }
  if(length(mod_cols)!=length(plot_models)) stop('number of colors needs to be equal to number of models to be plotted')
  
  par(bty='n', cex.lab=1.3, cex.axis=1.3,family='serif',
      mar=c(5,left_space,2,2),las=1, mfrow=c(1,1))
  
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
            rep(j+index_p[i],2),col=mod_cols[i],lwd=1,lty=mylty[i])
      lines(c(effect_ests[[plot_models[i]]][j,'10%'],
              effect_ests[[plot_models[i]]][j,'90%']),
            rep(j+index_p[i],2),col=mod_cols[i],lwd=3,lty=mylty[i])
    }
  }
  
  legend(place, col=mod_cols,lwd = 2,lty=mylty,
         title = 'Model',cex=1,bty='n',
         legend = names_plot_models)
}



plot_baseline_data = function(platcov_dat){
  
  baseline_ind = platcov_dat$Timepoint_ID==0
  bvl = aggregate(log10_viral_load ~ ID, platcov_dat[baseline_ind, ], median)
  
  hist(bvl$log10_viral_load,
       breaks = seq(1,8.5,by=.5),
       xlab='Baseline viral load (RNA copies per mL)',
       ylab ='Number of patients',xlim=c(1,8.5),
       main='', xaxt ='n')
  axis(1, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  grid(); 
  hist(bvl$log10_viral_load,breaks = seq(1,8.5,by=.5),add=T)
}

plot_serial_data = function(platcov_dat,trt_cols){
  
  PCR_dat = aggregate(log10_viral_load ~ ID + Timepoint_ID +
                        Trt_number + Trt, 
                      data = platcov_dat, mean)
  trt_smmry = aggregate(formula = log10_viral_load ~ Timepoint_ID+Trt_number, 
                        data = PCR_dat, FUN = median)
  PCR_dat$Timepoint_ID=jitter(PCR_dat$Timepoint_ID)
  
  gap.plot(PCR_dat$Timepoint_ID, PCR_dat$log10_viral_load,
           ylab = 'RNA copies per mL', panel.first=grid(),
           xlab = 'Time since randomisation (days)',
           gap = c(7.5,13.5), gap.axis = 'x',
           yticlab = '',ytics = 2, xtics = c(0,3,6,14),
           xlim = c(0,14), type='n', yaxt='n',
           col = trt_cols[PCR_dat$Trt_number],
           ylim = c(1, max(PCR_dat$log10_viral_load)))
  axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  IDs = unique(platcov_dat$ID)
  writeLines(sprintf('Plotting data for %s individuals', length(IDs)))
  for(id in IDs){
    ind = PCR_dat$ID==id
    gap.plot(PCR_dat$Timepoint_ID[ind],
             PCR_dat$log10_viral_load[ind],
             gap = c(7.5,13.5), gap.axis = 'x',add = T,
             col=adjustcolor(trt_cols[PCR_dat$Trt_number[ind]],.3))
  }
  gap.plot(trt_smmry$Timepoint_ID, trt_smmry$log10_viral_load,
           col= trt_cols[trt_smmry$Trt_number],
           gap = c(7.5,13.5), gap.axis = 'x',add = T,
           pch = 14+trt_smmry$Trt_number, cex=1.5)
  for(tt in unique(trt_smmry$Trt_number)){
    ind = trt_smmry$Trt_number==tt
    gap.plot(trt_smmry$Timepoint_ID[ind],
             trt_smmry$log10_viral_load[ind], 
             gap = c(7.5,13.5), gap.axis = 'x',add = T,
             type='l',
             col = trt_cols[tt],lwd=3)
  }
  legend('topright', col=trt_cols, 
         legend = unique(PCR_dat$Trt),title = 'Median',
         lwd=2,pch=14+unique(PCR_dat$Trt_number),
         cex=1.2, inset = 0.03)
}

bayes_R2 = function(mod_preds, mod_residuals) {
  var_pred = apply(mod_preds, 1, var)
  var_res = apply(mod_residuals, 1, var)
  var_pred / (var_pred + var_res)
}


make_stan_inputs = function(platcov_dat_fit, 
                            int_covs_base,
                            int_covs_full,
                            slope_covs_base,
                            slope_covs_full,
                            trt_col,
                            Dmax
){
  
  ## re-arrange so that censored values come last
  platcov_dat_fit = dplyr::arrange(platcov_dat_fit,
                                   log10_viral_load==log10_cens_vl)
  ind_dup = !duplicated(platcov_dat_fit$ID)
  
  platcov_dat_fit$RnaseP_scaled = t(scale(40 - platcov_dat_fit$CT_RNaseP, 
                                          scale = F))[1,] 
  platcov_dat_fit$Age_scaled = (platcov_dat_fit$Age-mean(platcov_dat_fit$Age[ind_dup]))/sd(platcov_dat_fit$Age[ind_dup])
  
  # make the covariate matrix
  # check no missing data
  print(all(!apply(platcov_dat_fit[, intersect(int_covs_full,slope_covs_full)], 2, function(x) any(is.na(x)))))
  
  # make the covariate matrix
  print(all(!apply(platcov_dat_fit[, c('Epoch','Age_scaled','Antibody_test','Symptom_onset','N_dose')], 2, function(x) any(is.na(x)))))
  
  X_intcpt_1 = model.matrix( ~ ., 
                             data = platcov_dat_fit[, int_covs_base])[, -1, drop=F]
  X_intcpt_2 = model.matrix( ~ ., 
                             data = platcov_dat_fit[, int_covs_full])[, -1, drop=F]
  X_slope_1 = model.matrix( ~ ., 
                            data = platcov_dat_fit[, slope_covs_base])[, -1, drop=F]
  X_slope_2 = model.matrix( ~ ., 
                            data = platcov_dat_fit[, slope_covs_full])[, -1, drop=F]
  nrow(X_intcpt_1) == nrow(platcov_dat_fit)
  nrow(X_slope_1) == nrow(platcov_dat_fit)
  
  cov_matrices = list(X_int=list(X_intcpt_1, X_intcpt_2),
                      X_slope=list(X_slope_1, X_slope_2))
  
  ID_map = data.frame(ID_key = platcov_dat_fit$ID,
                      ID_stan = as.numeric(as.factor(platcov_dat_fit$ID)))
  
  ind_cens = !platcov_dat_fit$log10_viral_load>
    platcov_dat_fit$log10_cens_vl
  
  writeLines(sprintf('%s%% of samples are below LOD',
                     round(100*mean(ind_cens),digits = 2)))
  
  analysis_data_stan = list(Ntot = nrow(platcov_dat_fit),
                            N_obs = sum(!ind_cens),
                            n_id = max(ID_map$ID_stan),
                            id = ID_map$ID_stan,
                            ind_start = which(!duplicated(ID_map$ID_stan)),
                            obs_day = platcov_dat_fit$Time,
                            log_10_vl = platcov_dat_fit$log10_viral_load,
                            log10_cens_vl = platcov_dat_fit$log10_cens_vl,
                            RNaseP = platcov_dat_fit$RnaseP_scaled,
                            Time_max = Dmax)
  
  writeLines('check stan data formatting:')
  all(analysis_data_stan$log_10_vl[1:analysis_data_stan$N_obs]>
        analysis_data_stan$log10_cens_vl[1:analysis_data_stan$N_obs]) &
    all(analysis_data_stan$log_10_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot] ==
          analysis_data_stan$log10_cens_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot])
  
  Trt_matrix = model.matrix( ~ ., data = platcov_dat_fit[, trt_col, drop=F])
  Trt_matrix[,1]=0 # first column is dummy
  
  analysis_inputs = list(cov_matrices=cov_matrices,
                         analysis_data_stan=analysis_data_stan,
                         Trt_matrix=Trt_matrix,
                         ID_map=ID_map)
  return(analysis_inputs)
}



# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}


checkStrict(make_stan_inputs)
checkStrict(plot_serial_data)
checkStrict(plot_effect_estimates)
