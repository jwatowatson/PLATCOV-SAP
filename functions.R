# list of functions for plotting data

plot_effect_estimates = function(effect_ests, #list of stan outputs
                                 plot_models, # indices of models to plot in list
                                 my_pch=1,
                                 mod_cols = NULL,
                                 study_threshold){
  
  if (length(my_pch)==1) my_pch = (1:length(plot_models))+15
  if (length(my_pch)!=length(plot_models)) stop('length of my_pch needs to be the same as the number of input models')
  if (is.null(mod_cols)){
    mod_cols = brewer.pal(n = length(plot_models), name = 'Set1')
    if(length(plot_models)>9) writeLines('too many models - supply user defined colors')
  }
  if(length(mod_cols)!=length(plot_models)) stop('number of colors needs to be equal to number of models to be plotted')
  
  K_treatments = nrow(effect_ests[[plot_models[1]]])
  
  xlims = (exp(range(c(0, range(sapply(effect_ests[plot_models], rbind)))) )-1)*100
  x_points = pretty((xlims),6)
  plot(NA, NA, xlim = range(x_points),
       ylim = c(0.75,K_treatments+.25),
       panel.first=grid(), ylab='', yaxt='n', type='n',
       xlab = 'Change in rate of clearance (%)',
       xaxt = 'n')
  
  index_p = rev(seq(-.2,.2, length.out = length(plot_models)))
  abline(v=0,lwd=2)
  polygon(c(-1000, 100*(study_threshold-1), 100*(study_threshold-1), -1000),
          c(-100, -100, 100, 100), border = NA,
          col = adjustcolor('grey',.4))
  
  sort_ind = order(effect_ests[[plot_models[1]]][,'50%'])
  for(i in 1:length(plot_models)){
    points((exp(effect_ests[[plot_models[i]]][sort_ind,'50%'])-1)*100,
           1:K_treatments+index_p[i],pch=my_pch[i],
           col=mod_cols[i],cex=1.5)
    for(j in 1:length(sort_ind)){
      kk = sort_ind[j]
      lines((exp(c(effect_ests[[plot_models[i]]][kk,'2.5%'],
                   effect_ests[[plot_models[i]]][kk,'97.5%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=1)
      lines((exp(c(effect_ests[[plot_models[i]]][kk,'10%'],
                   effect_ests[[plot_models[i]]][kk,'90%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=3)
    }
  }
  axis(1, at = x_points)
  axis(2, at = 1:K_treatments, labels = rownames(effect_ests[[plot_models[1]]])[sort_ind])
  
}



plot_baseline_data = function(input_data){
  
  baseline_ind = input_data$Timepoint_ID==0
  bvl = aggregate(log10_viral_load ~ ID, input_data[baseline_ind, ], median)
  
  hist(bvl$log10_viral_load,
       breaks = seq(1,9,by=.5),
       xlab='Baseline viral load (SARS CoV2 genomes/mL)',
       ylab ='Number of patients',xlim=c(1,9),
       main='', xaxt ='n')
  axis(1, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  grid(); 
}


plot_serial_data = function(xx, xlims=c(0,7)){
  
  daily_VL_data = xx %>% group_by(ID, Time) %>%
    mutate(daily_VL = mean(log10_viral_load,na.rm = T))
  
  summary_VL_data = daily_VL_data %>% group_by(Timepoint_ID, Trt) %>%
    summarise(daily_VL = mean(daily_VL,na.rm = T),
              trt_color=unique(trt_color))
  
  summary_dat = daily_VL_data %>% ungroup() %>% distinct(ID, .keep_all = T) %>%
    group_by(Trt) %>%
    summarise(n = n(),
              trt_color = unique(trt_color)) %>% ungroup() %>%
    mutate(legend = paste0(Trt, ' (n=',n,')'))
  
  par(las=1)
  plot(summary_VL_data$Timepoint_ID, summary_VL_data$daily_VL,
       ylab = 'SARS CoV2 genomes/mL', panel.first=grid(),
       xlab = 'Time since randomization (days)',
       xlim = xlims, yaxt='n',type='n',
       ylim = c(0.7, 8))
  axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  
  points(daily_VL_data$Time, daily_VL_data$daily_VL,
         col= adjustcolor(daily_VL_data$trt_color,.3),
         pch = 1)
  for(tt in unique(summary_VL_data$Trt)){
    ind = summary_VL_data$Trt==tt
    lines(summary_VL_data$Timepoint_ID[ind],
          summary_VL_data$daily_VL[ind], pch=15,
          type='b',col = summary_VL_data$trt_color[ind],lwd=3)
  }
  
  legend('topright', col = summary_dat$trt_color, 
         legend = summary_dat$legend,
         lwd=2, pch=15, cex=1, inset = 0.03)
}


bayes_R2 = function(mod_preds, mod_residuals) {
  var_pred = apply(mod_preds, 1, var)
  var_res = apply(mod_residuals, 1, var)
  var_pred / (var_pred + var_res)
}



make_stan_inputs = function(input_data_fit, 
                            int_covs_base,
                            int_covs_full,
                            slope_covs_base,
                            slope_covs_full,
                            trt_frmla,
                            epoch=NA,
                            Dmax
){
  
  ## check censored values come last
  if(!all(diff(input_data_fit$log10_viral_load == input_data_fit$log10_cens_vl)>=0)) stop()
  ind_dup = !duplicated(input_data_fit$ID)
  
  for(ll in unique(input_data_fit$Lab)){
    ind = input_data_fit$Lab==ll
    med_val = median(input_data_fit$CT_RNaseP[ind], na.rm = T)
    input_data_fit$CT_RNaseP[ind & is.na(input_data_fit$CT_RNaseP)]=med_val
    input_data_fit$CT_RNaseP[ind] = input_data_fit$CT_RNaseP[ind]-med_val
  }
  input_data_fit$RnaseP_scaled = input_data_fit$CT_RNaseP
  
  input_data_fit$Age_scaled = (input_data_fit$Age-mean(input_data_fit$Age[ind_dup]))/sd(input_data_fit$Age[ind_dup])
  
  # make the covariate matrix
  # check no missing data
  if(!all(!apply(input_data_fit[, union(int_covs_full,slope_covs_full), drop=F], 2, function(x) any(is.na(x))))){
    stop('Missing data in covariate matrix!')
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_base,drop=F], 2, 
                          function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_1 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_2 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_base,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_1 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_2 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  if(!nrow(X_intcpt_1) == nrow(input_data_fit)) stop()
  if(!nrow(X_slope_1) == nrow(input_data_fit)) stop()
  
  cov_matrices = list(X_int=list(X_intcpt_1, X_intcpt_2),
                      X_slope=list(X_slope_1, X_slope_2))
  
  ID_map = data.frame(ID_key = input_data_fit$ID,
                      ID_stan = as.numeric(as.factor(input_data_fit$ID)))
  writeLines(sprintf('There are a total of %s patients in the database with a total of %s PCRs analysable',
                     max(ID_map$ID_stan),
                     nrow(input_data_fit)))
  
  ind_cens = !input_data_fit$log10_viral_load>
    input_data_fit$log10_cens_vl
  
  writeLines(sprintf('%s%% of samples are below LOD',
                     round(100*mean(ind_cens),digits = 2)))
  
  analysis_data_stan = list(Ntot = nrow(input_data_fit),
                            N_obs = sum(!ind_cens),
                            n_id = max(ID_map$ID_stan),
                            id = ID_map$ID_stan,
                            ind_start = which(!duplicated(ID_map$ID_stan)),
                            obs_day = input_data_fit$Time,
                            log_10_vl = input_data_fit$log10_viral_load,
                            log10_cens_vl = input_data_fit$log10_cens_vl,
                            RNaseP = input_data_fit$RnaseP_scaled,
                            Time_max = Dmax)
  if(!is.na(epoch)){
    analysis_data_stan$epoch = as.numeric(as.factor(input_data_fit$Epoch))
    analysis_data_stan$K_epoch = max(analysis_data_stan$epoch)
  }
  ID_map = ID_map[!duplicated(ID_map$ID_key), ]
  
  writeLines('check stan data formatting:')
  all(analysis_data_stan$log_10_vl[1:analysis_data_stan$N_obs]>
        analysis_data_stan$log10_cens_vl[1:analysis_data_stan$N_obs]) &
    all(analysis_data_stan$log_10_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot] ==
          analysis_data_stan$log10_cens_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot])
  
  Trt_matrix = model.matrix(trt_frmla, data = input_data_fit)
  Trt_matrix[,1]=0 # first column is dummy
  
  analysis_inputs = list(cov_matrices=cov_matrices,
                         analysis_data_stan=analysis_data_stan,
                         Trt_matrix=Trt_matrix,
                         ID_map=ID_map)
  return(analysis_inputs)
}






get_rates_linear_mod = function(mod_out, # single model fit - not a list
                                analysis_data_stan){
  
  # get the indices of first datapoint for each individual
  ind_id = which(!duplicated(analysis_data_stan$id))
  thetas = 
    rstan::extract(mod_out, 
                   pars = c('beta_0','beta_cov','theta_rand_id','trt_effect'))
  beta_cov = x_slope*slope_coefs;
  
  beta_0*exp(trt_slope[i]+theta_rand_id[id[i]][2]+beta_cov[i])
}

# plot_data_model_fits = function(mod_out, # model fits
#                              models_plot, # which models to plot
#                              K_plots,
#                              mod_cols,
#                              ID_map,
#                              analysis_data_stan
# ){
#   
#   # extract posterior parameters and outputs
#   thetas = list()
#   for(mm in 1:length(mod_out)){
#     thetas[[mm]] = rstan::extract(mod_out[[mm]])
#   }
#   counter = 1
#   
#   ID_map$Trt = gsub(pattern = '\n',
#                     replacement = '',
#                     x = ID_map$Trt,fixed = T)
#   while(counter <= nrow(ID_map)){
#     
#     # draw individual model fit with data
#     ind = analysis_data_stan$id==ID_map$ID_stan[counter]
#     plot(analysis_data_stan$obs_day[ind],
#          analysis_data_stan$log_10_vl[ind],
#          xlab='', ylab='', 
#          xaxt='n', yaxt='n',
#          panel.first=grid(), xlim=c(0,7),
#          ylim = range(analysis_data_stan$log_10_vl))
#     # if(counter %% sqrt(K_plots) == 1){
#     #   mtext(text = 'RNA copies per mL',side = 2,
#     #         line = 3,las = 3)
#     # }
#     axis(1, at = c(0,3,7))
#     axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
#                                         expression(10^4),
#                                         expression(10^6),
#                                         expression(10^8)))
#     # if((counter%%K_plots) >= K_plots - sqrt(K_plots)){
#     #   mtext(text = 'Days',side = 1,line = 2)
#     # }
#     for(mm in models_plot){
#       ix = order(analysis_data_stan$obs_day[ind])
#       my_xs = analysis_data_stan$obs_day[ind][ix]
#       polygon(x = c(my_xs, rev(my_xs)),
#               y = c(apply(thetas[[mm]]$preds[,ind],2,
#                           quantile,probs=0.025)[ix],
#                     rev(apply(thetas[[mm]]$preds[,ind],2,
#                               quantile,probs=0.975)[ix])),
#               border = NA, 
#               col = adjustcolor(mod_cols[mm],alpha.f = .3))
#       lines(my_xs,
#             colMeans(thetas[[mm]]$preds[,ind])[ix],
#             col = mod_cols[mm],lwd=2)
#     }
#     points(analysis_data_stan$obs_day[ind],
#            analysis_data_stan$log_10_vl[ind],pch=16)
#     
#     mtext(text = paste0(ID_map$ID[counter],
#                         '\n',
#                         ID_map$Trt[counter]),
#           side = 3, line = -0.5, cex=0.8)
#     counter=counter+1
#   }
#   
# }

plot_data_model_fits = 
  function(model_list, # list of model fits
           models_to_plot, # which models to plot
           K_plots,
           mod_cols,
           ID_map,
           analysis_data_stan
  ){
    
    # extract posterior parameters and outputs
    thetas = list()
    for(mm in 1:length(model_list)){
      thetas[[mm]] = rstan::extract(model_list[[mm]], pars='preds')
    }
    counter = 1
    
    ID_map$Trt = gsub(pattern = '\n',
                      replacement = '',
                      x = ID_map$Trt,fixed = T)
    while(counter <= nrow(ID_map)){
      
      # draw individual model fit with data
      ind = analysis_data_stan$id==ID_map$ID_stan[counter]
      plot(analysis_data_stan$obs_day[ind],
           analysis_data_stan$log_10_vl[ind],
           xlab='', ylab='', 
           xaxt='n', yaxt='n',
           panel.first=grid(), xlim=c(0,7),
           ylim = range(analysis_data_stan$log_10_vl))
      
      axis(1, at = c(0,3,7))
      axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                          expression(10^4),
                                          expression(10^6),
                                          expression(10^8)))
      
      for(mm in models_to_plot){
        ix = order(analysis_data_stan$obs_day[ind])
        my_xs = analysis_data_stan$obs_day[ind][ix]
        polygon(x = c(my_xs, rev(my_xs)),
                y = c(apply(thetas[[mm]]$preds[,ind],2,
                            quantile,probs=0.025)[ix],
                      rev(apply(thetas[[mm]]$preds[,ind],2,
                                quantile,probs=0.975)[ix])),
                border = NA, 
                col = adjustcolor(mod_cols[mm],alpha.f = .3))
        lines(my_xs,
              colMeans(thetas[[mm]]$preds[,ind])[ix],
              col = mod_cols[mm],lwd=2)
      }
      points(analysis_data_stan$obs_day[ind],
             analysis_data_stan$log_10_vl[ind],pch=16)
      
      mtext(text = paste0(ID_map$ID_key[counter],
                          ' ',
                          ID_map$Trt[counter]),
            side = 3, line = 0.5, cex=0.7)
      counter=counter+1
    }
    
  }



plot_individ_curves = function(platcov_dat, IDs, xlims){
  
  if(!all(IDs %in% platcov_dat$ID)) stop('missing IDs!!')
  ylims = range(platcov_dat$log10_viral_load)
  platcov_dat$Day = floor(platcov_dat$Time)
  
  for(id in unique(IDs)){
    
    ind = platcov_dat$ID==id
    plot(platcov_dat$Time[ind],
         platcov_dat$log10_viral_load[ind],
         pch = as.numeric(platcov_dat$log10_cens_vl[ind]==platcov_dat$log10_viral_load[ind])+1,
         xlab='', ylab='', 
         xaxt='n', yaxt='n',
         panel.first=grid(), xlim=xlims,
         ylim = ylims)
    title(id)
    points(platcov_dat$Time[ind],
           platcov_dat$log10_viral_load[ind],
           pch = as.numeric(platcov_dat$log10_cens_vl[ind]==platcov_dat$log10_viral_load[ind])+16)
    axis(1, at = c(0,7,14,21))
    axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                        expression(10^4),
                                        expression(10^6),
                                        expression(10^8)))
    
  }
  legend('topright', legend = c('>LLOQ','<LLOQ'), pch = 16:17, inset = 0.03)
}


plot_coef_effects = function(stan_out, cov_mat, stan_inputs){
  
  thetas = rstan::extract(stan_out)
  alpha_coefs = apply(thetas$intercept_coefs,2,
                      quantile,probs=c(0.025,.1,.5,.9,0.975))
  xlims=range(alpha_coefs)
  
  cov_names_intercept =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_int[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  cov_names_slope =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_slope[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  plot(alpha_coefs['50%', ], 1:ncol(alpha_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Intercept (fold change)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(alpha_coefs)){
    lines(c(alpha_coefs['10%',i], alpha_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(alpha_coefs['2.5%',i], alpha_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(alpha_coefs), labels = cov_names_intercept,tick = F)
  
  x_points = signif(10^seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = log10(x_points), labels = x_points)
  
  beta_coefs = apply(thetas$slope_coefs,2,quantile,
                     probs=c(0.025,.1,.5,.9,0.975))
  xlims=range(beta_coefs)
  plot(beta_coefs['50%', ], 1:ncol(beta_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Slope (multiplicative effect)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(beta_coefs)){
    lines(c(beta_coefs['10%',i], beta_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(beta_coefs['2.5%',i], beta_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(beta_coefs), labels = cov_names_slope, tick = F)
  x_points = signif(10^seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = log10(x_points), labels = x_points)
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

calculate_fever_clearance = function(temp_dat,
                                     window_clear = 24/24, # look ahead window to define "fever clearance"
                                     threshold=37){
  
  if(!'ax_temperature' %in% colnames(temp_dat)) stop('needs to contain a ax_temperature column')
  
  temp_dat$clearance_time = NA
  # For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. 
  temp_dat$clearance_time_cens = 1
  
  temp_dat$fever_binary = temp_dat$ax_temperature>threshold
  temp_dat = dplyr::arrange(temp_dat, ID, Time) 
  temp_dat = temp_dat[!is.na(temp_dat$ax_temperature), ]
  
  for(id in unique(temp_dat$ID)){
    ind = temp_dat$ID==id
    if(all(!temp_dat$fever_binary[ind])){ # never fever
      temp_dat$clearance_time[ind]=0
    } else if(all(temp_dat$fever_binary[ind])){ # always fever
      writeLines(sprintf('all fever for %s with %s FUP points',id,sum(ind)))
      temp_dat$clearance_time[ind] = max(temp_dat$Time[ind])
      temp_dat$clearance_time_cens[ind] = 0 #censored obs
    } else { # fever cleared
      j_cleared = which(ind & !temp_dat$fever_binary)
      check_ahead=F
      for(j in j_cleared){
        if(!check_ahead){
          ind_check = 
            which(ind & 
                    temp_dat$Time>temp_dat$Time[j] &
                    temp_dat$Time<temp_dat$Time[j] + window_clear)
          if(length(ind_check)>0 & all(!temp_dat$fever_binary[ind_check])){
            temp_dat$clearance_time[ind]=temp_dat$Time[j]
            check_ahead=T
          }
        }
      }
      if(!check_ahead){
        temp_dat$clearance_time[ind]=tail(temp_dat$Time[ind],1)
        temp_dat$clearance_time_cens[ind]=0
      }
    }
  }
  
  return(temp_dat[!duplicated(temp_dat$ID), ])
}



make_slopes_plot = function(stan_out, 
                            analysis_data_stan,
                            ID_map, 
                            data_summary,
                            my_lims = c(5, 72), # hours
                            my_vals = c(7,24,48,72),
                            qq_show = 0.5){
  
  slopes = rstan::extract(stan_out, pars='slope')$slope
  
  t12_output = data.frame(t_12_med = 24*log10(.5)/apply(slopes,2,mean),
                          t_12_up = 24*log10(.5)/apply(slopes,2,quantile,.9),
                          t_12_low = 24*log10(.5)/apply(slopes,2,quantile,.1),
                          slope_median = apply(slopes,2,median),
                          ID_stan = analysis_data_stan$id[analysis_data_stan$ind_start])
  t12_output = merge(t12_output, ID_map, by = 'ID_stan')
  data_summary = merge(data_summary, t12_output, by.x = 'ID', by.y = 'ID_key')
  
  data_summary = dplyr::arrange(data_summary, Trt, t_12_med)
  
  par(mar=c(5,2,2,2))
  plot(data_summary$t_12_med, 1:nrow(data_summary),
       yaxt='n', xaxt='n',
       xlim=my_lims, panel.first=grid(), xlab='', 
       pch=15, ylab='', col=data_summary$trt_color)
  mtext(text = expression('t'[1/2] ~ ' (hours)'),side = 1,line=3)
  axis(1, at = my_vals,labels = my_vals)
  
  for(i in 1:nrow(data_summary)){
    lines(c(data_summary$t_12_low[i],
            data_summary$t_12_up[i]),
          rep(i,2),
          col=adjustcolor(data_summary$trt_color[i],alpha.f = .5))
  }
  
  for(kk in which(!duplicated(data_summary$Trt))){
    ind = data_summary$Trt==data_summary$Trt[kk]
    writeLines(sprintf('In %s the median clearance half life was %s (range %s to %s)',
                       data_summary$Trt[kk],
                       round(median(data_summary$t_12_med[ind]),1),
                       round(min(data_summary$t_12_med[ind]),1),
                       round(max(data_summary$t_12_med[ind]),1)))
    abline(v = quantile(data_summary$t_12_med[ind], probs=qq_show), 
           col=data_summary$trt_color[kk],lty=2,lwd=2)
  }
  
  legend('bottomright', legend = unique(data_summary$Trt),
         col = unique(data_summary$trt_color),
         pch=15,lwd=2,inset=0.03)
  return(data_summary)
}

get_itt_population = function(prefix_drop_rand){
  
  require(tidyverse)
  # TH58 and TH57 used envelopes so allocation is not recorded in the data-XXX.csv files
  
  rand.TH58 <- read.csv(paste0(prefix_drop_rand, "/rand-TH58.csv"))[0:9, ]
  rand.TH57 <- read.csv(paste0(prefix_drop_rand, "/rand-TH57.csv"))[0:10, ]
  rand.TH58$ID = paste('PLT-TH58-',rand.TH58$RandomisationID,sep='')
  rand.TH57$ID = paste('PLT-TH57-',rand.TH57$RandomisationID,sep='')
  
  ff_names = list.files(path = prefix_drop_rand, pattern = 'data',full.names = T)
  
  ff_names <- ff_names[-3] # Move back when NP has some data
  
  data_list = list()
  for(i in 1:length(ff_names)){
    data_list[[i]] = read.csv(ff_names[i])
    data_list[[i]]$Date = as.POSIXct(data_list[[i]]$Date,format='%a %b %d %H:%M:%S %Y')
    my_prefix=gsub(x = gsub(x = strsplit(ff_names[i], split = 'data-')[[1]][2], pattern = '.csv',replacement = ''),pattern = '0',replacement = '')
    data_list[[i]]$ID = paste('PLT-', my_prefix, '-', data_list[[i]]$randomizationID, sep='')
    data_list[[i]] = data_list[[i]][, c('ID', 'Treatment')]
  }
  data_list[[length(ff_names)+1]]=rand.TH57[, c('ID', 'Treatment')]
  data_list[[length(ff_names)+2]]=rand.TH58[, c('ID', 'Treatment')]
  
  xx = bind_rows(data_list)
  
  library(stringr)
  for(i in 1:nrow(xx)){
    id = unlist(strsplit(xx$ID[i],split = '-'))
    id[3] = str_pad(id[3], 3, pad = "0")
    id = paste(id, collapse = '-')
    xx$ID[i]=id
  }
  
  return(xx)
}




get_trt_colors = function(){
  trt_cols = array(dim = 15)
  names(trt_cols) = 
    c("Ivermectin",
      "Regeneron",
      'No study drug',
      "Remdesivir",
      "Favipiravir",
      "Nitazoxanide",           
      "Fluoxetine",
      "Molnupiravir",
      "Nirmatrelvir+Molnupiravir",
      "Evusheld",
      'Ensitrelvir',
      "Nirmatrelvir",
      "Hydroxychloroquine",
      "Metformin",
      "Atilotrelvir/ritonavir")
  trt_cols['No study drug'] = viridis::viridis(n = 10)[8]
  trt_cols['Fluoxetine'] = viridis::viridis(n = 10)[5]
  trt_cols['Nitazoxanide'] = viridis::magma(n = 10)[8]
  trt_cols['Evusheld'] = viridis::magma(n = 10)[1]
  trt_cols['Favipiravir'] = viridis::plasma(n = 100)[92]
  trt_cols['Ivermectin'] = viridis::plasma(n = 10)[4]
  trt_cols['Nirmatrelvir+Molnupiravir'] = viridis::plasma(n = 10)[3]
  trt_cols['Nirmatrelvir'] = viridis::plasma(n = 10)[1]
  trt_cols['Regeneron'] = viridis::inferno(n = 10)[5]
  trt_cols['Molnupiravir'] = viridis::inferno(n = 10)[7]
  trt_cols['Remdesivir'] = RColorBrewer::brewer.pal('Dark2',n=8)[8]
  trt_cols['Ensitrelvir'] = RColorBrewer::brewer.pal('Set1',n=8)[1]
  trt_cols['Hydroxychloroquine'] = RColorBrewer::brewer.pal('Paired',n=8)[1]
  trt_cols['Metformin'] = "#FF9F00"
  trt_cols['Atilotrelvir/ritonavir'] =  "#CD5656" 
  
  return(trt_cols)
}

plot_trt_colors = function(){
  trt_cols= get_trt_colors()
  plot(NA, NA, xaxt='n', xlab='', yaxt='n', ylab='',xlim=0:1, ylim=0:1,bty='n')
  legend('right', legend = names(trt_cols), fill = trt_cols,border = NA,cex=1.5)
}


plot_vl <- function(dataplot, trt_colors){
  dataplot <- dataplot %>%
    distinct(ID, Timepoint_ID, daily_VL, .keep_all = T)
  
  dataplot_median <- dataplot %>%
    group_by(Trt, Timepoint_ID) %>%
    summarise(median_VL = median(daily_VL), .groups = 'drop')
  
  colors <- trt_colors[names(trt_colors) %in% unique(dataplot$Trt)]
  colors <- colors[levels(dataplot$Trt)]
  labels <- names(colors)
  labels[labels == "Nirmatrelvir"] <- "Ritonavir-boosted nirmatrelvir"
  
  f_tab <- dataplot %>%
    distinct(ID, Trt) %>%
    group_by(Trt) %>%
    summarise(n = n()) %>%
    as.data.frame()
  
  f_tab$lab <- paste0(f_tab$Trt, ": ", f_tab$n)
  f_text <- paste(f_tab$lab, collapse = '\n')
  
  
  G <- ggplot(dataplot, aes(x = Timepoint_ID, y = log10_viral_load, col = Trt)) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.75, aes(shape = censor) ) +
    geom_line(data = dataplot_median, aes(x =  Timepoint_ID, y = median_VL, group = Trt, col = Trt), linewidth = 1, linetype = 1) +
    geom_point(data = dataplot_median, aes(x = Timepoint_ID, y = median_VL, fill = Trt)
               , size = 3.5, shape = 24, col = "black") +
    scale_shape_manual(values = c(25, 21), guide = NULL) +
    scale_color_manual(label = labels, values = colors, name = "") +
    scale_fill_manual(label = labels, values = colors, name = "") +
    annotate("text", x = 2.5, y = 9, label = f_text, hjust = 0, vjust = 1)  +
    theme_bw() +
    scale_x_continuous(breaks = 0:14) +
    scale_y_continuous(labels=label_math(), breaks = seq(0,10,2), limits = c(0,9)) +
    xlab("Time since randomisation (days)") +
    ylab("SARS-CoV-2 genomes/mL") + 
    theme(axis.title  = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text = element_text(size = 10)) +
    ggtitle("\nA) Viral load dynamics")
  
  G
}

plot_vl_box <- function(dataplot, trt_colors){
  dataplot$Timepoint_ID <- as.factor(dataplot$Timepoint_ID)
  dataplot <- dataplot %>%
    distinct(ID, Timepoint_ID, daily_VL, .keep_all = T)
  
  dataplot_median <- dataplot %>%
    group_by(Trt, Timepoint_ID) %>%
    summarise(median_VL = median(daily_VL), .groups = 'drop')
  
  f_tab <- dataplot %>%
    distinct(ID, Trt) %>%
    group_by(Trt) %>%
    summarise(n = n()) %>%
    as.data.frame()
  f_tab$lab <- paste0("n = ", f_tab$n)
  
  colors <- trt_colors[names(trt_colors) %in% unique(dataplot$Trt)]
  colors <- colors[levels(dataplot$Trt)]
  labels <- names(colors)
  labels[labels == "Nirmatrelvir"] <- "Ritonavir-boosted nirmatrelvir"
  
  G <- ggplot(dataplot, aes(x = Timepoint_ID, y = daily_VL, fill = Trt)) +
    geom_jitter(width = 0.2, alpha = 0.25, size = 1.5, aes(shape = censor, fill = Trt),
                stroke = 0.5) +
    geom_boxplot(width=0.65, size = 0.5, outlier.shape = NA,  coef = 0, aes(fill = Trt), alpha = 0.3) +
    geom_line(data = dataplot_median, aes(x =  Timepoint_ID, y = median_VL, group = Trt, col = Trt), 
              linewidth = 1.2, linetype = 1) +
    geom_point(data = dataplot_median, aes(x = Timepoint_ID, y = median_VL, fill = Trt)
               , size = 3.5, shape = 24, col = "black") +
    scale_shape_manual(values = c(25, 21), guide = NULL) +
    scale_color_manual(label = labels, values = colors, name = "") +
    scale_fill_manual(label = labels, values = colors, name = "") +
    facet_grid(.~Trt) +
    theme_bw(base_size = 13) +
    scale_y_continuous(labels=label_math(), breaks = seq(0,10,1), limits = c(0,9)) +
    xlab("Time since randomisation (days)") +
    ylab("SARS-CoV-2 genomes/mL") + 
    theme(axis.title  = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "none",
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13, face = 'bold')) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.75) +
    geom_text(data = f_tab, x = 4, y = 9, aes(label = lab),
              hjust = 0, vjust = 1) +
    ggtitle("Viral density dynamics")
  G
}

slope_to_hl  <- function(slope){
  24*log10(2)/(-(slope)) 
}

formatter <- function(x){  
  (x-1)*100 
}

plot_trt_effs <- function(effect_ests, model_cols, study_threshold, unblinded_arms, intervention){
  effect_ests_plots <- NULL
    for(i in 1:length(effect_ests)){
    effect_ests_plot <- as.data.frame(effect_ests[[i]])
    effect_ests_plot <- exp(effect_ests_plot)
    colnames(effect_ests_plot)[1:5] <- c("L95", "L80", "med", "U80", "U95")
    effect_ests_plot <- effect_ests_plot %>% arrange(med)
    effect_ests_plot$arm <- row.names(effect_ests_plot)
    if(intervention == "interim_all"){effect_ests_plot <- effect_ests_plot %>% filter(!arm %in% unblinded_arms)}
    effect_ests_plot$arm <- as.factor(effect_ests_plot$arm)
    effect_ests_plot$arm <- factor(effect_ests_plot$arm, levels = as.character(effect_ests_plot$arm))
    effect_ests_plot$model <- names(effect_ests)[i]
    effect_ests_plots <- rbind(effect_ests_plots, effect_ests_plot)
  }
  
  effect_ests_plots$model <- as.factor(effect_ests_plots$model)

  #Labeling reference arm
  lab_ref <- ref_arm
  lab_ref[lab_ref == "Nirmatrelvir"] <- "Ritonavir-boosted nirmatrelvir"
  #Labeling intervention arm
  
  my.labs <- levels(effect_ests_plot$arm)
  my.labs[my.labs == "Nirmatrelvir+Molnupiravir"] <- "Nirmatrelvir/Molnupiravir"
  my.labs[my.labs == "Regeneron"] <- "Casirivimab/imdevimab"
  my.labs[my.labs == "Evusheld"] <- "Tixagevimab/cilgavimab"
  
  title <- paste0("B) Estimated treatment effects \nrelative to ", tolower(lab_ref), " arm")
  
  names(model_cols) <- levels(effect_ests_plots$model)
  
  G <- ggplot(effect_ests_plots, 
         aes(x = arm, y = med, col = model, group = model)) +
    geom_rect(aes(ymin = min(0.75, min(L95)-0.05), 
                  ymax = study_threshold, xmin = 0, xmax = length(my.labs)+1), fill = "#7D7C7C", alpha = 0.1, col = NA) +
    geom_point(position = position_dodge(width = 0.5), size = 4) +
    geom_errorbar(aes(x = arm, ymin = L95, ymax = U95),position = position_dodge(width = 0.5), width = 0, linewidth = 0.65) +
    geom_errorbar(aes(x = arm, ymin = L80, ymax = U80),position = position_dodge(width = 0.5), width = 0, linewidth = 1.5) +
    scale_color_manual(values = model_cols) +
    coord_flip() +
    theme_bw() +
    geom_hline(yintercept = 1, col = "red", linetype = "dashed") +
    scale_y_continuous(labels = formatter, 
                       limits = c(min(0.75, min(effect_ests_plots$L95)-0.05), max( effect_ests_plots$U95) + .25), expand = c(0,0),
                       breaks = seq(0.2,3.6, 0.2)) +
    scale_x_discrete(labels= my.labs) +
    ylab("Change in viral clearance rate (%)") +
    xlab("") +
    ggtitle(title)  + 
    theme(axis.title  = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text = element_text(size = 10))
  G
}



plot_hl <- function(Half_life, trt_colors){
  Half_life_med <- Half_life %>%
    group_by(Trt) %>%
    summarise(med_hl = median(t_12_med))  %>%
    as.data.frame()
  
  colors <- trt_colors[names(trt_colors) %in% unique(Half_life$Trt)]
  colors <- colors[levels(Half_life$Trt)]
  labels <- names(colors)
  labels[labels == "Nirmatrelvir"] <- "Ritonavir-boosted nirmatrelvir"
  
  f_tab <- Half_life %>%
    distinct(ID, Trt) %>%
    group_by(Trt) %>%
    summarise(n = n()) %>%
    as.data.frame()
  f_tab$med_hl <- Half_life_med$med_hl
  
  f_tab$lab <- paste0(f_tab$Trt,": n = ", f_tab$n, "; t_1/2 = ", round(f_tab$med_hl,1), " h")
  freq_lab <- paste(f_tab$lab, collapse = '\n')
  
  
  G <- ggplot(Half_life, aes(x = t_12_med, y = ID, col = Trt)) +
    geom_errorbar(aes(xmin = t_12_low, xmax = t_12_up),width = 0, alpha = 0.4) +
    geom_point(size = 2) +
    geom_vline(data = Half_life_med, aes(xintercept = med_hl, col = Trt),linewidth = 1) +
    theme_bw() +  
    scale_y_discrete(expand = c(0.01,0.01), breaks = NULL) +
    scale_color_manual(label = labels, values = colors, name = "") +
    scale_x_continuous(breaks = seq(0,40,5), expand = c(0,0)) +
    guides(color = guide_legend(override.aes=list(linetype = rep(0, length(unique(Half_life$Trt)))))) +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text( face = "bold"),
          legend.position = "bottom") +
    coord_cartesian(xlim=c(0, 35)) +
    xlab("Estimated viral clearance half-life (h)") +
    ylab("") +
    ggtitle("A) Individual viral clearance half-life\n") +
    annotate("text", x = 13, y = nrow(Half_life)/6, label = freq_lab, hjust = 0, vjust = 1, size = 3) 
  G
  
  
}

plot_randomisation <- function(Baseline_data){
  Baseline_data$Rand_date <- as.Date(Baseline_data$Rand_date)
  
  ggplot(Baseline_data, aes(x = Rand_date, y = Trt, col = Country)) +
    geom_jitter(size = 3, alpha = 0.5, height = 0.1) +
    theme_bw() +
    xlab("Randomisation date") +
    ylab("") +
    scale_x_date(date_breaks = "3 months", date_minor_breaks = "1 months",
                 date_labels = "%b %y") +
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 11))
  
  
}


find_rebounds = function(platcov_dat, 
                         lower_bound=2,  # lower level such that VL is defined as non-detectable
                         upper_bound=4,  # upper level such that VL is defined as "high"
                         t_window=2  ,    # time window during which it has to be undetectable
                         day_min = 4
){
  platcov_dat$rebound=NA
  platcov_dat$Timepoint_ID_rebound=NA
  platcov_dat$Time_rebound=NA
  platcov_dat$logVL_rebound=NA
  for(id in unique(platcov_dat$ID)){
    ind=platcov_dat$ID==id
    res = assess_rebound(platcov_dat[ind,],
                         lower_bound = lower_bound,
                         upper_bound = upper_bound,
                         t_window = t_window,
                         day_min = day_min)
    platcov_dat$rebound[ind] = res$rebound
    platcov_dat$Timepoint_ID_rebound[ind] = res$details$Timepoint_ID_rebound
    platcov_dat$Time_rebound[ind] = res$details$Time_rebound
    platcov_dat$logVL_rebound[ind] = res$details$logVL_rebound
    
  }
  return(platcov_dat)
}

assess_rebound = function(patient_dat,
                          lower_bound=2,  # lower level such that VL is defined as non-detectable
                          upper_bound=3,  # upper level such that VL is defined as "high"
                          t_window=1.5,   # time window during which it has to be undetectable
                          day_min = 5
){
  xx = patient_dat %>% arrange(Time) %>% distinct(Timepoint_ID, .keep_all = T)
  rebound = virus_cleared = F
  details = data.frame("Timepoint_ID_rebound" = NA, "Time_rebound" = NA,
                       "logVL_rebound" = NA)
  if(nrow(xx)>3){
    for(i in 2:nrow(xx)){
      ind = which(xx$Time <= xx$Time[i] & (xx$Time >= (xx$Time[i]-t_window)))
      if(all(xx$daily_VL[ind] <= lower_bound)){
        virus_cleared=T
      }
      #print(virus_cleared)
      if((virus_cleared & xx$daily_VL[i] >= upper_bound & xx$Time[i]>day_min & xx$daily_VL[i-1] < upper_bound)|
         (virus_cleared & xx$Time[i] >= 8 & xx$daily_VL[i]  >= upper_bound)){
        rebound = T
        writeLines(sprintf('patient %s treated with %s had a rebound identified on day %s (%s days)', 
                           xx$ID[1], xx$Trt[1],xx$Timepoint_ID[i], round(xx$Time[i],1)))
        details = data.frame("Timepoint_ID_rebound" = xx$Timepoint_ID[i], "Time_rebound" = xx$Time[i],
                             "logVL_rebound" = xx$daily_VL[i])
      }
      # print(rebound)
    }
  }
  return(list("rebound" = rebound, "details" = details))
}


plot_rebound <- function(platcov_dat_analysis, platcov_dat_rebound, trt_colors){
  
  lab <- platcov_dat_analysis %>%
    distinct(ID, .keep_all = T) %>%
    group_by(Trt) %>%
    summarise(n = n()) 
  
  lab$n_rbound <- platcov_dat_rebound %>%
    filter(rebound) %>%
    distinct(ID, .keep_all = T) %>%
    group_by(Trt) %>%
    summarise(n = n()) %>%
    pull(n)
  
  lab$lab <- paste0(lab$Trt, " (n=", lab$n_rbound, "/", lab$n, ")")
  
  lab_labeller <- function(variable,value){
    return(lab$lab[value])
  }
  
  
  rebound_summary <-  platcov_dat_rebound  %>% filter(rebound) %>%
    group_by(Trt) %>%
    summarise(med_time = median(Time_rebound) %>% round(1) %>% format(nsmall = 1),
              Q1_time = quantile(Time_rebound, 0.25)%>% round(1) %>% format(nsmall = 1),
              Q3_time = quantile(Time_rebound, 0.75)%>% round(1) %>% format(nsmall = 1),
              med_vl_rebound = median(logVL_rebound)%>% round(1) %>% format(nsmall = 1),
              Q1_vl_rebound = quantile(logVL_rebound, 0.25)%>% round(1) %>% format(nsmall = 1),
              Q3_vl_rebound = quantile(logVL_rebound, 0.75)%>% round(1) %>% format(nsmall = 1),
              lab1 = paste0("Time-to-rebound (days): ",  med_time, " [IQR: ", Q1_time, " to ",  Q3_time, "]"),
              lab2 = paste0("Rebound viral loads (log genomes/mL): ",  med_vl_rebound, " [IQR: ", Q1_vl_rebound, " to ",  Q3_vl_rebound, "]"),
    )
  
  
  ggplot(platcov_dat_rebound  %>% filter(rebound), aes(x = Time, y = daily_VL)) +
    geom_point(shape = 21, aes(group = ID, col = Trt)) +
    geom_line(alpha = 0.5,  aes(group = ID, col = Trt)) +
    facet_wrap(Trt~., ncol = 1,
               labeller  = lab_labeller) +
    theme_bw(base_size = 13)+
    scale_color_manual(values = trt_colors) +
    geom_hline(yintercept = c(2,3), linetype = "dashed", linewidth = 0.4, alpha = 0.75) +
    scale_y_continuous(labels=label_math(), breaks = seq(0,8,2), limits = c(0,8)) +
    scale_x_continuous(limits = c(0,15), breaks = seq(0,15,3)) +
    xlab('Time since randomisation (days)') +
    ylab('SARS-CoV-2 genomes/mL') +
    theme(axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "none",
          panel.spacing = unit(1, "lines")) +
    geom_text(data = rebound_summary, x = 5, y = 7.5, aes(label = lab1), hjust = 0, size = 3) +
    geom_text(data = rebound_summary, x = 5, y = 6.5, aes(label = lab2), hjust = 0, size = 3)
  
  
}


checkStrict(make_stan_inputs)
checkStrict(plot_serial_data)
checkStrict(plot_effect_estimates)
checkStrict(plot_data_model_fits)
checkStrict(plot_coef_effects)
checkStrict(calculate_fever_clearance)
