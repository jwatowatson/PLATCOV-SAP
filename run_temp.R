rm(list=ls())
library(rstan)
options(mc.cores = parallel::detectCores())

platcov_dat = rbind(
  read.csv('REGN_analysis.csv'),
  read.csv('Remdesivir_analysis.csv'),
  read.csv('Paxlovid_Molnupiravir_analysis.csv'))

platcov_dat = platcov_dat[!duplicated(platcov_dat$BARCODE),]
range(platcov_dat$Time)
platcov_dat$Rand_date = as.POSIXct(platcov_dat$Rand_date)
platcov_dat = platcov_dat[platcov_dat$Rand_date < '2022-08-25', ]
table(platcov_dat$Site)

# make per protocol summary for all patients
PP=merge(aggregate(Timepoint_ID ~ ID+Trt,
                   platcov_dat[platcov_dat$Per_protocol_sample==1, ], max),
         aggregate(Per_protocol_sample ~ ID, platcov_dat, sum),by = 'ID')
PP$Include_mITT = PP$Timepoint_ID>=2


# We remove patients who only have undetectable virus
xx_undetectble = table(platcov_dat$ID, platcov_dat$CT_NS==40)
ids_neg = names(which(xx_undetectble[,1]==0))
writeLines(sprintf('All negative samples for id: %s', ids_neg))

# Exclude from mITT pop
PP$Include_mITT[PP$ID %in% ids_neg] = F

platcov_dat = platcov_dat[platcov_dat$ID %in% PP$ID[PP$Include_mITT], ]

platcov_dat$Trt[platcov_dat$Trt=='Regeneron'&platcov_dat$Variant=='Delta']='Regeneron_Delta'
platcov_dat$Trt[platcov_dat$Trt=='Regeneron'&platcov_dat$Variant!='Delta']='Regeneron_Omicron'

platcov_dat$Trt = factor(platcov_dat$Trt, levels=c('No study drug',
                                                   'Regeneron_Delta','Regeneron_Omicron',
                                                   'Nirmatrelvir + Ritonavir',
                                                   'Remdesivir',
                                                   'Molnupiravir'))
ind_dup = !duplicated(platcov_dat$ID)
table(platcov_dat$Trt[ind_dup])

platcov_dat$Variant = factor(platcov_dat$Variant, levels = c('Delta','BA.1','BA.2','BA.4','BA.5'))
table(platcov_dat$Variant[ind_dup])


covs_base = c('Variant','Site')
covs_full=c(covs_base, 'Age_scaled')


length(unique(platcov_dat$ID))

table(platcov_dat$Trt[ind_dup])

Dmax=8
source('functions.R')
platcov_dat = dplyr::arrange(platcov_dat, log10_viral_load==log10_cens_vl)
ind_fitting = platcov_dat$Time < Dmax
stan_inputs = 
  make_stan_inputs(input_data_fit = platcov_dat[ind_fitting,],
                   int_covs_base = covs_base,
                   int_covs_full = covs_full,
                   slope_covs_base = covs_base,
                   slope_covs_full = covs_full,
                   trt_frmla = formula('~ Trt'),
                   Dmax = Dmax)

colnames(stan_inputs$Trt_matrix)

analysis_data_stan = stan_inputs$analysis_data_stan
analysis_data_stan$trt_mat = stan_inputs$Trt_matrix
analysis_data_stan$K_trt = ncol(analysis_data_stan$trt_mat)

x_intercept = stan_inputs$cov_matrices$X_int[[1]]
if(ncol(x_intercept)==0) x_intercept = array(0, dim=c(nrow(x_intercept),1))
analysis_data_stan$x_intercept = x_intercept
analysis_data_stan$K_cov_intercept= ncol(x_intercept)


x_slope = stan_inputs$cov_matrices$X_slope[[1]]
if(ncol(x_slope)==0) x_slope = array(0, dim=c(nrow(x_slope),1))
analysis_data_stan$x_slope = x_slope
analysis_data_stan$K_cov_slope=ncol(x_slope)

mod = stan_model(file = 'Stan_models/Linear_model_RNaseP.stan') # compile 
source('priors.R')
out = sampling(mod, 
               data=c(analysis_data_stan,
                      all_priors[[1]]),
               iter=2000,
               chain=4,
               thin=4,
               warmup=1000,
               save_warmup = FALSE,
               seed=77,
               pars=c('L_Omega','theta_rand_id'), # we don't save these as it takes up vast memory!
               include=FALSE)

plot(out, pars='trt_effect')
thetas_trt = extract(out, pars='trt_effect')$trt_effect
colnames(stan_inputs$Trt_matrix)

names_xx = c('Casirivimab/imdevimab\nDelta',
             'Casirivimab/imdevimab\nOmicron',
             'Nirmatrelvir',
             'Remdesivir',
             'Molnupiravir'
             )
trt_matrix = apply(thetas_trt, 2, quantile, probs=c(0.025,.1,.5,.9,.975))

par(las=1, mar=c(5, 14, 2,2), family='serif',cex.lab=1.3, cex.axis=1.3, bty='n')
plot(trt_matrix[3, ], 1:ncol(trt_matrix), 
     panel.first=grid(nx = NA, ny = NULL),
     xlab='Change in rate of clearance (%)', ylab='',xaxt='n',
     xlim = range(trt_matrix, na.rm = T), yaxt='n', pch=16, cex=2)
abline(v=log(c(0.8, 1, 1.25, 1.5, 2)), 
       lty='dotted',  col='lightgray')
axis(1, at = log(c(0.85, 1, 1.25, 1.5, 2)), 
     labels = c(-15, 0, 25, 50, 100))
axis(2, at = 1:ncol(trt_matrix), labels = names_xx, tick = F)
for(j in 1:ncol(trt_matrix)){
  lines(trt_matrix[ c(1,5), j], c(j,j))
  lines(trt_matrix[ c(2,4), j], c(j,j), lwd=3)
}
abline(v=0,lwd=2,lty=2)
abline(v=log(1.125),lwd=2,lty=2)

