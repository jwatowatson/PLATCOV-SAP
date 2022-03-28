platcov_dat = read.csv('interim_dat.csv')
swabber = read.csv('../Data/Participant swabs by day and clinician.csv')
swabber$Swabber[swabber$Swabber=='-']=NA
platcov_dat$t_id = apply(platcov_dat[, c('ID','Timepoint_ID')],
                         1,paste0,collapse='_')
swabber$t_id = apply(swabber[, 1:2],1,paste0,collapse='_')
xx=aggregate(formula=CT_NS ~ t_id,
             data = platcov_dat[platcov_dat$Timepoint_ID>0,],
             FUN = function(x) abs(x[1]-x[2]))

xx=merge(xx, swabber, by = 't_id')

xx$Swabber = as.character(xx$Swabber)
table(xx$Swabber)

mod=lm(CT_NS ~ Swabber, data = xx)
summary(mod)
