SC = read.csv('Analysis_Data/interim_control_dat.csv')
table(SC$Lab)

control_dat = dplyr::arrange(SC, CT_NS, Plate)
control_dat$CT = control_dat$CT_NS
control_dat$CT[control_dat$CT_NS==40]=NA
control_dat$batch = as.factor(control_dat$Plate)
control_dat$Lab = as.factor(control_dat$Lab)
hist(table(control_dat$Plate))
library(lme4)
# library(rstanarm)
# options(mc.cores = parallel::detectCores())
conv_mod0 = lmer(CT ~ 1 + log10_true_density + (1+log10_true_density|batch),
                      data = control_dat[control_dat$Lab=='Thailand',])

summary(conv_mod0)

conv_mod2 = lmer(CT ~ 1 + log10_true_density + (1+log10_true_density|batch),
                 data = control_dat[control_dat$Lab!='Thailand',])

summary(conv_mod2)

# control = lmerControl(optimizer ="Nelder_Mead"))
# xx=posterior_predict(conv_mod0, re.form=NA, newdata=data.frame(Lab=c('Thailand','Brazil'),
#                                                   log10_true_density=4))
# colMeans(xx)
conv_mod = lmer(log10_true_density ~ 1 + CT*Lab + (1+CT|batch), 
                data = control_dat,
                control = lmerControl(optimizer ="Nelder_Mead"))
summary(conv_mod)

preds = predict(conv_mod)
par(las=1,cex.lab=1.3, cex.axis=1.3,family='serif')
my_cols = adjustcolor(c('lightgrey','pink'),.7)
plot(control_dat$CT, jitter(control_dat$log10_true_density), xlim=c(20,40),
     col = my_cols[as.numeric(control_dat$Lab=='Brazil')+1],
     xlab = 'CT value', ylab = 'Control density', panel.first=grid(), yaxt='n')
axis(2, at =2:7, labels = c(expression(10^2),
                            expression(10^3),
                            expression(10^4),
                            expression(10^5),
                            expression(10^6),
                            expression(10^7)))
for(bb in levels(control_dat$batch)){
  ind = control_dat$batch==bb
  lines(control_dat$CT[ind], preds[ind], 
        col= my_cols[as.numeric(control_dat$Lab[ind]=='Brazil')+1])
}
lines(20:40, predict(conv_mod, re.form=NA,data.frame(CT=20:40,Lab='Thailand',batch= -1)),
      lwd=3,col='brown')
lines(20:40, predict(conv_mod, re.form=NA,data.frame(CT=20:40,Lab='Brazil',batch= -1)),
      lwd=3,col='purple')

legend('topright', col=my_cols, legend = c('Thailand','Brazil'),
       cex=2,lwd=3,inset=0.03)


platcov_dat = read.csv('interim_dat.csv')
my_probs = c(0.05,.95)
par(mfrow=c(2,1),las=1)
ind=platcov_dat$Lab=='Brazil'
hist(platcov_dat$CT_RNaseP[ind], 
     breaks = 19:40,main='Brazil',xlab='RNaseP CT value')
qqs= quantile(platcov_dat$CT_RNaseP[ind],probs=my_probs,na.rm = T)
diff(qqs)
print(sd(platcov_dat$CT_RNaseP[ind &
                                 platcov_dat$CT_RNaseP<qqs[2]&
                                 platcov_dat$CT_RNaseP>qqs[1]],na.rm=T))
abline(v=median(platcov_dat$CT_RNaseP[ind],na.rm = T),lwd=3,col='red')

ind=platcov_dat$Lab=='Thailand'
hist(platcov_dat$CT_RNaseP[ind], 
     breaks = 19:40,main='Thailand',xlab='RNaseP CT value')
abline(v=median(platcov_dat$CT_RNaseP[ind],na.rm = T),lwd=3,col='red')
qqs= quantile(platcov_dat$CT_RNaseP[ind],probs=my_probs,na.rm = T)
diff(qqs)
print(sd(platcov_dat$CT_RNaseP[ind &
                                 platcov_dat$CT_RNaseP<qqs[2]&
                                 platcov_dat$CT_RNaseP>qqs[1]],na.rm=T))


ind = is.na(platcov_dat$Time) | platcov_dat$Time<0
table(ind)
platcov_dat$Time[ind]=platcov_dat$Timepoint_ID[ind]
# pdf('Brazil_profiles.pdf')
par(mfrow=c(4,4),cex.lab=1.1, cex.axis=1.1, mar=c(5,4,1,1),las=1)
for(id in unique(platcov_dat$ID[platcov_dat$Site=='br003'])){
  ind = platcov_dat$ID==id
  plot(platcov_dat$Time[ind], platcov_dat$log10_viral_load[ind], ylim = c(0.5, 9),
       panel.first=grid(), xlab="Time since randomisation",pch=16,
       ylab = 'RNA copies per mL',yaxt='n',cex=1.5,xlim=c(0,7))
  axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
}


for(id in sample(unique(platcov_dat$ID[platcov_dat$Site=='th001']),replace = F,size = 16)){
  ind = platcov_dat$ID==id
  plot(platcov_dat$Time[ind], platcov_dat$log10_viral_load[ind], ylim = c(0.5, 9),
       panel.first=grid(), xlab="Time since randomisation",pch=16,
       ylab = 'RNA copies per mL',yaxt='n',cex=1.5,xlim=c(0,7))
  axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
}

xx = aggregate(log10_viral_load ~ Timepoint_ID + ID + Site,
               data = platcov_dat,
               FUN = function(x) c(x[1]-x[2], mean(x)))

plot(xx$log10_viral_load[,2], xx$log10_viral_load[,1],
     col=as.numeric(as.factor(xx$Site)))

par(mfrow=c(2,2),cex.lab=1.3, cex.axis=1.3, las=1)
for(ss in unique(platcov_dat$Site)){
  ind = xx$Site==ss
  hist(xx$log10_viral_load[ind],main=ss,xlab='Difference between duplicates')
  print(ss)
  print(median(xx$log10_viral_load[ind]))
}

xx = aggregate(log10_viral_load ~ Timepoint_ID + ID + Site,
               data = platcov_dat,
               FUN = mean)

xx_pop = aggregate(log10_viral_load ~ Timepoint_ID + Site,
                   data = xx,
                   FUN = quantile, probs=c(0.05, 0.5, 0.95))
par(mfrow=c(2,2),
    cex.lab=1.3, cex.axis=1.3, las=1)
for(ss in unique(platcov_dat$Site)){
  ind = xx_pop$Site==ss
  plot(xx_pop$Timepoint_ID[ind], xx_pop$log10_viral_load[ind,2],
       ylim = c(0.5, 7.5),type='l',lwd=3,panel.first=grid(),
       xlim =c(0,7),main=ss,
       xlab='Time', ylab='log10 copies per ml')
  lines(xx_pop$Timepoint_ID[ind], xx_pop$log10_viral_load[ind,1])
  lines(xx_pop$Timepoint_ID[ind], xx_pop$log10_viral_load[ind,3])
}
legend('topright',lwd=c(3,1),lty=1,legend = c('Median','10/90th percentile'))



rm(list=ls())

