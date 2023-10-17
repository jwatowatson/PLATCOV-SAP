library(tidyverse)
library(mgcv)

f_name = 'Analysis_Data/interim_all_analysis.csv'
platcov_dat = read.csv(f_name)
platcov_dat$Lab = factor(platcov_dat$Lab, levels = c( "Thailand",
                                                      "Brazil",
                                                      "Laos",
                                                      "Pakistan"))
my_cols = RColorBrewer::brewer.pal(n = 4, name = 'Dark2')
my_cols[1] = adjustcolor(my_cols[1],alpha.f = .3)
names(my_cols)=levels(platcov_dat$Lab)
platcov_dat$Lab_col = my_cols[as.numeric(platcov_dat$Lab)]

platcov_dat = platcov_dat%>% filter(Timepoint_ID==0) %>%
  group_by(ID) %>%
  mutate(sd_RNaseP = sd(CT_RNaseP),
         mean_RNaseP= mean(CT_RNaseP),
         sd_viral = sd(CT_NS),
         mean_viral = mean(CT_NS)) %>%
  ungroup() %>%
  distinct(ID, .keep_all = T)
  
par(las=1, family='serif')
plot(platcov_dat$mean_RNaseP, platcov_dat$sd_RNaseP, col=platcov_dat$Lab_col,
     ylim = c(0,3),panel.first = grid(),
     xlab='Mean quartet RNaseP (baseline swabs)',
     ylab='SD of quartet RNaseP',
     pch = as.numeric(platcov_dat$Lab))
legend('topright', col=my_cols, legend = levels(platcov_dat$Lab),pch=1:4)

summary(lm(sd_RNaseP ~ Lab + mean_RNaseP, data = platcov_dat))
cor.test(platcov_dat$mean_RNaseP, platcov_dat$sd_RNaseP)
boxplot(sd_RNaseP ~ Lab, data = platcov_dat); grid()


plot(platcov_dat$mean_viral, platcov_dat$sd_viral, col=platcov_dat$Lab_col,
     ylim = c(0,8),panel.first = grid(),
     xlab='Mean quartet VL CT (baseline swabs)',
     ylab='SD of quartet VL CT',
     pch = as.numeric(platcov_dat$Lab))
mod = gam(sd_viral ~ s(mean_viral, k=4) + Lab, data = platcov_dat)
for(ss in levels(platcov_dat$Lab)){
lines(15:40, predict(mod, newdata = data.frame(mean_viral=15:40, Lab=ss)), col=my_cols[ss],lwd=2)
}
legend('topright', col=my_cols, legend = levels(platcov_dat$Lab),pch=1:4)

