library(rstan)
library(brms)
library(readr)
options(mc.cores = parallel::detectCores())
control_dat=read_csv(file = "../Analysis_Data/interim_control_dat.csv")
control_dat$Plate = control_dat$`Lot no.`
my_cols = RColorBrewer::brewer.pal(n = 4, name = 'Dark2')
control_dat$Lab = factor(control_dat$Lab, levels= c( "Thailand",
                                                     "Brazil",
                                                     "Laos",
                                                     "Pakistan"))
names(my_cols)=levels(control_dat$Lab)
control_dat$Lab_col = my_cols[as.numeric(control_dat$Lab)]

sort(table(control_dat$`Lot no.`))
# random slope and intercept
conv_mod = brm(log10_true_density ~ 1 + CT*Lab + (1+CT| Plate),
               data = control_dat)
summary(conv_mod)


preds=posterior_predict(object = conv_mod) 
newdat = expand.grid(CT=20:40, Lab=unique(control_dat$Lab),Plate=NA)

par(las=1)
plot(range(newdat$CT), range(colMeans(preds)), col=my_cols[newdat$Lab],type='n',
     xlab='CT value' , ylab='Control log10 density', panel.first=grid())
for(ss in unique(control_dat$Plate)){
  ind = control_dat$Plate==ss & !is.na(control_dat$CT)
  lab_ss = control_dat$Lab[ind][1]
  lines(control_dat$CT[ind], colMeans(preds)[ind], col = adjustcolor(my_cols[lab_ss],.3),lwd=.5)
}

preds=posterior_predict(object = conv_mod,newdata = newdat) 

for(ss in unique(control_dat$Lab)){
  ind = newdat$Lab==ss
  lines(newdat$CT[ind], colMeans(preds)[ind], col = my_cols[ss],lwd=3)
}
legend('topright', col=my_cols, legend = levels(control_dat$Lab),lwd=2)



