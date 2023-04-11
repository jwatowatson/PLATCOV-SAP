dat=read.csv('~/Dropbox/MORU/Adaptive Trials/PLATCOV/Data/Result Validate half vol.csv')


par(mfrow=c(2,2),las=1)
plot(dat$ct.N.S.Full.vol, dat$ct.N.S.Full.vol-dat$ct.N.S.Half.Vol, 
     xlab='CT value (full volume)', ylab='Difference in CT values')
title('Viral load experiment 1')
abline(h=0)
plot(dat$ct.N.S.Full.vol, dat$ct.N.S.Half.Vol,xlab='CT value (full volume)',ylab='CT value (half volume)')
title('Viral load experiment 1')
lines(0:40, 0:40)

plot(dat$ct.N.S.Full.vol, dat$ct.N.S.Full.vol-dat$ct.N.S.Half.Vol.1, 
     xlab='CT value (full volume)', ylab='Difference in CT values')
title('Viral load experiment 2')
abline(h=0)
plot(dat$ct.N.S.Full.vol, dat$ct.N.S.Half.Vol.1,xlab='CT value (full volume)',ylab='CT value (half volume)')
lines(0:40, 0:40)
title('Viral load experiment 2')




par(mfrow=c(2,2))
plot(dat$ct.Rnase.Full.vol, dat$ct.Rnase.Full.vol-dat$ct.Rnase.Half.Vol, 
     xlab='CT value (full volume)', ylab='Difference in CT values')
title('RNasep experiment 1')
plot(dat$ct.Rnase.Full.vol, dat$ct.Rnase.Half.Vol,xlab='CT value (full volume)',ylab='CT value (half volume)')
title('RNasep experiment 1')
lines(0:40, 0:40)

plot(dat$ct.Rnase.Full.vol, dat$ct.Rnase.Full.vol-dat$ct.Rnase.Half.Vol.1, 
     xlab='CT value (full volume)', ylab='Difference in CT values')
title('RNasep experiment 2')
abline(h=0)
plot(dat$ct.Rnase.Full.vol, dat$ct.Rnase.Half.Vol.1,xlab='CT value (full volume)',ylab='CT value (half volume)')
title('RNasep experiment 2')
lines(0:40, 0:40)
