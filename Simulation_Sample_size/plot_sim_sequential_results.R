library(waffle)
library(sjPlot)

f_sims = list.files('sims_out_sequential3',full.names = T, pattern = '*.csv')
names(f_sims) = as.numeric(gsub(gsub(f_sims, pattern = 'sims_out_sequential3/sim_sequential_',replacement = ''),
                                pattern = '.csv',replacement = ''))
f_sims = f_sims[as.character(sort(as.numeric(names(f_sims))))]

sim_res = lapply(f_sims, read.csv)
sim_seq_all = dplyr::bind_rows(sim_res)

Nmax = unique(sim_seq_all$Nmax)[1]

gg_out_list = gg_out_list2 = list(); i=1
ind_success = sim_seq_all$success==1
ind_futilty = sim_seq_all$futility==1

outcome_names = c('No result','Futility','Success')
my_cols = RColorBrewer::brewer.pal(n = 3,name = 'Set2')
names(my_cols) = outcome_names

for(tt in unique(sim_seq_all$intervention_effect)){
  for(threshold in unique(sim_seq_all$Futility_delta)){
    ind=sim_seq_all$intervention_effect==tt & sim_seq_all$Futility_delta==threshold
    result = sim_seq_all$futility[ind] + 2*sim_seq_all$success[ind]
    result = plyr::mapvalues(x = result,from = 0:2,to = outcome_names)

    result = factor(result, levels = outcome_names)
    mydat = data.frame(N=sim_seq_all$N_futility_success[ind],
                       result=result)
    mydat$N[is.na(mydat$N)]=Nmax
    result = table(result)
    writeLines(sprintf('The median sample size is %s for intervention effect of %s%% and NI threshold of %s%%',
                       median(mydat$N),
                       round(100*(tt-1)),
                       round(100*(exp(threshold)-1),1)))
    print(result/sum(result))
    plot_title = paste0('Effect=',round(100*(tt-1)),'%; threshold=',round(100*(exp(threshold)-1),1),'%')
    gg_out_list[[i]] = waffle(parts = result,rows = 20,
                              colors = my_cols,
                              title = plot_title)


    gg_out_list2[[i]] = ggplot(data = mydat, aes(x=N,fill = result)) +
      ggtitle(plot_title) +
      geom_histogram(breaks=seq(10,Nmax,by=10))+
      scale_fill_manual(values=my_cols[names(result[result>0])])

    i = i+1
  }
}
ggsave(plot_grid(gg_out_list, tags = rep('',9),margin = c(2,0,3,3)/10),filename = 'proportion_outcomes_SF.png',device = 'png')
ggsave(plot_grid(gg_out_list2, tags = rep('',9),margin = c(2,0,3,3)/10),filename = 'Sample_size_byoutcome_SF.png',device = 'png')


outcome_names = c('No result','Inferiority','Non-inferiority')
my_cols = RColorBrewer::brewer.pal(n = 3,name = 'Set2')
names(my_cols) = outcome_names

gg_out_list = gg_out_list2 = list(); i=1
for(tt in unique(sim_seq_all$intervention_effect)){
  for(threshold in unique(sim_seq_all$NI_delta)){
    ind=sim_seq_all$intervention_effect==tt &
      sim_seq_all$NI_delta==threshold &
      sim_seq_all$success==1

    result = sim_seq_all$inferiority[ind] + 2*sim_seq_all$non_inferiority[ind]
    result = plyr::mapvalues(x = result,from = 0:2,to = outcome_names)

    result = factor(result, levels = outcome_names)
    mydat = data.frame(N=ifelse(is.na(sim_seq_all$N_inferiority[ind]),
                                sim_seq_all$N_non_inferiority[ind],
                                sim_seq_all$N_inferiority[ind]),
                       result=result)
    mydat$N[is.na(mydat$N)]=Nmax
    result = table(result)
    
    writeLines(sprintf('The median sample size is %s for intervention effect of %s%% and NI threshold of %s%%',
                       median(mydat$N),
                       round(100*(tt-1)),
                       round(100*(exp(-threshold)-1),1)))
    print(result/sum(result))
    
    plot_title = paste0('Effect=',round(100*(tt-1)),'%; NI=',round(100*(exp(-threshold)-1),1),'%')
    gg_out_list[[i]] = waffle(parts = result,rows = 20,
                              colors = my_cols,
                              title = plot_title)


    gg_out_list2[[i]] = ggplot(data = mydat, aes(x=N,fill = result)) + ggtitle(plot_title) +
      geom_histogram(breaks=seq(10,Nmax,by=10))+
      scale_fill_manual(values=my_cols[names(result[result>0])])

    i = i+1
  }
}
ggsave(plot_grid(gg_out_list, tags = rep('',9),margin = c(2,0,3,3)/10),filename = 'proportion_outcomes_NI.png',device = 'png')
ggsave(plot_grid(gg_out_list2, tags = rep('',9),margin = c(2,0,3,3)/10),filename = 'Sample_size_byoutcome_NI.png',device = 'png')
