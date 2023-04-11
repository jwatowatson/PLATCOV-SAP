library(readr)
library(tidyverse)

dat = read_csv('~/Dropbox/MORU/Adaptive Trials/PLATCOV/Data/Serology/PLATCOV_batch1.csv')
clin_data = haven::read_dta('~/Dropbox/MORU/Adaptive Trials/PLATCOV/Data/InterimEnrolment.dta')

dat = dat[, c('SUBJECT\nNO','Time Point','Spike IgG - Original (ng/mL)_400',
              'Spike IgG - Original (ng/mL)_3200',
              'Spike IgG - Original (ng/mL)_3200_v2')] %>% 
  mutate(IgG_ngml = case_when(`Spike IgG - Original (ng/mL)_400` < 60000 ~ `Spike IgG - Original (ng/mL)_400`,
                              TRUE ~ `Spike IgG - Original (ng/mL)_3200`),
         IgG_ngml = ifelse(IgG_ngml<1,50,IgG_ngml))

dat = merge(clin_data, dat[, c('SUBJECT\nNO','Time Point','IgG_ngml')],
            all.y=T, by.y = 'SUBJECT\nNO', by.x = 'Label')

dat = dat %>% mutate(
  log10_IgG = log10(IgG_ngml),
  Day = gsub(x = `Time Point`, pattern = 'D', replacement=''),
  Day = as.numeric(gsub(x = Day, pattern = 'H0', replacement='')),
  ID = Label
)
boxplot(log10_IgG ~ Day, data = dat)

table(is.na(dat$log10_IgG))


write_csv(x = dat[, c('ID','Day','log10_IgG')], file = '../Analysis_Data/Serology_IgG.csv')
