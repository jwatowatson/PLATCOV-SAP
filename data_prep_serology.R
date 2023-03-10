library(readr)
library(tidyverse)

dat = read_csv('../Data/Serology/PLATCOV_batch1.csv')
clin_data = haven::read_dta('../Data/InterimEnrolment.dta')

dat = merge(clin_data, dat, all.y=T, by.y = 'SUBJECT\nNO', by.x = 'Label')

dat = dat %>% mutate(
  IgG = ifelse(`Spike IgG - Original (ng/mL)` <= 0, 1, `Spike IgG - Original (ng/mL)`),
  log10_IgG = log10(IgG),
  Day = gsub(x = `Time Point`, pattern = 'D', replacement=''),
  Day = as.numeric(gsub(x = Day, pattern = 'H0', replacement='')),
  ID = Label
)
boxplot(log10_IgG ~ Day, data = dat)

table(is.na(dat$log10_IgG))


write_csv(x = dat[, c('ID','Day','log10_IgG')], file = 'Analysis_Data/Serology_IgG.csv')
