library(readr)
library(tidyverse)

dat = read_csv('~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/Data/Serology/PLATCOV Results Master List 181023.csv')
clin_data = haven::read_dta('~/Dropbox/MORU/Adaptive Trials/PLATCOV_Analysis/Data/InterimEnrolment.dta')

colnames(dat)
dat = dat[, c('SUBJECT\nNO','Time Point',
              "Spike IgG - Original (ng/mL) 1/400 Dilution",
              "Spike IgG - Original (ng/mL) 1/50 Dilution Repeat",
              "Spike IgG - Original (ng/mL) 1/3200 Dilution",
              "Spike IgG - Original (ng/mL) 1/3200 Dilution Repeat",
              "Spike IgG - Original (ng/mL) 1/12800 Dilution Repeat"
              )] %>% ungroup() %>%
  mutate(
    IgG_ngml_50 = `Spike IgG - Original (ng/mL) 1/50 Dilution Repeat`,
    IgG_ngml_400 = `Spike IgG - Original (ng/mL) 1/400 Dilution`,
    IgG_ngml_3200 = `Spike IgG - Original (ng/mL) 1/3200 Dilution`,
    IgG_ngml_12800 = `Spike IgG - Original (ng/mL) 1/12800 Dilution Repeat`,
    IgG_ngml = case_when(IgG_ngml_400 <= 60000 & IgG_ngml_400 > 20000 ~ IgG_ngml_400,
                         IgG_ngml_400 <= 20000 ~ IgG_ngml_50,
                         (IgG_ngml_400 > 60000 & IgG_ngml_3200 < 450000) | (is.na(IgG_ngml_12800) & IgG_ngml_3200 < 600000) ~ IgG_ngml_3200,
                         T ~ IgG_ngml_12800),
    IgG_ngml = ifelse(IgG_ngml<1, 1, IgG_ngml))

dat = merge(clin_data, dat[, c('SUBJECT\nNO','Time Point','IgG_ngml')],
            all.y=T, by.y = 'SUBJECT\nNO', by.x = 'Label')

dat = dat %>% mutate(
  log10_IgG = log10(IgG_ngml),
  Day = gsub(x = `Time Point`, pattern = 'D', replacement=''),
  Day = as.numeric(gsub(x = Day, pattern = 'H0', replacement='')),
  ID = Label
) %>% select(ID, Day, log10_IgG)
boxplot(log10_IgG ~ Day, data = dat)

table(is.na(dat$log10_IgG))


write_csv(x = dat[, c('ID','Day','log10_IgG')], file = '../Analysis_Data/Serology_IgG.csv')
