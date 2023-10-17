prefix_dropbox <- "C:/Users/Phrutsamon/Dropbox/PLATCOV_Analysis"

fnames = list.files(paste0(prefix_dropbox, "/Data/CSV files"),full.names = T,recursive = T)

# column type specification for each csv file
my_specs = cols(
  `SUBJECT ID` = col_character(),
  BARCODE = col_character(),
  INITIALS = col_character(),
  Location = col_character(),
  `TIME-POINT` = col_character(),
  `Time Collected` = col_time(format = ""),
  `COLLECTION DATE` = col_character(),
  VOLUME = col_double(),
  `UNIT (mL)` = col_character(),
  `STORAGE DATE` = col_character(),
  `STORAGE TIME` = col_character(),
  BOX = col_double(),
  POSITION = col_character(),
  NOTE = col_character(),
  `Date Received` = col_character(),
  `Sample ID` = col_character(),
  `N/S Gene` = col_character(),
  RNaseP = col_character(),
  `Target conc. c/mL` = col_number(),
  `Lot no.` = col_character()
)

for(i in 1:length(fnames)){
  # writeLines(sprintf('Loading data from file:\n %s \n*********************************************************', fnames[i]))
  
  # read in file
  temp = readr::read_csv(fnames[i],col_types = my_specs)
  
  # check for duplicates
  if(any(duplicated(temp$BARCODE[!is.na(temp$BARCODE)]))){
    writeLines(sprintf('in file %s there are duplicate barcodes',fnames[i]))
  }
  # make sure that PCR plates have unique codes across sites
  if(length(grep(pattern = 'Brazil', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Brazil',x,sep='_'))
  }
  if(length(grep(pattern = 'Thailand', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Thailand',x,sep='_'))
  }
  if(length(grep(pattern = 'Laos', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Laos',x,sep='_'))
  }
  if(length(grep(pattern = 'Pakistan', x = fnames[i], ignore.case = F))>0){
    temp$`Lot no.` = apply(temp[, 'Lot no.', drop=F], 1, function(x) paste('Pakistan',x,sep='_'))
  }
  if(i==1){
    Res=temp
  } else {
    Res = rbind(Res,temp)
  }
}

# take out the summary PCR columns
ind_rm = Res$`Sample ID` %in% c('PC','R-squared','Efficiency (%)','Slope (M)')|
  (is.na(Res$`SUBJECT ID`) & is.na(Res$`Sample ID`))
sum(ind_rm)
Res = Res[!ind_rm, ]

# take out duplicate rows
ind_rm = !is.na(Res$BARCODE) & duplicated(Res$BARCODE)
sum(ind_rm)
Res = Res[!ind_rm, ]

sort(unique(Res$`Lot no.`))
Res$`Lot no.`[Res$`Lot no.`=='Thailand_D10 lot1']="Thailand_D10 Lot 1"
Res$`Lot no.`[Res$`Lot no.`=='Thailand_D10 lot2']="Thailand_D10 Lot 2"

## 
Res$`SUBJECT ID` = gsub(pattern = '_', replacement = '-', x = Res$`SUBJECT ID`,fixed = T)
Res$`SUBJECT ID` = gsub(pattern = 'PK01', replacement = 'PK1', x = Res$`SUBJECT ID`,fixed = T)

# make the plate/lab variables
Res$Lab = NA
Res$Lab[grep(pattern = 'Thailand',x = Res$`Lot no.`, ignore.case = T)]='Thailand'
Res$Lab[grep(pattern = 'Brazil',x = Res$`Lot no.`, ignore.case = T)]='Brazil'
Res$Lab[grep(pattern = 'Lao',x = Res$`Lot no.`, ignore.case = T)]='Laos'
Res$Lab[grep(pattern = 'Paki',x = Res$`Lot no.`, ignore.case = T)]='Pakistan'

Res$`Lot no.` = tolower(Res$`Lot no.`)
sort(unique(Res$`Lot no.`))
Res$Plate = as.numeric(as.factor(Res$`Lot no.`))
Res$Lab = as.factor(Res$Lab)

## Missing samples
ind_missing = grep('not',Res$`Sample ID`,ignore.case = T)
writeLines('Missing data for the following samples:')
Res$`TIME-POINT`[ind_missing]
Res = Res[-ind_missing, ]

## Extract standard curve data by plate
ind = grep('std', Res$`Sample ID`)
Res_new = Res[-ind, ]

ind_undetermined <- which(Res_new$RNaseP == "Undetermined")
ind_na <- which(is.na(Res_new$RNaseP))

writeLines('RNaseP of the following samples are undetermined:')
Res_new$`Sample ID`[c(ind_undetermined)]

Res_new<- Res_new[-union(ind_undetermined,ind_na),]

Res_new$RNaseP <- as.numeric(Res_new$RNaseP)

library(dplyr)
summarize_RNaseP <- Res_new %>% 
  group_by(Lab) %>% 
  summarise(Variance=var(RNaseP), Mean = mean(RNaseP),
            Median = median(RNaseP),
            Count = n())

Res_new$RNaseP

summarize_RNaseP$lab_var <- paste0("Var = ", round(summarize_RNaseP$Variance, 2))
summarize_RNaseP$lab_mean <- paste0("Mean = ", round(summarize_RNaseP$Mean, 2))
summarize_RNaseP$lab_n <- paste0("n = ", summarize_RNaseP$Count)


G <- ggplot(Res_new, aes(x = RNaseP)) + geom_histogram(aes(y = ..density..), bins = 50) +
  facet_grid(Lab~.) +
  theme_bw()  +
  xlab("CT values for RNaseP") +
  ylab("Density") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.3, "lines")) +
  ggtitle("Distributions of CT values for RNaseP") +
  geom_text(data = summarize_RNaseP, aes(x = 30, y = 0.25, label = lab_var), hjust = 0) +
  #geom_text(data = summarize_RNaseP, aes(x = 30, y = 0.3, label = lab_mean), hjust = 0) +
  geom_text(data = summarize_RNaseP, aes(x = 30, y = 0.3, label = lab_n), hjust = 0) +
   geom_vline(data = summarize_RNaseP, aes(xintercept = Median), col = "red", linetype = "dashed", linewidth = 0.75) +
   geom_hline(yintercept = 0, linetype = "dashed")

png(paste0("../Plots/RNaseP_histogram.png"), width = 6, height = 8, units = "in", res = 350)
G
dev.off()


?var

