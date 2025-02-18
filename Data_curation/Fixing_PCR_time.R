timepoint_id_not_matched %>%
  filter(Timepoint_ID <= 5)

ID_i <- "PLT-TH58-009"

timepoint_id_not_matched %>% filter(ID == ID_i) %>% arrange(Timepoint_ID)

Res %>% filter(ID == ID_i) %>%
  select(ID, BARCODE, Timepoint_ID, Time, Rand_date_time, `Lot no.`) %>%
  arrange(Timepoint_ID)

#clin_data %>% filter(Label == ID_i)

#log_data %>% filter(sl_barc %in% c("24AK019", "24AK016")) 

#Manual corrections
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-BR3-006' & Timepoint_ID == 0 & Time > 1, Time - 1, Time))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-LA8-001' & Timepoint_ID %in% c(5,6) & Time > 360, Time - 365, Time))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-LA8-027' & Time > 20, Time - 30, Time))
Res <- Res %>% mutate(Timepoint_ID = if_else(ID == 'PLT-TH1-1116' & BARCODE %in% c("23MD734", "23MD737"), 3, Timepoint_ID))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH1-1661' & Time > 20, Time - 30, Time))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH1-1743' & Time > 20, Time - 152, Time))
Res <- Res %>% mutate(Timepoint_ID = if_else(ID == 'PLT-TH1-1851' & BARCODE %in% c("24DR403", "24DR409"), 2, Timepoint_ID))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH1-1887' & Timepoint_ID == 5, Time + 1, Time))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH1-1960' & Timepoint_ID == 5, Time + 1, Time))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH1-1980' & Timepoint_ID == 5, Time + 1, Time))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH1-276', Time - 1, Time))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH1-682'& Timepoint_ID == 4, Time + 1, Time))
Res <- Res %>% mutate(Timepoint_ID = if_else(ID == 'PLT-TH1-781' & BARCODE %in% c("20HK232"), 4, Timepoint_ID))
Res <- Res %>% mutate(Timepoint_ID = if_else(ID == 'PLT-TH1-798' & BARCODE %in% c("20HK832"), 3, Timepoint_ID))
Res <- Res %>% mutate(Time = if_else(ID == 'PLT-TH58-009' & Timepoint_ID %in% 0:7, Time - 1, Time))






aa <- Res %>% filter(`SUBJECT ID` == ID_i) %>%
  select(`SUBJECT ID`, BARCODE, `TIME-POINT`, Location, `Lot no.`)
