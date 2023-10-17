library(rstan)
library(RColorBrewer)
library(dplyr)

data <- read.csv("Analysis_Data/interim_all_analysis.csv")
str(data)

data$Rand_date <- as.Date(data$Rand_date)

library(ggplot2)

data$Trt <- as.factor(data$Trt)
data$Trt <- factor(data$Trt, levels = c("No study drug", levels(data$Trt)[levels(data$Trt) != "No study drug"]) )

data$Variant <- as.factor(data$Variant)
data$Variant <- factor(data$Variant, levels = as.character(data$Variant[!duplicated(data$Variant)]))

nb.cols <- length(levels(data$Variant))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

data$Site <- as.factor(data$Site)

G <- ggplot(data = data, aes(x=Rand_date, y = Trt, col = Variant)) +
  geom_point(size = 3) +
  facet_grid(.~ Site) +
  theme_bw() +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = seq.Date(as.Date("2022-01-01"), as.Date("2024-01-01"), "6 months"),
             col = "red", linetype = "dashed") +
  scale_x_date(date_labels =  "%b %y", breaks = seq.Date(as.Date("2022-01-01"), as.Date("2024-01-01"), "3 months")) +
  theme(#axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) + 
  scale_color_manual(values = mycolors) 
G

pdf("Plots/overall_enrollments_2309.pdf", width = 12, height = 5) 
G
dev.off()
####################################################################################
list_files <- c("data-TH1", "data-PK01", "data-LA08", "data-BR3")
data_recruitment <- do.call(rbind, lapply(paste0("C:/Users/Phrutsamon/Dropbox/PLATCOV/",list_files,".csv"), read.csv))
data_recruitment$Date <- gsub("  ", " ", data_recruitment$Date)
date <- do.call("rbind", str_split(data_recruitment$Date, " "))
data_recruitment$Date_new <- paste(date[,5], date[,2], date[,3], sep = "-")
data_recruitment$Date_new <- as.Date(data_recruitment$Date_new, format = "%Y-%b-%d")

data_recruitment <- data_recruitment[order(data_recruitment$Date_new),]
data_recruitment$count <- 1


data_recruitment2 <- aggregate(list("count" =data_recruitment$count), by = list("Rand_date" = data_recruitment$Date_new), FUN = sum)


G2 <- ggplot(data_recruitment2, aes(x = Rand_date, y = cumsum(count))) +
  geom_line(linewidth = 1) +
  theme_bw() +
  xlab("") +
  ylab("Cumulative recruitment") +
  scale_x_date(date_labels =  "%b %y", breaks = seq.Date(as.Date("2021-01-01"), as.Date("2024-01-01"), "3 months")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) +
  geom_vline(xintercept = seq.Date(as.Date("2021-01-01"), as.Date("2024-01-01"), "12 months"),
             col = "red", linetype = "dashed") +
  scale_y_continuous(breaks = seq(0,1200,200)) +
  geom_hline(yintercept = 0, linetype = "dashed")

png("Plots/cumulative_recruitment_2310.png", width = 8, height = 5, units = "in", res = 350) 
G2
dev.off()
