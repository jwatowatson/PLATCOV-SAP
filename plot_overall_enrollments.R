library(rstan)
getwd()

mainDir <- "D:/PLATCOV-SAP"
setwd(mainDir)

data <- read.csv("Analysis_Data/interim_all_analysis.csv")
str(data)

data$Rand_date <- as.Date(data$Rand_date)

library(ggplot2)

data$Trt <- as.factor(data$Trt)
data$Trt <- factor(data$Trt, levels = c("No study drug", levels(data$Trt)[levels(data$Trt) != "No study drug"]) )

data$Variant <- as.factor(data$Variant)
data$Variant <- factor(data$Variant, levels = c("Delta", levels(data$Variant)[levels(data$Variant) != "Delta"]))

data$Site <- as.factor(data$Site)

G <- ggplot(data = data, aes(x=Rand_date, y = Trt, col = Variant)) +
  geom_point(size = 3) +
  facet_grid(.~ Site) +
  theme_bw() +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = seq.Date(as.Date("2022-01-01"), as.Date("2024-01-01"), "6 months"),
             col = "red", linetype = "dashed") +
  scale_x_date(date_labels =  "%b %y") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),
        strip.text = element_text(size = 12, face = "bold")) +
  scale_color_brewer(palette="Set1")

pdf("Plots/overall_enrollments.pdf", width = 12, height = 5) 
G
dev.off()


