library(rstan)
library(RColorBrewer)


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


