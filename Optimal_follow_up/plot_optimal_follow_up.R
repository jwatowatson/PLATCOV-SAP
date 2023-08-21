library(rstan)
library(ggplot2)
############################################################################
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
############################################################################
load('model_settings.RData')

# model_settings

Output <- NULL
tracePLOTs <- list()
for(i in 1:nrow(model_settings)){
  load(paste0('output/Rout/model_fits_',i, '.RData'))
  
  trt <- model_settings$intervention[i]
  ref <- model_settings$ref_arm[i]
  
  tracePLOTs[[i]] <- traceplot(out, par = c("trt_effect")) + 
    ggtitle(paste0(trt, " vs ", ref, ": ", model_settings$Dmax[i], " days")) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    ylim(-1.5,1.5)
  
  posterior <- unlist(extract(out, pars = c("trt_effect")))
  Dmax <- model_settings$Dmax[i]
  logeffects <- quantile(posterior, c(0.025, 0.5, 0.975))
  effects <- quantile(exp(posterior), c(0.025, 0.5, 0.975))
  
  
  output <- data.frame("Dmax" = Dmax, "log_eff_low" = logeffects[1], "log_eff_med" = logeffects[2], "log_eff_up" = logeffects[3],
             "eff_low" = effects[1], "eff_med" = effects[2], "eff_up" = effects[3], "itv_width" = logeffects[3] - logeffects[1],
             "trt" = trt, "ref" = ref)
  
  Output <- rbind(Output, output)
}
rownames(Output) <- NULL
############################################################################
Output$Dmax <- as.factor(Output$Dmax)
Output$Dmax <- factor(Output$Dmax, levels = as.character(1:14))

Output$trt[Output$trt == "Nirmatrelvir + Ritonavir"] <- "Paxlovid"
Output$ref[Output$ref == "Nirmatrelvir + Ritonavir"] <- "Paxlovid"
Output$trt[model_settings$data_ID == 3 & Output$trt == "Paxlovid"] <- "Recent Paxlovid"
Output$pairs <- paste0(Output$trt, " vs ", Output$ref)
Output$pairs <- as.factor(Output$pairs)
Output$pairs <- factor(Output$pairs, levels = c("Remdesivir vs No study drug", "Molnupiravir vs No study drug",
                                                "Molnupiravir vs Paxlovid", "Paxlovid vs No study drug", "Recent Paxlovid vs No study drug"))


G1 <- ggplot(Output, aes(x = Dmax, y = eff_med)) + 
  geom_point(size = 2) +
  facet_grid(.~pairs) +
  geom_errorbar(aes(ymin = eff_low, ymax = eff_up), width = 0.1, linewidth = 0.75) +
  theme_bw() +
  geom_hline(yintercept = 1, col = "red") +
  geom_hline(yintercept = 1.2, col = "red", linetype = "dashed") +
  xlab("Follow-up day") +
  ylab("Treatment effect size") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"))  +
  scale_x_discrete(drop=F) 


G2 <- ggplot(Output, aes(x = Dmax, y = itv_width)) + 
  geom_point(size = 3, col = "red") +
  facet_grid(.~pairs) +
 # geom_errorbar(aes(ymin = eff_low, ymax = eff_up), width = 0.1, linewidth = 0.75) +
  theme_bw() +
  #geom_hline(yintercept = 1, col = "red") +
  #geom_hline(yintercept = 1.2, col = "red", linetype = "dashed") +
  xlab("Follow-up day") +
  ylab("Credible interval width") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 14, face = "bold")) +
  scale_x_discrete(drop=F) 

ggarrange(G1, G2, ncol = 1, nrow = 2, align = "v")

?ggarrange
library(ggpubr)
ggarrange(plotlist=tracePLOTs, ncol = 4, nrow = 2, common.legend = T)

a















B <- ggplot(Output, aes(x = Dmax, y = itv_width)) + 
  geom_point(size = 3, col = "red") +
  theme_bw() +
  ylim(0, 1) +
  xlab("Follow-up day") +
  ylab("Credible interval width") +
  ggtitle("Remdesivir vs No study drug")

library(ggpubr)
ggarrange(A, B, ncol = 2, labels = "AUTO")







