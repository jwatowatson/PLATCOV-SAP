sim_vl2 <- aggregate(log10_viral_load~ID+Time+Trt, 
                       data = sim_vl, FUN = median)
sim_vl3 <- aggregate(log10_viral_load~Time+Trt, 
                     data = sim_vl2, FUN = median)

ggplot() +
  geom_point(data = sim_vl2, aes(x = Time, y = log10_viral_load, col = Trt), alpha = 0.4, size = 1.75, shape = 16) +
  geom_line(data = sim_vl2, aes(x =  Time, y = log10_viral_load, group = ID, col = Trt), 
            alpha = 0.5, linewidth = 0.5, linetype = 1) +
  scale_color_manual(values = c("#FA9494", "#87CEEB"), name = "") +
  ggnewscale::new_scale_color() +
  geom_point(data = sim_vl3, aes(x = Time, y = log10_viral_load, col = Trt), size = 3, shape = 17) +
  geom_line(data = sim_vl3, aes(x =  Time, y = log10_viral_load, group = Trt, col = Trt), linewidth = 1, linetype = 1) +
  scale_color_manual(values = c("#B31312", "#332FD0"), name = "") +
  theme_bw() +
  xlab("") +
  ylab("Log viral loads") + 
  theme(axis.title  = element_text(face = "bold"))


summary(exp(out_trt_effect))
