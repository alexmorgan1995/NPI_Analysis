library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/Vaccine Model/Alltogether")


# Trajectory Plot ---------------------------------------------------------
trajplots <- do.call(rbind,
                     lapply(list.files(path = "C:/Users/amorg/Documents/PhD/nCoV Work/Models/Vaccine Model/Alltogether", pattern = ".csv"), read.csv))

trajplots_i <- subset(trajplots, variable == "I_i")
trajplots_i$group <- as.factor(trajplots_i$group)

p1 <- ggplot(trajplots_i, aes(x = time, y = I, col = group)) + geom_line(position = "identity", size = 1.02) +
  theme_bw() + scale_y_continuous(limits = c(0, 0.1), expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + 
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_blank(), legend.position='none',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  labs(x ="\nTime (Days)", y = "Prevalence of Naturally Infected (i)\n") 

# Attack Rate Histogram ---------------------------------------------------

attack_rates <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/Vaccine Model/Alltogether/attackrate/attackrates.csv", header = TRUE)
colnames(attack_rates) <- c("attack_rate", "group")

p2 <- ggplot(attack_rates, aes(attack_rate)) + geom_histogram(bins = 20, fill = "grey", col = "black") +
  theme_bw() + scale_y_continuous(limits = c(0, 40), expand = c(0,0)) + 
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=12), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=12), legend.title = element_text(size=18), legend.position='bottom',
        axis.text.x=element_text(size=12),axis.text.y=element_blank(), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  labs(x ="Attack Rate", y = "Frequency") 

# Final Plotting ----------------------------------------------------------


p_iplot <- p1 + annotation_custom(ggplotGrob(p2), xmin = 0, xmax = 200, 
                                   ymin = 0.05, ymax = 0.1)

ggsave(p_iplot, filename = paste0("combplot.png"), 
       dpi = 300, type = "cairo", width = 12, height = 10, units = "in")