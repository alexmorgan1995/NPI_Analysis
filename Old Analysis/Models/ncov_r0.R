test <- data.frame(matrix(ncol = 0, nrow = 7))
test$paper <- c("imai2020","riou2020", "read2020", "liu2020", "zhao2020", "li2020","wu2020" )
test$r0 <- c(2.6, 2.2, 3.1, 2.9, 3.3, 2.2, 2.7)
test$lowci <- c(1.5, 1.4, 2.4, 2.3, 2.7, 1.4, 2.5) 
test$highci <- c(3.5, 3.8, 4.1, 3.6, 4,  3.9,  2.9) 

ggplot(data = test, aes(x = paper, y = r0)) + geom_point(size = 2) +
  geom_pointrange(aes(ymin = lowci, ymax = highci)) +
  labs(x ="Paper", y = "Basic Reproduction Number") + 
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) + 
  geom_hline(yintercept = 1, col = "red", lty=2, lwd = 1)