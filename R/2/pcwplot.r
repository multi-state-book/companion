rateCyA <- c(8.1,13.1,9.6)
ratePbo <- c(9.4,12.5,8.5)
pcwtime <- c(0,2,4,5)

# Collect data
plotdata <- data.frame(rates = c(rateCyA, ratePbo),
                       tment = c(rep("CyA", length(rateCyA)), rep("Placebo", length(ratePbo))),
                       times_s = rep(pcwtime[-4], 2),
                       times = rep(pcwtime[-1], 2))

# Create Figure 2.3
ggplot(aes(x = time, y = rates, linetype = tment),
                data = plotdata) +
  geom_segment(aes(x = times_s, y = rates, xend = times, yend = rates), size = 1) +
  scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA"))  +
  xlab("Time since randomization (years)") +
  ylab("Estimated hazard function (per 100 years)") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0,14), breaks = seq(0, 14, 2)) +
  theme_general
