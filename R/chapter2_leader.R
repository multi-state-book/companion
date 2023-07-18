# Set working directory - location of data
setwd("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB")

leader <- read.csv("data/leader_mi_3p.csv")
# One data set per endpoint type
leader_mi <- subset(leader, type == "recurrent_mi")
leader_3p <- subset(leader, type == "recurrent_comb_str_mi_cvdth")

#------------------------------------------------------------------#
# Table 2.14 
#------------------------------------------------------------------#

# Cox model first events
library(survival)
coxph(Surv(start, stop, status == 1) ~ treat,
      data = subset(leader_mi, eventno == 1), method = "breslow",)
    
# AG Cox type
coxph(Surv(start, stop, status == 1) ~ treat, data = leader_mi,
                   method = "breslow")

# PWP 2. event
coxph(Surv(start, stop, status == 1) ~ treat,
                     method = "breslow", subset = (eventno == 2),
                     data = leader_mi)

# PWP 3. event
coxph(Surv(start, stop, status == 1) ~ treat,
                   method = "breslow", subset = (eventno == 3),
                   data = leader_mi)

# PWP 4. event
coxph(Surv(start, stop, status == 1) ~ treat,
                   method = "breslow", subset = (eventno == 4),
                   data = leader_mi)

# PWP 5. event
coxph(Surv(start, stop, status == 1) ~ treat,
                   method = "breslow", subset = (eventno == 5),
                   data = leader_mi)

# PWP all
coxph(Surv(start, stop, status == 1) ~ treat + strata(eventno),
                   method = "breslow",
                   data = leader_mi)

# AG piece-wise constant
# Calculating cuts -> 5
alltimes <- seq(0,max(leader_mi$stop),length=99)
FunctionIntervalM <- function(a,b) {
  seq(from=min(a), to = max(a), by = (max(a)-min(a))/b)
}
cuts <- FunctionIntervalM(a = alltimes, b = 5)

# AG model, piece-wise constant hazards
cut_data <- survSplit(Surv(start, stop, status == 1) ~ .,
                      leader_mi,
                      cut = cuts[-1],
                      episode = "timegroup")

coxph(Surv(start, stop, event) ~ treat + strata(timegroup), data = cut_data)

#------------------------------------------------------------------#
# Figure 2.14
#------------------------------------------------------------------#

#------------------------------------------------------------------#
# Plotting style
#------------------------------------------------------------------#

library(ggplot2)
theme_general <- theme_bw() +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2,"line"),
        text = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 

# Nelson-Aalen estimate of the cumulative hazard of MI
fit <- survfit(Surv(start, stop, status == 1) ~ treat, data = leader_mi,
               id = id, ctype = 1)

plotdata <- data.frame(cumhaz = fit$cumhaz, time = fit$time,
                       treat = c(rep("Placebo", fit$strata[[1]]),
                                 rep("Liraglutide", fit$strata[[2]])))

fig2.15 <- ggplot(aes(x = time * 1 / (365.25 / 12), y = cumhaz), data = plotdata) +
  geom_step(aes(linetype = treat), size = 1) +
  scale_linetype_discrete("Treatment") +
  ylab("Cumulative hazard") +
  xlab("Time since randomization (months)") +
  theme_bw() +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)),
                     limits = c(0, 65), breaks = seq(0, 60, by = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)),
                     limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.02)) +
  theme_general

fig2.15
