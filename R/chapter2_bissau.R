bissau <- data.frame(read.csv("data/bissau.csv"))

#------------------------------------------------------------------#
# General plotting style 
#------------------------------------------------------------------#
library(ggplot2)
theme_general <- theme_bw() +
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        legend.position = "bottom", 
        legend.title=element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2,"line")
  ) 

#------------------------------------------------------------------#
# Table 2.12 and in-text estimate
#------------------------------------------------------------------#

# Cox model in column 1
bissau$agem <- as.integer(bissau$age/30.44)
coxfit0 <- coxph(Surv(fuptime, dead) ~ bcg, 
                 data = bissau, method = "breslow")
tidy(coxfit0)

# Cox model in column 1
coxfit1 <- coxph(Surv(fuptime, dead) ~ bcg + agem, 
                 data = bissau, method = "breslow")
tidy(coxfit1)

# Make age the time variable instead
bissau$agein <- bissau$age/(365.24/12)
bissau$ageout <- bissau$agein + bissau$fuptime/(365.24/12)

# Cox model in column 2
coxfit2 <- coxph(Surv(agein, ageout, dead) ~ bcg, 
                 data = bissau, method = "breslow")
tidy(coxfit2)

#------------------------------------------------------------------#
#---------------- Figure 2.12 -------------------------------------#
#------------------------------------------------------------------#

# Extract cumulative baseline hazard
coxcumhaz <- basehaz(coxfit1, centered = FALSE)

# Collect data for plot
coxdata <- data.frame(cumhaz = coxcumhaz$hazard, 
                      time = coxcumhaz$time, 
                      bcg = rep("1", nrow(coxcumhaz)), 
                      type = rep("Breslow estimate", nrow(coxcumhaz)))

# Create figure
fig2.12 <- ggplot(aes(x = time / (365.24/12), y = cumhaz), data = coxdata) + 
  geom_step(size = 1) + 
  xlab("Follow-up time (months)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), breaks = seq(0, 7, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  theme_general

fig2.12


#------------------------------------------------------------------#
# Figure 2.13 
#------------------------------------------------------------------#

# Extract cumulative baseline hazard
coxcumhaz <- basehaz(coxfit2, centered = FALSE)
# Collect data for plot
coxdata <- data.frame(cumhaz = coxcumhaz$hazard, 
                      time = coxcumhaz$time, 
                      bcg = rep("1", nrow(coxcumhaz)), 
                      type = rep("Breslow estimate", nrow(coxcumhaz)))

# Create figure
fig2.13 <- ggplot(aes(x = time, y = cumhaz), data = coxdata) + 
  geom_step(size = 1) + 
  xlab("Age (months)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), breaks = seq(0, 15, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1)), breaks = seq(0, 0.15, by = 0.05), 
                     labels = c("0", "0.05", "0.10", "0.15")) +
  theme_general

fig2.13





