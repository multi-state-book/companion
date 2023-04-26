#------------------------------------------------------------------#
#------- Chapter 2, R code, Bissau data ---------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

bissau <- read.csv("data/bissau.csv")
head(bissau)

bissau <- data.frame(bissau)


#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------# 
#------------------------------------------------------------------#

theme_general <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 

theme_general

#------------------------------------------------------------------#
#---------------- Table 2.12 --------------------------------------#
#------------------------------------------------------------------#

# Cox model in column 1
bissau$agem <- as.integer(bissau$age/30.44)
coxfit1 <- coxph(Surv(fuptime, dead) ~ as.factor(bcg) + agem, 
                 data = bissau, method = "breslow")

summary(coxfit1)

# Make age the time variable instead
bissau$age12 <- bissau$age/(365.24/12)
bissau$ageout <- bissau$age12 + bissau$fuptime/(365.24/12)

# Cox model in column 2
coxfit2 <- coxph(Surv(age12, ageout, dead) ~ as.factor(bcg), 
                 data = bissau, method = "breslow")

summary(coxfit2)

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
fig212 <- ggplot(aes(x = time / (365.24/12), y = cumhaz), data = coxdata) + 
  geom_step(size = 1) + 
  xlab("Time since randomization (months)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), breaks = seq(0, 7, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  theme_general

fig212


ggsave("figures/j_bissaubasefupHENRIK.pdf", plot = fig212, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 2.13 -------------------------------------#
#------------------------------------------------------------------#

# Extract cumulative baseline hazard
coxcumhaz <- basehaz(coxfit2, centered = FALSE)

# Collect data for plot
coxdata <- data.frame(cumhaz = coxcumhaz$hazard, 
                      time = coxcumhaz$time, 
                      bcg = rep("1", nrow(coxcumhaz)), 
                      type = rep("Breslow estimate", nrow(coxcumhaz)))

# Create figure
fig213 <- ggplot(aes(x = time, y = cumhaz), data = coxdata) + 
  geom_step(size = 1) + 
  xlab("Age (months)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), breaks = seq(0, 15, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1)), breaks = seq(0, 0.15, by = 0.05), 
                     labels = c("0", "0.05", "0.10", "0.15")) +
  theme_general

fig213


ggsave("figures/j_bissaubaseageHENRIK.pdf", plot = fig213, 
       width = 29.7, height = 21, units = "cm")






