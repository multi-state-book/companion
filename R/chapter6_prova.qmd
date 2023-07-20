#------------------------------------------------------------------#
#------- Chapter 6, R code, PROVA data ----------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")


# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

prova <- read.csv("data/prova.csv", na.strings = c("."))
head(prova)

prova <- data.frame(prova)
prova <- prova %>% mutate(timebleed = ifelse(bleed == 1, timebleed, timedeath),
                          outof0 = ifelse(bleed ==1, 1, death),
                          wait = ifelse(bleed ==1, timedeath - timebleed, NA))

# Faster processing
require(data.table)
prova <- setDT(prova)

# Setting start times, s = 1 year and s = 2 years
s1 = 365.25
s2 = 2*365.25

#Censoring times
cens_outof0 <- 1509
cens_wait <- 1363

#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------#
#------------------------------------------------------------------#

theme_general <- theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

#------------------------------------------------------------------#
#---------------- Figure 6.10 -------------------------------------#
#------------------------------------------------------------------#
library(survival) #Cox models, Kaplan-Meier, ect. 
library(tidyverse) #Plots
library(data.table) #Faster aggregation of data
library(mstate)
library(xtable)

# Creating data frame
# Converting data to long format
# Transition matrix for irreversible illness-death model
tmat <- trans.illdeath(names = c("Non-bleeding", "Bleeding", "Dead"))
long_format_scle0 <- msprep(time = c(NA, "timebleed", "timedeath"), status = c(NA, "bleed", "death"), data = subset(as.data.frame(prova), scle ==0), trans = tmat, keep = "scle")
long_format_scle1 <- msprep(time = c(NA, "timebleed", "timedeath"), status = c(NA, "bleed", "death"), data = subset(as.data.frame(prova), scle ==1), trans = tmat, keep = "scle")

# Estimation of transition probabilities with the Landmark Aalen-Johansen estimator
# Warning is not related to the estimate of pstate2
(LMAaJ_scle0 <- LMAJ(long_format_scle0, s = s1, from = 1)) #P(V(t) = 1 | V(1) = 0)
(LMAaJ_scle1 <- LMAJ(long_format_scle1, s = s1, from = 1)) #P(V(t) = 1 | V(1) = 0)

LMAJ_scle0_data <- as.data.frame(cbind(LMAaJ_scle0$time, LMAaJ_scle0$pstate2, "No"))
colnames(LMAJ_scle0_data) <- c("time", "p01", "Sclerotherapy")
LMAJ_scle1_data <- as.data.frame(cbind(LMAaJ_scle1$time, LMAaJ_scle1$pstate2, "Yes"))
colnames(LMAJ_scle1_data) <- c("time", "p01", "Sclerotherapy")

LMAJ_scle <- as.data.frame(rbind(LMAJ_scle0_data, LMAJ_scle1_data))


LMAJ_scle <- read.csv("data/datafig610.csv")

fig610 <- ggplot(LMAJ_scle,aes(x = as.numeric(time)/365.25, y = as.numeric(p01), group = Sclerotherapy)) + 
  geom_step(size = 1) +
  xlab("Time since randomization (years)") +
  ylab("Probability") + 
  theme_bw() +
  scale_x_continuous(expand = expansion(),limits = c(0.94,4.2)) +  
  aes(linetype=Sclerotherapy) + theme_general + 
  theme(legend.key.width = unit(1,"cm"),legend.text = element_text(size = 20)) 
fig610

ggsave("figures/h_provascle.pdf", plot = fig610,
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Table 6.7 ---------------------------------------#
#------------------------------------------------------------------#

## SEE REFERENCE PAPER




