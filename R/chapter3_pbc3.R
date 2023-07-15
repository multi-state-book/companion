#------------------------------------------------------------------#
#------- Chapter 3, R code, PBC3 data  ----------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")


# Load required packages (should be installed if not already)
require(haven)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

# Read pcb data and create a data frame
pbc3 <- read.csv("data/pbc3.csv")
pbc3 <- data.frame(pbc3)

# See the head of the data
head(pbc3)

# Summary of the data
summary(pbc3)

# tment = 0 (placebo), tment = 1 (cyA)
with(pbc3, table(status, tment))

pbc3$tment_char <- ifelse(pbc3$tment == 0, "Placebo", "CyA")

# Add transformations of covariates 
pbc3$albnorm <- with(pbc3, (alb-35)*(alb>35))
pbc3$alb10 <- with(pbc3, alb/10)
pbc3$alb2 <- with(pbc3, alb10*alb10)

pbc3$bilihigh <- with(pbc3, (bili-17.1)*(bili>17.1))
pbc3$bilitoohigh <- with(pbc3, (bili-34.2)*(bili>34.2))
pbc3$bilimuchtoohigh <- with(pbc3, (bili-51.3)*(bili>51.3))
pbc3$bili100 <- with(pbc3, bili/100)
pbc3$bili2 <- with(pbc3, bili100*bili100)

pbc3$log2bili <- with(pbc3, log2(bili))
pbc3$logbilihigh <- with(pbc3, (log2bili-log2(17.1))*(bili>17.1))
pbc3$logbilitoohigh <- with(pbc3, (log2bili-log2(34.2))*(bili>34.2))
pbc3$logbilimuchtoohigh <- with(pbc3, (log2bili-log2(51.3))*(bili>51.3))
pbc3$log2bili2 <- with(pbc3, log2bili*log2bili)



#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------# 
#------------------------------------------------------------------#

# General theme
theme_general <- theme_bw() +
                 theme(legend.position = "bottom", 
                       text = element_text(size = 20), 
                       axis.text.x = element_text(size = 20), 
                       axis.text.y = element_text(size = 20)) 
                 
theme_general

#------------------------------------------------------------------#
#---------------- Table 3.10 --------------------------------------#
#------------------------------------------------------------------#


# Treatment 
summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(tment), 
              data = pbc3, tt = function(x,t, ...) (x==1)*t, 
              method = "breslow"))


summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(tment), 
              data = pbc3, tt = function(x,t, ...) (x==1)*log(t), 
              method = "breslow"))


summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(tment), 
              data = pbc3, tt = function(x,t, ...) (x==1)*(t > 2 * 365.25), 
              method = "breslow"))


# Log bili
summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(log2bili), 
              data = pbc3, tt = function(x,t, ...) x*t, 
              method = "breslow"))


summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(log2bili), 
              data = pbc3, tt = function(x,t, ...) x*log(t), 
              method = "breslow"))


summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(log2bili), 
              data = pbc3,  tt = function(x,t, ...) (x)*(t > 2 * 365.25), 
              method = "breslow"))



# Alb 
summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(alb), 
              data = pbc3, tt = function(x,t, ...) x*t, 
              method = "breslow"))


summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(alb), 
              data = pbc3, tt = function(x,t, ...) x*log(t), 
              method = "breslow"))


summary(coxph(Surv(days, status != 0) ~ as.factor(tment) + alb + log2bili + tt(alb), 
              data = pbc3,  tt = function(x,t, ...) (x)*(t > 2 * 365.25), 
              method = "breslow"))



