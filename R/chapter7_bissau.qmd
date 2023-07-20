#------------------------------------------------------------------#
#------- Chapter 8, R code, Bissau data ---------------------------#
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
bissau$agem <- as.integer(bissau$age/30.44)

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
#---------------- Table 8.1 ---------------------------------------#
#------------------------------------------------------------------#

# Cox model - full cohort
coxph(Surv(fuptime, dead != 0) ~ factor(bcg) + factor(agem), data = bissau)

# 12.5 pct sub-cohort
bissau$s <- rbinom(nrow(bissau), 1, prob = 0.125)

with(bissau, table(dead, s))

# Case-cohort data
bis2 <- subset(bissau, s == 0 & dead != 0 | s == 1)
ccfit <- cch(Surv(fuptime, dead!= 0) ~ factor(bcg) + factor(agem), 
             data = bis2, 
             subcoh=bis2$s, id=bis2$id, cohort.size=nrow(bissau))
summary(ccfit)


# Nested-case control
require(multipleNCC)
require(Epi)


# Add noise to remove ties
bis2$fuptime_noise <- jitter(bis2$fuptime, factor=1, amount = NULL)

nccdata <- Epi::ccwc(exit=fuptime_noise, 
                     fail=dead!=0, 
                     data=bis2, 
                     include=list(bcg, agem), controls=3, silent=TRUE)
head(nccdata)

summary(clogit(Fail ~ factor(bcg) + factor(agem), data=nccdata))

