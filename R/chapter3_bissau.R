#------------------------------------------------------------------#
#------- Chapter 3, R code, Bissau data ---------------------------#
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
#---------------- Table 3.1 ---------------------------------------#
#------------------------------------------------------------------#

table(bissau$bcg, bissau$dtp)

table(bissau$bcg, bissau$dtp) / rowSums(table(bissau$bcg, bissau$dtp))

#------------------------------------------------------------------#
#---------------- Table 3.2 ---------------------------------------#
#------------------------------------------------------------------#

# Add extra variables
bissau$age12 <- with(bissau, age/(365.24/12))
bissau$ageout <- with(bissau, age12+fuptime/(365.24/12))
bissau$dtpany <- with(bissau, dtp>0)
bissau$bcg <- as.factor(bissau$bcg)
bissau$dtpany <- as.factor(bissau$dtpany)

bissau$bcg <- as.factor(bissau$bcg)

mod1 <- coxph(Surv(age12, ageout, dead != 0) ~ bcg, data = bissau, 
              method = "breslow", timefix = FALSE)
summary(mod1)


mod2 <- coxph(Surv(age12, ageout, dead != 0) ~ dtpany, data = bissau, 
              method = "breslow", timefix = FALSE)
summary(mod2)


mod3 <- coxph(Surv(age12, ageout, dead != 0) ~ bcg + dtpany, data = bissau, 
              method = "breslow", timefix = FALSE)
summary(mod3)


mod4 <- coxph(Surv(age12, ageout, dead != 0) ~ bcg * dtpany, data = bissau, 
              method = "breslow", timefix = FALSE)
summary(mod4)
