#------------------------------------------------------------------#
#------- Chapter 4, R code, PROVA data ----------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")


# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

provany <- read.csv("data/prova.csv", na.strings = c("."))
head(provany)


provany <- data.frame(provany)

# Structure
str(provany)
summary(provany)

provany$beh <- with(provany, scle*2 + beta)
provany$log2bili <- with(provany, log2(bili))

#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------#
#------------------------------------------------------------------#

theme_general <- theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

#------------------------------------------------------------------#
#---------------- Figure 4.25 -------------------------------------#
#------------------------------------------------------------------#

# Make KM estimate of censoring distribution
provany$time <- with(provany, ifelse(!is.na(timebleed), timebleed + timedeath, timedeath))
censdist <- survfit(Surv(time, death == 0) ~ 1, data = provany)
censdist

# Make data ready for plotting
pdata <- data.frame(time = censdist$time,
                    surv = censdist$surv)


# Create Figure 4.24
fig424 <- ggplot(aes(x = time / 365.25, y = surv),
                 data = pdata) +
  geom_step(size = 1) +
  xlab("Time since randomization (years)") +
  ylab('Probability of no censoring') +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005)),
                     limits = c(0, 7),
                     breaks = seq(0, 7, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.005)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  theme(legend.box = "vertical",
        legend.key.size = unit(1.5, 'cm'))


fig424


ggsave("figures/j_provacens.pdf", plot = fig424,
       width = 29.7, height = 21, units = "cm")

#------------------------------------------------------------------#
#---------------- Table 4.13 --------------------------------------#
#------------------------------------------------------------------#

# null model
censcoxnull <- coxph(Surv(time, death == 0) ~ 1, data = provany)
summary(censcoxnull)

# treat
censcoxtreat <- coxph(Surv(time, death == 0) ~ beh, data = provany)
summary(censcoxtreat)

require(lmtest)
lrtest(censcoxtreat, censcoxnull)


# size 
censcoxvapr <- coxph(Surv(time, death == 0) ~ varsize, data = provany)
summary(censcoxvapr)

lrtest(censcoxvapr, censcoxnull)

# sex
censcoxsex <- coxph(Surv(time, death == 0) ~ sex, data = provany)
summary(censcoxsex)

# coag
censcoxkoag <- coxph(Surv(time, death == 0) ~ coag, data = provany)
summary(censcoxkoag)

# bili
censcoxbili <- coxph(Surv(time, death == 0) ~ log2bili, data = provany)
summary(censcoxbili)

# age
censcoxage <- coxph(Surv(time, death == 0) ~ age, data = provany)
summary(censcoxage)


