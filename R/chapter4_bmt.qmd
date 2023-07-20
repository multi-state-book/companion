#------------------------------------------------------------------#
#------- Chapter 4, R code, bmt data ------------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

# Load required packages (should be installed if not already)
require(haven)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

# Read data and create a data frame
bmt <- read.csv("data/bmt.csv")
bmt <- data.frame(bmt)

# See the head of the data
head(bmt)

# Summary of the data
summary(bmt)

# Reformatting and adding extra variables
bmt$intxsurv <- bmt$timedeath
bmt$dead <- bmt$death
bmt$intxrel <- ifelse(bmt$rel == 1, bmt$timerel, bmt$timedeath)
bmt$trm <- ifelse(bmt$rel == 0 & bmt$death == 1, 1, 0)
bmt$tgvhd <- ifelse(bmt$gvhd == 1, bmt$timegvhd, bmt$intxrel)
bmt$state0 <- bmt$rel + 2*bmt$trm

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
#---------------- Table 4.8 ---------------------------------------#
#------------------------------------------------------------------#

# Relapse
coxph(Surv(intxrel, rel == 1) ~ bmonly + all + age, data = bmt, 
      ties = "breslow")

coxph(Surv(intxrel, rel == 1) ~ bmonly + all + age + cluster(team), 
      data = bmt, 
      ties = "breslow")

# Relapse-free survival
coxph(Surv(intxrel, state0 != 0) ~ bmonly + all + age, data = bmt, 
      ties = "breslow")

coxph(Surv(intxrel, state0 != 0) ~ bmonly + all + age + cluster(team), 
      data = bmt, 
      ties = "breslow")

# overall survival
coxph(Surv(intxsurv, dead != 0) ~ bmonly + all + age, data = bmt, 
      ties = "breslow")

coxph(Surv(intxsurv, dead != 0) ~ bmonly + all + age + cluster(team), 
      data = bmt, 
      ties = "breslow")


#------------------------------------------------------------------#
#---------------- Figure 4.17 -------------------------------------#
#------------------------------------------------------------------#

# Relapse-free survival 
fit1 <- survfit(Surv(intxrel, state0 != 0) ~ 1, data = bmt)

# relapse
require(mets)
fit2 <- cif(Event(intxrel, state0) ~ 1, data = bmt, cause = 1)

# death in remission
fit3 <- cif(Event(intxrel, state0) ~ 1, data = bmt, cause = 2)

# overall survival
fit4 <- survfit(Surv(intxsurv, dead == 1) ~ 1, data = bmt)


# We need the same time for all probabilities
require(dplyr)
require(tidyr)
m1 <- stepfun(x = fit1$time, y = c(1, fit1$surv)) 
m2 <- stepfun(x = fit2$times, y = c(0, fit2$mu))

m3 <- stepfun(x = fit3$times, y = c(0, fit3$mu))
m4 <- stepfun(x = fit4$time, y = c(0, 1-fit4$surv))



unitimes <- sort(unique(c(fit1$time, fit2$times, fit3$times, fit4$time)))
m <- data.frame(time = unitimes, 
                q0 = m1(unitimes),
                c1 = m2(unitimes), 
                c2 = m3(unitimes), 
                c23 = m4(unitimes))

# tail(m %>% group_by(time) %>%
#   fill(q0, .direction = "down") %>%
#   ungroup())


m$q2 <-m$c2
m$q3 <- m$c23 - m$c2
m$q1 <- m$c1 - m$q3
m$sum <- with(m, q0+q1+q2+q3)

m$prev <- with(m, q1 / (q0 + q1))


# Overview
head(m); tail(m)

# Prepare data for plotting
plotdata <- with(m, 
                 data.frame(time = c(time, time), 
                            prob = c(prev, q1),
                            type = c(rep("Prevalence of relapse", length(time)), 
                                     rep("Probability of being alive with relapse", length(time)))
                 ))


# Create Figure 4.17
fig417 <- ggplot(aes(x = time, y = prob, linetype = type), data = plotdata) + 
  geom_step(linewidth = 1) + 
  scale_linetype_discrete("Type") + 
  xlab("Time since bone marrow transplantation (months)") + 
  ylab("Probability") + 
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)), 
                     limits = c(0, 156), breaks = seq(0, 156, by = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)), 
                     limits = c(0, 0.05), 
                     breaks = seq(0, 0.05, 0.01)) +
  theme_general + 
  theme(legend.box = "vertical",
        text = element_text(size=21), 
        legend.key.size = unit(1, 'cm'))

fig417
fig417<-fig417+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig417

ggsave("figures/h_bmtrelapse.pdf", plot = fig417, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Table 4.11 --------------------------------------#
#------------------------------------------------------------------#

# Relapse-free survival 
summary(coxph(Surv(intxrel, state0 > 0) ~ bmonly + all + age, 
      data = bmt, 
      ties = "breslow"))

# Overall survival
summary(coxph(Surv(intxsurv, dead != 0) ~ bmonly + all + age, 
              data = bmt, 
              ties = "breslow"))

bmt1 <- bmt
bmt1$gvhdny <- bmt$gvhd
bmt1$nytgvhd <- ifelse(bmt1$gvhdny == 1, bmt1$tgvhd, bmt1$intxsurv)
bmt1$gvhdny <- ifelse(bmt1$dead == 1, 1, bmt1$gvhdny)

# gvhd free survival
summary(coxph(Surv(nytgvhd, gvhdny != 0) ~ bmonly + all + age, 
              data = bmt1, 
              ties = "breslow"))

