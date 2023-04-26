#------------------------------------------------------------------#
#------- Chapter 2, R code, affective data ------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")


# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)


affective <- read.csv("data/affective.csv")
head(affective)

affective <- data.frame(affective)

affective$wait <- with(affective, stop - start)

#------------------------------------------------------------------#
#---------------- Table 2.13 --------------------------------------#
#------------------------------------------------------------------#

# Cox model for 1., 2., 3., 4. episode 'Markov': Column 1
coxfit1_markov <- coxph(Surv(start, stop, status == 1) ~ as.factor(bip), 
                        data = subset(affective, episode == 1 & state == 0), 
                        method = "breslow")

summary(coxfit1_markov)


coxfit2_markov <- coxph(Surv(start, stop, status == 1) ~ as.factor(bip), 
                        data = subset(affective, episode == 2 & state == 0), 
                        method = "breslow")

summary(coxfit2_markov)

coxfit3_markov <- coxph(Surv(start, stop, status == 1) ~ as.factor(bip), 
                        data = subset(affective, episode == 3 & state == 0), 
                        method = "breslow")

summary(coxfit3_markov)

coxfit4_markov <- coxph(Surv(start, stop, status == 1) ~ as.factor(bip), 
                        data = subset(affective, episode == 4 & state == 0), 
                        method = "breslow")

summary(coxfit4_markov)

# Cox model for 1., 2., 3., 4. episode 'Gap time': Column 2
coxfit1_gap <- coxph(Surv(wait, status == 1) ~ as.factor(bip), 
                        data = subset(affective, episode == 1 & state == 0), 
                        method = "breslow")

summary(coxfit1_gap)

coxfit2_gap <- coxph(Surv(wait, status == 1) ~ as.factor(bip), 
                     data = subset(affective, episode == 2 & state == 0), 
                     method = "breslow")

summary(coxfit2_gap)

coxfit3_gap <- coxph(Surv(wait, status == 1) ~ as.factor(bip), 
                     data = subset(affective, episode == 3 & state == 0), 
                     method = "breslow")

summary(coxfit3_gap)

coxfit4_gap <- coxph(Surv(wait, status == 1) ~ as.factor(bip), 
                     data = subset(affective, episode == 4 & state == 0), 
                     method = "breslow")

summary(coxfit4_gap)


# AG cox model, total time
coxfit_ag <- coxph(Surv(start, stop, status == 1) ~ as.factor(bip), 
                        data = subset(affective, state == 0), 
                        method = "breslow")

summary(coxfit_ag)


# AG cox model, gap time
coxfit_ag_gap <- coxph(Surv(wait, status == 1) ~ as.factor(bip), 
                   data = subset(affective, state == 0), 
                   method = "breslow")

summary(coxfit_ag_gap)


# PWP cox model, total time
coxfit_pwp <- coxph(Surv(start, stop, status == 1) ~ strata(episode) + as.factor(bip), 
                   data = subset(affective, state == 0), 
                   method = "breslow")

summary(coxfit_pwp)

# PWP cox model, gap time
coxfit_pwp_gap <- coxph(Surv(wait, status == 1) ~ strata(episode) + as.factor(bip), 
                    data = subset(affective, state == 0), 
                    method = "breslow")

summary(coxfit_pwp_gap)












