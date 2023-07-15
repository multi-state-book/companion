#------------------------------------------------------------------#
# Table 2.13 
#------------------------------------------------------------------#

affective <- data.frame(read.csv("data/affective.csv"))
affective$wait <- with(affective, stop - start)


# Cox model for 1., 2., 3., 4. episode 'Markov': Column 1
library(survival)
tidy(coxph(Surv(start, stop, status == 1) ~ bip, 
                        data = subset(affective, episode == 1 & state == 0), 
                        method = "breslow"))

tidy(coxph(Surv(start, stop, status == 1) ~ bip, 
                        data = subset(affective, episode == 2 & state == 0), 
                        method = "breslow"))

tidy(coxph(Surv(start, stop, status == 1) ~ bip, 
                        data = subset(affective, episode == 3 & state == 0), 
                        method = "breslow"))

tidy(coxph(Surv(start, stop, status == 1) ~ bip, 
                        data = subset(affective, episode == 4 & state == 0), 
                        method = "breslow"))

# Cox model for 1., 2., 3., 4. episode 'Gap time': Column 2
tidy(coxph(Surv(wait, status == 1) ~ bip, 
                        data = subset(affective, episode == 1 & state == 0), 
                        method = "breslow"))

tidy(coxph(Surv(wait, status == 1) ~ bip, 
                     data = subset(affective, episode == 2 & state == 0), 
                     method = "breslow"))

tidy(coxph(Surv(wait, status == 1) ~ bip, 
                     data = subset(affective, episode == 3 & state == 0), 
                     method = "breslow"))

tidy(coxph(Surv(wait, status == 1) ~ bip, 
                     data = subset(affective, episode == 4 & state == 0), 
                     method = "breslow"))


# AG cox model, total time
tidy(coxph(Surv(start, stop, status == 1) ~ bip, 
                        data = subset(affective, state == 0), 
                        method = "breslow"))

# AG cox model, gap time
tidy(coxph(Surv(wait, status == 1) ~ bip, 
                   data = subset(affective, state == 0), 
                   method = "breslow"))

# PWP cox model, total time
tidy(coxph(Surv(start, stop, status == 1) ~ strata(episode) + bip, 
                   data = subset(affective, state == 0), 
                   method = "breslow"))

# PWP cox model, gap time
tidy(coxph(Surv(wait, status == 1) ~ strata(episode) + bip, 
                    data = subset(affective, state == 0), 
                    method = "breslow"))
