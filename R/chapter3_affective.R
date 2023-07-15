#------------------------------------------------------------------#
#------- Chapter 3, R code, affective data ------------------------#
#------------------------------------------------------------------#

# Load required packages (should be installed if not already)
library(survival)

affective <- data.frame(read.csv("data/affective.csv"))
affective$wait <- with(affective, stop - start)

#------------------------------------------------------------------#
#---------------- Table 3.6 ---------------------------------------#
#------------------------------------------------------------------#

coxph(Surv(start, stop, status == 1) ~ bip,
      data = subset(affective, state == 0), ties = "breslow")

coxph(Surv(start, stop, status == 1) ~ bip + episode,
      data = subset(affective, state == 0), ties = "breslow")

coxph(Surv(start, stop, status == 1) ~ bip + episode + I(episode*episode),
      data = subset(affective, state == 0), ties = "breslow")

