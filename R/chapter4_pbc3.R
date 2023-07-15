#------------------------------------------------------------------#
#------- Chapter 4, R code, PBC3 data  ----------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

# Load required packages (should be installed if not already)
require(haven)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)
require(mets)
require(cmprsk)

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

# 
theme_general <- theme_bw() +
                 theme(legend.position = "bottom", 
                       text = element_text(size = 20), 
                       axis.text.x = element_text(size = 20), 
                       axis.text.y = element_text(size = 20)) 
                 
theme_general

#------------------------------------------------------------------#
#---------------- Table 4.1 ---------------------------------------#
#------------------------------------------------------------------#

#### 3-year restricted means - ESTIMATES ####

# Non-parametric, Kaplan Meier
np_km <- survfit(Surv(days, status != 0) ~ tment, data = pbc3)

kmdata <- data.frame(surv = np_km$surv, 
                      time = np_km$time, 
                     strata = c(rep(names(np_km$strata[1]), np_km$strata[[1]]), 
                              rep(names(np_km$strata[2]), np_km$strata[[2]])))
# Restrict to each treat
kmdata0 <- subset(kmdata, strata == "tment=0")
kmdata1 <- subset(kmdata, strata == "tment=1")

# rmst
require(RISCA)
rmst0 <- rmst(times = kmdata0$time, surv.rates = kmdata0$surv, max.time = 3 * 365.25, type = "s")
rmst1 <- rmst(times = kmdata1$time, surv.rates = kmdata1$surv, max.time = 3 * 365.25, type = "s")
rmst0; rmst1; 

# Cox model, alb = 38, bili = 45
# Cox model fit with covariates tment, alb and log2bili
coxfit <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili, data = pbc3, 
                method = "breslow")

# Unique followup times 
fu <- sort(unique(pbc3$days))

# Data for prediction
preddata1 <- data.frame(tment = c(rep(0, length(fu)), rep(1, length(fu))), 
                       alb = rep(38, length(fu) * 2),
                       log2bili = rep(log2(45), length(fu) * 2), 
                       days = c(fu, fu),
                       status = rep(1, length(fu) * 2)
)

# Linear predictor
preds1 <- predict(coxfit, newdata = preddata1, type = "survival")

preddata1$preds <- preds1

cox10 <- subset(preddata1, tment == "0")
cox11 <- subset(preddata1, tment == "1")

# Rmst 
rmstcox10 <- rmst(times = cox10$days, surv.rates = cox10$preds, max.time = 3 * 365.25, type = "s")
rmstcox11 <- rmst(times = cox11$days, surv.rates = cox11$preds, max.time = 3 * 365.25, type = "s")
rmstcox10; rmstcox11; 

# Cox for alb = 20 and bili = 90 

# Data for prediction
preddata2 <- data.frame(tment = c(rep(0, length(fu)), rep(1, length(fu))), 
                        alb = rep(20, length(fu) * 2),
                        log2bili = rep(log2(90), length(fu) * 2), 
                        days = c(fu, fu),
                        status = rep(1, length(fu) * 2)
)

# Linear predictor
preds2 <- predict(coxfit, newdata = preddata2, type = "survival")

preddata2$preds <- preds2

cox20 <- subset(preddata2, tment == "0")
cox21 <- subset(preddata2, tment == "1")

# Rmst 
rmstcox20 <- rmst(times = cox20$days, surv.rates = cox20$preds, max.time = 3 * 365.25, type = "s")
rmstcox21 <- rmst(times = cox21$days, surv.rates = cox21$preds, max.time = 3 * 365.25, type = "s")
rmstcox20; rmstcox21; 


# Cox model, g-formula

# We want to predict responses if all subjects had been given both CyA and placebo (tment = 0,1)
# We need to have n survival curves * 2
# make a double data set with extra Z
pbc3_counterfact <- pbc3
pbc3_counterfact$tment <- ifelse(pbc3_counterfact$tment == 1, 0, 1) # opposite treatment
pbc3_double <- rbind(pbc3, pbc3_counterfact)

# Check - one of each treat!
with(pbc3_double, table(tment, id))

# Baseline survival
pred <- survfit(coxfit, newdata = data.frame(tment = 0, alb = 0, log2bili = 0))

allsurv <-
  lapply(1:nrow(pbc3_double),
         function(i)
           pred$surv ^
           exp(coef(coxfit)[1] * pbc3_double$tment[i] +
                 coef(coxfit)[2] * pbc3_double$alb[i] +
                 coef(coxfit)[3] * pbc3_double$log2bili[i]))

potout <-
  data.frame(surv = unlist(allsurv),
             tment = rep(pbc3_double$tment, each = length(pred$time)),
             time = rep(pred$time, times = nrow(pbc3)*2)
  )

# Average over values
require(dplyr)
sumdata <- potout %>%
  group_by(tment, time) %>%
  summarise(average_pred = mean(surv, na.rm = TRUE),
            .groups = c("keep"))
sumdata <- as.data.frame(sumdata)


# Split data per group
coxg0 <- subset(sumdata, tment == "0")
coxg1 <- subset(sumdata, tment == "1")

# Rmst 
rmstcoxg0 <- rmst(times = coxg0$time, surv.rates = coxg0$average_pred, max.time = 3 * 365.25, type = "s")
rmstcoxg1 <- rmst(times = coxg1$time, surv.rates = coxg1$average_pred, max.time = 3 * 365.25, type = "s")
rmstcoxg0; rmstcoxg1;

##### BOOTSTRAP FOR SE'S #####

# Resample data sets 

B <- 200
bootdata <- list()
kmres <- cox1res <- cox2res <- coxgres <- list()

colnames(kmres) <- colnames(cox1res) <- colnames(cox2res) <- colnames(coxgres) <- c("rmst0", "rmst1")


for (b in 1:B){
  bootdata[[b]] <- pbc3[sample(1:nrow(pbc3), size = nrow(pbc3), replace = T),]
  
  ###### KM #######
  np_km <- survfit(Surv(days, status != 0) ~ tment, data = bootdata[[b]])
  
  kmdata <- data.frame(surv = np_km$surv, 
                       time = np_km$time, 
                       strata = c(rep(names(np_km$strata[1]), np_km$strata[[1]]), 
                                  rep(names(np_km$strata[2]), np_km$strata[[2]])))
  # Restrict to each treat
  kmdata0 <- subset(kmdata, strata == "tment=0")
  kmdata1 <- subset(kmdata, strata == "tment=1")
  
  # rmst
  rmst0 <- rmst(times = kmdata0$time, surv.rates = kmdata0$surv, max.time = 3 * 365.25, type = "s")
  rmst1 <- rmst(times = kmdata1$time, surv.rates = kmdata1$surv, max.time = 3 * 365.25, type = "s")

  kmres[[b]] <- c(rmst0, rmst1)
  
  
  ###### Cox 1 #######
  
  coxfit <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili, data = bootdata[[b]], 
                  method = "breslow")
  
  # Unique followup times 
  fu <- sort(unique(pbc3$days))
  
  # Data for prediction
  preddata1 <- data.frame(tment = c(rep(0, length(fu)), rep(1, length(fu))), 
                          alb = rep(38, length(fu) * 2),
                          log2bili = rep(log2(45), length(fu) * 2), 
                          days = c(fu, fu),
                          status = rep(1, length(fu) * 2)
  )
  
  # Linear predictor
  preds1 <- predict(coxfit, newdata = preddata1, type = "survival")
  
  preddata1$preds <- preds1
  
  cox10 <- subset(preddata1, tment == "0")
  cox11 <- subset(preddata1, tment == "1")
  
  # Rmst 
  rmstcox10 <- rmst(times = cox10$days, surv.rates = cox10$preds, max.time = 3 * 365.25, type = "s")
  rmstcox11 <- rmst(times = cox11$days, surv.rates = cox11$preds, max.time = 3 * 365.25, type = "s")
  
  cox1res[[b]] <- c(rmstcox10, rmstcox11)
  
  ###### Cox 2 #######

  # Cox for alb = 20 and bili = 90   
  # Data for prediction
  preddata2 <- data.frame(tment = c(rep(0, length(fu)), rep(1, length(fu))), 
                          alb = rep(20, length(fu) * 2),
                          log2bili = rep(log2(90), length(fu) * 2), 
                          days = c(fu, fu),
                          status = rep(1, length(fu) * 2)
  )
  
  # Linear predictor
  preds2 <- predict(coxfit, newdata = preddata2, type = "survival")
  
  preddata2$preds <- preds2
  
  cox20 <- subset(preddata2, tment == "0")
  cox21 <- subset(preddata2, tment == "1")
  
  # Rmst 
  rmstcox20 <- rmst(times = cox20$days, surv.rates = cox20$preds, max.time = 3 * 365.25, type = "s")
  rmstcox21 <- rmst(times = cox21$days, surv.rates = cox21$preds, max.time = 3 * 365.25, type = "s")
  
  
  cox2res[[b]] <- c(rmstcox20, rmstcox21) 
  
  
  ##### Cox - g formula #######
  
  # We want to predict responses if all subjects had been given both CyA and placebo (tment = 0,1)
  # We need to have n survival curves * 2
  # make a double data set with extra Z
  pbc3_counterfact <- bootdata[[b]]
  pbc3_counterfact$tment <- ifelse(pbc3_counterfact$tment == 1, 0, 1) # opposite treatment
  pbc3_double <- rbind(bootdata[[b]], pbc3_counterfact)
  
  # Baseline survival
  pred <- survfit(coxfit, newdata = data.frame(tment = 0, alb = 0, log2bili = 0))
  
  allsurv <-
    lapply(1:nrow(pbc3_double),
           function(i)
             pred$surv ^
             exp(coef(coxfit)[1] * pbc3_double$tment[i] +
                   coef(coxfit)[2] * pbc3_double$alb[i] +
                   coef(coxfit)[3] * pbc3_double$log2bili[i]))
  
  potout <-
    data.frame(surv = unlist(allsurv),
               tment = rep(pbc3_double$tment, each = length(pred$time)),
               time = rep(pred$time, times = nrow(pbc3)*2)
    )
  
  # Average over values
  require(dplyr)
  sumdata <- potout %>%
    group_by(tment, time) %>%
    summarise(average_pred = mean(surv, na.rm = TRUE),
              .groups = c("keep"))
  sumdata <- as.data.frame(sumdata)
  
  
  # Split data per group
  coxg0 <- subset(sumdata, tment == "0")
  coxg1 <- subset(sumdata, tment == "1")
  
  # Rmst 
  rmstcoxg0 <- rmst(times = coxg0$time, surv.rates = coxg0$average_pred, max.time = 3 * 365.25, type = "s")
  rmstcoxg1 <- rmst(times = coxg1$time, surv.rates = coxg1$average_pred, max.time = 3 * 365.25, type = "s")
  rmstcoxg0; rmstcoxg1;
  
  coxgres[[b]] <- c(rmstcoxg0, rmstcoxg1)
  
}

kmreso <- do.call("rbind", kmres)
cox1reso <- do.call("rbind", cox1res)
cox2reso <- do.call("rbind", cox2res)
coxgreso <- do.call("rbind", coxgres)

# For non-parametric, KM
apply(kmreso, 2, mean)
apply(kmreso, 2, sd)

# For Cox 1
apply(cox1reso, 2, mean)
apply(cox1reso, 2, sd)

# For Cox 2
apply(cox2reso, 2, mean)
apply(cox2reso, 2, sd)

# For Cox g formula
apply(coxgreso, 2, mean)
apply(coxgreso, 2, sd)


#------------------------------------------------------------------#
#---------------- Table 4.2 ---------------------------------------#
#------------------------------------------------------------------#

# Reference levels 
pbc3$sex <- relevel(as.factor(pbc3$sex), ref = "1")

# Death without transplantation 
summary(coxph(Surv(days, status == 2) ~ tment + alb + log2bili + sex + age, data = pbc3, 
                method = "breslow"))

# Transplantation 
summary(coxph(Surv(days, status == 1) ~ tment + alb + log2bili + sex + age, data = pbc3, 
              method = "breslow"))


# Failure of medical treatment 
summary(coxph(Surv(days, status %in% c(1, 2)) ~ tment + alb + log2bili + sex + age, data = pbc3, 
              method = "breslow"))



#------------------------------------------------------------------#
#---------------- Table 4.3 ---------------------------------------#
#------------------------------------------------------------------#

# No adjustment 
surv <- coxph(Surv(days, status != 0) ~ strata(tment), data = pbc3, 
                      method = "breslow")

# death w/o transplant
cox1 <- coxph(Surv(days, status == 2) ~ strata(tment), data = pbc3, 
              method = "breslow")

# transplant
cox2 <- coxph(Surv(days, status == 1) ~ strata(tment), data = pbc3, 
              method = "breslow")

time <- basehaz(cox1, center = F)$time
surv <- exp(-basehaz(surv, center = F)$hazard)
lamd1 <- basehaz(cox1, center = F)$hazard
lamd2 <- basehaz(cox2, center = F)$hazard
strat <- basehaz(cox1, center = F)$strata
  
# Cumulative incidence
F01_plac <- cumsum(surv[strat == "tment=0"] * diff(c(0,lamd1[strat == "tment=0"])))
F01_cya <- cumsum(surv[strat == "tment=1"] * diff(c(0,lamd1[strat == "tment=1"])))

F02_plac <- cumsum(surv[strat == "tment=0"] * diff(c(0,lamd2[strat == "tment=0"])))
F02_cya <- cumsum(surv[strat == "tment=1"] * diff(c(0,lamd2[strat == "tment=1"])))

timep <- time[strat == "tment=0"]
timec <- time[strat == "tment=1"]

# Years lost, AUC

# F01 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F01_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F01_cya[timec <= 3 * 365.25])

# F02 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F02_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F02_cya[timec <= 3 * 365.25])


# Sex = F, age = 40, alb = 38, bili = 45
surv <- coxph(Surv(days, status != 0) ~ tment + sex + age + alb + log2bili, data = pbc3, 
              method = "breslow")

# death w/o transplant
cox1 <- coxph(Surv(days, status == 2) ~ tment + sex + age + alb + log2bili, data = pbc3, 
              method = "breslow")

# transplant
cox2 <- coxph(Surv(days, status == 1) ~ tment + sex + age + alb + log2bili, data = pbc3, 
              method = "breslow")


# Predictions 
newd <- data.frame(tment = c(0, 1),
                   sex = c(0, 0),
                   alb = c(38, 38), 
                   age = c(40, 40),
                   log2bili = log2(c(45, 45)))

# predictions
lpsurv <- predict(surv, newd, type = "lp", reference = "zero")
lpcox1 <- predict(cox1, newd, type = "lp", reference = "zero")
lpcox2 <- predict(cox2, newd, type = "lp", reference = "zero")

time <- basehaz(cox1, center = F)$time
surv <- exp(-basehaz(surv, center = F)$hazard)
lamd1 <- basehaz(cox1, center = F)$hazard
lamd2 <- basehaz(cox2, center = F)$hazard

# Linear predictors
surv_plac <- surv^exp(lpsurv[1])
surv_cya <- surv^exp(lpsurv[2])

lamd1_plac <- lamd1 * exp(lpcox1[1])
lamd1_cya <- lamd1 * exp(lpcox1[2])

lamd2_plac <- lamd2 * exp(lpcox2[1])
lamd2_cya <- lamd2 * exp(lpcox2[2])

# Cumulative incidence
F01_plac <- cumsum(surv_plac * diff(c(0,lamd1_plac)))
F01_cya <- cumsum(surv_cya * diff(c(0,lamd1_cya)))

F02_plac <- cumsum(surv_plac * diff(c(0,lamd2_plac)))
F02_cya <- cumsum(surv_cya * diff(c(0,lamd2_cya)))

timep <- timec <- time

# Years lost, AUC

# F01 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F01_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F01_cya[timec <= 3 * 365.25])

# F02 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F02_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F02_cya[timec <= 3 * 365.25])


##### Sex = F, age = 40, alb = 20, bili = 90 #####

# Predictions 
newd <- data.frame(tment = c(0, 1),
                   sex = c(0, 0),
                   alb = c(20, 20), 
                   age = c(40, 40),
                   log2bili = log2(c(90, 90)))

# predictions
lpsurv <- predict(surv, newd, type = "lp", reference = "zero")
lpcox1 <- predict(cox1, newd, type = "lp", reference = "zero")
lpcox2 <- predict(cox2, newd, type = "lp", reference = "zero")

time <- basehaz(cox1, center = F)$time
surv <- exp(-basehaz(surv, center = F)$hazard)
lamd1 <- basehaz(cox1, center = F)$hazard
lamd2 <- basehaz(cox2, center = F)$hazard

# Linear predictors
surv_plac <- surv^exp(lpsurv[1])
surv_cya <- surv^exp(lpsurv[2])

lamd1_plac <- lamd1 * exp(lpcox1[1])
lamd1_cya <- lamd1 * exp(lpcox1[2])

lamd2_plac <- lamd2 * exp(lpcox2[1])
lamd2_cya <- lamd2 * exp(lpcox2[2])

# Cumulative incidence
F01_plac <- cumsum(surv_plac * diff(c(0,lamd1_plac)))
F01_cya <- cumsum(surv_cya * diff(c(0,lamd1_cya)))

F02_plac <- cumsum(surv_plac * diff(c(0,lamd2_plac)))
F02_cya <- cumsum(surv_cya * diff(c(0,lamd2_cya)))

timep <- timec <- time

# Years lost, AUC

# F01 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F01_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F01_cya[timec <= 3 * 365.25])

# F02 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F02_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F02_cya[timec <= 3 * 365.25])




##### Sex = F, age = 60, alb = 38, bili = 45 #####

# Predictions 
newd <- data.frame(tment = c(0, 1),
                   sex = c(0, 0),
                   alb = c(38, 38), 
                   age = c(60, 60),
                   log2bili = log2(c(45, 45)))

# predictions
lpsurv <- predict(surv, newd, type = "lp", reference = "zero")
lpcox1 <- predict(cox1, newd, type = "lp", reference = "zero")
lpcox2 <- predict(cox2, newd, type = "lp", reference = "zero")

time <- basehaz(cox1, center = F)$time
surv <- exp(-basehaz(surv, center = F)$hazard)
lamd1 <- basehaz(cox1, center = F)$hazard
lamd2 <- basehaz(cox2, center = F)$hazard

# Linear predictors
surv_plac <- surv^exp(lpsurv[1])
surv_cya <- surv^exp(lpsurv[2])

lamd1_plac <- lamd1 * exp(lpcox1[1])
lamd1_cya <- lamd1 * exp(lpcox1[2])

lamd2_plac <- lamd2 * exp(lpcox2[1])
lamd2_cya <- lamd2 * exp(lpcox2[2])

# Cumulative incidence
F01_plac <- cumsum(surv_plac * diff(c(0,lamd1_plac)))
F01_cya <- cumsum(surv_cya * diff(c(0,lamd1_cya)))

F02_plac <- cumsum(surv_plac * diff(c(0,lamd2_plac)))
F02_cya <- cumsum(surv_cya * diff(c(0,lamd2_cya)))

timep <- timec <- time

# Years lost, AUC

# F01 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F01_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F01_cya[timec <= 3 * 365.25])

# F02 : death w/p transplant
# Placebo
sum(diff(c(0, timep[timep <= 3 * 365.25])) * F02_plac[timep <= 3 * 365.25])
# CyA
sum(diff(c(0, timec[timec <= 3 * 365.25])) * F02_cya[timec <= 3 * 365.25])





#------------------------------------------------------------------#
#---------------- Table 4.5 ---------------------------------------#
#------------------------------------------------------------------#

library("survRM2")


pbcny <- subset(pbc3,!is.na(alb))
time <- pbcny$days
status <- pbcny$status!=0
arm <- pbcny$tment
alb <- pbcny$alb

logbili <- log2(pbcny$bili)
x <- cbind(alb, logbili)
rmst2(time, status, arm, tau= 3 * 365.25, covariates=x)


#------------------------------------------------------------------#
#---------------- Figure 4.2 --------------------------------------#
#------------------------------------------------------------------#

# Kaplan-Meier estimate per treatment
kmfit <- survfit(Surv(days, status != 0) ~ tment, data = pbc3)

# Collect data for plot
# Note that the standard errors produced by survfit are for the cumulative hazard
kmdata <- data.frame(surv = kmfit$surv, 
                     time = kmfit$time, 
                     tment = c(rep(names(kmfit$strata)[1], kmfit$strata[1]), 
                               rep(names(kmfit$strata)[2], kmfit$strata[2])))

# Create Figure 4.2
fig42 <- ggplot(aes(x = time / 365.25, y = surv, linetype = tment), data = kmdata) + 
  geom_step(size = 1) + 
  scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA")) + 
  xlab("Time since randomization (years)") + 
  ylab("Survival probability") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6), breaks = seq(0, 6, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), limits = c(0,1)) +
  theme_general

fig42
fig42<-fig42+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"))
fig42
ggsave("figures/h_pbc3kmtreat.pdf", plot = fig42, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.4 --------------------------------------#
#------------------------------------------------------------------#

# Poisson model fit (like in Chapter 2)

# Cuts
cuts <- c(0, 2, 4) * 365.25

# event/failure indicator
pbc3$fail <- ifelse(pbc3$status != 0, 1, 0)

# Make the data ready using survSplit
pbc3mult <- survSplit(Surv(days, fail) ~ ., 
                      pbc3,
                      cut = cuts[-1], 
                      episode = "timegroup")

# Number of observations per interval (timegroup)
with(pbc3mult, table(tment, timegroup))

# Risk time
pbc3mult$risktime <- pbc3mult$days - cuts[pbc3mult$timegroup] 


# Summarize
require(dplyr)
sumdata <- pbc3mult %>% 
  group_by(tment, timegroup) %>% 
  summarise(fail = sum(fail), 
            risktime = sum(days - cuts[timegroup])#,
            #logrisktime = sum(log(days))
  )
sumdata <- as.data.frame(sumdata)
sumdata



# KM data from figure 4.2 model fit
tment1 <- subset(kmdata, tment == "tment=1")

# Estimated hazard per time group
sumdata$hazard_timegroup <- sumdata$fail / sumdata$risktime

# Add a numeric version of the treatment to the NA estimates
kmdata$tmentnum <- ifelse(kmdata$tment == "tment=0", 0, 1)

# Add piecewise constant hazard to data
kmdata$pwch <- NULL

# Between time 0 and 2
kmdata$pwch[kmdata$time <= 2 * 365.25] <- kmdata$time[kmdata$time <= 2  * 365.25] * 
  (sumdata$hazard_timegroup[1] * (1-kmdata$tmentnum[kmdata$time <= 2  * 365.25]) + 
     sumdata$hazard_timegroup[4] * (kmdata$tmentnum[kmdata$time <= 2 * 365.25]))

# Between time 2 and 4
kmdata$pwch[kmdata$time > 2  * 365.25 & kmdata$time <= 4 * 365.25] <- 2  * 365.25 * 
  (sumdata$hazard_timegroup[1] * (1-kmdata$tmentnum[kmdata$time > 2  * 365.25& kmdata$time <= 4 * 365.25]) + 
     sumdata$hazard_timegroup[4] * (kmdata$tmentnum[kmdata$time > 2  * 365.25& kmdata$time <= 4 * 365.25])) + 
  (kmdata$time[kmdata$time > 2 * 365.25 & kmdata$time <= 4 * 365.25] - 2 * 365.25) * 
  (sumdata$hazard_timegroup[2] * (1-kmdata$tmentnum[kmdata$time > 2  * 365.25& kmdata$time <= 4 * 365.25]) + 
     sumdata$hazard_timegroup[5] * (kmdata$tmentnum[kmdata$time > 2 * 365.25 & kmdata$time <= 4 * 365.25]))

# After time 4
kmdata$pwch[kmdata$time > 4 * 365.25] <- 2 * 365.25 * 
  (sumdata$hazard_timegroup[1] * (1-kmdata$tmentnum[kmdata$time > 4 * 365.25]) + 
     sumdata$hazard_timegroup[4] * (kmdata$tmentnum[kmdata$time > 4 * 365.25])) + 
  2 * 365.25 *
  (sumdata$hazard_timegroup[2] * (1-kmdata$tmentnum[kmdata$time > 4 * 365.25]) + 
     sumdata$hazard_timegroup[5] * (kmdata$tmentnum[kmdata$time > 4 * 365.25])) + 
  (kmdata$time[kmdata$time > 4 * 365.25] - 4 * 365.25) * 
  (sumdata$hazard_timegroup[3] * (1-kmdata$tmentnum[kmdata$time > 4 * 365.25]) + 
     sumdata$hazard_timegroup[6] * (kmdata$tmentnum[kmdata$time > 4 * 365.25]))


# Change to estimated survival (plug-in formula)
kmdata$pwcs <- exp(-kmdata$pwch)

# View the data 
head(kmdata)

# Reformat for plot
piecepdata <- data.frame(surv = c(kmdata$surv, kmdata$pwcs), 
                         time = rep(kmdata$time, 2),
                         tmentnum = rep(kmdata$tmentnum, 2),
                         type = c(rep("Kaplan-Meier", length(kmdata$time)), 
                                  rep("Piece-wise exponential", length(kmdata$time))))
# Only for treatment 1
piecepdata1 <- subset(piecepdata, tmentnum == 1)


# Create Figure 4.4
fig44 <- ggplot(aes(x = time / 365.25, y = surv, linetype = type), 
                data = subset(piecepdata1, type == "Kaplan-Meier")) + 
  geom_step(size = 1) + 
  geom_line(aes(x = time/ 365.25, y = surv, linetype = type), size = 1,
            data = subset(piecepdata1, type == "Piece-wise exponential")) + 
  labs(linetype = "Type") + 
  xlab("Time since randomization (years)") + 
  ylab("Survival probability") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general + 
  theme(legend.position="bottom",
        legend.box="vertical",
        text = element_text(size=21))
fig44
fig44<-fig44+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"))
fig44

ggsave("figures/h_pbc3survpwchplac.pdf", plot = fig44, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.5 --------------------------------------#
#------------------------------------------------------------------#


# Cox model fit with covariates tment, alb and log2bili
coxfit <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili, data = pbc3, 
                 method = "breslow")

# Unique days times 
fu <- sort(unique(pbc3$days))


# Data for prediction
preddata <- data.frame(tment = c(rep(0, length(fu)), rep(1, length(fu))), 
                       alb = rep(38, length(fu) * 2),
                       log2bili = rep(log2(45), length(fu) * 2), 
                       days = c(fu, fu),
                       status = rep(1, length(fu) * 2)
                       )

# Linear predictor
preds <- predict(coxfit, newdata = preddata, type = "survival")

preddata$preds <- preds

# Look at the predictions 
head(preddata)

# Create Figure 4.5
fig45 <- ggplot(aes(x = days / 365.25, y = preds, linetype = as.factor(tment)), 
                data = preddata) + 
  geom_step(size = 1) + 
  scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA")) + 
  xlab("Time since randomization (years)") + 
  ylab("Estimated survival function") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general 
fig45


fig45<-fig45+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"))
fig45

ggsave("figures/h_pbc3predsurvcox.pdf", plot = fig45, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 4.6 --------------------------------------#
#------------------------------------------------------------------#

# Add log(-log(S(t)))
preddata$logminlogsurv <- with(preddata, log(-log(preds)))

# Look at the predictions 
head(preddata)

# Create Figure 4.6
fig46 <- ggplot(aes(x = days / 365.25, y = logminlogsurv, linetype = as.factor(tment)), 
                data = preddata) + 
  geom_step(size = 1) + 
  scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA")) + 
  xlab("Time since randomization (years)") + 
  ylab("log(-log(survival function))") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05))) +
  theme_general

fig46

fig46<-fig46+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"))
fig46

ggsave("figures/h_pbc3cloglogsurv.pdf", plot = fig46, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.7 --------------------------------------#
#------------------------------------------------------------------#

# Fit Cox model with tment, alb and log2bili as covariates
summary(coxfit)

# We want to predict responses if all subjects had been given both CyA and placebo (tment = 0,1)
# We need to have n survival curves * 2
# make a double data set with extra Z
pbc3_counterfact <- pbc3
pbc3_counterfact$tment <- ifelse(pbc3_counterfact$tment == 1, 0, 1) # opposite treatment
pbc3_double <- rbind(pbc3, pbc3_counterfact)

# Check - one of each treat!
with(pbc3_double, table(tment, id))

# Baseline survival
pred <- survfit(coxfit, newdata = data.frame(tment = 0, alb = 0, log2bili = 0))

allsurv <-
  lapply(1:nrow(pbc3_double),
       function(i)
         pred$surv ^
         exp(coef(coxfit)[1] * pbc3_double$tment[i] +
               coef(coxfit)[2] * pbc3_double$alb[i] +
               coef(coxfit)[3] * pbc3_double$log2bili[i]))

potout <-
  data.frame(surv = unlist(allsurv),
           tment = rep(pbc3_double$tment, each = length(pred$time)),
           time = rep(pred$time, times = nrow(pbc3)*2)
)

# Average over values
require(dplyr)
sumdata <- potout %>%
  group_by(tment, time) %>%
  summarise(average_pred = mean(surv, na.rm = TRUE),
            .groups = c("keep"))
sumdata <- as.data.frame(sumdata)

# Create Figure 4.7
fig47 <- ggplot(aes(x = time / 365.25, y = average_pred, linetype = as.factor(tment)),
                data = sumdata) +
  geom_step(size = 1) +
  scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA")) +
  xlab("Time since randomization (years)") +
  ylab("Estimated survival function (g-formula)") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general 

fig47

fig47<-fig47+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"))
fig47

ggsave("figures/h_pbc3gformulabreslow.pdf", plot = fig47, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.8 --------------------------------------#
#------------------------------------------------------------------#

# Survival prediction
preds <- basehaz(coxfit, centered = F)

# For tment=0, alb=38, log2bili=log2(45)
dA_tment0 <- diff(c(0, preds$hazard * exp(coef(coxfit)[1]*0 + coef(coxfit)[2]*38 + coef(coxfit)[3] * log2(45) )))
surv_tment0 <- cumprod(1 - dA_tment0)

# For tment=1, alb=38, log2bili=log2(45)
dA_tment1 <- diff(c(0, preds$hazard * exp(coef(coxfit)[1]*1 + coef(coxfit)[2]*38 + coef(coxfit)[3] * log2(45) )))
surv_tment1 <- cumprod(1 - dA_tment1)


pdata <- data.frame(surv = c(surv_tment0, surv_tment1),
                    time = c(preds$time, preds$time), 
                    tment = c(rep("0", length(preds$time)),
                              rep("1", length(preds$time))))


# Create Figure 4.8
fig48 <- ggplot(aes(x = time / 365.25, y = surv, linetype = as.factor(tment)), 
                data = pdata) + 
  geom_step(size = 1) + 
  scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA")) + 
  xlab("Time since randomization (years)") + 
  ylab("Estimated survival function") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general

fig48

fig48<-fig48+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"))
fig48

ggsave("figures/h_pbc3predsurvcoxpl.pdf", plot = fig48, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.10 -------------------------------------#
#------------------------------------------------------------------#

# Cumulative incidences are estimated using Aalen-Johansen

# Overall survival 
overall_surv <- survfit(Surv(days, status != 0) ~ 1, 
                        data = subset(pbc3, tment_char == "Placebo"))

# Cause 1 survival: transplantation
cause1_cif <- cif(Event(days, status) ~ 1, 
                  data = subset(pbc3, tment_char == "Placebo"), 
                  cause = 1)

# Cause 2 survival: death w/o transplantation
cause2_cif <- cif(Event(days, status) ~ 1, 
                  data = subset(pbc3, tment_char == "Placebo"), 
                  cause = 2)

# Get them on the same time scale - book keeping

alltimes <- overall_surv$time
cause1st <- stepfun(x = cause1_cif$cumhaz[,1], y = c(0, cause1_cif$cumhaz[,2]))
cause1times <- cause1st(v = alltimes)

cause2st <- stepfun(x = cause2_cif$cumhaz[,1], y = c(0, cause2_cif$cumhaz[,2]))
cause2times <- cause2st(v = alltimes)


# Collect the data 
data_comb <- data.frame(cif = c(overall_surv$surv + cause1times + cause2times, 
                                cause1times + cause2times, 
                                cause1times),
                        time = c(alltimes, alltimes, alltimes),
                        type = c(rep("Overall", length(alltimes)), 
                                 rep("Transplantation + death without transplantation", length(alltimes)), 
                                 rep("Transplantation", length(alltimes)))
                        )

subset(data_comb, type == "Transplantation")

# Create Figure 4.10
fig410 <- ggplot(aes(x = time / 365.25, y = cif, linetype = type), 
                data = data_comb) + 
  geom_step(size = 1) + 
  scale_linetype_discrete("Type") + 
  xlab("Time since randomization (years)") + 
  ylab("Stacked cumulative incidence and survival") + 
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)), 
                     limits = c(0, 6), 
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)), 
                     limits = c(0, 1.0), 
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general + 
  guides(linetype = guide_legend(nrow = 3, byrow = TRUE)) + 
  theme(legend.position="bottom",
        legend.box="vertical",
        text = element_text(size=21), 
        legend.key.width = unit(1, "cm"))

fig410
fig410<-fig410+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig410

ggsave("figures/h_pbc3stackplacok.pdf", plot = fig410, 
       width = 29.7, height = 21, units = "cm")




#------------------------------------------------------------------#
#---------------- Figure 4.11 -------------------------------------#
#------------------------------------------------------------------#

# Cumulative incidences are (wrongly) estimated using Kaplan-Meier

# Overall survival 
overall_surv <- survfit(Surv(days, status != 0) ~ 1, 
                        data = subset(pbc3, tment_char == "Placebo"))

# Cause 1 survival: transplantation, KM
cause1_cif_w <- survfit(Surv(days, status == 1) ~ 1, 
                        data = subset(pbc3, tment_char == "Placebo"))

# Cause 2 survival: death w/o transplantation, KM
cause2_cif_w <- survfit(Surv(days, status == 2) ~ 1, 
                      data = subset(pbc3, tment_char == "Placebo"))

# Get them on the same time scale - book keeping
alltimes <- overall_surv$time
cause1st <- stepfun(x = cause1_cif_w$time, y = c(0, cause1_cif_w$surv))
cause1times <- cause1st(v = alltimes)

cause2st <- stepfun(x = cause2_cif_w$time, y = c(0, cause2_cif_w$surv))
cause2times <- cause2st(v = alltimes)

# Collect the data 
data_comb <- data.frame(cif_w = c(overall_surv$surv + 1-cause1times + 1-cause2times, 
                                1-cause1times + 1-cause2times, 
                                1-cause1times),
                        time = c(alltimes, alltimes, alltimes),
                        type = c(rep("Overall", length(alltimes)), 
                                 rep("Transplantation + death without transplantation", length(alltimes)), 
                                 rep("Transplantation", length(alltimes)))
)



# Create Figure 4.11
fig411 <- ggplot(aes(x = time / 365.25, y = cif_w, linetype = type), 
                 data = data_comb) + 
  geom_step(size = 1) + 
  scale_linetype_discrete("Type") + 
  xlab("Time since randomization (years)") + 
  ylab("Stacked cumulative incidence and survival") + 
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)), 
                     limits = c(0, 6), 
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)), 
                     limits = c(0, 1.1), 
                     breaks = seq(0, 1.1, 0.1)) +
  theme_general + 
  guides(linetype = guide_legend(nrow = 3, byrow = TRUE)) + 
  theme(legend.position="bottom",
        legend.box="vertical",
        text = element_text(size=21), 
        legend.key.width = unit(1, "cm"))
fig411
fig411<-fig411+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig411
ggsave("figures/h_pbc3stackplacwrong.pdf", plot = fig411, 
       width = 29.7, height = 21, units = "cm")




#------------------------------------------------------------------#
#---------------- Figure 4.12 -------------------------------------#
#------------------------------------------------------------------#

# Overall
overall_cox <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili + sex + age,
                        data = subset(pbc3, !is.na(alb)),
                     method = "breslow", eps=1e-9)

# Cause 1: transplantation
cause1_cox <- coxph(Surv(days, status == 1) ~ tment + alb + log2bili + sex + age,
                    data = subset(pbc3, !is.na(alb)), 
                    method = "breslow", eps=1e-9)

# Cause 2: death w/o transplantation
cause2_cox <- coxph(Surv(days, status == 2) ~ tment + alb + log2bili + sex + age,
                    data = subset(pbc3, !is.na(alb)), 
                    method = "breslow", eps=1e-9)


# for tment = placebo, age = 40, alb = 38, bili = log2(45) and sex = F
newd <- data.frame(tment = 0, age = 40, 
                   alb = 38, log2bili = log2(45), 
                   sex = 0)
est_cause1 <- survfit(cause1_cox, ctype = 2, newdata = newd)$cumhaz 

est_cause2 <- survfit(cause2_cox, ctype = 2, newdata = newd)$cumhaz 

est_overall <- survfit(overall_cox, ctype = 2, newdata = newd)$surv

# Calculate S, and F1 and F2
alltimes <- basehaz(overall_cox, centered = F)$time

A <- cumsum(diff(c(0, est_cause1 + est_cause2)))
dA <- diff(c(0, A))
S <- cumprod(1-dA) #exp(-cumsum(dA)) 
A1 <- est_cause1
A2 <- est_cause2

dat <- data.frame(cbind(alltimes, A, dA, S, A1, A2))
dat2 <- subset(dat, dA > 0)
dat2$lagS <- with(dat2, c(0, S[-length(S)])) 

dat2$F1 <- with(dat2, cumsum(lagS * diff(c(0, A1))))
dat2$F2 <- with(dat2, cumsum(lagS * diff(c(0, A2))))


# Collect the data
data_comb <- with(dat2, 
                  data.frame(cif_w = c(rep(1, length(alltimes)),
                                       F1 + F2,
                                       F1),
                             time = c(alltimes, alltimes, alltimes),
                             type = c(rep("Overall", length(alltimes)),
                                      rep("Transplantation + death without transplantation", length(alltimes)),
                                      rep("Transplantation", length(alltimes)))
                  ))


# Create Figure 4.12
fig412 <- ggplot(aes(x = time / 365.25, y = cif_w, linetype = type),
                 data = data_comb) +
  geom_step(size = 1) +
  scale_linetype_discrete("Type") +
  xlab("Time since randomization (years)") +
  ylab("Stacked cumulative incidence and survival") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  guides(linetype = guide_legend(nrow = 3, byrow = TRUE)) + 
  theme(legend.position="bottom",
        legend.box="vertical",
        text = element_text(size=22), 
        legend.key.width = unit(1, "cm"))


fig412

fig412<-fig412+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig412


ggsave("figures/h_pbc3predcumincPlacex1.pdf", plot = fig412,
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 4.13 -------------------------------------#
#------------------------------------------------------------------#

# for tment = placebo, age = 40, alb = 20, bili = log2(90) and sex = F
newd <- data.frame(tment = 0, age = 40, 
                    alb = 20, log2bili = log2(90), 
                    sex = 0)

est_cause1 <- survfit(cause1_cox, ctype = 2, newdata = newd)$cumhaz

est_cause2 <- survfit(cause2_cox, ctype = 2, newdata = newd)$cumhaz

est_overall <- survfit(overall_cox, ctype = 2, newdata = newd)$surv


# Calculate S, and F1 and F2
alltimes <- basehaz(overall_cox, centered = F)$time

A <- cumsum(diff(c(0, est_cause1 + est_cause2)))
dA <- diff(c(0, A))
S <- cumprod(1-dA) #exp(-cumsum(dA)) 
A1 <- est_cause1
A2 <- est_cause2

dat <- data.frame(cbind(alltimes, A, dA, S, A1, A2))
dat2 <- subset(dat, dA > 0)
dat2$lagS <- with(dat2, c(0, S[-length(S)])) 

dat2$F1 <- with(dat2, cumsum(lagS * diff(c(0, A1))))
dat2$F2 <- with(dat2, cumsum(lagS * diff(c(0, A2))))


# Collect the data
data_comb <- with(dat2, 
                  data.frame(cif_w = c(rep(1, length(alltimes)),
                                  F1 + F2,
                                  F1),
                        time = c(alltimes, alltimes, alltimes),
                        type = c(rep("Overall", length(alltimes)),
                                 rep("Transplantation + death without transplantation", length(alltimes)),
                                 rep("Transplantation", length(alltimes)))
))


# Create Figure 4.13
fig413 <- ggplot(aes(x = time / 365.25, y = cif_w, linetype = type),
                 data = data_comb) +
  geom_step(size = 1) +
  scale_linetype_discrete("Type") +
  xlab("Time since randomization (years)") +
  ylab("Stacked cumulative incidence and survival") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  guides(linetype = guide_legend(nrow = 3, byrow = TRUE)) + 
  theme(legend.position="bottom",
        legend.box="vertical",
        text = element_text(size=22), 
        legend.key.width = unit(1, "cm"))


fig413

fig413<-fig413+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig413

ggsave("figures/h_pbc3predcumincPlacex2.pdf", plot = fig413,
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 4.14 -------------------------------------#
#------------------------------------------------------------------#

# for tment = placebo, age = 60, alb = 38, bili = log2(45) and sex = F
newd <- data.frame(tment = 0, age = 60, 
                   alb = 38, log2bili = log2(45), 
                   sex = 0)

est_cause1 <- survfit(cause1_cox, ctype = 2, newdata = newd)$cumhaz

est_cause2 <- survfit(cause2_cox, ctype = 2, newdata = newd)$cumhaz

est_overall <- survfit(overall_cox, ctype = 2, newdata = newd)$surv


# Calculate S, and F1 and F2
alltimes <- basehaz(overall_cox, centered = F)$time

A <- cumsum(diff(c(0, est_cause1 + est_cause2)))
dA <- diff(c(0, A))
S <- cumprod(1-dA) #exp(-cumsum(dA)) 
A1 <- est_cause1
A2 <- est_cause2

dat <- data.frame(cbind(alltimes, A, dA, S, A1, A2))
dat2 <- subset(dat, dA > 0)
dat2$lagS <- with(dat2, c(0, S[-length(S)])) 

dat2$F1 <- with(dat2, cumsum(lagS * diff(c(0, A1))))
dat2$F2 <- with(dat2, cumsum(lagS * diff(c(0, A2))))


# Collect the data
data_comb <- with(dat2, 
                  data.frame(cif_w = c(rep(1, length(alltimes)),
                                       F1 + F2,
                                       F1),
                             time = c(alltimes, alltimes, alltimes),
                             type = c(rep("Overall", length(alltimes)),
                                      rep("Transplantation + death without transplantation", length(alltimes)),
                                      rep("Transplantation", length(alltimes)))
                  ))


# Create Figure 4.14
fig414 <- ggplot(aes(x = time / 365.25, y = cif_w, linetype = type),
                 data = data_comb) +
  geom_step(size = 1) +
  scale_linetype_discrete("Type") +
  xlab("Time since randomization (years)") +
  ylab("Stacked cumulative incidence and survival") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general + 
  guides(linetype = guide_legend(nrow = 3, byrow = TRUE)) + 
  theme(legend.position="bottom",
        legend.box="vertical",
        text = element_text(size=22), 
        legend.key.width = unit(1, "cm"))


fig414

fig414<-fig414+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig414

ggsave("figures/h_pbc3predcumincPlacex3.pdf", plot = fig414,
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Table 4.5 ---------------------------------------#
#------------------------------------------------------------------#

library("survRM2")

pbcny <- subset(pbc3,!is.na(alb))
time <- pbcny$days
status <- pbcny$status!=0
arm <- pbcny$tment
alb <- pbcny$alb
logbili <- log2(pbcny$bili)
x <- cbind(alb, logbili)

rmst2(time, status, arm, tau=3, covariates=x)

#------------------------------------------------------------------#
#---------------- Table 4.6 ---------------------------------------#
#------------------------------------------------------------------#


# Fine-Gray model for transplant
# Two step
fg_c1 <- finegray(Surv(days, as.factor(status)) ~ tment + alb + log2bili + sex + age,
                  etype = 1,
                  data = subset(pbc3, !is.na(alb)))

fg_cox1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ tment + alb + log2bili + sex + age, 
                weight = fgwt, data = fg_c1, ties = "breslow", eps=1e-9)
summary(fg_cox1)


# Fine-Gray model for death w/o transplant
# Two step
fg_c2 <- finegray(Surv(days, as.factor(status)) ~ tment + alb + log2bili + sex + age,
                  etype = 2,
                  data = subset(pbc3, !is.na(alb)))

fg_cox2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ tment + alb + log2bili + sex + age, 
                weight = fgwt, data = fg_c2, ties = "breslow", eps=1e-9)
summary(fg_cox2)

#------------------------------------------------------------------#
#---------------- Figure 4.19 -------------------------------------#
#------------------------------------------------------------------#


# Fine-Gray model for death w/o transplant
# Predictions 
# Placebo
newd_t0 <- data.frame(tment = 0, age = 40, 
                      alb = 38, log2bili = log2(45), 
                      sex = 0)
C2_t0 <- survfit(fg_cox2, ctype = 2, newdata = newd_t0)$cumhaz 

# CyA
newd_t1 <- data.frame(tment = 1, age = 40, 
                      alb = 38, log2bili = log2(45), 
                      sex = 0)
C2_t1 <- survfit(fg_cox2, ctype = 2, newdata = newd_t1)$cumhaz 

time <- survfit(fg_cox2, ctype = 2, newdata = newd_t1)$time

# Make data ready for plotting 
pdata <- data.frame(time = c(time, time), 
                    cif = c(C2_t0, C2_t1), 
                    tment = c(rep("Placebo", length(time)), 
                              rep("CyA", length(time))))


# Create Figure 4.19
fig419 <- ggplot(aes(x = time / 365.25, y = cif, linetype = tment), 
                 data = pdata) + 
  geom_step(size = 1) + 
  scale_linetype_manual("Treatment", values = c("dashed", "solid"),guide = guide_legend(reverse = TRUE)) + 
  xlab("Time since randomization (years)") + 
  ylab('Cumulative incidence for death w/o transplantation') + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6), 
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), 
                     limits = c(0, 1.0), 
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general + 
  theme(legend.position="bottom", legend.box = "vertical")
      #  legend.key.size = unit(1.5, 'cm'))

fig419
fig419<-fig419+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(3,"line"))
fig419


ggsave("figures/h_pbc3deathFG.pdf", plot = fig419, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.20 -------------------------------------#
#------------------------------------------------------------------#


# Fine-Gray model for transplant
# Predictions 
# Placebo
C1_t0 <- survfit(fg_cox1, ctype = 2, newdata = newd_t0)$cumhaz 

# CyA
C1_t1 <- survfit(fg_cox1, ctype = 2, newdata = newd_t1)$cumhaz 

time <- survfit(fg_cox1, ctype = 2, newdata = newd_t1)$time

# Make data ready for plotting 
pdata <- data.frame(time = c(time, time), 
                    cif = c(C1_t0, C1_t1), 
                    tment = c(rep("Placebo", length(time)), 
                              rep("CyA", length(time))))


# Create Figure 4.20
fig420 <- ggplot(aes(x = time / 365.25, y = cif, linetype = tment), 
                 data = pdata) + 
  geom_step(size = 1) + 
  scale_linetype_manual("Treatment", values = c("dashed", "solid"),guide = guide_legend(reverse = TRUE)) + 
  xlab("Time since randomization (years)") + 
  ylab('Cumulative incidence for transplantation') + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6), 
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), 
                     limits = c(0, 1.0), 
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general + 
  theme(legend.box = "vertical",
        legend.key.size = unit(1.5, 'cm'))


fig420

fig420<-fig420+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig420

ggsave("figures/h_pbc3transplFG.pdf", plot = fig420, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Table 4.7 ---------------------------------------#
#------------------------------------------------------------------#


library(survival)

rmtl.ipcw <- function(times, event, eoi=1, tau, cov=NULL, strata=FALSE, group=NULL){
  
  if(is.null(group) & strata==TRUE){stop('Please specify a factor variable to statify weights.')}
  if(is.null(cov)){print('Warning: Fitting intercept-only model.')}
  
  # Round event times to avoid issues with survival() package rounding differently
  y <- round(times,4)
  id <- 1:length(y)
  
  # Recode so delta1 reflects event of interest, delta2 reflects all competing events. Assumes 0=censoring.
  delta1 <- ifelse(event==eoi, 1, 0)
  delta2 <- ifelse(event!=0 & event!=eoi, 1, 0)
  
  # Overall quantities
  x <- cbind(int=rep(1, length(y)), cov)
  p <- length(x[1,])
  if(is.null(group)){group <- as.factor(rep(1, length(y)))}
  
  # Recode event indicators to reflect status at chosen tau
  delta1[y>tau] <- 0
  delta2[y>tau] <- 0
  
  y <- pmin(y, tau)
  y1 <- y*delta1
  
  d0 <- 1 - (delta1 + delta2) # censoring indicator
  d0[y==tau] <- 0  # If follow-up lasts til tau, the event will not count as 'censored' in IPCW weights
  weights <- NULL
  
  ## Calculate IPCW weights (option to stratify by group) ## 
  
  if(strata==TRUE){
    for(aa in 1:length(unique(group))){
      # Subset the group
      a <- unique(group)[aa]
      d0.a <- d0[group==a]
      delta1.a <- delta1[group==a]
      y.a <- y[group==a]
      x.a <- x[group==a,]
      n.a <- length(d0.a)
      orig.id.a0 <- orig.id.a <- id[group==a]
      
      # Order the event times
      id.a <- order(y.a)
      y.a <- y.a[id.a]
      d0.a <- d0.a[id.a]
      delta1.a <- delta1.a[id.a]
      x.a <- x.a[id.a,]
      orig.id.a <- orig.id.a[id.a]
      
      # Derive IPCW
      fit <- survfit(Surv(y.a, d0.a) ~ 1)
      weights.a <- (1-d0.a)/rep(fit$surv, table(y.a))
      
      # Need to assign weights accordig to original ID, not ordered by event time
      linked.weights.a <- cbind(orig.id.a, weights.a, delta1.a, d0.a, y.a)
      weights <- rbind(weights, linked.weights.a)
    }
  } else {
    
    # Order the event times
    id.a <- order(y)
    y.a <- y[id.a]
    d0.a <- d0[id.a]
    delta1.a <- delta1[id.a]
    x.a <- x[id.a,]
    orig.id.a <- id[id.a]
    
    # Derive IPCW
    fit <- survfit(Surv(y.a, d0.a) ~ 1)
    weights.a <- (1-d0.a)/rep(fit$surv, table(y.a))
    
    # Need to assign weights accordig to original ID, not ordered by event time
    linked.weights.a <- cbind(orig.id.a, weights.a, delta1.a, d0.a, y.a)
    weights <- rbind(weights, linked.weights.a)
  }
  
  
  ## Fit linear model ## 
  
  # Link weights to original data frame
  #colnames(weights) <- c('id', 'weights')
  #data <- merge(data0, weights, by='id')
  #summary(lm(tau-y ~ x-1, weights=weights, data=data))
  
  # Or, sort weights and use vectors
  w <- weights[order(weights[, 1]),2]
  lm.fit <- lm(delta1*(tau-y) ~ x-1, weights=w)
  
  
  ## Derive SE ##
  
  beta0 <- lm.fit$coef
  error <- tau - y - as.vector(x %*% beta0)
  score <- x * w * error
  
  # Kappa (sandwich variance components) stratified by group
  kappa <- NULL
  
  for(aa in 1:length(unique(group))){
    
    # Subset the group
    a <- unique(group)[aa]
    d0.a <- d0[group==a]
    delta1.a <- delta1[group==a]
    y.a <- y[group==a]
    x.a <- x[group==a,]
    n.a <- length(d0.a)
    orig.id.a0 <- orig.id.a <- id[group==a]
    score.a <- score[group==a,]
    
    # Kappa calculations for sandwich variance
    kappa.a <- matrix(0, n.a, p)
    
    for(i in 1:n.a){
      kappa1 <- score.a[i,]
      
      kappa2 <- apply(score.a[y.a>=y.a[i],,drop=F], 2, sum)*(d0.a[i])/sum(y.a>=y.a[i])
      
      kappa3 <- rep(0, p)
      
      for(k in 1:n.a){
        if(y.a[k]<=y.a[i]){
          kappa3 <- kappa3+apply(score.a[y.a>=y.a[k],,drop=F], 2, sum)*(d0.a[k])/(sum(y.a>=y.a[k]))^2
        }
      }
      
      kappa.a[i,] <- kappa1+kappa2-kappa3
    }
    kappa <- rbind(kappa, kappa.a)
  }
  
  # Transpose the kappas rbinded from each group gives pxp matrix
  gamma <- t(kappa) %*% kappa
  
  A <- t(x) %*% x
  varbeta <- solve(A) %*% gamma %*% solve(A)
  se <- sqrt(diag(varbeta))
  
  
  #--- Return results ---
  
  res <- cbind(beta=lm.fit$coef, se=se, cil=lm.fit$coef-(1.96*se), ciu=lm.fit$coef+(1.96*se), 
               z=lm.fit$coef/se, p=2*(1-pnorm(abs(lm.fit$coef/se))))
  #rownames(res) <- c("Intercept", colnames(x[,-1])) FJERNET! PKA 030621
  
  allres <- list(res=res, varbeta=varbeta)
  print(round(res, 3))
  invisible(allres)
  return(res[,1]) #tilf?jet 040322
}  


rmtl.ipcw(pbc3$days,pbc3$status,eoi=1,tau=3,cov=pbc3$tment)

rmtl.ipcw(pbc3$days,pbc3$status,eoi=2,tau=3,cov=pbc3$tment)

pbcny <- subset(pbc3,!is.na(alb))
time <- pbcny$days
status <- pbcny$status
arm <- pbcny$tment
sex <- pbcny$sex
age <- pbcny$age
alb <- pbcny$alb
logbili <- log2(pbcny$bili)

x <- cbind(arm, sex, age, alb, logbili)

rmtl.ipcw(time, status, eoi=1, tau=3, x)
rmtl.ipcw(time, status, eoi=2, tau=3, x)

x1 <- cbind(arm, alb, logbili)

rmtl.ipcw(time, status, eoi=1, tau=3, x1)
rmtl.ipcw(time, status, eoi=2, tau=3, x1)



#------------------------------------------------------------------#
#---------------- Figure 4.23 -------------------------------------#
#------------------------------------------------------------------#

# Overall survival 
overall_surv <- survfit(Surv(days, status != 0) ~ 1, data = pbc3)

# Censoring
cens_surv <- survfit(Surv(days, status == 0) ~ 1, data = pbc3)


# Make data ready for plotting 
pdata <- data.frame(time = c(overall_surv$time, cens_surv$time), 
                    surv = c(overall_surv$surv, cens_surv$surv), 
                    type = c(rep("Treatment failure", length(overall_surv$time)), 
                             rep("Censoring", length(overall_surv$time))))


# Create Figure 4.23
fig423 <- ggplot(aes(x = time / 365.25, y = surv, linetype = type), 
                 data = pdata) + 
  geom_step(size = 1) + 
  scale_linetype_manual("Type", values = c("dashed", "solid")) + 
  xlab("Time since randomization (years)") + 
  ylab('Probability') + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 6), 
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), 
                     limits = c(0, 1.0), 
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general + 
  theme(legend.box = "vertical",
        legend.key.size = unit(1.5, 'cm'))


fig423

fig423<-fig423+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig423


ggsave("figures/h_pbc3censsurv.pdf", plot = fig423, 
       width = 29.7, height = 21, units = "cm")

