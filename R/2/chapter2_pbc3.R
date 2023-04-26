#------------------------------------------------------------------#
#------- Chapter 2, R code, PBC3 data  ----------------------------#
#------------------------------------------------------------------#

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


# Note that status represents status at exit 
# (0: censored, 1: liver transpl, 2 : dead)
# Per example 1.1.1, the main event of interest is "failure of medical treatment" 
# this is either death or liver transplant (status %in% c(1,2))


# tment = 0 (placebo), tment = 1 (cyA)
with(pbc3, table(status, tment))

pbc3$tment_char <- ifelse(pbc3$tment == 0, "Placebo", "CyA")

#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------# 
#------------------------------------------------------------------#

# # Empty data example
# pdata <- data.frame()
# 
# playout <- ggplot(pdata) + 
#            scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
#            scale_y_continuous(expand = expansion(mult = c(0, 0))) 
# playout

# Extra wishes from PKA & HNRV to be added here
theme_general <- theme_bw() +
                 theme(legend.position = "bottom", 
                       text = element_text(size = 20), 
                       axis.text.x = element_text(size = 20), 
                       axis.text.y = element_text(size = 20)) 
                 
theme_general

# Question: Legend position?




#------------------------------------------------------------------#
#---------------- Figure 2.1 --------------------------------------#
#------------------------------------------------------------------#

# Fit a Nelson Aalen estimate of the cumulative hazard using the pbc3 data
# Stratifying on treatment (i.e. doing it per treatment)
nafit <- survfit(Surv(days, status != 0) ~ tment, data = pbc3)

# Summary of fit
summary(nafit)

# Collect data for plot
# Note that the standard errors produced by survfit are for the cumulative hazard
nadata <- data.frame(cumhaz = nafit$cumhaz, 
                     cumhaz_se = nafit$std.err, 
                     time = nafit$time, 
                     tment = c(rep(names(nafit$strata)[1], nafit$strata[1]), 
                               rep(names(nafit$strata)[2], nafit$strata[2])))

# Create Figure 2.1
fig21 <- ggplot(aes(x = time / 365.25, y = cumhaz, linetype = tment), data = nadata) + 
          geom_step(size = 1) + 
          scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA")) + 
          xlab("Time since randomization (years)") + 
          ylab("Cumulative hazard") + 
          scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                             limits = c(0, 6),
                             breaks = seq(0, 6, 1))  + 
          scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
          theme_general + 
      theme(legend.position="bottom",
            legend.box="vertical")

fig21

ggsave("figures/j_pbc3NAatreat.pdf", plot = fig21, 
       width = 29.7, height = 21, units = "cm")

# Estimates after two years
# Treatment 0, observation after two years
tail(subset(nadata, nadata$time <= 2 & nadata$tment == "tment=0"),1)

# Treatment 1, observation after two years
tail(subset(nadata, nadata$time <= 2 & nadata$tment == "tment=1"),1)


# Log rank test for treatment 
survdiff(Surv(days, status != 0) ~ tment, data = pbc3)


#------------------------------------------------------------------#
#---------------- Table 2.1 ---------------------------------------#
#------------------------------------------------------------------#

# Descriptive statistics
with(pbc3, table(tment, status))

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


# Poisson model for treatment = 1
require(stats)
poismod_trt1 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) - 1, 
                    data = subset(pbc3mult, tment == 1), 
                    family = poisson)

summary(poismod_trt1)

# Transform back to original scale (delta method)
exp(poismod_trt1$coefficients)

# Standard errors
sqrt(diag(vcov(poismod_trt1))) * exp(poismod_trt1$coefficients)


# Poisson model for treatment = 0
poismod_trt0 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) - 1, 
                    data = subset(pbc3mult, tment == 0), 
                    family = poisson)

summary(poismod_trt0)

# Transform back to original scale (delta method)
exp(poismod_trt0$coefficients)

# Standard errors
sqrt(diag(vcov(poismod_trt0))) * exp(poismod_trt0$coefficients)



#------------------------------------------------------------------#
#---------------- Figure 2.3 --------------------------------------#
#------------------------------------------------------------------#

# Rates + cuts from Table 2.1
rateCyA <- c(8.1,13.1,9.6)
ratePbo <- c(9.4,12.5,8.5)
pcwtime <- c(0,2,4,5)

# Collect data
plotdata <- data.frame(rates = c(rateCyA, ratePbo),
                       tment = c(rep("CyA", length(rateCyA)), rep("Placebo", length(ratePbo))),
                       times_s = rep(pcwtime[-4], 2),
                       times = rep(pcwtime[-1], 2))


# Create Figure 2.3
fig23 <- ggplot(aes(x = time, y = rates, linetype = tment),
                data = plotdata) +
  geom_segment(aes(x = times_s, y = rates, xend = times, yend = rates), size = 1) +
  scale_linetype_discrete("Treatment", labels = c("Placebo", "CyA"))  +
  xlab("Time since randomization (years)") +
  ylab("Estimated hazard function (per 100 years)") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 5),
                     breaks = seq(0, 5, 1))  +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0,14), breaks = seq(0, 14, 2)) +
  theme_general

fig23
 
ggsave("figures/j_pbc3pwchestimates.pdf", plot = fig23, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 2.4 --------------------------------------#
#------------------------------------------------------------------#

# NA data from figure 2.1 model fit
tment1 <- subset(nadata, tment == "tment=1")

# Estimated hazard per time group
sumdata$hazard_timegroup <- sumdata$fail / sumdata$risktime

# View data
head(sumdata)

# Add a numeric version of the treatment to the NA estimates
nadata$tmentnum <- ifelse(nadata$tment == "tment=1", 1, 0)

head(nadata)

# Add piecewise constant hazard to data
nadata$pwch <- NULL

# Between time 0 and 2
nadata$pwch[nadata$time <= 2 * 365.25] <- nadata$time[nadata$time <= 2 * 365.25] * 
                                 (sumdata$hazard_timegroup[1] * (1-nadata$tmentnum[nadata$time <= 2 * 365.25]) + 
                                  sumdata$hazard_timegroup[4] * (nadata$tmentnum[nadata$time <= 2 * 365.25]))
# Between time 2 and 4
nadata$pwch[nadata$time > 2 * 365.25 & nadata$time <= 4* 365.25 ] <- 2 * 365.25 * 
  (sumdata$hazard_timegroup[1] * (1-nadata$tmentnum[nadata$time > 2 * 365.25 & nadata$time <= 4 * 365.25]) + 
     sumdata$hazard_timegroup[4] * (nadata$tmentnum[nadata$time > 2 * 365.25 & nadata$time <= 4 * 365.25])) + 
  (nadata$time[nadata$time > 2 * 365.25 & nadata$time <= 4 * 365.25] - 2 * 365.25) * 
  (sumdata$hazard_timegroup[2] * (1-nadata$tmentnum[nadata$time > 2 * 365.25 & nadata$time <= 4 * 365.25]) + 
     sumdata$hazard_timegroup[5] * (nadata$tmentnum[nadata$time > 2 * 365.25 & nadata$time <= 4 * 365.25]))
  
# After time 4
nadata$pwch[nadata$time > 4 * 365.25] <- 2 * 365.25 * 
  (sumdata$hazard_timegroup[1] * (1-nadata$tmentnum[nadata$time > 4 * 365.25]) + 
     sumdata$hazard_timegroup[4] * (nadata$tmentnum[nadata$time > 4 * 365.25])) + 
   2 * 365.25 * 
  (sumdata$hazard_timegroup[2] * (1-nadata$tmentnum[nadata$time > 4 * 365.25]) + 
     sumdata$hazard_timegroup[5] * (nadata$tmentnum[nadata$time > 4 * 365.25])) + 
  (nadata$time[nadata$time > 4 * 365.25] - 4 * 365.25) * 
  (sumdata$hazard_timegroup[3] * (1-nadata$tmentnum[nadata$time > 4 * 365.25]) + 
     sumdata$hazard_timegroup[6] * (nadata$tmentnum[nadata$time > 4 * 365.25]))

# View the data 
head(nadata)


# Reformat for plot
piecepdata <- data.frame(cumhaz = c(nadata$cumhaz, nadata$pwch), 
                         time = rep(nadata$time, 2),
                         tmentnum = rep(nadata$tmentnum, 2),
                         type = c(rep("Nelson-Aalen", length(nadata$time)), 
                                  rep("Piece-wise exponential", length(nadata$time))))
# Only for treatment 1
piecepdata1 <- subset(piecepdata, tmentnum == 1)


# Create Figure 2.4
fig24 <- ggplot(aes(x = time / 365.25, y = cumhaz, linetype = type), 
                data = subset(piecepdata1, type == "Nelson-Aalen")) + 
  geom_step(size = 1) + 
  geom_line(aes(x = time / 365.25, y = cumhaz, linetype = type), size = 1,
            data = subset(piecepdata1, type == "Piece-wise exponential")) + 
  labs(linetype = "Type") + 
  xlab("Time since randomization (years)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1))  + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  theme_general

fig24

ggsave("figures/j_pbc3pwchplac.pdf", plot = fig24, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 2.5 --------------------------------------#
#------------------------------------------------------------------#

# Fit a Cox model using the pbc3 data with treatment as a covariate
coxfit <- coxph(Surv(days, status != 0) ~ tment, data = pbc3, method = "breslow")

# Summary of fit
summary(coxfit)

# Extract cumulative baseline hazard
coxcumhaz <- basehaz(coxfit, centered = FALSE)

# Collect data for plot
coxdata <- data.frame(cumhaz = coxcumhaz$hazard, 
                      time = coxcumhaz$time, 
                      tment = rep("0", nrow(coxcumhaz)), 
                      type = rep("Breslow estimate", nrow(coxcumhaz)))

# Create Figure 2.5
fig25 <- ggplot(aes(x = time / 365.25, y = cumhaz, linetype = tment), data = coxdata) + 
  geom_step(size = 1) + 
  xlab("Time since randomization (years)") + 
  ylab("Cumulative baseline hazard") + 
  scale_linetype_discrete("Treatment", labels = c("Placebo")) + 
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  theme_general

fig25

ggsave("figures/j_pbc3breslow.pdf", plot = fig25, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Table 2.3 ---------------------------------------#
#------------------------------------------------------------------#

# Calculate average covariate values of albumin and bilirubin per treatment group

# Summarize
require(dplyr)
sumdata <- pbc3 %>% 
  group_by(tment, tment_char) %>% 
  summarise(n = sum(id != 0),  
            average_albumin = mean(alb, na.rm = TRUE), # NOTE: Removes missing observations from mean computation
            average_bilirubin = mean(bili, na.rm = TRUE), 
            )
sumdata <- as.data.frame(sumdata)
sumdata


#------------------------------------------------------------------#
#---------------- Table 2.4 ---------------------------------------#
#------------------------------------------------------------------#

# Cox model with treatment, albumin and bilirubin as covariates 

coxfit <- coxph(Surv(days, status != 0) ~ tment + alb + bili,
                data = pbc3)

summary(coxfit)

#------------------------------------------------------------------#
#---------------- Table 2.5 ---------------------------------------#
#------------------------------------------------------------------#

# Poisson model with treatment, albumin and bilirubin as covariates

pbc3mult$tment_char <- as.factor(pbc3mult$tment_char)
pbc3mult$tment_char <- relevel(pbc3mult$tment_char, ref = "Placebo")

poismod_t25 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + bili,
                   data = pbc3mult,
                   family = poisson)

summary(poismod_t25)

#------------------------------------------------------------------#
#---------------- Table 2.6 ---------------------------------------#
#------------------------------------------------------------------#

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



# Make the data ready using survSplit
pbc3mult <- survSplit(Surv(days, fail) ~ ., 
                      pbc3,
                      cut = cuts[-1], 
                      episode = "timegroup")

# Risk time
pbc3mult$risktime <- pbc3mult$days - cuts[pbc3mult$timegroup] 


## Poisson models ##

# For comparison, LR test
base_poismod <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + bili,
                data = pbc3mult,
                family = poisson)

base_poismod_log <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili,
                    data = pbc3mult,
                    family = poisson)

## Linear effects 
# Poisson model 1
poismod_t26_l1 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + bili + albnorm,
                   data = pbc3mult,
                   family = poisson)

summary(poismod_t26_l1)

# LR test
require(lmtest)
lrtest(poismod_t26_l1, 
       base_poismod)
# Small error in table, test statistic = 0.44 (not 0.43)


# Poisson model 2
poismod_t26_l2 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + 
                       bili + bilihigh + bilitoohigh + bilimuchtoohigh,
                     data = pbc3mult,
                     family = poisson)

summary(poismod_t26_l2)

# LR test
lrtest(poismod_t26_l2, 
       base_poismod)

# Poisson model 3
poismod_t26_l3 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + 
                       log2bili + logbilihigh + logbilitoohigh + logbilimuchtoohigh,
                     data = pbc3mult,
                     family = poisson)

summary(poismod_t26_l3)

# LR test
lrtest(poismod_t26_l3, 
       base_poismod_log)


## Quadratic effects

# For comparison, LR test
base_poismod2 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb10 + bili100,
                    data = pbc3mult,
                    family = poisson)

base_poismod2_log <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb10 + log2bili,
                        data = pbc3mult,
                        family = poisson)



# Poisson model 1
poismod_t26_q1 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb10 + alb2 + bili100,
                     data = pbc3mult,
                     family = poisson)

summary(poismod_t26_q1)

# LR test
lrtest(poismod_t26_q1, 
       base_poismod2)

# Poisson model 2
poismod_t26_q2 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb10 + bili2 + bili100,
                      data = pbc3mult,
                      family = poisson)

summary(poismod_t26_q2)

# LR test
lrtest(poismod_t26_q2, 
       base_poismod2)


# Poisson model 3
poismod_t26_q3 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb10 + log2bili + log2bili2,
                      data = pbc3mult,
                      family = poisson)

summary(poismod_t26_q3)

# LR test
lrtest(poismod_t26_q3, 
       base_poismod2_log)


## Cox models ##

# Models for LR tests

base_coxmod <- coxph(Surv(days, status != 0) ~ tment + alb + bili, eps = 1e-8, method = "breslow",
                       data = pbc3)

base_coxmod_log <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili, eps = 1e-8, method = "breslow",
                     data = pbc3)


## Linear effects 
# Cox model 1
coxmod_t26_l1 <- coxph(Surv(days, status != 0) ~ tment + alb + bili + albnorm, eps = 1e-8, method = "breslow",
                 data = pbc3)

summary(coxmod_t26_l1)

# LR test
lrtest(coxmod_t26_l1, base_coxmod)

# or alternatively anova(coxmod_t26_l1)

# Cox model 2
coxmod_t26_l2 <- coxph(Surv(days, status != 0) ~ tment + alb + bili + bilihigh + bilitoohigh + bilimuchtoohigh,
                       eps = 1e-8, method = "breslow",
                       data = pbc3)

summary(coxmod_t26_l2)

# LR test
lrtest(coxmod_t26_l2, base_coxmod)

# Cox model 3
coxmod_t26_l3 <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili + logbilihigh + logbilitoohigh + logbilimuchtoohigh,
                       eps = 1e-8, method = "breslow",
                       data = pbc3)

summary(coxmod_t26_l3)

# LR test
lrtest(coxmod_t26_l3, base_coxmod_log)


## Quadratic effects 
# Cox model 1
coxmod_t26_q1 <- coxph(Surv(days, status != 0) ~ tment + alb10 + bili100 + alb2,
                       data = pbc3)

summary(coxmod_t26_q1)

# Cox model 2
coxmod_t26_q2 <- coxph(Surv(days, status != 0) ~ tment + alb10 + bili100 + bili2,
                       data = pbc3)

summary(coxmod_t26_q2)

# Cox model 3
coxmod_t26_q3 <- coxph(Surv(days, status != 0) ~ tment + alb10 + log2bili + log2bili2,
                       data = pbc3)

summary(coxmod_t26_q3)



#------------------------------------------------------------------#
#---------------- Figure 2.6 --------------------------------------#
#------------------------------------------------------------------#

# The below linear predictors include estimates from the following models, 


pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "3")
pbc3mult$tment_char <- as.factor(pbc3mult$tment_char)
pbc3mult$tment_char <- relevel(pbc3mult$tment_char, ref = "Placebo")

base_poismod <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + bili,
                    data = pbc3mult,
                    family = poisson)
summary(base_poismod)


poismod_t26_l1 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + bili + albnorm,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t26_l1) 

poismod_t26_l3 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + 
                        bili + bilihigh + bilitoohigh + bilimuchtoohigh,
                      data = pbc3mult,
                      family = poisson)

summary(poismod_t26_l3) 


# Make a dataset with linear predictor information
lin2 <- pbc3

lin2$lp1 <- with(pbc3, 1.6076-0.1123*alb)
lin2$lp2 <- with(pbc3, 0.7828-0.0864*alb-0.0474*albnorm)
lin2$lp3 <- with(pbc3, 1.6076+0.0085*bili-38.7*0.1123)
lin2$lp4 <- with(pbc3, -0.5534+0.0617*bili-0.0168*bilihigh+0.0027*bilitoohigh
                 -0.0428*bilimuchtoohigh-0.0865*38.7)


lin2row <- data.frame(effect = c(rep("Linear effect", nrow(lin2)),
                                 rep("Effect as linear spline", nrow(lin2))), 
                      lp = c(lin2$lp1, lin2$lp2),
                      alb = c(lin2$alb, lin2$alb))

require(ggplot2)

fig26 <- ggplot(aes(x = alb, y = lp, linetype = effect), data = lin2row) + 
         geom_line(size = 1) + 
         xlab("Albumin") + 
         ylab("Linear predictor") + 
         scale_linetype_discrete("Effect") + 
         scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) + 
         scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
         theme_general

fig26

ggsave("figures/j_pbc3linalb.pdf", plot = fig26, 
       width = 29.7, height = 21, units = "cm")




#------------------------------------------------------------------#
#---------------- Figure 2.7 --------------------------------------#
#------------------------------------------------------------------#

# bilirubin and the two last linear predictors
lin2row2 <- data.frame(effect = c(rep("Linear effect", nrow(lin2)),
                                  rep("Effect as linear spline", nrow(lin2))), 
                      lp = c(lin2$lp3, lin2$lp4),
                      bili = c(lin2$bili, lin2$bili))

require(ggplot2)

fig27 <- ggplot(aes(x = bili, y = lp, linetype = effect), data = lin2row2) + 
  geom_line(size = 1) + 
  xlab("Bilirubin") + 
  ylab("Linear predictor") + 
  scale_linetype_discrete("Effect") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig27

ggsave("figures/j_pbc3linbili.pdf", plot = fig27, 
       width = 29.7, height = 21, units = "cm")




#------------------------------------------------------------------#
#---------------- Figure 2.8 --------------------------------------#
#------------------------------------------------------------------#


# The below linear predictors include estimates from the following models, 

base_poismod2_log <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb10 + log2bili,
                         data = pbc3mult,
                         family = poisson)
summary(base_poismod2_log)

poismod_t26_l3 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + 
                        log2bili + logbilihigh + logbilitoohigh + logbilimuchtoohigh,
                        data = pbc3mult,
                        family = poisson)

summary(poismod_t26_l3)


# Make a dataset with linear predictor information
log2 <- pbc3

log2$lp3 <- with(pbc3, -2.0162+0.6469*log2bili-38.7*0.087)
log2$lp4 <- with(pbc3, -0.6194+0.198*log2bili+0.8815*logbilihigh-0.2336*logbilitoohigh
                 -0.3139*logbilimuchtoohigh-0.0844*38.7)


log2row <- data.frame(effect = c(rep("Linear effect", nrow(log2)),
                                 rep("Effect as linear spline", nrow(log2))), 
                      lp = c(log2$lp3, log2$lp4),
                      log2bili = c(log2$log2bili, log2$log2bili))

# Make plot
fig28 <- ggplot(aes(x = log2bili, y = lp, linetype = effect), data = log2row) + 
  geom_line(size = 1) + 
  xlab(expression(log[2] * "(bilirubin)")) +
  ylab("Linear predictor") + 
  scale_linetype_discrete("Effect") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig28

ggsave("figures/j_pbc3linlogbili.pdf", plot = fig28, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Table 2.7 ---------------------------------------#
#------------------------------------------------------------------#

# Cox model
coxmod_t27 <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili,
                    eps = 1e-8, method = "breslow",
                    data = pbc3)

summary(coxmod_t27)


# Poisson model
poismod_t27 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili,
                   data = pbc3mult,
                   family = poisson)
summary(poismod_t27) 



#------------------------------------------------------------------#
#---------------- Table 2.8 ---------------------------------------#
#------------------------------------------------------------------#


## Cox models

# Model for LR comparison
coxmod_t28_base <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili,
                      eps = 1e-8, method = "breslow",
                      data = pbc3)


# Cox model 1
coxmod_t28_1 <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili + tment*alb,
                      eps = 1e-8, method = "breslow",
                      data = pbc3)

summary(coxmod_t28_1)

# LR test
anova(coxmod_t28_1, 
      coxmod_t28_base)


# Cox model 2
coxmod_t28_2 <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili + tment*log2bili,
                      eps = 1e-8, method = "breslow",
                      data = pbc3)

summary(coxmod_t28_2)

# LR test
anova(coxmod_t28_2, 
      coxmod_t28_base)




## Poisson models

# Model for LR comparison - no interaction
poismod_t28_base <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili,
                     data = pbc3mult,
                     family = poisson)

# Poisson model 1
pbc3mult$tment_char <- relevel(pbc3mult$tment_char, ref = "Placebo")
poismod_t28_1 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + tment_char*alb,
                    data = pbc3mult,
                    family = poisson)
summary(poismod_t28_1) 


pbc3mult$tment_char <- relevel(pbc3mult$tment_char, ref = "CyA")
poismod_t28_11 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + tment_char*alb,
                     data = pbc3mult,
                     family = poisson)
summary(poismod_t28_11) 

# LR test
lrtest(poismod_t28_base, 
       poismod_t28_1)


# Poisson model 2
pbc3mult$tment_char <- relevel(pbc3mult$tment_char, ref = "Placebo")
poismod_t28_2 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + tment_char*log2bili,
                     data = pbc3mult,
                     family = poisson)
summary(poismod_t28_2) 


pbc3mult$tment_char <- relevel(pbc3mult$tment_char, ref = "CyA")
poismod_t28_21 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + tment_char*log2bili,
                     data = pbc3mult,
                     family = poisson)
summary(poismod_t28_21) 


# LR test
lrtest(poismod_t28_base, 
       poismod_t28_2)


#------------------------------------------------------------------#
#---------------- Table 2.9 ---------------------------------------#
#------------------------------------------------------------------#


# Column 1; 
pbc3mult$tment_char <- relevel(pbc3mult$tment_char, ref = "Placebo")
pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "1")

poismod_t29_11 <- glm(fail ~ offset(log(risktime)) + timegroup + tment_char + alb + log2bili + tment_char*timegroup,
                     data = pbc3mult,
                     family = poisson)
summary(poismod_t29_11) 

pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "2")
poismod_t29_12 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + tment_char*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_12) 


pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "3")
poismod_t29_13 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + tment_char*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_13) 



# Column 2; 
pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "1")
poismod_t29_21 <- glm(fail ~ offset(log(risktime)) + timegroup + tment_char + alb + log2bili + alb*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_21) 

pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "2")
poismod_t29_22 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + alb*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_22) 


pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "3")
poismod_t29_23 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + alb*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_23) 


# Column 3; 
pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "1")
poismod_t29_31 <- glm(fail ~ offset(log(risktime)) + timegroup + tment_char + alb + log2bili + log2bili*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_31) 

pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "2")
poismod_t29_32 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + log2bili*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_32) 


pbc3mult$timegroup <- relevel(as.factor(pbc3mult$timegroup), ref = "3")
poismod_t29_33 <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili + log2bili*timegroup,
                      data = pbc3mult,
                      family = poisson)
summary(poismod_t29_33) 


## LR tests 

# Base model; no interaction

poismod_t29_base <- glm(fail ~ offset(log(risktime)) + as.factor(timegroup) + tment_char + alb + log2bili,
                      data = pbc3mult,
                      family = poisson)
# Column 1
lrtest(poismod_t29_base, 
       poismod_t29_11)

# Column 2
lrtest(poismod_t29_base, 
       poismod_t29_21)

# Column 3
lrtest(poismod_t29_base, 
       poismod_t29_31)



#------------------------------------------------------------------#
#---------------- Figure 2.9 --------------------------------------#
#------------------------------------------------------------------#

coxmod_f29 <- coxph(Surv(days, status != 0) ~ strata(tment) + alb + log2bili,
                      eps = 1e-8, method = "breslow",
                      data = pbc3)

summary(coxmod_f29)

# Extracting cumulative baseline hazards per treatment
cumhaz_treat <- basehaz(coxmod_f29, centered = FALSE)
cumhaz_treat <- as.data.frame(cumhaz_treat)

# Per treatment
hazard_t0 <- cumhaz_treat[cumhaz_treat$strata == "tment=0",]
hazard_t0[1,] <- c(0, 0, "tment=0")
hazard_t0$time<- as.numeric(hazard_t0$time)
hazard_t0$hazard <- as.numeric(hazard_t0$hazard)
hazard_t1 <- cumhaz_treat[cumhaz_treat$strata == "tment=1",]

# Match times
alltimes <- sort(unique(cumhaz_treat$time))

hazard_t0_allt <- as.numeric(sapply(1:length(alltimes), function(k) tail(hazard_t0$hazard[hazard_t0$time <= alltimes[k]], 1)))
hazard_t1_allt <- as.numeric(sapply(1:length(alltimes), function(k) tail(hazard_t1$hazard[hazard_t1$time <= alltimes[k]], 1))) 

hazards <- data.frame(hazard_t0_allt, 
                      hazard_t1_allt)

# Extract coefficient
coxmod_f29_t <- coxph(Surv(days, status != 0) ~ tment + alb + log2bili,
                      eps = 1e-8, method = "breslow",
                      data = pbc3)

summary(coxmod_f29_t)




# Make plot
fig29 <- ggplot(aes(x = hazard_t0_allt, y = hazard_t1_allt), data = hazards) + 
  geom_step(size = 1) + 
  geom_abline(intercept = 0, slope = exp(coef(coxmod_f29_t)[["tment"]]), linetype = "dashed", size = 1) + 
  xlab("Cumulative baseline hazard: placebo") +
  ylab("Cumulative baseline hazard: CyA") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig29

ggsave("figures/j_pbc3gofplot.pdf", plot = fig29, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 2.10 -------------------------------------#
#------------------------------------------------------------------#

# Additive Aalen models - available with timereg
require(timereg)

# Fit model
nonparmod <- aalen(Surv(days, status != 0) ~ tment, data = pbc3)


cumhazdata <- data.frame(eventtimes = nonparmod$cum[,1],
                         basecumhaz = nonparmod$cum[,2], 
                         cumhaztreat = nonparmod$cum[,3], 
                         cumhaztreat_ll = nonparmod$cum[,3]-1.96*sqrt(nonparmod$var.cum[,3]),
                         cumhaztreat_ul = nonparmod$cum[,3]+1.96*sqrt(nonparmod$var.cum[,3])
                         )
cumhazdata


# Extend lines to last observed time
cumhazdata[nrow(cumhazdata)+1,] <- c(max(pbc3$days), tail(cumhazdata, 1)[-1])



# Figure LHS
cumhazplot1 <- ggplot(aes(x = eventtimes / 365.25, y = basecumhaz), data = cumhazdata) + 
  geom_step(size = 1) + 
  xlab("Time since randomization (years)") +
  ylab("Cumulative baseline hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

cumhazplot1


# Figure RHS
cumhazplot2 <- ggplot(aes(x = eventtimes / 365.25, y = cumhaztreat), data = cumhazdata) + 
  geom_step(size = 1) + 
  geom_step(size = 1, aes(x = eventtimes / 365.25, y = cumhaztreat_ll), linetype = "dashed") +  
  geom_step(size = 1, aes(x = eventtimes / 365.25, y = cumhaztreat_ul), linetype = "dashed") +
  xlab("Time since randomization (years)") +
  ylab("Cumulative treatment effect") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

cumhazplot2

# Combine plots
require(gridExtra)
fig210 <- grid.arrange(cumhazplot1, cumhazplot2, ncol = 2)
fig210

ggsave("figures/j_pbc3aalen1.pdf", plot = fig210, 
       width = 29.7, height = 12, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 2.11 -------------------------------------#
#------------------------------------------------------------------#

# Make Aalen model fit
nonparmod2 <- aalen(Surv(days, status != 0) ~ tment + alb + bili, data = pbc3) 

summary(nonparmod2)


cumhazdata <- data.frame(eventtimes = nonparmod2$cum[,1],
                         basecumhaz = nonparmod2$cum[,2], 
                         cumhaztreat = nonparmod2$cum[,3], 
                         cumhazalb = nonparmod2$cum[,4], 
                         cumhazbili= nonparmod2$cum[,5]
)
cumhazdata


# Extend lines to last observed time
cumhazdata[nrow(cumhazdata)+1,] <- c(max(pbc3$days), tail(cumhazdata, 1)[-1])


# Figure treat
cumhazplot1 <- ggplot(aes(x = eventtimes / 365.25, y = cumhaztreat), data = cumhazdata) + 
  geom_step(size = 1) + 
  xlab("Time since randomization (years)") +
  ylab("Cumulative treatment effect") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

cumhazplot1


# Figure albumin
cumhazplot2 <- ggplot(aes(x = eventtimes / 365.25, y = cumhazalb), data = cumhazdata) + 
  geom_step(size = 1) + 
  xlab("Time since randomization (years)") +
  ylab("Cumulative albumin effect") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05)), 
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

cumhazplot2


# Figure bilirubin
cumhazplot3 <- ggplot(aes(x = eventtimes / 365.25, y = cumhazbili), data = cumhazdata) + 
  geom_step(size = 1) + 
  xlab("Time since randomization (years)") +
  ylab("Cumulative bilirubin effect") + 
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05)),
                     limits = c(0, 6),
                     breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

cumhazplot3


# Combine plots
require(gridExtra)
fig211 <- grid.arrange(cumhazplot1, cumhazplot2, cumhazplot3, ncol = 1)
fig211

ggsave("figures/j_pbc3aalen3.pdf", plot = fig211, 
       width = 29.7, height = 35, units = "cm")




#------------------------------------------------------------------#
#---------------- Table 2.10 --------------------------------------#
#------------------------------------------------------------------#

# -- Only treatment -- #
# In text p-values prior to Table 2.10
nonparmod0 <- aalen(Surv(days, status != 0) ~ tment, data = pbc3) 
summary(nonparmod0)
# Hypothesis of a time-constant hazard difference for treatment

# Constant effect of treatment
nonparmod01 <- aalen(Surv(days, status != 0) ~ const(tment), data = pbc3) 
summary(nonparmod01)

# -- Treatment, albumin, bilirubin -- #
# Table 2.10, first two columns
nonparmod1 <- aalen(Surv(days, status != 0) ~ tment + alb + bili, data = pbc3) 
summary(nonparmod1)


# Table 2.10, last columns
nonparmod2 <- aalen(Surv(days, status != 0) ~ const(tment) + const(alb) + const(bili), data = pbc3) 
summary(nonparmod2)


# In-text
nonparmod3 <- aalen(Surv(days, status != 0) ~ const(tment) + alb + bili, data = pbc3) 
summary(nonparmod3)

# -- Quadratic effects -- #
nonparmod41 <- aalen(Surv(days, status != 0) ~ tment + I(alb/10) + I(bili/100) + I((bili/100)^2), data = pbc3) 
summary(nonparmod41)


nonparmod42 <- aalen(Surv(days, status != 0) ~ tment + I(alb/10) + I(bili/100) + I((alb/10)^2), data = pbc3) 
summary(nonparmod42)

nonparmod43 <- aalen(Surv(days, status != 0) ~ const(tment) + I(alb/10) + I(bili/100) + I((bili/100)^2), data = pbc3) 
summary(nonparmod43)

nonparmod44 <- aalen(Surv(days, status != 0) ~ const(tment) + I(alb/10) + I(bili/100) + I((alb/10)^2), data = pbc3) 
summary(nonparmod44)


# -- Interactions -- #
nonparmod51 <- aalen(Surv(days, status != 0) ~ const(tment) + const(alb) + const(bili) + const(tment * bili), data = pbc3) 
summary(nonparmod51)


nonparmod52 <- aalen(Surv(days, status != 0) ~ const(tment) + const(alb) + const(bili) + const(tment * alb), data = pbc3) 
summary(nonparmod52)



#------------------------------------------------------------------#
#---------------- Table 2.11 --------------------------------------#
#------------------------------------------------------------------#


# -- In-text -- #

# Additive hazards model with piecewise constant baseline hazards
# Model with only treatment as covariate

head(pbc3mult)

# update data set
pbc3add <- pbc3mult
pbc3add$time1 <- with(pbc3add, (timegroup == 1)*risktime)
pbc3add$time2 <- with(pbc3add, (timegroup == 2)*risktime)
pbc3add$time3 <- with(pbc3add, (timegroup == 3)*risktime)
pbc3add$tment0 <- with(pbc3add, (tment == 0)*risktime)
pbc3add$tment1 <- with(pbc3add, (tment == 1)*risktime)
pbc3add$albny <- with(pbc3add, (alb-35)/100*risktime)
pbc3add$biliny <- with(pbc3add, (bili-50)/1000*risktime)


# In-text
piece_aalen1 <- glm(fail ~ time1 + time2 + time3 + tment1 - 1,
                    data = pbc3add, start = c(0.1, 0.1, 0.1, 0),
                    family = poisson(link = "identity"))
summary(piece_aalen1)


# Table 2.11 - questionable fit
piece_aalen2 <- glm(fail ~ time1 + time2 + time3 + tment1 + albny + biliny -1,
                    data = pbc3add, start = c(0.3, 0.35, 0.4, -0.05, -0.8, 2),
                    family = poisson(link = "identity"))
summary(piece_aalen2)




