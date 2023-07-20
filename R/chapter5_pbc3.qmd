#------------------------------------------------------------------#
#------- Chapter 5, R code, PBC3 data  ----------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

# Load required packages (should be installed if not already)
require(haven)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)
require(mets)
require(cmprsk)
require(pseudo)

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
#---------------- Table 5.6 ---------------------------------------#
#------------------------------------------------------------------#

library(mets)

pbcny <- subset(pbc3, !is.na(alb))

# COLUMN 1: t0 = 2
dirbin2tment <- binreg(Event(days, status) ~ tment, data = pbc3,
                     cause = 2, time = 2 * 365.25, cens.code = 0)
summary(dirbin2tment)


dirbin2 <- binreg(Event(days, status) ~ tment + I(alb - 40) + I(log2(bili) - 4.6), 
                  data = pbcny, cause = 2, time = 2 * 365.25, cens.code = 0)
summary(dirbin2)


# COLUMN 2: t1,t2,t3 = 1,2,3

library(foreign)


# trytimereg1 <- comp.risk(Event(days, status) ~ const(tment) + const(I(alb - 40)) +  
#                          const(I(log2(bili) - 4.6)), data = pbcny, cause = 2, times = 2 * 365.25,
#                        model = "logistic",resample.iid = 1, n.sim = 100, monotone = 1)
# 
# summary(trytimereg1)$coef

trytimereg3 <- comp.risk(Event(days, status) ~ const(tment), data = pbc3, 
                         cause = 2, times = c(1,2,3) * 365.25,
                         model = "logistic", resample.iid = 1, n.sim = 100, monotone = 1)

summary(trytimereg3)$coef


trytimereg3 <- comp.risk(Event(days, status) ~ const(tment) + const(I(alb - 40)) +  
                         const(I(log2(bili) - 4.6)), data = pbcny, cause = 2, times = c(1,2,3) * 365.25,
                         model = "logistic", resample.iid = 1, n.sim = 100, monotone = 1)

summary(trytimereg3)$coef


#------------------------------------------------------------------#
#---------------- Figure 5.12 -------------------------------------#
#------------------------------------------------------------------#


# Make Cox model fit
pbc3$followup <- pbc3$days / 365.25
fit <- cox.aalen(Surv(followup, status != 0) ~ prop(tment) + prop(alb) + prop(bili),
                 data = pbc3, n.sim = 0,
                 residuals = 1)
#summary(fit)

# Cumulative martingale residuals
set.seed(061166)
cum_res <- cum.residuals(fit, pbc3, cum.resid=1, max.point.func = 50, n.sim = 1000)
summary(cum_res)


cumresdata_alb <- data.frame(alb = unname(cum_res$proc.cumz[[1]][,1]),
                             cum_mg_res = unname(cum_res$proc.cumz[[1]][,2]))

cumresdata_albsim <- data.frame(alb = rep(unname(cum_res$proc.cumz[[1]][,1]), times = 50),
                                cum_mg_res = c(cum_res$sim.test.proccumz[[1]]),
                                sim = rep(1:50, each = length(unname(cum_res$proc.cumz[[1]][,1])))
)

# Figure 5.12
fig512 <- ggplot(aes(x = alb, y = cum_mg_res), data = cumresdata_alb) +
  geom_step(aes(x = alb, y = cum_mg_res, group = sim), color = "grey", size = 0.8, data = cumresdata_albsim) +
  geom_step(size = 1) +
  xlab("Albumin") +
  ylab("Cumulative martingale residuals") +
  #geom_text(x = 52, y = -5, label = paste0("P-value, \n supremum test: \n", cum_res$pval.test[1]), size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig512
cum_res$pval.test[1]

ggsave("figures/h_CumulativeResidualsalb.pdf", plot = fig512, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 5.13 -------------------------------------#
#------------------------------------------------------------------#


# Figure bili
cumresdata_bili <- data.frame(bili = unname(cum_res$proc.cumz[[2]][,1]),
                              cum_mg_res = unname(cum_res$proc.cumz[[2]][,2]))

cumresdata_bilisim <- data.frame(bili= rep(unname(cum_res$proc.cumz[[2]][,1]), times = 50),
                                 cum_mg_res = c(cum_res$sim.test.proccumz[[2]]),
                                 sim = rep(1:50, each = length(unname(cum_res$proc.cumz[[2]][,1])))
)

# Figure 5.13
fig513 <- ggplot(aes(x = bili, y = cum_mg_res), data = cumresdata_bili) +
  geom_step(aes(x = bili, y = cum_mg_res, group = sim), color = "grey", size = 0.8, data = cumresdata_bilisim) +
  geom_step(size = 1) +
  xlab("Bilirubin") +
  ylab("Cumulative martingale residuals") + 
  #geom_text(x = 375, y = -12, label = paste0("P-value, \n supremum test: \n", cum_res$pval.test[2]), size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig513


ggsave("figures/h_CumulativeResidualsbili.pdf", plot = fig513, 
       width = 29.7, height = 21, units = "cm")

#------------------------------------------------------------------#
#---------------- Figure 5.14 -------------------------------------#
#------------------------------------------------------------------#


# Make Cox model fit with log2bili
fit2 <- cox.aalen(Surv(followup, status != 0) ~ prop(tment) + prop(alb) + prop(log2bili),
                  data = pbc3, n.sim = 0,
                  residuals = 1)

# Cumulative martingale residuals
set.seed(061166)
cum_res <- cum.residuals(fit2, pbc3, cum.resid=1, max.point.func = 50, n.sim = 1000)
summary(cum_res)


cumresdata_log2bili <- data.frame(log2bili = unname(cum_res$proc.cumz[[2]][,1]),
                                  cum_mg_res = unname(cum_res$proc.cumz[[2]][,2]))

cumresdata_log2bilisim <- data.frame(log2bili = rep(unname(cum_res$proc.cumz[[2]][,1]), times = 50),
                                     cum_mg_res = c(cum_res$sim.test.proccumz[[2]]),
                                     sim = rep(1:50, each = length(unname(cum_res$proc.cumz[[2]][,1])))
)

# Figure 5.14
fig514 <- ggplot(aes(x = log2bili, y = cum_mg_res), data = cumresdata_log2bili) +
  geom_step(aes(x = log2bili, y = cum_mg_res, group = sim), color = "grey", size = 0.8, data = cumresdata_log2bilisim) +
  geom_step(size = 1) +
  xlab(expression("log" [2] * "(bilirubin)"))+
  ylab("Cumulative martingale residuals") +
#  geom_text(x = 8, y = -6, label = paste0("P-value, \n supremum test: \n", cum_res$pval.test[2]), size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig514
cum_res$pval.test[2]


ggsave("figures/h_CumulativeResidualslogbili.pdf", plot = fig514, 
       width = 29.7, height = 21, units = "cm")

#------------------------------------------------------------------#
#---------------- Figure 5.15 -------------------------------------#
#------------------------------------------------------------------#


# Compute Schoenfeld residuals (standardized)
set.seed(130966)
fit3 <- cox.aalen(Surv(followup, status != 0) ~ prop(tment) + prop(alb) + prop(log2bili),
                  data = pbc3, n.sim = 1000, residuals = 1,
                  weighted.test = 1)
summary(fit3)
par(mfrow = c(1,1))
plot(fit3, score = T)

time <- fit3$residuals$time
sim <- fit3$sim.test.procProp
obs <- fit3$test.procProp

data_tment_obs <- data.frame(time = time, 
                             res = obs[,2])

data_tment_sim <- data.frame(res = do.call("rbind", sim)[,1],
                             time = rep(time, times = 50), 
                             sim = rep(1:50, each = length(time))
)


# Figure 5.15
fig515 <- ggplot(aes(x = time, y = res), data = data_tment_obs) +
  geom_step(aes(x = time, y = res, group = sim), 
            color = "grey", size = 0.8, data = data_tment_sim) +
  geom_step(size = 1) +
  xlab(expression("Time since randomization (years)"))+
  ylab("Standardized score process") +
#  geom_text(x = 4.5, y = -2, label = paste0("P-value, \n supremum test: \n", fit3$pval.Prop[1]), size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig515
fit3$pval.Prop[1]

ggsave("figures/h_ScoreProcesstment.pdf", plot = fig515, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 5.16 -------------------------------------#
#------------------------------------------------------------------#

# Same but for albumin
data_alb_obs <- data.frame(time = time, 
                           res = obs[,3])

data_alb_sim <- data.frame(res = do.call("rbind", sim)[,2],
                           time = rep(time, times = 50), 
                           sim = rep(1:50, each = length(time))
)


# Figure 5.16
fig516 <- ggplot(aes(x = time, y = res), data = data_alb_obs) +
  geom_step(aes(x = time, y = res, group = sim), 
            color = "grey", size = 0.8, data = data_alb_sim) +
  geom_step(size = 1) +
  xlab(expression("Time since randomization (years)"))+
  ylab("Standardized score process") +
  #geom_text(x = 4.5, y = -2.5, label = paste0("P-value, \n supremum test: \n", fit3$pval.Prop[2]), size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig516
fit3$pval.Prop[2]


ggsave("figures/h_ScoreProcessalb.pdf", plot = fig516, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 5.17 -------------------------------------#
#------------------------------------------------------------------#


# Same but for log2bili
data_log2bili_obs <- data.frame(time = time, 
                                res = obs[,4])

data_log2bili_sim <- data.frame(res = do.call("rbind", sim)[,3],
                                time = rep(time, times = 50), 
                                sim = rep(1:50, each = length(time))
)


# Figure 5.17
fig517 <- ggplot(aes(x = time, y = res), data = data_log2bili_obs) +
  geom_step(aes(x = time, y = res, group = sim), 
            color = "grey", size = 0.8, data = data_log2bili_sim) +
  geom_step(size = 1) +
  xlab(expression("Time since randomization (years)"))+
  ylab("Standardized score process") +
#  geom_text(x = 4.5, y = -2, label = paste0("P-value, \n supremum test: \n", fit3$pval.Prop[3]), size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_general

fig517
fit3$pval.Prop[3]


ggsave("figures/h_ScoreProcesslogbili.pdf", plot = fig517, 
       width = 29.7, height = 21, units = "cm")


