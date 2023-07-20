#------------------------------------------------------------------#
#------- Chapter 6, R code, PBC3 data  ----------------------------#
#------------------------------------------------------------------#

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


#------------------------------------------------------------------#
#---------------- Figure 6.1 --------------------------------------#
#------------------------------------------------------------------#


# Add failure variable 
pbc3$fail <- as.numeric(with(pbc3, status > 0))


# Calculate pseudo-observations of the survival function for subjects
# The subjects with ptno=415 and ptno=458 are selected
pseudo_allt <- pseudosurv(pbc3$days, pbc3$fail)

# Re-arrange the data into a long data set
b <- NULL
for(it in 1:length(pseudo_allt$time)){
  b <- rbind(b,cbind(pbc3,
                     pseudo = pseudo_allt$pseudo[,it],
                     tpseudo = pseudo_allt$time[it],
                     id = 1:nrow(pbc3)))
}
b <- b[order(b$id),]


pseudo_alltid <- b

head(pseudo_alltid)

# Subset the two subjects 
subdat <- subset(pseudo_alltid, id %in% c("305", "325"))
head(subdat)

# Collect data for plot
pseudodata <- data.frame(tpseudo = subdat$tpseudo / 365.25, 
                         pseudo = subdat$pseudo, 
                         id = as.factor(subdat$id))

# Create Figure 6.1

fig61 <- ggplot(aes(x = tpseudo, y = pseudo, linetype = id), data = pseudodata) + 
  geom_hline(yintercept = c(0, 1), color = "darkgrey", size = 1) + 
  geom_step(size = 1) + 
  scale_linetype_manual("Patient number", values = c("dashed", "dotted")) + 
  xlab("Time since randomization (years)") + 
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme_general +
  theme(legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(1, 'cm'))

fig61
ggsave("figures/h_pbc3pseudo1.pdf", plot = fig61, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 6.2a --------------------------------------#
#------------------------------------------------------------------#

theme_general_1 <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 26), 
        axis.text.x = element_text(size = 26), 
        axis.text.y = element_text(size = 26)) 

# Select the pseudo-observations at times 
times <- c(366, 743, 1105)

# At 1 year
pseudo_t1 <- subset(pseudo_alltid, tpseudo == times[1])

# Collect data for plot
pseudodata <- data.frame(tpseudo = pseudo_t1$tpseudo / 365.25, 
                         days = pseudo_t1$days,
                         pseudo = pseudo_t1$pseudo, 
                         failtype = as.factor(pseudo_t1$fail))

# Create Figure 6.2
fig62 <- ggplot(aes(x = days  / 365.25, y = pseudo, shape = failtype), data = pseudodata) + 
  geom_point(size = 6) + 
  scale_shape_manual("Fail", values = c(4, 1)) + 
  xlab(expression("X"[i]*" (years)")) + 
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-0.9, 1.1)) +
  theme_general_1 +   
  theme(legend.position="none")

fig62

ggsave("figures/h_pbc3pseudo1yr.pdf", plot = fig62, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 6.2b --------------------------------------#
#------------------------------------------------------------------#

# At 2 years
pseudo_t2 <- subset(pseudo_alltid, tpseudo == times[2])

# Collect data for plot
pseudodata <- data.frame(tpseudo = pseudo_t2$tpseudo / 365.25, 
                         days = pseudo_t2$days,
                         pseudo = pseudo_t2$pseudo, 
                         failtype = as.factor(pseudo_t2$fail))

# Create Figure 6.3
fig63 <- ggplot(aes(x = days / 365.25, y = pseudo, shape = failtype), data = pseudodata) + 
  geom_point(size = 6) + 
  scale_shape_manual("Fail", values = c(4, 1)) + 
  xlab(expression("X"[i]*" (years)")) +  
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-0.9, 1.1)) +
  theme_general_1 +
  theme(legend.position="none")

fig63

ggsave("figures/h_pbc3pseudo2yr.pdf", plot = fig63, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 6.2c -------------------------------------#
#------------------------------------------------------------------#

# At 3 years
pseudo_t3 <- subset(pseudo_alltid, tpseudo == times[3])

# Collect data for plot
pseudodata <- data.frame(tpseudo = pseudo_t3$tpseudo / 365.25, 
                         days = pseudo_t3$days,
                         pseudo = pseudo_t3$pseudo, 
                         failtype = as.factor(pseudo_t3$fail))

# Create Figure 6.4
fig64 <- ggplot(aes(x = days/365.25, y = pseudo, shape = failtype), data = pseudodata) + 
  geom_point(size = 6) + 
  scale_shape_manual("Fail", values = c(4, 1)) + 
  xlab(expression("X"[i]*" (years)")) + 
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(0, 6, 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-0.9, 1.1)) +
  theme_general_1 +
  theme(legend.position="none")

fig64

ggsave("figures/h_pbc3pseudo3yr.pdf", plot = fig64, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 6.3 --------------------------------------#
#------------------------------------------------------------------#


# Collect data for plot
pseudodata <- data.frame(tpseudo = pseudo_t2$tpseudo / 365.25, 
                         days = pseudo_t2$days,
                         pseudo = pseudo_t2$pseudo, 
                         log2bili =  pseudo_t2$log2bili, 
                         bili = pseudo_t2$bili,
                         alb = pseudo_t2$alb,
                         tment = pseudo_t2$tment,
                         id = pseudo_t2$id,
                         failtype = as.factor(pseudo_t2$fail))

bili_loess <- loess(pseudo ~ bili, data = pseudodata, span = 0.8, degree = 1)
pseudodata$loesspred <- predict(bili_loess)

pseudodata$t_loesspred <- with(pseudodata, 
                               ifelse(loesspred > 0 , log(-log(loesspred)), NA))

# Create Figure 6.3 (a)
fig63_a <- ggplot(aes(x = bili, y = pseudo), data = pseudodata) + 
  geom_point(size = 2, shape = 4) + 
  geom_line(aes(x = bili, y = loesspred), size = 1) +
  xlab("Bilirubin") +  
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(-0.4, 1.2, by = 0.2), 
                     limits = c(-0.4, 1.2)) +
  theme_general 

fig63_a


# Create Figure 6.3 (b)

fig63_b <- ggplot(aes(x = bili, y = t_loesspred), data = pseudodata) + 
  geom_line(size = 1) + 
  xlab("Bilirubin") +  
  ylab("log(-log(predicted pseudo-values))") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(-5, 2, by = 1), 
                     limits = c(-5, 2)) +
  theme_general 

fig63_b

require(gridExtra)

fig63 <- grid.arrange(fig63_a, fig63_b, ncol = 2)
fig63


ggsave("figures/j_pbc3pseudoscatter.pdf", plot = fig63, 
       width = 29.7, height = 14, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 6.4 --------------------------------------#
#------------------------------------------------------------------#


# Loess for log2 bili
log2bili_loess <- loess(pseudo ~ log2bili, data = pseudodata, span = 0.8, degree = 1)
pseudodata$log2bili_loesspred <- predict(log2bili_loess)

pseudodata$t_log2bili_loesspred <- with(pseudodata, 
                               ifelse(log2bili_loesspred > 0 , log(-log(log2bili_loesspred)), NA))

# Create Figure 6.4 (a)
fig64_a <- ggplot(aes(x = log2bili, y = pseudo), data = pseudodata) + 
  geom_point(size = 2, shape = 4) + 
  geom_line(aes(x = log2bili, y = log2bili_loesspred), size = 1) +
  xlab(expression("log" [2] * "(bilirubin)")) +  
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(1, 9), 
                     breaks = seq(1, 9, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(-0.4, 1.2, by = 0.2), 
                     limits = c(-0.4, 1.2)) +
  theme_general 

fig64_a


# Create Figure 6.4 (b)
fig64_b <- ggplot(aes(x = log2bili, y = t_log2bili_loesspred), data = pseudodata) + 
  geom_line(size = 1) + 
  xlab(expression("log" [2] * "(bilirubin)")) +  
  ylab("log(-log(predicted pseudo-values))") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(1, 9), 
                     breaks = seq(1, 9, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(-5, 2, by = 1), 
                     limits = c(-5, 2)) +
  theme_general 

fig64_b

require(gridExtra)

fig64 <- grid.arrange(fig64_a, fig64_b, ncol = 2)
fig64


ggsave("figures/j_pbc3pseudoscatterlog.pdf", plot = fig64, 
       width = 29.7, height = 14, units = "cm")



#------------------------------------------------------------------#
#----------------Table 6.1 ----------------------------------------#
#------------------------------------------------------------------#

require(geepack)


# At time year 2
fit1 <- geese(formula = I(1 - pseudo) ~ log2bili + alb + tment, 
              data = subset(pseudodata, !is.na(alb)), 
              id = id, 
              mean.link = "cloglog",
              variance = "gaussian")
summary(fit1)


# At time year c(1, 2, 3) 

# Select the pseudo-observations at times 
times <- c(366, 743, 1105)

# At 1 year
pseudo_t1 <- subset(pseudo_alltid, tpseudo == times[1])
pseudo_t2 <- subset(pseudo_alltid, tpseudo == times[2])
pseudo_t3 <- subset(pseudo_alltid, tpseudo == times[3])

pseudo_t <- rbind(pseudo_t1, pseudo_t2, pseudo_t3)


pseudodata_t <- data.frame(tpseudo = pseudo_t$tpseudo / 365.25, 
                         days = pseudo_t$days,
                         pseudo = pseudo_t$pseudo, 
                         log2bili =  pseudo_t$log2bili, 
                         bili = pseudo_t$bili,
                         alb = pseudo_t$alb,
                         tment = pseudo_t$tment,
                         id = pseudo_t$id,
                         failtype = as.factor(pseudo_t$fail))


# At time year c(1, 2, 3)
fit2 <- geese(formula = I(1 - pseudo) ~ as.factor(tpseudo) + log2bili + alb + tment - 1, 
              data = subset(pseudodata_t, !is.na(alb)), 
              id = id, 
              mean.link = "cloglog",
              variance = "gaussian",
              scale.value = 0, 
              scale.fix = 1)
summary(fit2)



#------------------------------------------------------------------#
#---------------- Figure 6.5 --------------------------------------#
#------------------------------------------------------------------#

# Use regression estimates from fit2 


pseudodata_t$pred <- rep(NA, nrow(pseudodata_t))

# At time 1
pseudodata_t[pseudodata_t$tpseudo == times[1] / 365.25,]$pred <- 
  with(pseudodata_t[pseudodata_t$tpseudo == times[1] / 365.25,],
        fit2$beta[1] + fit2$beta[4]*log2bili + fit2$beta[5]*alb + fit2$beta[6]*tment)

# At time 2
pseudodata_t[pseudodata_t$tpseudo == times[2] / 365.25,]$pred <- 
  with(pseudodata_t[pseudodata_t$tpseudo == times[2] / 365.25,],
       fit2$beta[2] + fit2$beta[4]*log2bili + fit2$beta[5]*alb + fit2$beta[6]*tment)

# At time 3
pseudodata_t[pseudodata_t$tpseudo == times[3] / 365.25,]$pred <- 
  with(pseudodata_t[pseudodata_t$tpseudo == times[3] / 365.25,],
       fit2$beta[3] + fit2$beta[4]*log2bili + fit2$beta[5]*alb + fit2$beta[6]*tment)

# Transform back
pseudodata_t$res <- with(pseudodata_t, pseudo - exp(-exp(pred)))
pseudodata_t$time <- as.factor(pseudodata_t$tpseudo * 365.25)

# Make loess smooth per time
log2bili_res_loess_1 <- loess(res ~ log2bili, 
                              data = subset(pseudodata_t, time == 366), span = 0.7, degree = 1)
log2bili_res_loess_2 <- loess(res ~ log2bili, 
                              data = subset(pseudodata_t, time == 743), span = 0.7, degree = 1)
log2bili_res_loess_3 <- loess(res ~ log2bili, 
                              data = subset(pseudodata_t, time == 1105), span = 0.7, degree = 1)

log2bili_res_pred_1 <- predict(log2bili_res_loess_1, newdata = subset(pseudodata_t, time == 366))
log2bili_res_pred_2 <- predict(log2bili_res_loess_2, newdata = subset(pseudodata_t, time == 743))
log2bili_res_pred_3 <- predict(log2bili_res_loess_3, newdata = subset(pseudodata_t, time == 1105))

pseudodata_t$log2bili_res_pred <- c(log2bili_res_pred_1, log2bili_res_pred_2, log2bili_res_pred_3)


# Create Figure 6.5
fig65 <- ggplot(aes(x = log2bili, y = res, shape = time), data = pseudodata_t) + 
  geom_line(aes(x = log2bili, y = log2bili_res_pred , linetype = time), size = 1) +
  geom_point(size = 3) + 
  scale_shape_manual("Year", values = c(4, 2, 0), labels = c("1", "2", "3")) +
  scale_linetype_manual("Year", values = c("dashed", "solid", "dotted"), labels = c("1", "2", "3")) +
  xlab(expression("log" [2] * "(bilirubin)")) +  
  ylab("Pseudo-residuals") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(1, 9), 
                     breaks = seq(1, 9, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(-2, 1, by = 1), 
                     limits = c(-2, 1)) +
  theme_general + 
  theme(legend.box = "vertical",
        text = element_text(size=20), 
        legend.key.size = unit(1.5, 'cm'))


fig65

fig65<-fig65+theme(legend.text = element_text(size = 20))
fig65
ggsave("figures/h_pbc3pseudoresvslogbili.pdf", plot = fig65, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Table 6.2 ---------------------------------------#
#------------------------------------------------------------------#

# At time year c(1, 2, 3) using log-link instead
fit2 <- geese(formula = pseudo ~ as.factor(tpseudo) + bili + alb + tment - 1,
              data = subset(pseudodata_t, !is.na(alb)),
              id = id,
              mean.link = "log",
              variance = "gaussian",
              scale.value = 0,
              scale.fix = 1, maxit = 50)
summary(fit2)

## - log scale used as link in Table 6.2 (add "-" in front of all estimates)
# No default "-log(x)" link function in R

#------------------------------------------------------------------#
#---------------- Figure 6.6 --------------------------------------#
#------------------------------------------------------------------#


# Create Figure 6.6 (a) - done earlier
fig66_a <- fig63_a

fig66_a

# Same as earlier: 
bili_loess <- loess(pseudo ~ bili, data = pseudodata, span = 0.8, degree = 1)
pseudodata$loesspred <- predict(bili_loess)

# Use log-trans instead of cloglog
pseudodata$log_loesspred <- with(pseudodata, 
                               ifelse(loesspred > 0 , -log(loesspred), NA))
# Create Figure 6.6 (b)
fig66_b <- ggplot(aes(x = bili, y = log_loesspred), data = pseudodata) + 
  geom_line(size = 1) + 
  xlab("Bilirubin") +  
  ylab("-log(predicted pseudo-values)") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(-1, 3, by = 1), 
                     limits = c(-1, 3)) +
  theme_general 

fig66_b

fig66 <- grid.arrange(fig66_a, fig66_b, ncol = 2)
fig66


ggsave("figures/j_pbc3pseudoscatteradd.pdf", plot = fig66, 
       width = 29.7, height = 14, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 6.7 --------------------------------------#
#------------------------------------------------------------------#


# Use regression estimates from fit2 

pseudodata_t$pred <- rep(NA, nrow(pseudodata_t))

# At time 1
pseudodata_t[pseudodata_t$tpseudo == times[1] / 365.25,]$pred <- 
  with(pseudodata_t[pseudodata_t$tpseudo == times[1] / 365.25,],
       fit2$beta[1] + fit2$beta[4]*bili + fit2$beta[5]*alb + fit2$beta[6]*tment)

# At time 2
pseudodata_t[pseudodata_t$tpseudo == times[2] / 365.25,]$pred <- 
  with(pseudodata_t[pseudodata_t$tpseudo == times[2] / 365.25,],
       fit2$beta[2] + fit2$beta[4]*bili + fit2$beta[5]*alb + fit2$beta[6]*tment)

# At time 3
pseudodata_t[pseudodata_t$tpseudo == times[3] / 365.25,]$pred <- 
  with(pseudodata_t[pseudodata_t$tpseudo == times[3] / 365.25,],
       fit2$beta[3] + fit2$beta[4]*bili + fit2$beta[5]*alb + fit2$beta[6]*tment)

# Transform back
pseudodata_t$res <- with(pseudodata_t, pseudo - exp(pred))
pseudodata_t$time <- as.factor(pseudodata_t$tpseudo * 365.25)

# Make loess smooth per time
bili_res_loess_1 <- loess(res ~ bili, 
                              data = subset(pseudodata_t, time == 366), span = 0.7, degree = 1)
bili_res_loess_2 <- loess(res ~ bili, 
                              data = subset(pseudodata_t, time == 743), span = 0.7, degree = 1)
bili_res_loess_3 <- loess(res ~ bili, 
                              data = subset(pseudodata_t, time == 1105), span = 0.7, degree = 1)

bili_res_pred_1 <- predict(bili_res_loess_1, newdata = subset(pseudodata_t, time == 366))
bili_res_pred_2 <- predict(bili_res_loess_2, newdata = subset(pseudodata_t, time == 743))
bili_res_pred_3 <- predict(bili_res_loess_3, newdata = subset(pseudodata_t, time == 1105))

pseudodata_t$bili_res_pred <- c(bili_res_pred_1, bili_res_pred_2, bili_res_pred_3)



# Create Figure 6.7
fig67 <- ggplot(aes(x = bili, y = res, shape = time), data = pseudodata_t) + 
  geom_line(aes(x = bili, y = bili_res_pred , linetype = time), size = 1) +
  geom_point(size = 3) + 
  scale_shape_manual("Year", values = c(4, 2, 0), labels = c("1", "2", "3")) +
  scale_linetype_manual("Year", values = c("dashed", "solid", "dotted"), labels = c("1", "2", "3")) +  xlab(expression("Bilirubin")) +  
  ylab("Pseudo-residuals") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(0, 500), 
                     breaks = seq(0, 500, by = 100)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = seq(-2, 1, by = 1), 
                     limits = c(-2, 1)) +
  theme_general + 
  theme(legend.box = "vertical",
        text = element_text(size=20), 
        legend.key.size = unit(2, 'cm'))


fig67
fig67<-fig67+theme(legend.text = element_text(size = 20))
fig67

ggsave("figures/h_pbc3pseudoresvsbili.pdf", plot = fig67, 
       width = 29.7, height = 21, units = "cm")

#------------------------------------------------------------------#
#---------------- Table 6.3 ---------------------------------------#
#------------------------------------------------------------------#

# Calculate pseudo-observations of the restricted mean for subjects
pseudo_allt <- pseudomean(pbc3$days / 365.25, pbc3$fail, tmax = 3)

# Re-arrange the data into a long data set
a <- cbind(pbc3, pseudo = pseudo_allt)

outmean3 <- a

head(outmean3)

# Fit GEE model
fit_resm <- geese(pseudo ~ tment + log2bili + alb, id = id, 
                  data = subset(outmean3, !is.na(alb)),
                  mean.link = "identity", corstr = "independence", 
                  scale.fix = FALSE,
                  family = "gaussian")
summary(fit_resm)



#------------------------------------------------------------------#
#---------------- Figure 6.8 --------------------------------------#
#------------------------------------------------------------------#


# Collect data for plot
pseudodata <- data.frame(days = outmean3$days,
                         pseudo = outmean3$pseudo, 
                         failtype = as.factor(outmean3$fail))

# Create Figure 6.8
fig68 <- ggplot(aes(x = days / 365.25, y = pseudo, shape = failtype), data = pseudodata) + 
  geom_point(size = 6) + 
  scale_shape_manual("Fail", values = c(4, 1)) + 
  xlab(expression("X"[i]*" (years)")) + 
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(0, 6), 
                     breaks = seq(0, 6, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(0, 4)) +
  theme_general +
  theme(legend.position="none")

fig68

ggsave("figures/h_pbc3pseudomean3vstime.pdf", plot = fig68, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 6.9 --------------------------------------#
#------------------------------------------------------------------#


# Loess for log2 bili
log2bili_loess <- loess(pseudo ~ log2bili, data = outmean3, span = 0.8, degree = 1)
outmean3$log2bili_loesspred <- predict(log2bili_loess)

#outmean3$t_log2bili_loesspred <- with(outmean3, 
#                                        ifelse(log2bili_loesspred > 0 , 
#                                               log(log2bili_loesspred/(1-log2bili_loesspred)), NA))

# Create Figure 6.9
fig69 <- ggplot(aes(x = log2bili, y = pseudo), data = outmean3) + 
  geom_point(size = 6, shape = 4) + 
  geom_line(aes(x = log2bili, y = log2bili_loesspred), size = 1) +
  xlab(expression("log" [2] * "(bilirubin)")) +  
  ylab("Pseudo-values") + 
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)), limits = c(1, 9), 
                     breaks = seq(1, 9, by = 1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)), breaks = seq(0, 4, by = 1), 
                     limits = c(0, 4)) +
  theme_general 

fig69



ggsave("figures/j_pbc3pseudomean3vslogbili.pdf", plot = fig69, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Table 6.4 ---------------------------------------#
#------------------------------------------------------------------#

# Pseudo-obs
pseudo2 <- pseudoci(time = pbc3$days / 365.25,
                    event = pbc3$status, 
                    tmax = 2)

b2 <- NULL
for(it in 1:length(pseudo2$time)){
  b2 <- rbind(b2,cbind(pbc3,pseudo = pseudo2$pseudo[[2]][,it],
                       tpseudo = pseudo2$time[it]))
}
b2 <- b2[order(b2$id),]

# GEE fits
fit2 <- geese(pseudo ~ tment + I(alb-40) + I(log2(bili) - 4.6), data = b2, 
              id = id, jack = TRUE, scale.fix = TRUE, family = gaussian,
              mean.link = "logit", corstr = "independence")
summary(fit2)


fit2cloglog <- geese(pseudo ~ tment + I(alb-40) + I(log2(bili) - 4.6), data = b2, 
                     id = id, jack = TRUE, scale.fix = TRUE, family = gaussian,
                     mean.link = "cloglog", corstr = "independence")
summary(fit2cloglog)

# NU ALLE ptt (kun TMENT):
pseudo2tot <- pseudoci(time = pbc3$days / 365.25,
                       event = pbc3$status,
                       tmax = 2)

b2tot <- NULL
for(it in 1:length(pseudo2tot$time)){
  b2tot <- rbind(b2tot,cbind(pbc3,pseudo = pseudo2tot$pseudo[[2]][,it],
                             tpseudo = pseudo2tot$time[it],id=1:nrow(pbc3)))
}
b2tot <- b2tot[order(b2tot$id),]


fit2tment <- geese(pseudo ~ tment, data = b2tot, 
                   id = id, jack = TRUE, scale.fix = TRUE, family = gaussian,
                   mean.link = "logit", corstr = "independence")
summary(fit2tment)

fit2tmentcloglog <- geese(pseudo ~ tment, data = b2tot, 
                          id = id, jack = TRUE, scale.fix = TRUE, family = gaussian,
                          mean.link = "cloglog", corstr = "independence")
summary(fit2tmentcloglog)




#------------------------------------------------------------------#
#---------------- Table 6.5 ---------------------------------------#
#------------------------------------------------------------------#


pseudo123 <- pseudoci(time = pbc3$days / 365.25,
                      event = pbc3$status,
                      tmax = c(1,2,3))

b123 <- NULL
for(it in 1:length(pseudo123$time)){
  b123 <- rbind(b123,cbind(pbc3,pseudo = pseudo123$pseudo[[2]][,it],
                           tpseudo = pseudo123$time[it]))
}
b123 <- b123[order(b123$id),]

# fit123 <- geese(pseudo ~ as.factor(tpseudo) + tment + I(alb-40) + I(log2(bili) - 4.6), 
#                 data = b123, 
#                 id = id, jack = TRUE, scale.fix = TRUE, family = gaussian,
#                 mean.link = "cloglog", corstr = "independence")
# summary(fit123)

b123$sex <- relevel(as.factor(b123$sex), ref = "1")

fit123ny <- geese(pseudo ~ as.factor(tpseudo) + tment + I(alb-40) + I(log2(bili) - 4.6)
                  + as.factor(sex) + age, data = b123, 
                  id = id, jack = TRUE, scale.fix = TRUE, family = gaussian,
                  mean.link =  "cloglog", corstr = "independence")
summary(fit123ny)

pseudo10 <- pseudoci(time = pbc3$days / 365.25,
                     event = pbc3$status,
                     tmax = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))

b10 <- NULL
for(it in 1:length(pseudo10$time)){
  b10 <- rbind(b10,cbind(pbc3,pseudo = pseudo10$pseudo[[2]][,it],
                         tpseudo = pseudo10$time[it],id=1:nrow(pbc3)))
}
b10 <- b10[order(b10$id),]
b10$sex <- relevel(as.factor(b10$sex), ref = "1")

fit10 <- geese(pseudo ~ as.factor(tpseudo) + tment + I(alb-40) + I(log2(bili) - 4.6)
               + as.factor(sex) + age, data = b10, 
               id = id, jack = TRUE, scale.fix = TRUE, family = gaussian,
               mean.link = "cloglog", corstr = "independence")
summary(fit10)



#------------------------------------------------------------------#
#---------------- Table 6.6 ---------------------------------------#
#------------------------------------------------------------------#

# Lost years
pseudolost <- pseudoyl(time = pbc3$days / 365.25, 
                       event = as.integer(pbc3$status),
                       tmax = 3)

a1 <- cbind(pbc3,pseudo = pseudolost$pseudo[[1]],id=1:nrow(pbc3))

a2 <- cbind(pbc3,pseudo = pseudolost$pseudo[[2]],id=1:nrow(pbc3))

summary(fityl1 <- geese(pseudo ~ tment,
                        data = a1, id = id, jack = TRUE, family = gaussian,
                        corstr = "independence", scale.fix = FALSE))


summary(fityl2 <- geese(pseudo ~ tment,
                        data = a2, id = id, jack = TRUE, family = gaussian,
                        corstr = "independence", scale.fix = FALSE))

summary(fityl1 <- geese(pseudo ~ tment + alb + log2(bili),
                        data = subset(a1,!is.na(alb)), id = id, jack = TRUE, 
                        family = gaussian,
                        corstr = "independence", scale.fix = FALSE))

summary(fityl2 <- geese(pseudo ~ tment + alb + log2(bili),
                        data = subset(a2,!is.na(alb)), id = id, jack = TRUE, 
                        family = gaussian,
                        corstr = "independence", scale.fix = FALSE))




#------------------------------------------------------------------#
#---------------- Figure 6.10 -------------------------------------#
#------------------------------------------------------------------#


#------------------------------------------------------------------#
#---------------- Table 6.7 ---------------------------------------#
#------------------------------------------------------------------#


#------------------------------------------------------------------#
#---------------- Figure 6.11 -------------------------------------#
#------------------------------------------------------------------#

# comp. end point t=2 smlg med IJ
pbc3$nystatus<-ifelse(pbc3$status > 0, 1, 0)

pseudo <- pseudosurv(time = pbc3$days,
                     event = pbc3$nystatus,
                     tmax=2 * 365.25)

fit <- survfit(Surv(days, status > 0) ~ 1, data = pbc3)

IJpseudo <- pseudo(fit, times = 2 * 365.25, addNA = TRUE, data.frame = TRUE, minus1 = TRUE)



pdata <- data.frame(IJpseudo = IJpseudo$pseudo, 
                    pseudo = pseudo$pseudo, 
                    diff = IJpseudo$pseudo - pseudo$pseudo)



# Create Figure 6.11
fig611a <- ggplot(aes(x = pseudo, y = IJpseudo), data = pdata) + 
  geom_point(size = 4, shape = 1) + 
  geom_line(size = 1) +
  ylab("IJ pseudo-value") +  
  xlab("Pseudo-value") + 
  scale_x_continuous(limits = c(-0.2, 1), 
                     breaks = seq(-0.2, 1, by = 0.2)) + 
  scale_y_continuous(limits = c(-0.2, 1), 
                     breaks = seq(-0.2, 1, by = 0.2)) +
  theme_general 

fig611a

fig611b <- ggplot(aes(x = pseudo, y = IJpseudo - pseudo), data = pdata) + 
  geom_point(size = 4, shape = 1) + 
  ylab("Difference") +  
  xlab("Pseudo-value") + 
  scale_x_continuous(limits = c(-0.2, 1), 
                     breaks = seq(-0.2, 1, by = 0.2)) + 
  scale_y_continuous(limits = c(-0.01, 0.01), 
                     breaks = seq(-0.01, 0.01, by = 0.005)) +
  theme_general 

fig611b

require(gridExtra)
fig611 <- grid.arrange(fig611a, fig611b, ncol = 2)
fig611


ggsave("figures/j_IJpseudo.pdf", plot = fig611, 
       width = 29.7, height = 15, units = "cm")










