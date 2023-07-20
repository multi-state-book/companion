#------------------------------------------------------------------#
#------- Chapter 4, R code, affective data-------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")


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
# -------- General plotting style ---------------------------------# 
#------------------------------------------------------------------#

# General theme
theme_general <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 

#------------------------------------------------------------------#
#---------------- Table 4.4 ---------------------------------------#
#------------------------------------------------------------------#

# Read mstate data from sas file
#require(haven)
#affectivemstate <- read_sas(data_file = "C:/Users/jukf/University of Copenhagen/UCPH_MSB - General/Julie programs/data/affectivemstate.sas7bdat")
#affectivemstate <- as.data.frame(affectivemstate)


# Make dataset ready for mstate 

# From Out -> In, trans = 1

# From Out -> Dead, trans = 2


# From In -> Out, trans = 3

# From In -> Dead, trans = 4

# + update status variable


require(dplyr)

affectivemstate__ <- affective %>% mutate(statusnew = ifelse(status == 3, 0, 1), 
                                  trans = case_when(state == 0 & status == 1 ~ 1, 
                                                    state == 0 & status == 2 ~ 2, 
                                                    state == 1 & status == 0 ~ 3,
                                                    state == 1 & status == 2 ~ 4, 
                                                    state == 0 & status == 3 ~ 1, 
                                                    state == 1 & status == 3 ~ 3)
)

head(affectivemstate__)

# For each transition, we should have a censoring for the trans to the other state
affectivemstate_ <- affectivemstate__ %>% mutate(statusnew = 0, 
                                               trans = case_when(trans == 1 ~ 2, 
                                                                 trans == 2 ~ 1, 
                                                                 trans == 3 ~ 4, 
                                                                 trans == 4 ~ 3))


affectivemstate <- rbind(affectivemstate__, affectivemstate_) %>% arrange(id, start)


affectivemstate <- affectivemstate %>% mutate(from = case_when(trans == 1 ~ 1, 
                                                               trans == 2 ~ 1, 
                                                               trans == 3 ~ 2, 
                                                               trans == 4 ~ 2), 
                                              to =   case_when(trans == 1 ~ 2, 
                                                               trans == 2 ~ 3, 
                                                               trans == 3 ~ 1, 
                                                               trans == 4 ~ 3),
                                              starty = start /12, 
                                              stopy = stop / 12)

head(affectivemstate)


require(mstate)

# Subset data
affective0 <- subset(affectivemstate, bip == 0)
affective1 <- subset(affectivemstate, bip == 1)


# Set-up transition matrix
tmat <- matrix(NA, 3, 3)
tmat[1, 2:3] <- 1:2
tmat[2, c(1,3)] <- 3:4
statenames <- c("Out of hospital", "In hospital", "Dead")
dimnames(tmat) <- list(from = statenames, to = statenames)
tmat

# For unipolar (bip = 0)
attr(affective0, 'class') <- c("msdata","data.frame")
attr(affective0, 'trans') <- tmat
head(affective0)

# Transitions
# affective0$trans  <- NA
# for(i in 1:nrow(affective0)) { affective0$trans[i] <- tmat[affective0$from[i],affective0$to[i]] }
# events(affective0)

# Fit cox models per trans
c0 <- coxph(Surv(starty, stopy, statusnew) ~ strata(trans), 
            data = affective0)
summary(c0)

# Make a mstate object
msf0 <- msfit(object=c0, trans=tmat)

summary(msf0)
pt0 <- probtrans(msf0, predt=0)
summary(pt0, from=2)


# For bipolar (bip = 1)
attr(affective1, 'class') <- c("msdata","data.frame")
attr(affective1, 'trans') <- tmat
head(affective1)

# affective1$trans  <- NA
# for(i in 1:NROW(affective1)) { affective1$trans[i] <- tmat[affective1$from[i],affective1$to[i]] }
# events(affective1)

c1 <- coxph(Surv(starty, stopy, statusnew) ~ strata(trans), data=affective1)
summary(c1)
msf1 <- msfit(object=c1, trans=tmat)

summary(msf1)
pt1 <- probtrans(msf1, predt=0)
summary(pt1, from=2)



# Boostrap

regcoefvec <- function(data, tmat, tau) {
  cx <- coxph(Surv(starty, stopy, statusnew) ~ strata(trans), data=data)
  msf0 <- msfit(object = cx, trans = tmat)
  pt0 <- probtrans(msf0, predt=0)
  mat <- ELOS(pt0, tau=tau)
  return(mat[2,])
}

# Unipolar
regcoefvec(affective0, tmat, 15)
set.seed(1234)
res <- msboot(theta=regcoefvec, data=affective0, B=1000, id="id", tmat=tmat, tau=15)
mean(res[1,]); sqrt(var(res[1,]))
mean(res[2,]); sqrt(var(res[2,]))
mean(res[3,]); sqrt(var(res[3,]))


# Bipolar
regcoefvec(affective1, tmat, 15)
set.seed(1234)
res <- msboot(theta=regcoefvec, data=affective1, B=1000, id="id", tmat=tmat, tau=15)
mean(res[1,]); sqrt(var(res[1,]))
mean(res[2,]); sqrt(var(res[2,]))
mean(res[3,]); sqrt(var(res[3,]))


#------------------------------------------------------------------#
#---------------- Figure 4.18 -------------------------------------#
#------------------------------------------------------------------#

## Make data set with probabilities - predictions from state 2 (in hospital)

# bip = 0 
t0 <- data.frame(
  time = pt0[[2]]$time, 
  pstate1 = pt0[[2]]$pstate1, 
  pstate2 = pt0[[2]]$pstate2, 
  pstate3 = pt0[[2]]$pstate3, 
  bip = rep("No", nrow(pt0[[2]]))
)

# bip = 1
t1 <- data.frame(
  time = pt1[[2]]$time, 
  pstate1 = pt1[[2]]$pstate1, 
  pstate2 = pt1[[2]]$pstate2, 
  pstate3 = pt1[[2]]$pstate3, 
  bip = rep("Yes", nrow(pt1[[2]]))
)

head(t0); tail(t0)
head(t1); tail(t1)

pstate1 <- data.frame(type = rep("Out of hospital", nrow(pt0[[2]]) + nrow(pt1[[2]])), 
                      bip = c(rep("No", nrow(pt0[[2]])), rep("Yes", nrow(pt1[[2]]))),
                      pstate = c(pt0[[2]]$pstate1, pt1[[2]]$pstate1), 
                      time = c(pt0[[2]]$time, pt1[[2]]$time))

pstate2 <- data.frame(type = rep("In hospital", nrow(pt0[[2]]) + nrow(pt1[[2]])), 
                      bip = c(rep("No", nrow(pt0[[2]])), rep("Yes", nrow(pt1[[2]]))),
                      pstate = c(pt0[[2]]$pstate2, pt1[[2]]$pstate2), 
                      time = c(pt0[[2]]$time, pt1[[2]]$time))

pstate3 <- data.frame(type = rep("Dead", nrow(pt0[[2]]) + nrow(pt1[[2]])), 
                      bip = c(rep("No", nrow(pt0[[2]])), rep("Yes", nrow(pt1[[2]]))),
                      pstate = c(pt0[[2]]$pstate3, pt1[[2]]$pstate3), 
                      time = c(pt0[[2]]$time, pt1[[2]]$time))

# pstate1$pstate1 <- pstate1$pstate
# pstate2$pstate2 <- pstate2$pstate
# pstate3$pstate3 <- pstate3$pstate
# 
# probs_wide <- cbind(pstate1, pstate2, pstate3)
# probs_wide$psum <- with(probs_wide, pstate1 + pstate2 + pstate3)


probs <- rbind(pstate1, pstate2, pstate3)



# A couple of places, a very small negative probability is predicted
# We fix it + sum here

#probs[probs$pstate <0,]$pstate <-0
t <- probs[probs$pstate <0,]$time #<-0

probs[probs$bip == "No" & probs$time == t[1], "pstate"]  <- 
  c(probs[probs$bip == "No" & probs$time == t[1], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "No" & probs$time == t[1], "pstate"] %*% c(0,0,1))

probs[probs$bip == "No" & probs$time == t[2], "pstate"]  <- 
  c(probs[probs$bip == "No" & probs$time == t[2], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "No" & probs$time == t[2], "pstate"] %*% c(0,0,1))

probs[probs$bip == "No" & probs$time == t[3], "pstate"]  <- 
  c(probs[probs$bip == "No" & probs$time == t[3], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "No" & probs$time == t[3], "pstate"] %*% c(0,0,1))

probs[probs$bip == "No" & probs$time == t[4], "pstate"]  <- 
  c(probs[probs$bip == "No" & probs$time == t[4], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "No" & probs$time == t[4], "pstate"] %*% c(0,0,1))


probs[probs$bip == "Yes" & probs$time == t[5], "pstate"]  <- 
  c(probs[probs$bip == "Yes" & probs$time == t[5], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "Yes" & probs$time == t[5], "pstate"] %*% c(0,0,1))


probs[probs$bip == "Yes" & probs$time == t[6], "pstate"]  <- 
  c(probs[probs$bip == "Yes" & probs$time == t[6], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "Yes" & probs$time == t[6], "pstate"] %*% c(0,0,1))


probs[probs$bip == "Yes" & probs$time == t[7], "pstate"]  <- 
  c(probs[probs$bip == "Yes" & probs$time == t[7], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "Yes" & probs$time == t[7], "pstate"] %*% c(0,0,1))

probs[probs$bip == "Yes" & probs$time == t[8], "pstate"]  <- 
  c(probs[probs$bip == "Yes" & probs$time == t[8], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "Yes" & probs$time == t[8], "pstate"] %*% c(0,0,1))

probs[probs$bip == "Yes" & probs$time == t[9], "pstate"]  <- 
  c(probs[probs$bip == "Yes" & probs$time == t[9], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "Yes" & probs$time == t[9], "pstate"] %*% c(0,0,1))

probs[probs$bip == "Yes" & probs$time == t[10], "pstate"]  <- 
  c(probs[probs$bip == "Yes" & probs$time == t[10], "pstate"] %*% c(1,1,0), 
    0, 
    probs[probs$bip == "Yes" & probs$time == t[10], "pstate"] %*% c(0,0,1))




probs2 <- probs[order(probs$bip, probs$type, probs$time, probs$pstate, decreasing = F),]
#View(probs2)

require(ggpattern)
p1 <- ggplot(aes(x = time), data = subset(probs2, bip == "No")) +
  geom_area_pattern(aes(y = pstate, 
                        pattern = type,
                        pattern_fill = type ), 
                    fill = 'white', 
                    colour = 'black', 
                    #pattern_density = 0.02, 
                    pattern_aspect_ratio = 1,
                    pattern_fill = 'darkgrey', 
                    pattern_color = 'black', 
                    pattern_spacing = 0.02,
                    linewidth = 0.7) + 
  ylab("Probability") + 
  xlab("Time since first admission (years)") + 
  scale_pattern_discrete(name = c("State"), 
                         choices = c("circle", "stripe", "crosshatch")) + 
  scale_pattern_fill_discrete(name = c("State")) + 
  scale_pattern_spacing_discrete(name = c("State")) + 
  theme_general + ggtitle("Unipolar") + 
  theme(legend.key.size = unit(1, 'cm'))
p1

p2 <- ggplot(aes(x = time), data = subset(probs2, bip == "Yes")) +
  geom_area_pattern(aes(y = pstate, 
                        pattern = type,
                        pattern_fill = type ), 
                    fill = 'white', 
                    colour = 'black', 
                    #pattern_density = 0.02, 
                    pattern_aspect_ratio = 1,
                    pattern_fill = 'darkgrey', 
                    pattern_color = 'black', 
                    pattern_spacing = 0.02,
                    linewidth = 0.7) + 
  ylab("Probability") + 
  xlab("Time since first admission (years)") + 
  scale_pattern_discrete(name = c("State"), choices = c("circle", "stripe", "crosshatch")) + 
  scale_pattern_fill_discrete(name = c("State")) + 
  scale_pattern_spacing_discrete(name = c("State")) + 
  theme_general + ggtitle("Bipolar") + 
  theme(legend.key.size = unit(1, 'cm'))
p2

# common legend
require(grid); require(gridExtra)
plots <- list(p1, p2)
g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
tmp <- arrangeGrob(p1 + theme(legend.position = "none"), 
                   p2 + theme(legend.position = "none"), 
                   layout_matrix = matrix(c(1, 2), nrow = 1))


fig418 <- grid.arrange(tmp, legend, ncol = 1, 
                       heights = unit.c(unit(1, "npc") - lheight, lheight))


# Save figure
ggsave("figures/j_angstprobtrans.pdf", plot = fig418, 
       width = 29.7, height = 15, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.18 -------------------------------------#
#------------------------------------------------------------------#

# In years
affective <- affective %>% mutate(starty = start / 12, stopy = stop / 12) %>% group_by(id) %>% 
                mutate(prevy1 = lag(starty, n = 1, default = 0), 
                       prevy2 = lag(stopy, n = 1, default = 0),
                       prevy = ifelse(state == 1, prevy2, prevy1))


# Use Nelson-Aalen to estimate marginal mean 
# Incorrectly censoring for death 
naa_est <- survfit(Surv(prevy, stopy, status == 1) ~ bip + cluster(id), 
                   data = subset(affective, stopy > prevy & (state == 0 | status %in% c(2,3))), 
                   ctype = 1)

# Collect data for plotting
plotdata <- data.frame(time = naa_est$time, 
                       mu = naa_est$cumhaz, 
                       bip = c(rep("No", naa_est$strata[[1]]), 
                               rep("Yes", naa_est$strata[[2]])))


# Create Figure 4.18
fig418 <- ggplot(aes(x = time, y = mu, linetype = bip), data = plotdata) + 
  geom_step(size = 1) + 
  xlab("Time since first admission (years)") + 
  ylab("Expected number of episodes") + 
  scale_linetype_manual("Bipolar", values = c("dashed", "solid"),labels=c("Unipolar","Bipolar")) + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 30), breaks = seq(0, 30, by = 5)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  theme_general + 
  theme(legend.box = "vertical",
        legend.key.size = unit(1.5, 'cm'))
fig418
fig418<-fig418+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig418
# Save figure
ggsave("figures/h_angstcmfgal.pdf", plot = fig418, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Table 4.8 ---------------------------------------#
#------------------------------------------------------------------#


# Mortality not taken into account 
fit1 <- coxph(Surv(prevy, stopy, status == 1) ~ bip + cluster(id), 
              data = subset(affective, state == 0 | status %in% c(2,3)), 
              ties = "breslow")

summary(fit1)


# Mortality taken into account
require(mets)

fit2 <- recreg(Event(prevy, stopy, status) ~ bip + cluster(id),
               data = subset(affective, state == 0 | status %in% c(2,3)), 
               cause = 1, cens.code = 3, death.code = 2)

summary(fit2)

#------------------------------------------------------------------#
#---------------- Table 4.10 --------------------------------------#
#------------------------------------------------------------------#

# Number of recurrences and death

with(subset(affective, episode < 5), table(episode, status))

apply(with(subset(affective, episode < 5), table(episode, status)), 2, cumsum)

#------------------------------------------------------------------#
#---------------- Table 4.11 --------------------------------------#
#------------------------------------------------------------------#

# Make WLW data ready
affectivewlw <- read.csv("data/affectivewlw.csv")
head(affectivewlw)

affectivewlw <- affectivewlw %>% mutate(bip1 = bip * (stratum == 1), 
                                        bip2 = bip * (stratum == 2), 
                                        bip3 = bip * (stratum == 3), 
                                        bip4 = bip * (stratum == 4))



# Composite endpoint
fit1 <- coxph(Surv(time, dc %in% c(1, 2)) ~ bip1 + bip2 + bip3 + bip4 + cluster(id) + strata(stratum), 
              data = affectivewlw, 
              ties = "breslow")

summary(fit1)

fit2 <- coxph(Surv(time, dc %in% c(1, 2)) ~ bip + cluster(id) + strata(stratum), 
              data = affectivewlw, 
              ties = "breslow")

summary(fit2)


# Cause-specific hazard of recurrence
fit3 <- coxph(Surv(time, dc %in% c(1)) ~ bip1 + bip2 + bip3 + bip4 + cluster(id) + strata(stratum), 
              data = affectivewlw, 
              ties = "breslow")

summary(fit3)

fit4 <- coxph(Surv(time, dc %in% c(1)) ~ bip + cluster(id) + strata(stratum), 
              data = affectivewlw, 
              ties = "breslow")

summary(fit4)




#------------------------------------------------------------------#
#---------------- Figure 4.19 -------------------------------------#
#------------------------------------------------------------------#

# Mortality taken into account - mean number are estimated using (4.13)
require(mets)

xr <- phreg(Surv(prevy, stopy, status == 1) ~ strata(bip) + cluster(id), data = subset(affective, state == 0 | status %in% c(2,3)))
xd <- phreg(Surv(prevy, stopy, status == 2) ~ strata(bip) + cluster(id), data = subset(affective, state == 0 | status %in% c(2,3)))

out <- recurrentMarginal(xr, xd)


pout <- data.frame(time = out$cumhaz[,1], 
                   mu = out$cumhaz[,2],
                   bip = as.factor(out$strata))

  
NAa_fit <- survfit(Surv(prevy, stopy, status == 1) ~ strata(bip),
                   data = subset(affective, state == 0 | status %in% c(2,3)), id = id,
                   ctype = 1, timefix = FALSE)
  
KM_fit <- survfit(Surv(prevy, stopy, status == 2) ~ strata(bip),
                  data = subset(affective, state == 0 | status %in% c(2,3)), id = id, timefix = FALSE)
  
# Adjust hat(mu)
lS0 <- dplyr::lag(KM_fit$surv[1:(KM_fit$strata[1])], default = 1)
dA0 <- diff(NAa_fit$cumhaz[1:NAa_fit$strata[1]])
mu_adj0 <- cumsum(lS0 * c(0, dA0))


lS1 <- dplyr::lag(KM_fit$surv[(KM_fit$strata[1]+1):(KM_fit$strata[1] + KM_fit$strata[2])], default = 1)
dA1 <- diff(NAa_fit$cumhaz[(KM_fit$strata[1]+1):(KM_fit$strata[1] + KM_fit$strata[2])])
mu_adj1 <- cumsum(lS1 * c(0, dA1))


plotdata2 <- data.frame(time = KM_fit$time, 
                       mu = c(mu_adj0, mu_adj1), 
                       bip = c(rep("No", length(mu_adj0)), 
                               rep("Yes", length(mu_adj1))))


# Create Figure 4.19
fig419 <- ggplot(aes(x = time, y = mu, linetype = bip), data = plotdata2) + 
  geom_step(linewidth = 1) + 
  xlab("Time since first admission (years)") + 
  ylab("Expected number of episodes") + 
  scale_linetype_manual("Bipolar", values = c("dashed", "solid"),labels=c("Unipolar","Bipolar") ) + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 30), breaks = seq(0, 30, by = 5)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  theme_general + 
  theme(legend.box = "vertical",
        legend.key.size = unit(1.5, 'cm'))

fig419
fig419<-fig419+theme(legend.title=element_blank(), legend.text = element_text(size = 20))
fig419

# Save figure
ggsave("figures/h_angstcmfok.pdf", plot = fig419, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.23 -------------------------------------#
#------------------------------------------------------------------#

# Last observations
cens <- affective %>% group_by(id) %>% slice(c(n()))

# Censoring dist, KM
censdist <- survfit(Surv(stop, status == 3) ~ 1, 
                    data = cens)
censdist

# Make data ready for plotting
pdata <- data.frame(time = censdist$time,
                    surv = censdist$surv)

# Create Figure 4.23
fig423 <- ggplot(aes(x = time / 12, y = surv), data = pdata) +
  geom_step(size = 1) +
  xlab("Time since first admission (years)") +
  ylab('Probability of no censoring') +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)),
                     limits = c(0, 30),
                     breaks = seq(0, 30, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.005)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  theme(legend.box = "vertical",
        legend.key.size = unit(1.5, 'cm'))

fig423
ggsave("figures/h_angstcens.pdf", plot = fig423,
       width = 29.7, height = 21, units = "cm")


