#------------------------------------------------------------------#
#------- Chapter 3, R code, bmt data ------------------------------#
#------------------------------------------------------------------#

# Make zeros print as "0" always
library(stringr)
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}

# Load required packages (should be installed if not already)
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

# Add extra variables
require(dplyr)

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
                 
#------------------------------------------------------------------#
#---------------- Figure 3.5 --------------------------------------#
#------------------------------------------------------------------#

# Model fit 
fit1 <- survfit(Surv(tgvhd, gvhd != 0) ~ 1, data = bmt)

# Plot data 
pdata1 <- data.frame(cumhaz = fit1$cumhaz, 
                     time = fit1$time) 
pdata1

# Create Figure 3.5
fig35 <- ggplot(aes(x = time, y = cumhaz), data = pdata1) + 
  geom_step(linewidth = 1) + 
  xlab("Time since bone marrow transplantation (months)") + 
  ylab("Cumulative GvHD hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 156), breaks = seq(0, 156, by = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.1)) +
  theme_general

fig35

# Save figure
ggsave("figures/h_bmtgvhdcumhaz.pdf", plot = fig35, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 3.6 --------------------------------------#
#------------------------------------------------------------------#

# Reformatting
bmt$nyreltime <- with(bmt, ifelse(tgvhd < intxrel, tgvhd, intxrel))
bmt$nyrel <- with(bmt, ifelse(tgvhd < intxrel, 0, rel))
bmt$nytrm <- with(bmt, ifelse(tgvhd < intxrel, 0, trm))

# Cumulative relapse rate without GvHD 
fit11 <- survfit(Surv(nyreltime, nyrel != 0) ~ 1, data = bmt)

# Cumulative relapse rate after GvHD
fit12 <- survfit(Surv(tgvhd, intxrel, rel != 0) ~ 1, data = subset(bmt, gvhd == 1 & tgvhd < intxrel))

# Collect plot data
pdata11 <- data.frame(cumhaz = fit11$cumhaz, 
                      time = fit11$time) 

head(pdata11)


pdata12 <- data.frame(cumhaz = fit12$cumhaz, 
                      time = fit12$time) 

head(pdata12)


# Create Figure 3.6
fig36 <- ggplot(aes(x = time, y = cumhaz), data = pdata11) + 
  geom_step(linewidth = 1) +
  geom_step(aes(x = time, y = cumhaz), data = pdata12, linewidth = 1, linetype = "dashed") + 
  xlab("Time since bone marrow transplantation (months)") + 
  ylab("Cumulative relapse hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 156), breaks = seq(0, 156, by = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 0.2), breaks = seq(0, 0.2, by = 0.1)) +
  theme_general

fig36

# Save figure
ggsave("figures/h_bmtrelcumhazgvhd.pdf", plot = fig36, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 3.7 --------------------------------------#
#------------------------------------------------------------------#


# Need to get the cumulative hazard estimates at the same time points
# All times
all_t <- sort(unique(c(unique(pdata11$time), unique(pdata12$time))))
all_t


# Evaluate each step function at these time points 
step11 <- stepfun(x = pdata11$time, y = c(0, pdata11$cumhaz))
step12 <- stepfun(x = pdata12$time, y = c(0, pdata12$cumhaz))

pdata11_a <- data.frame(time = all_t, 
                        cumhaz = step11(all_t))
pdata12_a <- data.frame(time = all_t, 
                        cumhaz = step12(all_t))

# Collect data
pdata_b <- data.frame(cumhaz1 = pdata11_a$cumhaz, 
                      cumhaz2 = pdata12_a$cumhaz, 
                      time = all_t)

# Create Figure 3.7
fig37 <- ggplot(aes(x = cumhaz1, y = cumhaz2), data = pdata_b) + 
  geom_step(linewidth = 1) +
  geom_abline(aes(intercept = 0, slope = 0.858), linewidth = 1, linetype = "dashed") + 
  xlab("Cumulative relapse hazard: no GvHD") + 
  ylab("Cumulative relapse hazard: GvHD") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 0.2), breaks = seq(0, 0.2, by = 0.05),labels = prettyZero) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 0.2), breaks = seq(0, 0.2, by = 0.05),labels = prettyZero) +
  theme_general

fig37

# Save figure
ggsave("figures/h_bmtrelgof.pdf", plot = fig37, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 3.8 --------------------------------------#
#------------------------------------------------------------------#


# Cumulative death rate without GvHD 
fit11 <- survfit(Surv(nyreltime, nytrm != 0) ~ 1, data = bmt)

# Cumulative death rate with GvHD
fit12 <- survfit(Surv(tgvhd, intxrel, trm != 0) ~ 1, data = subset(bmt, gvhd == 1 & tgvhd < intxrel))

# Collect plot data
pdata11 <- data.frame(cumhaz = fit11$cumhaz, 
                      time = fit11$time) 

pdata12 <- data.frame(cumhaz = fit12$cumhaz, 
                      time = fit12$time) 


# Need to get the cumulative hazard estimates at the same time points
# All times
all_t <- sort(unique(c(unique(pdata11$time), unique(pdata12$time))))
all_t

# Evaluate each step function at these time points 
step11 <- stepfun(x = pdata11$time, y = c(0, pdata11$cumhaz))
step12 <- stepfun(x = pdata12$time, y = c(0, pdata12$cumhaz))

pdata11_a <- data.frame(time = all_t, 
                        cumhaz = step11(all_t))
pdata12_a <- data.frame(time = all_t, 
                        cumhaz = step12(all_t))

# Collect data
pdata_b <- data.frame(cumhaz1 = pdata11_a$cumhaz, 
                      cumhaz2 = pdata12_a$cumhaz, 
                      time = all_t)

# Create Figure 3.8
fig38 <- ggplot(aes(x = cumhaz1, y = cumhaz2), data = pdata_b) + 
  geom_step(linewidth = 1) +
  geom_abline(aes(intercept = 0, slope = 3.113), 
              linewidth = 1, linetype = "dashed") + 
  xlab("Cumulative death hazard: no GvHD") + 
  ylab("Cumulative death hazard: GvHD") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 0.2), breaks = seq(0, 0.2, by = 0.05),labels = prettyZero) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2),labels = prettyZero) +
  theme_general

fig38

# Save figure
ggsave("figures/h_bmtdeadgof.pdf", plot = fig38, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 3.4 --------------------------------------#
#------------------------------------------------------------------#


# Cumulative death hazards wrt relapse; 

# Cumulative death rate no relapse or GvHD
fit11 <- survfit(Surv(intxrel, trm != 0) ~ 1, data = bmt)

# Cumulative death rate after relapse
bmt$intxsurv <- with(bmt, ifelse(intxrel == intxsurv & rel == 1, intxsurv + 0.01, intxsurv))
fit12 <- survfit(Surv(intxrel, intxsurv, dead != 0) ~ 1, data = subset(bmt, rel == 1))

# Cumulative death rate after GvHD
fit13 <- survfit(Surv(tgvhd, intxrel, trm != 0) ~ 1, data = subset(bmt, gvhd == 1))


# Collect plot data
pdata11 <- data.frame(cumhaz = fit11$cumhaz, 
                      time = fit11$time) 

pdata12 <- data.frame(cumhaz = fit12$cumhaz, 
                      time = fit12$time) 

pdata13 <- data.frame(cumhaz = fit13$cumhaz, 
                      time = fit13$time) 

tail(pdata12)

# Create Figure 3.4
fig34 <- ggplot(aes(x = time, y = cumhaz), data = pdata11) + 
  geom_step(linewidth = 1) +
  geom_step(aes(x = time, y = cumhaz), data = pdata12, linewidth = 1, linetype = "dashed") + 
  geom_step(aes(x = time, y = cumhaz), data = pdata13, linewidth = 1, linetype = "dotted") + 
  xlab("Time since bone marrow transplantation (months)") + 
  ylab("Cumulative hazard of death") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 156), breaks = seq(0, 156, by = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  theme_general

fig34

# Save figure
ggsave("figures/h_bmtcumhazdead3.pdf", plot = fig34, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Table 3.11 --------------------------------------#
#------------------------------------------------------------------#

# In SAS


#------------------------------------------------------------------#
#---------------- Table 3.13 --------------------------------------#
#------------------------------------------------------------------#


# Prepare data set for analysis - double
bmtext <- bmt %>% mutate(nygsource = bmonly, 
                         nydisease = all, 
                         age10 = age/10,
                         time = ifelse(gvhd == 1, tgvhd, intxrel), 
                         status = ifelse(gvhd == 1, 0, trm),
                         entry = 0,
                         ) 


double1 <- bmtext  %>% 
  mutate(age1 = age10, 
         age2 = 0,
         gsource1 = nygsource, 
         gsource2 = 0,
         disease1 = nydisease,
         disease2 = 0, 
         stratum = 1)

double2 <- bmtext %>% filter(gvhd == 1) %>% 
  mutate(time = intxrel, 
         status = trm,  
         entry = tgvhd,
         age1 = 0, 
         age2 = age10,
         gsource1 = 0, 
         gsource2 = nygsource,
         disease1 = 0,
         disease2 = nydisease, 
         stratum = 2)

doubledod <- as.data.frame(rbind(double1, double2))
head(doubledod)


# Row 1
r1 <- coxph(Surv(entry, time, status != 0) ~ strata(stratum) + disease1 + disease2 + age1 + age2, 
            data = doubledod, ties = "breslow")
summary(r1)

# Row 2
doubledod$type <- doubledod$stratum - 1
r2 <- coxph(Surv(entry, time, status != 0) ~ disease1 + disease2 + age1 + age2 + type, 
            data = doubledod, ties = "breslow")
summary(r2)

# Row 3
r3 <- coxph(Surv(entry, time, status != 0) ~ nydisease + age10 + type, 
            data = doubledod, ties = "breslow")
summary(r3)



#------------------------------------------------------------------#
#---------------- Figure 3.10 -------------------------------------#
#------------------------------------------------------------------#

# Extract cumulative hazard from r1 
survr1 <- basehaz(r1, center = F)
pcumhaz <- data.frame(cumhaz = survr1$hazard, 
                      time = survr1$time, 
                      strata = survr1$strata
) 

# Create Figure 3.10
fig310 <- ggplot(aes(x = time, y = cumhaz, linetype = strata), data = pcumhaz) + 
  geom_step(linewidth = 1) + 
  xlab("Time since bone marrow transplantation (months)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 156), breaks = seq(0, 156, by = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  scale_linetype_discrete("Stratum", labels = c("1", "2")) + 
  theme_general + theme(legend.position = "none")

fig310

ggsave("figures/h_bmtjointcumhaztrm.pdf", plot = fig310, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Table 3.14 --------------------------------------#
#------------------------------------------------------------------#

# Gamma
coxph(Surv(intxrel, state0 != 0) ~ bmonly + all + age + frailty(team), 
      data = bmt, ties = "breslow")
