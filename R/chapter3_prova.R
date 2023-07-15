#------------------------------------------------------------------#
#------- Chapter 3, R code, PROVA data ----------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

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
require(dplyr)

prova <- read.csv("data/prova.csv", na.strings = c("."))
head(prova)

prova <- data.frame(prova)

str(prova)
summary(prova)

# Treat
prova$beh <- with(prova, as.factor(scle + beta*2))
prova$log2bili <- with(prova, log2(bili))

# Extra variables
provany <- prova

provany$btime <- ifelse(provany$bleed == 1, provany$timebleed, provany$timedeath)
provany$d0time <- ifelse(provany$bleed == 1, provany$timebleed, provany$timedeath)
provany$dead0 <- ifelse(provany$bleed == 1, 0, provany$death)
provany$outof0 <- ifelse(provany$bleed == 1, 1, provany$death)

provany$bdtime <- ifelse(provany$bleed == 1, provany$timedeath, NA)
provany$deadb <- ifelse(provany$bleed == 1, provany$death, NA)
provany$wait <- ifelse(provany$bleed == 1, provany$bdtime - provany$timebleed, NA)

#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------# 
#------------------------------------------------------------------#

theme_general <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 


#------------------------------------------------------------------#
#---------------- Table 3.3 ---------------------------------------#
#------------------------------------------------------------------#


## Column 1
# Variceal bleeding
coxph(Surv(btime, bleed) ~ beh, data = provany, ties = "breslow")

# Death without bleeding
coxph(Surv(d0time, dead0) ~ beh, data = provany, ties = "breslow")


## Column 2
# Variceal bleeding
coxph(Surv(btime, bleed) ~ beh + sex + coag + log2bili + factor(varsize), data = provany, ties = "breslow")

# Death without bleeding
coxph(Surv(d0time, dead0) ~ beh + sex + coag + log2bili + factor(varsize), data = provany, ties = "breslow")



#------------------------------------------------------------------#
#---------------- Table 3.4 ---------------------------------------#
#------------------------------------------------------------------#

# Composite
coxph(Surv(btime, outof0) ~ beh + sex + coag + log2bili + factor(varsize), data = provany, ties = "breslow")


#------------------------------------------------------------------#
#---------------- Table 3.7 ---------------------------------------#
#------------------------------------------------------------------#

## Column 1

# Time since randomisation
fitrand <- coxph(Surv(btime, bdtime, deadb != 0) ~ beh + sex + log2bili, 
                 data = provany, ties = "breslow")
summary(fitrand)


# Duration
fitdur <- coxph(Surv(wait, deadb != 0) ~ beh + sex + log2bili, data = provany, ties = "breslow")
summary(fitdur)

## Column 2
# -> In SAS


#------------------------------------------------------------------#
#---------------- Figure 3.2 --------------------------------------#
#------------------------------------------------------------------#

# Extract cumulative baseline hazard
coxcumhaz <- survfit(fitrand, 
                     newdata = data.frame(sex = 0, 
                                          beh = "0", 
                                          log2bili = 0))

# Collect data for plot
coxdata <- data.frame(cumhaz = append(0,coxcumhaz$cumhaz), 
                      time = append(0,coxcumhaz$time), 
                      type = rep("Breslow estimate", 1+length(coxcumhaz$time)))

# Create Figure 3.2
fig32 <- ggplot(aes(x = time / 365.25, y = cumhaz), data = coxdata) + 
  geom_step(linewidth = 1) + 
  xlab("Time since randomization (years)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  theme_general

fig32

ggsave("figures/h_PROVA12time.pdf", plot = fig32, 
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 3.3 --------------------------------------#
#------------------------------------------------------------------#


# Extract cumulative baseline hazard
coxcumhaz <- survfit(fitdur, 
                     newdata = data.frame(sex = 0, 
                                          beh = "0", 
                                          log2bili = 0))

# Collect data for plot
coxdata <- data.frame(cumhaz = append(0,coxcumhaz$cumhaz), 
                      time = append(0,coxcumhaz$time), 
                      type = rep("Breslow estimate", 1+length(coxcumhaz$time)))

# Create Figure 3.3
fig33 <- ggplot(aes(x = time / 365.25, y = cumhaz), data = coxdata) + 
  geom_step(linewidth = 1) + 
  xlab("Duration (years)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), ) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)),labels = prettyZero) +
  theme_general

fig33

ggsave("figures/h_PROVA12dur.pdf", plot = fig33, 
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Table 3.8 ---------------------------------------#
#------------------------------------------------------------------#

# All bleeds
provasplit11 <- subset(provany, bleed == 1)

# Split by duration
provasplit1 <- survSplit(Surv(wait, deadb != 0) ~ ., data = provasplit11, 
                           cut = c(5, 10),
                           episode = "wint")

provasplit1$start <- with(provasplit1, 
                            btime + ifelse(wint == 1, 0, ifelse(wint == 2, 5 , 10))
                            )

provasplit1$stop <- with(provasplit1, btime + wait)

provasplit1$risktime <- with(provasplit1, wait - tstart)

provasplit1$logrisktime <- log(provasplit1$risktime)
provasplit1$fail <- provasplit1$event

head(provasplit1,10)

# Split by time since rand (t)
provasplit2 <- survSplit(Surv(start, stop, fail) ~ ., data = provasplit1, 
                           cut = c((1) * 365.25, (2) * 365.25),
                           episode = "tint")

provasplit2$risktime2 <- with(provasplit2, stop - start)
provasplit2$risktimeys2 <- provasplit2$risktime2 / 365.25
  
provasplit2$logrisktime2 <- log(provasplit2$risktime2)
provasplit2$fail2 <- provasplit2$fail


head(provasplit2,10)


# Summarize the data, Table 3.8 output
aggregate(cbind(fail2, risktimeys2) ~ tint + wint, provasplit2, 
          FUN = function(x) c(count = length(x),
                              #mean = mean(x), 
                              sum = sum(x)
                              )
)


#------------------------------------------------------------------#
#---------------- Table 3.9 ---------------------------------------#
#------------------------------------------------------------------#

# Column 1
summary(glm(fail ~ offset(log(risktime)) + beh + relevel(as.factor(wint), ref = "3") + sex + log2bili, 
                    data = provasplit1, 
                    family = poisson)
)


# Column 2
summary(glm(fail2 ~ offset(log(risktime2)) + beh + relevel(as.factor(tint), ref = "3") + sex + log2bili, 
            data = provasplit2, 
            family = poisson)
)


# Column 3
summary(glm(fail2 ~ offset(log(risktime2)) + beh + relevel(as.factor(tint), ref = "3") + 
            relevel(as.factor(wint), ref = "3") + sex + log2bili, 
            data = provasplit2, 
            family = poisson)
)



# Interaction model, in-text

summary(glm(fail2 ~ offset(log(risktime2)) + beh + as.factor(tint) * as.factor(wint) + sex + log2bili, 
            data = provasplit2, 
            family = poisson)
)



#------------------------------------------------------------------#
#---------------- Table 3.12 --------------------------------------#
#------------------------------------------------------------------#

# Prepare data set for analysis - double
double1 <- provany %>% mutate(time = d0time, 
                              status = dead0, 
                              entrytime = 0, 
                              sex1 = sex, 
                              sex2 = 0, 
                              age1 = age, 
                              age2 = 0, 
                              bili1 = log2bili, 
                              bili2 = log2bili * 0,
                              stratum = 1)

double2 <- provany %>% filter(bleed == 1) %>% 
                       mutate(time = bdtime, 
                              status = deadb, 
                              entrytime = btime, 
                              sex1 = 0, 
                              sex2 = sex, 
                              age1 = 0, 
                              age2 = age, 
                              bili1 = log2bili * 0, 
                              bili2 = log2bili,
                              stratum = 2)
double <- as.data.frame(rbind(double1, double2))
head(double)


# Row 1
r1 <- coxph(Surv(entrytime, time, status != 0) ~ strata(stratum) + sex1 + sex2 + bili1 + bili2, 
            data = double, ties = "breslow")
summary(r1)

# Row 2
r2 <- coxph(Surv(entrytime, time, status != 0) ~ strata(stratum) + sex + bili1 + bili2, 
            data = double, ties = "breslow")
summary(r2)

# Row 3
r3 <- coxph(Surv(entrytime, time, status != 0) ~ strata(stratum) + sex + bili1, 
            data = double, ties = "breslow")
summary(r3)



#------------------------------------------------------------------#
#---------------- Figure 3.9 --------------------------------------#
#------------------------------------------------------------------#

# Extract cumulative hazard from r1 
survr1 <- basehaz(r1, center = F)

#pcumhaz <- data.frame(cumhaz = survr1$hazard, 
#                      time = survr1$time, 
#                      strata = survr1$strata
#) 

pcumhaz <- data.frame(
  cumhaz = c(survr1$hazard[survr1$strata=="stratum=1"],0,survr1$hazard[survr1$strata=="stratum=2"]), 
  time   = c(survr1$time[survr1$strata=="stratum=1"],0,survr1$time[survr1$strata=="stratum=2"]), 
  strata = c(survr1$strata[survr1$strata=="stratum=1"],"2",survr1$strata[survr1$strata=="stratum=2"])
) 

# Create Figure 3.9
fig39 <- ggplot(aes(x = time / 365.25, y = cumhaz, linetype = strata), data = pcumhaz) + 
  geom_step(linewidth = 1) + 
  xlab("Time since randomization (years)") + 
  ylab("Cumulative hazard") + 
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  scale_linetype_discrete("Stratum", labels = c("1", "2")) + 
  theme_general + theme(legend.position = "none")

fig39

ggsave("figures/h_provamortalities.pdf", plot = fig39, 
       width = 29.7, height = 21, units = "cm")



