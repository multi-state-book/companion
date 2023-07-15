#------------------------------------------------------------------#
#------- Chapter 5, R code, bmt data ------------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

# Load required packages (should be installed if not already)
require(haven)
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

# Reformatting and adding extra variables
bmt$intxsurv <- bmt$timedeath
bmt$dead <- bmt$death
bmt$intxrel <- ifelse(bmt$rel == 1, bmt$timerel, bmt$timedeath)
bmt$trm <- ifelse(bmt$rel == 0 & bmt$death == 1, 1, 0)
bmt$tgvhd <- ifelse(bmt$gvhd == 1, bmt$timegvhd, bmt$intxrel)
bmt$state0 <- bmt$rel + 2*bmt$trm
bmt$dytxanc5 <- bmt$timeanc500 * 30

#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------#
#------------------------------------------------------------------#

# Make zeros print as "0" always
library(stringr)
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}

# General theme
theme_general <- theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        legend.text = element_text(size = 26))

#------------------------------------------------------------------#
#---------------- Table 5.1 ---------------------------------------#
#------------------------------------------------------------------#

# Make landmark data at 0.5, 1, 1.5, 2, 2.5
require(dplyr)
dat05 <- bmt
dat05 <- subset(bmt, intxrel >= 0.5) %>%
  mutate(time = pmin(intxrel, 6.5),
         status = ifelse(time < 6.5, state0, 0),
         landmark = 0.5,
         entry = 0.5,
         anc = ifelse(anc500 == 1 & timeanc500 <= 0.5, 1, 0),
         gvh = ifelse(gvhd == 1 & tgvhd <= 0.5, 1, 0))
# dat05

dat10 <- bmt
dat10 <- subset(bmt, intxrel >= 1.0) %>%
  mutate(time = pmin(intxrel, 7),
         status = ifelse(time < 7, state0, 0),
         landmark = 1.0,
         entry = 1.0,
         anc = ifelse(anc500 == 1 & timeanc500 <= 1, 1, 0),
         gvh = ifelse(gvhd == 1 & tgvhd <= 1, 1, 0))

dat15 <- bmt
dat15 <- subset(bmt, intxrel >= 1.5) %>%
  mutate(time = pmin(intxrel, 7.5),
         status = ifelse(time < 7.5, state0, 0),
         landmark = 1.5,
         entry = 1.5,
         anc = ifelse(anc500 == 1 & timeanc500 <= 1.5, 1, 0),
         gvh = ifelse(gvhd == 1 & tgvhd <= 1.5, 1, 0))

dat20 <- bmt
dat20 <- subset(bmt, intxrel >= 2.0) %>%
  mutate(time = pmin(intxrel, 8),
         status = ifelse(time < 8, state0, 0),
         landmark = 2,
         entry = 2,
         anc = ifelse(anc500 == 1 & timeanc500 <= 2, 1, 0),
         gvh = ifelse(gvhd == 1 & tgvhd <= 2, 1, 0))

dat25 <- bmt
dat25 <- subset(bmt, intxrel >= 2.5) %>%
  mutate(time = pmin(intxrel, 8.5),
         status = ifelse(time < 8.5, state0, 0),
         landmark = 2.5,
         entry = 2.5,
         anc = ifelse(anc500 == 1 & timeanc500 <= 2.5, 1, 0),
         gvh = ifelse(gvhd == 1 & tgvhd <= 2.5, 1, 0))


landmark <- rbind(dat05, dat10, dat15, dat20, dat25)

head(landmark)

# Table 5.1 output
rowSums(with(landmark, table(landmark, anc > 0 | gvh > 0))) # at risk
with(landmark, table(landmark, anc)) # anc
with(landmark, table(landmark, gvh)) # gvhd


# Drop observations with stoptime = starttime
landmark <- subset(landmark, time > entry)

#------------------------------------------------------------------#
#---------------- Table 5.2 ---------------------------------------#
#------------------------------------------------------------------#

# Add time-varying covariates
landmarkw <- landmark %>%
  mutate(anc05 = (landmark == 0.5) * anc,
         gvh05 = (landmark == 0.5) * gvh,
         anc10 = (landmark == 1.0) * anc,
         gvh10 = (landmark == 1.0) * gvh,
         anc15 = (landmark == 1.5) * anc,
         gvh15 = (landmark == 1.5) * gvh,
         anc20 = (landmark == 2.0) * anc,
         gvh20 = (landmark == 2.0) * gvh,
         anc25 = (landmark == 2.5) * anc,
         gvh25 = (landmark == 2.5) * gvh
         )

head(landmarkw)

# with(landmarkw, table(landmark, status))


cox_land <- coxph(Surv(entry, time, status != 0) ~ cluster(id) + strata(landmark) +
                    anc05 + anc10 + anc15 + anc20 + anc25 +
                    gvh05 + gvh10 + gvh15 + gvh20 + gvh25 ,
                  data = landmarkw,
                  method = "breslow", timefix = FALSE, eps = 1e-9)
summary(cox_land)

#------------------------------------------------------------------#
#---------------- Figure 5.1 --------------------------------------#
#------------------------------------------------------------------#

# Extract estimated survival probabilities
surv <- survfit(cox_land,
                newdata = data.frame(anc05 = 0, anc10 = 0, anc15 = 0, anc20 = 0, anc25 = 0,
                                     gvh05 = 0, gvh10 = 0, gvh15 = 0, gvh20 = 0, gvh25 = 0))

# Order data for plotting
pdata <- data.frame(surv = surv$surv,
                    time = surv$time,
                    landmark = c(rep("0.5", surv$strata[[1]]),
                                 rep("1.0", surv$strata[[2]]),
                                 rep("1.5", surv$strata[[3]]),
                                 rep("2.0", surv$strata[[4]]),
                                 rep("2.5", surv$strata[[5]]))
)

# Add prob 1 in beginning for all
pdata2 <- as.data.frame(pdata %>% group_by(landmark) %>%
  group_modify(~ add_row(.x, surv = 1, time = 0, .before=0)))


# Create Figure 5.1 (left)
fig51_a <- ggplot(aes(x = time, y = surv, linetype = landmark), data = pdata2) +
  geom_step(size = 1) +
  scale_linetype_discrete("Landmark") +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Conditional survival probability") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 9),
                     breaks = seq(0, 9, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  theme(legend.key.size = unit(1.5, 'cm')) +
  guides(linetype = guide_legend(nrow = 1, byrow = TRUE))

fig51_a

ggsave("figures/h_bmtLM10.pdf", plot = fig51_a,
       width = 29.7, height = 21, units = "cm")


# For anc = gvh = 1
# Add LP to data
pdata2$lp <- with(pdata2,
                  ifelse(landmark == '0.5', coef(cox_land)[1] + coef(cox_land)[6],
                    ifelse(landmark == '1.0', coef(cox_land)[2] + coef(cox_land)[7],
                      ifelse(landmark == '1.5', coef(cox_land)[3] + coef(cox_land)[8],
                        ifelse(landmark == '2.0', coef(cox_land)[4] + coef(cox_land)[9],
                          ifelse(landmark == '2.5', coef(cox_land)[5] + coef(cox_land)[10],
                    0))))))

pdata2_wc <- pdata2
pdata2_wc$survlp <- with(pdata2, surv ^ exp(lp))


# Create Figure 5.1 (right)
fig51_b <- ggplot(aes(x = time, y = survlp, linetype = landmark), data = pdata2_wc) +
  geom_step(size = 1) +
  scale_linetype_discrete("Landmark") +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Conditional survival probability") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 9),
                     breaks = seq(0, 9, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  theme(legend.key.size = unit(1.5, 'cm')) +
  guides(linetype = guide_legend(nrow = 1, byrow = TRUE))

fig51_b

ggsave("figures/h_bmtLM11.pdf", plot = fig51_b,
       width = 29.7, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Table 5.3 ---------------------------------------#
#------------------------------------------------------------------#
## COLUMN 1 
# Add extra variables
landmarkw2 <- landmark %>%
  mutate(anctime  = anc * (landmark - 0.5)/2,
         anctime2 = anc * ((landmark - 0.5)/2)^2,
         gvhtime  = gvh * (landmark - 0.5)/2,
         gvhtime2 = gvh * ((landmark - 0.5)/2)^2
        )

landmarkw2

# Fit model
cox_land2 <- coxph(Surv(entry, time, status != 0) ~ cluster(id) + strata(landmark) +
                    anc + anctime + anctime2 +
                    gvh + gvhtime + gvhtime2,
                  data = landmarkw2,
                  method = "breslow", timefix = FALSE, eps = 1e-9)
summary(cox_land2)


## COLUMN 2


# Add extra variables
landmarkw3 <- landmarkw2 %>%
  mutate(strtime = (landmark - 0.5)/2,
         strtime2 = strtime^2
  )

# Fit model
cox_land3 <- coxph(Surv(entry, time, status != 0) ~ cluster(id) +
                     anc + anctime + anctime2 +
                     gvh + gvhtime + gvhtime2 +
                     strtime + strtime2,
                   data = landmarkw3,
                   method = "breslow")
summary(cox_land3)


#------------------------------------------------------------------#
#---------------- Figure 5.2 --------------------------------------#
#------------------------------------------------------------------#


# Extract estimated survival probabilities
surv <- survfit(cox_land2,
                newdata = data.frame(anc = 0, anctime = 0, anctime2 = 0,
                                     gvh = 0, gvhtime = 0, gvhtime2 = 0))

# Order data for plotting
pdata <- data.frame(surv = surv$surv,
                    time = surv$time,
                    landmark = c(rep("0.5", surv$strata[[1]]),
                                 rep("1.0", surv$strata[[2]]),
                                 rep("1.5", surv$strata[[3]]),
                                 rep("2.0", surv$strata[[4]]),
                                 rep("2.5", surv$strata[[5]]))
)

# Add prob 1 in beginning for all
pdata2 <- as.data.frame(pdata %>% group_by(landmark) %>%
                          group_modify(~ add_row(.x, surv = 1, time = 0, .before=0)))


# Create Figure 5.2 (left)
fig52_a <- ggplot(aes(x = time, y = surv, linetype = landmark), data = pdata2) +
  geom_step(size = 1) +
  scale_linetype_discrete("Landmark") +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Conditional survival probability") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 9),
                     breaks = seq(0, 9, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  theme(legend.key.size = unit(1.5, 'cm')) +
  guides(linetype = guide_legend(nrow = 1, byrow = TRUE))

fig52_a

ggsave("figures/h_bmtLM20.pdf", plot = fig52_a,
       width = 29.7, height = 21, units = "cm")


# For anc = gvh = 1
# Add LP to data
pdata2$landmarknum <- as.numeric(pdata2$landmark)
pdata2$lp <- with(pdata2,
                  coef(cox_land2)[1] +
                    coef(cox_land2)[2] * (landmarknum - 0.5) / 2 +
                     coef(cox_land2)[3] * ((landmarknum - 0.5) / 2)^2 +
                      coef(cox_land2)[4] +
                       coef(cox_land2)[5] * (landmarknum - 0.5) / 2 +
                        coef(cox_land2)[6] * ((landmarknum - 0.5) / 2)^2
                        )



pdata2_wc <- pdata2
pdata2_wc$survlp <- with(pdata2, surv ^ exp(lp))


# Create Figure 5.2 (right)
fig52_b <- ggplot(aes(x = time, y = survlp, linetype = landmark), data = pdata2_wc) +
  geom_step(size = 1) +
  scale_linetype_discrete("Landmark") +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Conditional survival probability") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 9),
                     breaks = seq(0, 9, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  theme(legend.key.size = unit(1.5, 'cm')) +
  guides(linetype = guide_legend(nrow = 1, byrow = TRUE))

fig52_b

ggsave("figures/h_bmtLM21.pdf", plot = fig52_b,
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Figure 5.3 --------------------------------------#
#------------------------------------------------------------------#



# Extract estimated survival probabilities
surv <- survfit(cox_land3,
                newdata = data.frame(anc = 0, anctime = 0, anctime2 = 0,
                                     gvh = 0, gvhtime = 0, gvhtime2 = 0 ,
                                     strtime = 0, strtime2 = 0))

# Order data for plotting
pdata <- data.frame(surv = surv$surv,
                    time = surv$time)

require(dplyr)
base305 <- pdata %>% filter(time < 6.5) %>% 
  mutate(landmark = 0.5,
         km0 = surv, 
         lpt = 0,
         lpz = coef(cox_land3)[1] + coef(cox_land3)[4],
         km1 = km0^exp(lpz))


base310 <- pdata %>% filter(time > 1 & time < 7) %>%  
  mutate(landmark = 1,
         lpt = coef(cox_land3)[7] * 0.5 / 2 + coef(cox_land3)[8] * (0.5 / 2)^2,
         km0 = (surv / 0.9784098)^exp(lpt), 
         lpz = coef(cox_land3)[1] + coef(cox_land3)[4] + 
           coef(cox_land3)[2] * 0.5 / 2 + 
           coef(cox_land3)[5] * 0.5 / 2 +
           coef(cox_land3)[3] * (0.5 / 2)^2 + 
           coef(cox_land3)[6] * (0.5 / 2)^2 ,
         km1 = km0^exp(lpz))


base315 <- pdata %>% filter(time > 1.5 & time < 7.5) %>%  
  mutate(landmark = 1.5,
         lpt = coef(cox_land3)[7] * 1 / 2 + coef(cox_land3)[8] * (1 / 2)^2,
         km0 = (surv / 0.9554490)^exp(lpt), 
         lpz = coef(cox_land3)[1] + coef(cox_land3)[4] + 
           coef(cox_land3)[2] * 1 / 2 + 
           coef(cox_land3)[5] * 1 / 2 +
           coef(cox_land3)[3] * (1 / 2)^2 + 
           coef(cox_land3)[6] * (1 / 2)^2 ,
         km1 = km0^exp(lpz))

base320 <- pdata %>% filter(time > 2 & time < 8) %>%  
  mutate(landmark = 2,
         lpt = coef(cox_land3)[7] * 1.5 / 2 + coef(cox_land3)[8] * (1.5 / 2)^2,
         km0 = (surv / 0.9395795)^exp(lpt), 
         lpz = coef(cox_land3)[1] + coef(cox_land3)[4] + 
           coef(cox_land3)[2] * 1.5 / 2 + 
           coef(cox_land3)[5] * 1.5 / 2 +
           coef(cox_land3)[3] * (1.5 / 2)^2 + 
           coef(cox_land3)[6] * (1.5 / 2)^2 ,
         km1 = km0^exp(lpz))

base325 <- pdata %>% filter(time > 2.5 & time < 8.5) %>%  
  mutate(landmark = 2.5,
         lpt = coef(cox_land3)[7] * 2 / 2 + coef(cox_land3)[8] * (2 / 2)^2,
         km0 = (surv / 0.9126288)^exp(lpt), 
         lpz = coef(cox_land3)[1] + coef(cox_land3)[4] + 
           coef(cox_land3)[2] * 2 / 2 + 
           coef(cox_land3)[5] * 2 / 2 +
           coef(cox_land3)[3] * (2 / 2)^2 + 
           coef(cox_land3)[6] * (2 / 2)^2 ,
         km1 = km0^exp(lpz))
base3slut <- as.data.frame(rbind(base305, base310, base315, base320, base325))


# Add prob 1 in beginning for all
pdata3 <- as.data.frame(base3slut %>% group_by(landmark) %>%
                           group_modify(~ add_row(.x, km0 = 1, km1  = 1, time = 0, .before=0)))

# Create Figure 5.3 (left)
fig53_a <- ggplot(aes(x = time, y = km0, linetype = as.factor(landmark)), data = pdata3) +
  geom_step(size = 1) +
  scale_linetype_discrete("Landmark",labels=c("0.5","1.0","1.5", "2.0","2.5")) +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Conditional survival probability") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 9),
                     breaks = seq(0, 9, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
  theme_general +
  theme(legend.key.size = unit(1.5, 'cm')) + 
  guides(linetype = guide_legend(nrow = 1, byrow = TRUE))

fig53_a

ggsave("figures/h_bmtLM30.pdf", plot = fig53_a,
       width = 29.7, height = 21, units = "cm")


# Create Figure 5.3 (right)
fig53_b <- ggplot(aes(x = time, y = km1, linetype = as.factor(landmark)), data = pdata3) +
  geom_step(size = 1) +
  scale_linetype_discrete("Landmark",labels=c("0.5","1.0","1.5", "2.0","2.5")) +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Conditional survival probability") +
  scale_x_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 9),
                     breaks = seq(0, 9, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.05)),
                     limits = c(0, 1.0),
                     breaks = seq(0, 1.0, 0.1)) +
   theme_general +
   theme(legend.key.size = unit(1.5, 'cm'))  +
  guides(linetype = guide_legend(nrow = 1, byrow = TRUE))

fig53_b


ggsave("figures/h_bmtLM31.pdf", plot = fig53_b,
       width = 29.7, height = 21, units = "cm")



#------------------------------------------------------------------#
#---------------- Table 5.5 ---------------------------------------#
#------------------------------------------------------------------#
# and
#------------------------------------------------------------------#
#---------------- Figure 5.9+5.10 ---------------------------------#
#------------------------------------------------------------------#
#
# See Eva's great program bmt230217.R
