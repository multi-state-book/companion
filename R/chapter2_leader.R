#------------------------------------------------------------------#
#------- Chapter 2, R code, LEADER data ---------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
setwd("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB")

# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

leader <- read.csv("data/leader_mi_3p.csv")

# One data set per endpoint type
leader_mi <- subset(leader, type == "recurrent_mi")
leader_3p <- subset(leader, type == "recurrent_comb_str_mi_cvdth")


#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------# 
#------------------------------------------------------------------#

# General theme
theme_general <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 

theme_general


#------------------------------------------------------------------#
#---------------- Table 2.14 --------------------------------------#
#------------------------------------------------------------------#

# Cox model first events
coxmod <- coxph(Surv(start, stop, status == 1) ~ treat,
                method = "breslow",
                data = subset(leader_mi, eventno == 1))
summary(coxmod)

# AG Cox type
agcoxmod <- coxph(Surv(start, stop, status == 1) ~ treat,
                   method = "breslow",
                   data = leader_mi)
summary(agcoxmod)


# PWP 2. event
pwp_2_mod <- coxph(Surv(start, stop, status == 1) ~ treat,
                     method = "breslow", subset = (eventno == 2),
                     data = leader_mi)
summary(pwp_2_mod)

# PWP 3. event
pwp_3_mod <- coxph(Surv(start, stop, status == 1) ~ treat,
                   method = "breslow", subset = (eventno == 3),
                   data = leader_mi)
summary(pwp_3_mod)

# PWP 4. event
pwp_4_mod <- coxph(Surv(start, stop, status == 1) ~ treat,
                   method = "breslow", subset = (eventno == 4),
                   data = leader_mi)
summary(pwp_4_mod)

# PWP 5. event
pwp_5_mod <- coxph(Surv(start, stop, status == 1) ~ treat,
                   method = "breslow", subset = (eventno == 5),
                   data = leader_mi)
summary(pwp_5_mod)


# PWP all
pwpmod_a <- coxph(Surv(start, stop, status == 1) ~ treat + strata(eventno),
                   method = "breslow",
                   data = leader_mi)
summary(pwpmod_a)


# AG piecewise constant
# Calculating cuts -> 5
alltimes <- seq(0,max(leader_mi$stop),length=99)

FunctionIntervalM <- function(a,b) {
  seq(from=min(a), to = max(a), by = (max(a)-min(a))/b)
}

cuts <- FunctionIntervalM(a = alltimes, b = 5)
cuts

# AG model, piecewise constant hazards
cut_data <- survSplit(Surv(start, stop, status == 1) ~ .,
                      leader_mi,
                      cut = cuts[-1],
                      episode = "timegroup")

ag_piecewise_constant <- coxph(Surv(start, stop, event) ~
                                 treat + strata(timegroup),
                               data = cut_data)

summary(ag_piecewise_constant)

#------------------------------------------------------------------#
#---------------- Figure 2.14 -------------------------------------#
#------------------------------------------------------------------#

# Nelson-Aalen estimate of the cumulative hazard of MI
fit <- survfit(Surv(start, stop, status == 1) ~ treat,
                 data = leader_mi, id = id,
                 ctype = 1)

plotdata <- data.frame(cumhaz = fit$cumhaz,
                        time = fit$time,
                        treat = c(rep("Placebo", fit$strata[[1]]),
                                  rep("Liraglutide", fit$strata[[2]])))

# Create figure 2.14
fig214 <- ggplot(aes(x = time * 1 / (365.25 / 12), y = cumhaz), data = plotdata) +
  geom_step(aes(linetype = treat), size = 1) +
  scale_linetype_discrete("Treatment") +
  ylab("Cumulative hazard") +
  xlab("Time since randomization (months)") +
  theme_bw() +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)),
                     limits = c(0, 65), breaks = seq(0, 60, by = 12)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)),
                     limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.02)) +
  theme_general

fig214
fig214<-fig214+theme(legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"))
fig214

setwd("c:/Users/hnrv/OneDrive - Novo Nordisk/Book/Julie programs copied 19Jun23/R programs")
ggsave("figures/h_LEADER-NAa.pdf", plot = fig214, 
       width = 29.7, height = 21, units = "cm")


