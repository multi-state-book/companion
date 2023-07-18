#------------------------------------------------------------------#
#------- Chapter 3, R code, LEADER data ---------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
setwd("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB")

# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)
require(mets)


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
#---------------- Table 3.15 --------------------------------------#
#------------------------------------------------------------------#

# Cox model, frailty
coxfrail <- coxph(Surv(start, stop, status == 1) ~ treat + frailty(id), 
                  method = "breslow", 
                  data = leader_mi)
summary(coxfrail)

# PieceWise constant, frailty
# Make cuts 
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

pwch_frail <- coxph(Surv(start, stop, event) ~ 
                    treat + strata(timegroup) + frailty(id), 
                               data = cut_data)

summary(pwch_frail)

# Joint frailty model, piecewise constant hazards
require(frailtypack)
leader_mi$death <- ifelse(leader_mi$status == 2, 1, 0)
jointfrail_pc_eq_mi <- frailtyPenal(Surv(start, stop, status == 1) ~
                                      cluster(id) + treat + terminal(death),
                                    formula.terminalEvent = ~ treat,
                                    data = leader_mi,
                                    hazard = "Piecewise-equi", nb.int = c(5, 5),
                                    recurrentAG = TRUE)
jointfrail_pc_eq_mi
summary(jointfrail_pc_eq_mi)
