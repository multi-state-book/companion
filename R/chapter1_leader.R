#------------------------------------------------------------------#
#------- Chapter 1, R code, LEADER data ---------------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
setwd("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB")

# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

# Read data 
leader <- read.csv("data/leader_mi_3p.csv")

head(leader)

# One data set per endpoint type
leader_mi <- subset(leader, type == "recurrent_mi")
leader_3p <- subset(leader, type == "recurrent_comb")

#------------------------------------------------------------------#
#---------------- Table 1.3 ---------------------------------------#
#------------------------------------------------------------------#

# For recurrent MI
with(subset(leader_mi, status == 1), table(eventno, treat))

colSums(with(subset(leader_mi, status == 1), table(eventno, treat)))

with(leader_mi, table(status, eventno, treat))

colSums(with(leader_mi, table(status, eventno, treat)))

# For recurrent 3p MACE
with(subset(leader_3p, status == 1), table(eventno, treat))

colSums(with(subset(leader_3p, status == 1), table(eventno, treat)))

with(leader_3p, table(status, eventno, treat))

# colSums(with(leader_3p, table(status, eventno, treat)))
