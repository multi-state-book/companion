#------------------------------------------------------------------#
#------- Chapter 1, R code, Holter data ---------------------------#
#------------------------------------------------------------------#

# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

# Read data
holter <- read.csv("data/cphholter.csv")

head(holter)

#------------------------------------------------------------------#
#---------------- Table 1.5 ---------------------------------------#
#------------------------------------------------------------------#

with(holter, table(esvea))


with(holter, table(afib, esvea))

with(holter, table(stroke, esvea))

with(holter, table(death, esvea))


with(holter, table(afib, stroke, death, esvea))

with(subset(holter, timeafib < timestroke), table(afib, stroke, death, esvea))

with(subset(holter, timestroke < timeafib), table(afib, stroke, death, esvea))
