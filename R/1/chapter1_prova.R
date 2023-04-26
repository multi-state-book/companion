#------------------------------------------------------------------#
#------- Chapter 1, R code, PROVA data ----------------------------#
#------------------------------------------------------------------#

# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

prova <- read.csv("data/prova.csv", na.strings = c("."))
head(prova)

prova <- data.frame(prova)

str(prova)
summary(prova)

#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------#
#------------------------------------------------------------------#

theme_general <- theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

theme_general

#------------------------------------------------------------------#
#---------------- Table 1.2 ---------------------------------------#
#------------------------------------------------------------------#

prova$beh <- with(prova, scle*2 + beta)

with(prova, table(beh))

# Bleeding
with(prova, table(beh, bleed))

# Death
with(prova, table(beh, death))

# Death w/o bleed
with(subset(prova, bleed == 0), table(beh, death))

# Death after bleed
with(subset(prova, bleed == 1), table(beh, death))


