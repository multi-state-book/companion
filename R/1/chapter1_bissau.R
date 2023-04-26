#------------------------------------------------------------------#
#------- Chapter 1, R code, Bissau data ---------------------------#
#------------------------------------------------------------------#

# Load required packages (should be installed if not already)
require(survival)
require(ggplot2)
require(timereg)
require(gridExtra)

bissau <- read.csv("data/bissau.csv")
head(bissau)

bissau <- data.frame(bissau)

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
#---------------- Table 1.1 ---------------------------------------#
#------------------------------------------------------------------#

with(bissau, table(bcg, dead))

rowSums(with(bissau, table(bcg, dead)))
colSums(with(bissau, table(bcg, dead)))


with(bissau, table(bcg, dead)) / rowSums(with(bissau, table(bcg, dead)))


colSums(with(bissau, table(bcg, dead))) / (sum(colSums(with(bissau, table(bcg, dead)))))
