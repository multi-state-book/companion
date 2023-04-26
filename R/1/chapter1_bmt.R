#------------------------------------------------------------------#
#------- Chapter 1, R code, bmt data ------------------------------#
#------------------------------------------------------------------#

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
#---------------- Table 1.4 ---------------------------------------#
#------------------------------------------------------------------#

with(bmt, table(rel))
with(bmt, table(rel)) / sum(with(bmt, table(rel)))

with(bmt, table(death))
with(bmt, table(death)) / sum(with(bmt, table(death)))

with(bmt, table(rel,death))
with(bmt, table(rel,death)) / rowSums(with(bmt, table(rel,death)))


with(bmt, table(gvhd))
with(bmt, table(gvhd)) / sum(with(bmt, table(gvhd)))

with(bmt, table(gvhd, rel))
with(bmt, table(gvhd, rel)) / rowSums(with(bmt, table(gvhd, rel)))

with(bmt, table(gvhd, death))
with(bmt, table(gvhd, death)) / rowSums(with(bmt, table(gvhd, death)))


