# R options set globally
options(width = 60)

# SAS path
saspath<-"C:\\Program Files\\SASHome\\SASFoundation\\9.4\\sas.exe"

# chunk options set globally
knitr::opts_chunk$set(comment = NA, collapse = FALSE)

#
require(survival)
require(tidyverse)
require(timereg)
require(gridExtra)

# Plot setup
theme_general <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 15), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15)) 

