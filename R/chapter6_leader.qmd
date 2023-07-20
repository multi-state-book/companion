#------------------------------------------------------------------#
#------- Chapter 6, R code, LEADER data ---------------------------#
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
#---------------- Figure 6.11 -------------------------------------#
#------------------------------------------------------------------#

require(mets)
require(ggpubr)

np_ests <- function(endpointdat){
  # Fit NAa
  NAa_fit <- survfit(Surv(start, stop, status == 1) ~ treat, 
                     data = endpointdat, id = id,
                     ctype = 1)
  # Fit KM
  KM_fit <- survfit(Surv(start, stop, status == 2) ~ treat,
                    data = endpointdat, id = id,
                    ctype = 1)

    # Adjust hat(mu)
  mu_adj <- c(cumsum(KM_fit$surv[1:KM_fit$strata[[1]]] * c(0, diff(NAa_fit$cumhaz[1:NAa_fit$strata[[1]]]))),
              cumsum(KM_fit$surv[(KM_fit$strata[[1]] + 1):(KM_fit$strata[[1]] + KM_fit$strata[[2]])] * 
                       c(0, diff(NAa_fit$cumhaz[(NAa_fit$strata[[1]] + 1):(NAa_fit$strata[[1]] + NAa_fit$strata[[2]])])))
  )
  
  dat_adj <- data.frame(mu = mu_adj, 
                        time = NAa_fit$time, 
                        treat = c(rep("Liraglutide", NAa_fit$strata[[1]]), rep("Placebo", NAa_fit$strata[[2]])),
                        type = rep("Mortality taken into account (GL)", length(NAa_fit$time)))
  
  pdat <- rbind(dat_adj)
  
  mu <- ggplot(aes(x = time  * 1 / (365.25 / 12), y = mu), data = pdat) + 
    geom_step(aes(linetype = treat), size = 1.05) + 
    xlab("Time since randomization (months)") + 
    ylab(expression(hat(mu)* "(t)")) + 
    scale_linetype_discrete("Treatment") + 
    theme_general + 
    theme( legend.title=element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(2,"line"),
      legend.direction = "horizontal", 
      legend.box = "horizontal",
      legend.key.width = unit(1.5, "cm")) + 
    scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                       limits = c(0, 65), breaks = seq(0, 65, by = 12)) + 
    scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                       limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.02))
  
  mu
  
  ## S
  
  dat_S <- data.frame(surv = KM_fit$surv, 
                          time = KM_fit$time, 
                          treat = c(rep("Liraglutide", KM_fit$strata[[1]]), rep("Placebo", KM_fit$strata[[2]])))
  
  surv <- ggplot(aes(x = time  * 1 / (365.25 / 12), y = surv), data = dat_S) + 
    geom_step(aes(linetype = treat), size = 1.05) + 
    xlab("Time since randomization (months)") + 
    ylab(expression(hat(S)* "(t)")) + 
    scale_linetype_discrete("Treatment") + 
    theme_general + 
    theme( 
      legend.direction = "horizontal", 
      legend.box = "horizontal",
      legend.key.width = unit(1.5, "cm")) + 
    scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                       limits = c(0, 65), breaks = seq(0, 65, by = 12)) + 
    scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                       limits = c(0, 1), breaks = seq(0, 1, by = 0.1))
  
  surv
  
  both <- ggarrange(mu, surv, common.legend = T, legend = "bottom")
  both
  
}

fig611 <- np_ests(endpointdat = leader_mi)
fig611

setwd("c:/Users/hnrv/OneDrive - Novo Nordisk/Book/Julie programs copied 19Jun23/R programs")
ggsave("figures/h_LEADERest.pdf", plot = fig611, 
       width = 29.7, height = 15, units = "cm")
