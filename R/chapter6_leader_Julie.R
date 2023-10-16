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


# # MU 
# xr <- phreg(Surv(start, stop, status == 1) ~ strata(trt) + cluster(id), data = leader_mi)
# xd <- phreg(Surv(start, stop, status == 2) ~ strata(trt) + cluster(id), data = leader_mi)
# out <- recurrentMarginal(xr, xd)
# #out
# 
# pout <- data.frame(time = out$cumhaz[,1], 
#                    mu = out$cumhaz[,2],
#                    rx = as.factor(out$strata))
# 
# # extend to max obs 
# tmax <- max(leader_mi$stop)
# leader_mi[leader_mi$stop == tmax, ]
# 
# maxmu1 <- max(pout[pout$rx == 1,]$mu)
# maxmu0 <- max(pout[pout$rx == 0,]$mu)
# pout[nrow(pout) + 1,] <- c(tmax, maxmu1, 1)
# pout[nrow(pout) + 1,] <- c(tmax, maxmu0, 0)
# 
# mu_np <- ggplot(aes(x = time * 12/365.25, y = mu, linetype = rx), data = pout) + geom_step(size = 1) + 
#   ylab(expression(hat(mu)(t))) + 
#   xlab("Time since randomisation (months)") + 
#   scale_x_continuous(breaks = seq(0, 63, by = 20), limits = c(0, 63)) + 
#   theme_bw() + scale_linetype_discrete("Treatment", labels = c("Placebo", "Liraglutide")) + 
#   theme(legend.position="bottom",
#         legend.box="vertical")
# mu_np
# 
# # S - instead Lambda0d
# 
# pout_lamd <- data.frame(time = xd$cumhaz[,1], 
#                         cumhaz = xd$surv[,2],
#                         rx = as.factor(xd$strata.jumps))
# 
# LamD_np <- ggplot(aes(x = time* 12/365.25, y = cumhaz, color = rx), data = pout_lamd) + 
#   geom_step(size = 1) + 
#   ylab(expression(hat(Lambda)^D*(t))) + 
#   xlab("Time since randomisation (months)") + 
#   scale_x_continuous(breaks = seq(0, 63, by = 20), limits = c(0, 63)) + 
#   theme_bw() + scale_color_discrete("Treatment", labels = c("Placebo", "Liraglutide")) + 
#   theme(legend.position="bottom",
#         legend.box="vertical")
# LamD_np
# 
# leader_np <- ggarrange(mu_np, LamD_np, ncol = 2, common.legend = T, legend = "bottom")
# leader_np


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
    theme( 
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

setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

ggsave("figures/j_LEADERest.pdf", plot = fig611, 
       width = 29.7, height = 15, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 6.8 --------------------------------------#
#------------------------------------------------------------------#


## while it is on GitHub
# require(devtools)
# devtools::install_github("JulieKFurberg/recurrentpseudo", force = TRUE)
require(recurrentpseudo)

leader_pseudo <- pseudo.twodim(tstart = leader_mi$start,
                               tstop = leader_mi$stop,
                               status = leader_mi$status,
                               id = leader_mi$id,
                               covar_names = "treat",
                               tk = c(20, 30, 40) * (365.25 / 12),
                               data = leader_mi)

fit_leader <- pseudo.geefit(pseudodata = leader_pseudo,
                            covar_names = c("treat"))

fit_leader

# Treatment differences
xi_diff_2d <- as.matrix(c(fit_leader$xi[4],
                          fit_leader$xi[8]), ncol = 1)

mslabels <- c("treat, mu", "treat, surv")
rownames(xi_diff_2d) <- mslabels
colnames(xi_diff_2d) <- ""
xi_diff_2d

# Variance matrix for differences
sigma_diff_2d <- matrix(c(fit_leader$sigma[4,4],
                          fit_leader$sigma[4,8],
                          fit_leader$sigma[4,8],
                          fit_leader$sigma[8,8]),
                        ncol = 2, nrow = 2,
                        byrow = T)

rownames(sigma_diff_2d) <- colnames(sigma_diff_2d) <- mslabels
sigma_diff_2d

# Correlation matrix
cov2cor(sigma_diff_2d)
