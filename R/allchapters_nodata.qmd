#------------------------------------------------------------------#
#------- All chapters, R code, non-data plots ---------------------#
#------------------------------------------------------------------#

# Set working directory - location of data
#setwd("C:/Users/JUKF/University of Copenhagen/UCPH_MSB - General/Julie programs")

# Packages
require(survival)
require(ggplot2)


#------------------------------------------------------------------#
# -------- General plotting style ---------------------------------# 
#------------------------------------------------------------------#

# Make zeros print as "0" always
library(stringr)
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}
# General theme
theme_general <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 

theme_general_x <- theme_bw() +
  theme(legend.position = "bottom", 
        text = element_text(size = 26), 
        axis.text.x = element_text(size = 26), 
        axis.text.y = element_text(size = 26)) 

#------------------------------------------------------------------#
#---------------- Figure 1.8 --------------------------------------#
#------------------------------------------------------------------#

# Create example data
times <- c(5,6,7,8,9,12,13,15,16,20,22,23)
age <- c(12,0,0,10,6,6,9,3,8,0,0,2)
event <- c(1,0,1,1,0,0,1,1,1,0,0,1)
id <- 1:12
exitage <- age + times
simpledata <- as.data.frame(cbind(id, times, event, age, exitage))
head(simpledata)


# Create Figure 1.8, part 1
fig18_1 <- ggplot(aes(x = times, y = id),
                data = simpledata) +
  geom_segment(aes(x = 0, y = id, xend = times - 0.25, yend = id), size = 1) +
  geom_point(aes(x = times, y = id, shape = as.factor(event)), size = 6) +
  scale_shape_manual(values = c(1, 16)) + 
  xlab("Time (t)") +
  ylab("Subject number") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 25)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0.5, 12), breaks = seq(1, 12, 1)) +
  theme_general_x +
  theme(legend.position = "none")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

fig18_1

ggsave("figures/h_simpledata.pdf", plot = fig18_1, 
       width = 21, height = 21, units = "cm")


# Create Figure 1.8, part 2
fig18_2 <- ggplot(aes(x = exitafe, y = id),
                  data = simpledata) +
  geom_segment(aes(x = age, y = id, xend = exitage - 0.25, yend = id), size = 1) +
  geom_point(aes(x = exitage, y = id, shape = as.factor(event)), size = 6) +
  scale_shape_manual(values = c(1, 16)) + 
  xlab("Age") +
  ylab("Subject number") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 25)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0.5, 12), breaks = seq(1, 12, 1)) +
  theme_general_x + 
  theme(legend.position = "none")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )


fig18_2

ggsave("figures/h_simpledataAge.pdf", plot = fig18_2, 
       width = 21, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 1.9 --------------------------------------#
#------------------------------------------------------------------#

# Kaplan-Meier time as time variable
kmfit <- survfit(Surv(times, event) ~ 1)

simpledata_km <- data.frame(surv = c(1, kmfit$surv), 
                            time = c(0, kmfit$time))
simpledata_km

# Create Figure 1.9, part 1
fig19_1 <- ggplot(aes(x = time, y = surv), data = simpledata_km) +
  geom_step(size = 1) +
  xlab("Time (t)") +
  ylab("Estimate of S(t)") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 25)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_general_x + 
  theme(legend.position = "none")

fig19_1

ggsave("figures/h_simpledata_KM_time.pdf", plot = fig19_1, 
       width = 21, height = 21, units = "cm")


# Kaplan-Meier age as time variable
kmfit2 <- survfit(Surv(age, exitage, event) ~ 1)

simpledata_km_age <- data.frame(surv = c(1, kmfit2$surv), 
                                time = c(0, kmfit2$time))
simpledata_km_age

# Create Figure 1.9, part 2
fig19_2 <- ggplot(aes(x = time, y = surv), data = simpledata_km_age) +
  geom_step(size = 1) +
  xlab("Age") +
  ylab("Estimate of S(age)") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 25)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_general_x + 
  theme(legend.position = "none")

fig19_2


ggsave("figures/h_simpledata_KM_age.pdf", plot = fig19_2, 
       width = 21, height = 21, units = "cm")

#------------------------------------------------------------------#
#---------------- Figure 1.11 -------------------------------------#
#------------------------------------------------------------------#


# N(t)
sfit <- survfit(Surv(times, event) ~ 1)
p1 <- data.frame(time = c(0, sfit$time, 26), 
                 n = c(0, cumsum(sfit$n.event), tail(cumsum(sfit$n.event),1))
                 )
p2 <- data.frame(time = sfit$time[sfit$n.event==1], 
                 n = cumsum(sfit$n.event)[sfit$n.event==1]-1)
p3 <- data.frame(time = sfit$time[sfit$n.event==1], 
                 n = cumsum(sfit$n.event)[sfit$n.event==1])

# Create Figure 1.11, part 1
fig111_1 <- ggplot(aes(x = time, y = n), data = p1) +
  geom_step(size = 1) + 
  geom_point(aes(x = time, y = n), data = p2, shape = 21, size = 6, color='black', fill='white') +
  geom_point(aes(x = time, y = n), data = p3, shape = 16, size = 6) +
  xlab("Time (t)") +
  ylab("N(t)") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005)), 
                     limits = c(0, 26)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 12), breaks = seq(0, 12, 1)) +
  theme_general_x + 
  theme(legend.position = "none")

fig111_1

ggsave("figures/h_simpledata_N_time.pdf", plot = fig111_1, 
       width = 21, height = 21, units = "cm")



# Y(t)
sfit <- survfit(Surv(times, event) ~ 1)
p1 <- data.frame(time = c(0, sfit$time, 26), 
                 risk = c(sfit$n.risk, 0, 0)
)
p2 <- data.frame(time = sfit$time, 
                 risk = sfit$n.risk - 1)
p3 <- data.frame(time = sfit$time, 
                 risk = sfit$n.risk)

# Create Figure 1.11, part 2
fig111_2 <- ggplot(aes(x = time, y = risk), data = p1) +
  geom_step(size = 1) + 
  geom_point(aes(x = time, y = risk), data = p2, shape = 21, size = 6, color='black', fill='white') +
  geom_point(aes(x = time, y = risk), data = p3, shape = 16, size = 6) +
  xlab("Time (t)") +
  ylab("Y(t)") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005)), 
                     limits = c(0, 26)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 12), breaks = seq(0, 12, 1)) +
  theme_general_x + 
  theme(legend.position = "none")

fig111_2

ggsave("figures/h_simpledata_atrisk_time.pdf", plot = fig111_2, 
       width = 21, height = 21, units = "cm")

# N(age)
asfit <- survfit(Surv(age, exitage, event) ~ 1)
p1 <- data.frame(time = c(0, asfit$time, 26), 
                 n = c(0, cumsum(asfit$n.event), tail(cumsum(asfit$n.event), 1))
)
p2 <- data.frame(time = asfit$time[asfit$n.event > 0], 
                 n = cumsum(asfit$n.event)[asfit$n.event > 0])
p3 <- data.frame(time = asfit$time[asfit$n.event > 0], 
                 n = c(0, 1, 2, 4, 5, 6))

# Create Figure 1.11, part 3
fig111_3 <- ggplot(aes(x = time, y = n), data = p1) +
  geom_step(size = 1) + 
  geom_point(aes(x = time, y = n), data = p2, shape = 16, size = 6) +
  geom_point(aes(x = time, y = n), data = p3, shape = 21, size = 6, color='black', fill='white') +
  xlab("Age") +
  ylab("N(age)") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005)), 
                     limits = c(0, 26)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 12), breaks = seq(0, 12, 1)) +
  theme_general_x + 
  theme(legend.position = "none")

fig111_3

ggsave("figures/h_simpledata_N_age.pdf", plot = fig111_3, 
       width = 21, height = 21, units = "cm")


# Y(age)
atriskage <- vector(length=25)
for(i in 1:25){
  atriskage[i] <- sum((age<=i & exitage>i))
}

p1 <- data.frame(time = 0:26, 
                 risk = c(atriskage[1], atriskage, 0)
)
p2 <- data.frame(time = c(2,3,6,7,8,9,10,12,15,17,18,20,22,24,25), 
                 risk = c(5,6,7,6,7,8, 9,10, 9, 8, 5, 4, 2, 1,0))
p3 <- data.frame(time = c(2,3,6,7,8,9,10,12,15,17,18,20,22,24,25), 
                 risk = c(4,5,6,7,6,7, 8,9, 10, 9, 8, 5, 4, 2,1))

# Create Figure 1.11, part 4
fig111_4 <- ggplot(aes(x = time, y = risk), data = p1) +
  geom_step(size = 1) + 
  geom_point(aes(x = time, y = risk), data = p2, shape = 21, size = 6, color='black', fill='white') +
  geom_point(aes(x = time, y = risk), data = p3, shape = 16, size = 6) +
  xlab("Age") +
  ylab("Y(age)") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005)), 
                     limits = c(0, 26)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 12), breaks = seq(0, 12, 1)) +
  theme_general_x + 
  theme(legend.position = "none")

fig111_4

ggsave("figures/h_simpledata_atrisk_age.pdf", plot = fig111_4, 
       width = 21, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 4.3 --------------------------------------#
#------------------------------------------------------------------#

# Plot content
x <- seq(0, 1, by = 0.001)
y1 <- -log(1-x)

# Plot data 
pdata <- data.frame(x = x,
                    y2 = x,
                    y1 = y1)
# Look at data
head(pdata)

# Create the plot
fig4_3 <- ggplot(aes(x = x, y = y2), data = pdata) +
  geom_line(size = 1) +
  geom_line(aes(x = x, y = y1), data = pdata, size = 1) +
  ylab("-log(1-x)") +
  xlab("x") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), 
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_general + 
  theme(legend.position = "none")

fig4_3

ggsave("figures/j_approxlog.pdf", plot = fig4_3, 
       width = 21, height = 21, units = "cm")


#------------------------------------------------------------------#
#---------------- Figure 5.4 --------------------------------------#
#------------------------------------------------------------------#

# Create mock survival data
surv0 <- function(t) exp(-t*0.6)
surv1 <- function(t) exp(-t*0.2)

ts <- (1:10-1) + (1:100-1)/100


pdata <- data.frame(type = c(rep("0", length(ts)), rep("1", length(ts))),
                    surv = c(surv0(ts), surv1(ts)), 
                    ts = c(ts, ts)
                    )


# Create the plot
fig54 <- ggplot(aes(x = ts, y = surv, linetype = type), data = pdata) +
  geom_line(size = 1) +
  ylab("Survival probability, U") +
  xlab("Survival time, T") +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.005)), 
                     limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.015)), 
                     limits = c(0, 10), breaks = seq(0, 10, 1)) +
  theme_general + 
  theme(legend.position = "none") + 
  geom_segment(aes(x=0, xend=-log(0.4)/0.2, y=0.4, yend=0.4), color = "grey", size = 1) + 
  geom_segment(aes(x=0, xend=pmin(10, -log(0.1)/0.2), y=0.1, yend=0.1), color = "grey", size = 1) + 
  geom_segment(aes(x=-log(0.4)/0.2, xend=-log(0.4)/0.2, y=0.4, yend=0), color = "grey", size = 1) +
  geom_segment(aes(x=-log(0.4)/0.6, xend=-log(0.4)/0.6, y=0.4, yend=0), color = "grey", size = 1) +
  geom_segment(aes(x=-log(0.1)/0.6, xend=-log(0.1)/0.6, y=0.1, yend=0), color = "grey", size = 1) +
  geom_segment(aes(x=pmin(10,-log(0.1)/0.2), xend=pmin(10,-log(0.1)/0.2), y=0.1, yend=0), color = "grey", size = 1) 

fig54

ggsave("figures/j_microsimfig.pdf", plot = fig54, 
       width = 29.7, height = 21, units = "cm")




