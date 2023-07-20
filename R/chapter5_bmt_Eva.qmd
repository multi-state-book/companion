################################################ DOWNLOADING PACKAGES ###################################################################################
library(survival)
library(mstate)
library(tidyverse) #Mostly for ggplot2
library(ggpubr) #Extra functions to ggplot2
library(xtable) #Converts table to LaTeX code
library(data.table) #Faster than data.frame
library(compiler) #Improves speed of functions a tiny bit
library(gridExtra) #Graphics

# General theme
theme_general <- theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.text = element_text(size = 22))

################################################# DATA PREPARATIONS ######################################################################################### 

# Read data and create a data frame
bmt <- read.csv("data/bmt.csv")
bmt <- bmt %>% mutate(outof0 = ifelse(gvhd == 1 | rel == 1 | death == 1, 1, 0),
                              outof0and1 = ifelse(rel == 1 | death == 1, 1,0),
                              dead = death,
                              timerel = ifelse(rel == 1, timerel, timedeath),
                              timegvhd = ifelse(gvhd ==1, timegvhd, timerel),
                              wait2 = timedeath - timerel,
                              intxsurv = timedeath,
                              intxrel = timerel,
                              tgvhd = timegvhd,)

################################# PREPARATIONS FOR AALEN-JOHANSEN & LANDMARK AALEN-JOHANSEN ############################################################################

#Transition matrix
tmat <- matrix(NA, 4,4)
tmat[1,2:4] <- 1:3
tmat[2,3:4] <- 4:5
tmat[3, 4] <- 6
dimnames(tmat) <- list(from = c("BMT", "GvHD", "Relapse", "Dead"), to = c("BMT", "GvHD", "Relapse", "Dead"))
tmat

#Creating data in long format, i.e. row for every transition
long_format <- msprep(time = c(NA, "tgvhd", "intxrel", "intxsurv"), status = c(NA,"gvhd", "rel", "dead"), data = as.data.frame(bmt), trans = tmat)
#Warning is due to 4 observations where time of relapse and death are the same.

#Cox model - no covariates, only strata
cox_base <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = long_format, method = "breslow")
#Warning due to start time being equal to stop time for some of the observations

#Subject specific transition hazards
msfit <- msfit(cox_base, trans = tmat)
#Warning due to the same as above


############################################### S = 3 ############################################################################################

#Aalen-Johansen
AJ3 <- probtrans(msfit, predt = 3, direction = "forward", method = "aalen")[[1]]
AJ3_data <- as.data.frame(cbind(AJ3$time, AJ3$pstate3))
colnames(AJ3_data) <- c("time", "p02")

AJ3_max60 <- subset(AJ3, time <= 60)


#Landmark Aalen-Johansen
LMAJ3 <- LMAJ(long_format, s = 3, from = 1, method = "aalen") #Same warning as before
LMAJ3_data <- as.data.frame(cbind(LMAJ3$time, LMAJ3$pstate3))
colnames(LMAJ3_data) <- c("time", "p02")

LMAJ3_max60 <- subset(LMAJ3, time <= 60)

############################################## S = 9 #############################################################################################

#Aalen-Johansen
AJ9 <- probtrans(msfit, predt = 9, direction = "forward", method = "aalen")[[1]]
AJ9_data <- as.data.frame(cbind(AJ9$time, AJ9$pstate3))
colnames(AJ9_data) <- c("time", "p02")

AJ9_max60 <- subset(AJ9, time <= 60)

#Landmark Aalen-Johansen
LMAJ9 <- LMAJ(long_format, s = 9, from = 1, method = "aalen")
LMAJ9_data <- as.data.frame(cbind(LMAJ9$time, LMAJ9$pstate3))
colnames(LMAJ9_data) <- c("time", "p02")

LMAJ9_max60 <- subset(LMAJ9, time <= 60)

####################################################################################################################################################
####################################################################################################################################################

########################### PREPARATIONS FOR SEMI-MARKOV AND PLUG-IN ESTIMATORS #####################################################################################################

#Transforming data to data.table for faster speed
bmt_dt <- setDT(bmt)

#Number of patients in state 0 at x, Y0(x)
patients_in_0 <- function(x){
  return(nrow(bmt_dt[tgvhd > (x - 0.000001)]))}

#Number of patients in state 1 at x, Y1(x)
patients_in_1 <- function(x){
  return(nrow(bmt_dt[gvhd ==1 & tgvhd < x & intxrel > (x - 0.000001)]))}

#Number of patients in state 2 at x, Y2(x)
patients_in_2 <- function(x){
  return(nrow(bmt_dt[rel ==1 & intxrel < x & intxsurv > (x - 0.000001)]))}

#alpha01 = dN01(u)/Y0(u)
alpha01 <- function(u){ 
  dN01 <- nrow(bmt_dt[tgvhd == u & gvhd == 1])
  Y0 <- patients_in_0(u)
  return(dN01/Y0)}

#alpha02 = dN02(u)/Y0(u)
alpha02 <- function(u){ 
  dN02 <- nrow(bmt_dt[gvhd == 0 & rel ==1 & intxrel == u])
  Y0 <- patients_in_0(u)
  return(dN02/Y0)}

#alpha12 = dN12(u)/Y1(u)
alpha12 <- function(u){ 
  dN12 <- nrow(bmt_dt[gvhd == 1 & rel ==1 & intxrel == u])
  Y1 <- patients_in_1(u)
  return(dN12/Y1)}

#P_00(s,u -) = \Pi_{u \in (s,t]} (1 - dN0(u)/Y0(u))
p00 <- function(start, end){
  exits <- bmt_dt[outof0 == 1 & tgvhd > start & tgvhd <= (end - 0.000001)]$tgvhd #Exits from state 0 in (s,t]
  if (length(exits) == 0){prob <- 1}
  else{
    exits_unique <- data.frame(table(exits))
    Y0 <- sapply(as.numeric(as.character(exits_unique$exits)), patients_in_0)
    frac <- (1 - exits_unique$Freq/Y0)
    prob <- prod(frac)}
  return(prob)}

#P_11(s,u-) = \Pi_{u \in (s,t]} (1 - dN1(u)/Y1(u))
p11 <- function(start,end){
  exits <- bmt_dt[gvhd ==1 & intxrel <= (end- 0.000001) & intxrel > start & (rel ==1 | dead ==1)]$intxrel #Time of exits from state 1 in (start, end)
  if (length(exits) == 0){prob <- 1}
  else {
    exits_unique <- data.frame(table(exits))
    Y1 <- sapply(as.numeric(as.character(exits_unique$exits)), patients_in_1)
    frac <- (1 - exits_unique$Freq/Y1)
    prob <- prod(frac)}
  return(prob)}

#Times of 0 -> 1 transitions
time01 <- sort(unique(bmt_dt[gvhd == 1]$tgvhd))

#Times of 1 -> 2 transitions
time12 <- sort(unique(bmt_dt[gvhd ==1 & rel == 1]$intxrel)) 

#################### FOR REPRODUCING THE TABLE p11alpha12. Takes approximately 10 minutes  ###################################################
#p11alpha12_table_raw <- matrix(NA, nrow = length(time01), ncol = length(time12))
#colnames(p11alpha12_table_raw) <- time12
#rownames(p11alpha12_table_raw) <- time01
#for (i in 1:length(time01)){
#  for (j in 1:length(time12)){
#    if (time01[i] > time12[j]){p11alpha12_table_raw[i,j] <- NA}
#    else {p11alpha12_table_raw[i,j] <- p11(time01[i],time12[j])*alpha12(time12[j])}}}
#p11alpha12_table <- as.data.table(p11alpha12_table_raw, keep.rownames = TRUE)
#############################################################################################################################################

#Data.table of P11(u,x-)alpha12(x)
p11alpha12_table <- read.csv("mellemregninger/p11alpha12tabel.csv")
colnames(p11alpha12_table) <- c("rn",time12)
p11alpha12_table <- as.data.table(p11alpha12_table)

#P_11(u,x-)alpha12(x)
p11alpha12 <- function(start,end){
  start_point <- length(time01[time01 <= start])
  end_point <- length(time12[time12 <= end]) +1
  return(as.numeric(p11alpha12_table[start_point,end_point, with = FALSE]))}

#Compilation of functions for tiny increase in speed performance
alpha01 <- cmpfun(alpha01)
alpha02 <- cmpfun(alpha02)
alpha12 <- cmpfun(alpha12)
patients_in_0 <- cmpfun(patients_in_0)
patients_in_1 <- cmpfun(patients_in_1)
patients_in_2 <- cmpfun(patients_in_2)
p11alpha12 <- cmpfun(p11alpha12)

########################################## S = 3 ##################################################################################################

#Time of 0 -> 1 transitions for u > 3
time01_s3 <- sort(unique(bmt_dt[tgvhd > 3 & gvhd == 1]$tgvhd))

#Time of 0 -> 2 transitions for u > 3  
time02_s3 <- sort(unique(bmt_dt[gvhd == 0 & rel == 1 & intxrel > 3]$intxrel))

#Matrix where 1st column is times for 0 -> 1 transition, u,
#and 4th column is P00(3,u -)*alpha01(u)
p00alpha01_s3 <- matrix(NA, nrow = length(time01_s3), ncol = 4)
p00alpha01_s3[, 1] <- time01_s3 #u, times of 0->1 transitions after 3 months
p00alpha01_s3[, 2] <- sapply(p00alpha01_s3[,1], p00, start = 3) #P00(3,u -)
p00alpha01_s3[, 3] <- sapply(p00alpha01_s3[,1], alpha01) #alpha01(u)
p00alpha01_s3[, 4] <- p00alpha01_s3[,2]*p00alpha01_s3[,3] #P00(3,u -)*alpha01(u)

#Matrix where 1st column is times for 0 -> 2 transition, u,
#and 4th column is P00(3,u -)*alpha02(u)
p00alpha02_s3 <- matrix(NA, nrow = length(time02_s3), ncol = 4)
p00alpha02_s3[, 1] <- time02_s3 #u, times of 0->2 transitions after 3 months
p00alpha02_s3[, 2] <- sapply(p00alpha02_s3[,1], p00, start = 3) #P00(3,u -)
p00alpha02_s3[, 3] <- sapply(p00alpha02_s3[,1], alpha02) #alpha02(u) 
p00alpha02_s3[, 4] <- p00alpha02_s3[,2]*p00alpha02_s3[,3] #P00(3,u -)*alpha(u)

######################################## S = 9 #####################################################################################################

#Time of 0 -> 1 transitions for u > 9
time01_s9 <- sort(unique(bmt_dt[tgvhd > 9 & gvhd == 1]$tgvhd))

#Time of 0 -> 2 transitions for u > 9
time02_s9 <- sort(unique(bmt_dt[gvhd == 0 & rel == 1 & intxrel > 9]$intxrel))

#Matrix where 1st column is times for 0 -> 1 transition, u,
#and 4th column is P00(6,u -)*alpha01(u)
p00alpha01_s9 <- matrix(NA, nrow = length(time01_s9), ncol = 4)
p00alpha01_s9[, 1] <- time01_s9 #u, times of 0->1 transitions after 9 months
p00alpha01_s9[, 2] <- sapply(p00alpha01_s9[,1], p00, start = 9) #P00(9,u -)
p00alpha01_s9[, 3] <- sapply(p00alpha01_s9[,1], alpha01) #alpha01(u)
p00alpha01_s9[, 4] <- p00alpha01_s9[,2]*p00alpha01_s9[,3] #P00(9,u -)*alpha01(u)

#Matrix where 1st column is times for 0 -> 2 transition, u,
#and 4th column is P00(9,u -)*alpha02(u)
p00alpha02_s9 <- matrix(NA, nrow = length(time02_s9), ncol = 4)
p00alpha02_s9[, 1] <- time02_s9 #u, times of 0->2 transitions after 9 months
p00alpha02_s9[, 2] <- sapply(p00alpha02_s9[,1], p00, start = 9) #P00(9,u -)
p00alpha02_s9[, 3] <- sapply(p00alpha02_s9[,1], alpha02) #alpha02(u) 
p00alpha02_s9[, 4] <- p00alpha02_s9[,2]*p00alpha02_s9[,3] #P00(9,u -)*alpha(u)

###################################################################################################################################################
####################################################################################################################################################

############################################## SEMI-MARKOV #################################################################################################

#Kaplan-Meier for being in state 2 under semi-Markov assumption.
s2_semi <- survfit(Surv(wait2, dead) ~1, data = bmt_dt[rel ==1], type = "kaplan-meier") 
#Table containing time of observed durations in state 2, and estimates of S2_semi for each time.
s2_semi_table <- as.data.table(cbind(s2_semi$time, s2_semi$surv))

#P_22(s,t) = P_22(t-s) under semi-Markov assumption.
p22_semi <- function(start,end){
  duration <- end - start
  sub_tabel <- s2_semi_table[V1 <= duration]$V2
  p22 <- sub_tabel[length(sub_tabel)]
  return(p22)}

##################################### S = 3 ########################################################################################################### 

#P_{0->2}(3,t) = \int_{3}^{t} P_00(3,u-)alpha02(u)P_22(u,t)
p02_direct_semi_s3 <- function(t){
  time <- p00alpha02_s3[, 1][p00alpha02_s3[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha2 <- p00alpha02_s3[, 4][p00alpha02_s3[, 1] <= t]
    p22 <- sapply(time, p22_semi, end = t)
    prob02 <- sum(p00alpha2*p22)}
  return(prob02)}

#P_{0->1->2}(3,t) = \int_{3}^{t} (P_00(3,u-)alpha01(u) \int_{u}^{t} (P_11(u,x-)alpha12(x)P_22(x,t)))
p02_non_direct_semi_s3 <- function(t){
  time <- p00alpha01_s3[, 1][p00alpha01_s3[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha1 <- p00alpha01_s3[, 4][p00alpha01_s3[, 1] <= t]
    p12 <- rep(0, length(time))
    possible_time12 <- time12[time12 < t]
    for (i in 1:length(time)){
      time12 <- possible_time12[possible_time12 > time[i]]
      if (length(time12) == 0){p12[i] <- 0}
      else{
        p11alpha2 <- vapply(time12, p11alpha12, start = time[i], numeric(1))
        p22 <- sapply(time12, p22_semi, end = t)
        p12[i] <- sum(p11alpha2*p22)}}
    prob02 <- sum(p00alpha1*p12)}
  return(prob02)}

#P_02(s,t) = P_{0->2}(s,t) + P_{0 -> 1 -> 2}(s,t)
p02_semi_s3 <- function(t){
  return(p02_direct_semi_s3(t) + p02_non_direct_semi_s3(t))}

########################################### s = 9 ##############################################################################################

#P_{0->2}(9,t) = \int_{9}^{t} P_00(9,u-)alpha02(u)P_22(u,t)
p02_direct_semi_s9 <- function(t){
  time <- p00alpha02_s9[, 1][p00alpha02_s9[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha2 <- p00alpha02_s9[, 4][p00alpha02_s9[, 1] <= t]
    p22 <- sapply(time, p22_semi, end = t)
    prob02 <- sum(p00alpha2*p22)}
  return(prob02)}

#P_{0 -> 1 -> 2}(9,t) = \int_{9}^{t} (P_00(9,u-)alpha01(u) \int_{u}^{t} (P_11(u,x-)alpha12(x)P_22(x,t)))
p02_non_direct_semi_s9 <- function(t){
  time <- p00alpha01_s9[, 1][p00alpha01_s9[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha1 <- p00alpha01_s9[, 4][p00alpha01_s9[, 1] <= t]
    p12 <- rep(0, length(time))
    possible_time12 <- time12[time12 < t]
    for (i in 1:length(time)){
      time12 <- possible_time12[possible_time12 > time[i]]
      if (length(time12) == 0){p12[i] <- 0}
      else{
        p11alpha2 <- sapply(time12, p11alpha12, start = time[i])
        p22 <- sapply(time12, p22_semi, end = t)
        p12[i] <- sum(p11alpha2*p22)}}
    prob02 <- sum(p00alpha1*p12)}
  return(prob02)}

#P_02(s,t) = P_{0->2}(s,t) + P_{0 -> 1 -> 2}(s,t)
p02_semi_s9 <- function(t){
  return(p02_direct_semi_s9(t) + p02_non_direct_semi_s9(t))}

########################################################################################################################################################
########################################################################################################################################################

########################################### Landmark Pepe #####################################################################################################

############################################# s = 3 ###################################################################################################

#S_012, Kaplan-Meier for being in state 0,1 or 2 among patients in state 0 at time 3
Pepe_S012_s3 <- survfit(Surv(intxsurv, dead) ~1, data = bmt_dt[tgvhd > 3]) 
Pepe_S012_s3_table <- as.data.table(cbind(Pepe_S012_s3$time, Pepe_S012_s3$surv))
#S_01, Kaplen-Meier for being in state 0 or 1 among patients in state 0 at time 3
Pepe_S01_s3 <- survfit(Surv(intxrel, outof0and1)~1, data = bmt_dt[tgvhd > 3]) 
Pepe_S01_s3_table <- as.data.table(cbind(Pepe_S01_s3$time, Pepe_S01_s3$surv))

#S_012(t) - S_01(t) among subjects in state 0 at time 3
LMPepe_s3 <- function(t){ 
  S012t <- Pepe_S012_s3_table$V2[nrow(Pepe_S012_s3_table[V1 <= t])]
  S01t <- Pepe_S01_s3_table$V2[nrow(Pepe_S01_s3_table[V1 <= t])]
  if (length(S012t) == 0){S012t <- 1}
  if (length(S01t) == 0){S01t <-1}
  Pepe_est <- S012t - S01t
  return(Pepe_est)}

############################################# s = 9 ###################################################################################################

#S_012, Kaplan-Meier for being in state 0,1 or 2 among patients in state 0 at time 9
Pepe_S012_s9 <- survfit(Surv(intxsurv, dead) ~1, data = bmt_dt[tgvhd > 9]) 
Pepe_S012_s9_table <- as.data.table(cbind(Pepe_S012_s9$time, Pepe_S012_s9$surv))
#S_01, Kaplan-Meier for being in state 0 or 1 among patients in state 0 at time 9
Pepe_S01_s9 <- survfit(Surv(intxrel, outof0and1)~1, data = bmt_dt[tgvhd > 9]) 
Pepe_S01_s9_table <- as.data.table(cbind(Pepe_S01_s9$time, Pepe_S01_s9$surv))

#S_012(t) - S_01(t) among subjects in state 0 at time 9
LMPepe_s9 <- function(t){ 
  S012t <- Pepe_S012_s9_table$V2[nrow(Pepe_S012_s9_table[V1 <= t])]
  S01t <- Pepe_S01_s9_table$V2[nrow(Pepe_S01_s9_table[V1 <= t])]
  if (length(S012t) == 0){S012t <- 1}
  if (length(S01t) == 0){S01t <-1}
  Pepe_est <- S012t - S01t
  return(Pepe_est)}

########################################################################################################################################################
########################################################################################################################################################

######################################### PLUG-IN WITH DURATION OF GVHD AS LINEAR EFFECT ###############################################################

#Subset of data only containing subjects with GvHD = 1.
gvhd_subset <- subset(bmt, gvhd ==1)

#Cox model for alpha12 (GvHD -> Relapse) with a linear effect of duration in state 1 as covariate
cox_p12_tgvhd <- coxph(Surv(tgvhd, intxrel, rel) ~ tt(tgvhd), tt = function(x,t,...) t - x, data = gvhd_subset, method = "breslow") 
(beta12 <- cox_p12_tgvhd$coefficients)

#Baseline hazard for alpha12
basehaz12 <- function(t){
  dN12 <- nrow(bmt_dt[gvhd == 1 & rel == 1 & intxrel == t])
  Y1 <- bmt_dt[gvhd ==1 & tgvhd < t & intxrel > (t-0.0000001)]
  Y1beta <- sum(exp(beta12*(t - Y1$tgvhd)))
  base_alpha12 <- dN12/Y1beta
  return(base_alpha12)}
basehaz12_table <- as.data.table(cbind(time12,sapply(time12,basehaz12)))

#Hazard for alpha12 with linear effect of duration in GvHD after s.
haz12 <- function(s,x){
  if (x %in% time12){
    haz <- basehaz12_table[basehaz12_table$time12 == x]$V2*exp(beta12*(x-s))}
  else {haz <- 0}
  return(haz)}

#Time of 1 -> 3 transitions
time13 <- unique(sort(bmt_dt[gvhd ==1 & rel == 0 & dead == 1]$intxsurv))

#Cox model for alpha13 (GvHD -> Dead) with a linear effect of duration in state 1 as covariate
cox_p13_tgvhd <- coxph(Surv(tgvhd, intxsurv, dead) ~ tt(tgvhd), tt = function(x,t,...) t - x, data = gvhd_subset, method = "breslow") 
(beta13 <- cox_p13_tgvhd$coefficients)

#Baseline hazard for alpha13, dN13(x)/sum_i(Y_1i(x)*exp(beta(x-T01i)))
basehaz13 <- function(t){
  dN13 <- nrow(bmt_dt[gvhd == 1 & dead == 1 & rel == 0 & intxsurv == t])
  Y1 <- bmt_dt[gvhd ==1 & tgvhd < t & intxrel > (t-0.0000001)]
  Y1beta <- sum(exp(beta13*(t - Y1$tgvhd)))
  base_alpha13 <- dN13/Y1beta
  return(base_alpha13)}
basehaz13_table <- as.data.table(cbind(time13,sapply(time13,basehaz13)))

#Hazard for alpha12 with linear effect of duration in GvHD after s.
haz13 <- function(s,x){
  if (x %in% time13){
    haz <- basehaz13_table[basehaz13_table$time13 == x]$V2*exp(beta13*(x-s))}
  else {haz <- 0}
  return(haz)}

#P_11(s,t-) with linear effect of duration in state 1
p11_linear_gvhd <- function(start,end){
  if (start > (end - 0.000001)){p11 <- NA}
  else{
    time12_sub <- time12[start < time12 & time12 <= (end - 0.000001)]
    if (length(time12_sub) == 0){int12 <- 0}
    else {int12 <- sum(sapply(time12_sub, haz12, s = start))}
    time13_sub <- time13[start < time13 & time13 <= (end - 0.000001)]
    if (length(time13_sub) == 0){int13 <- 0}
    else{int13 <- sum(sapply(time13_sub, haz13, s = start))}
    p11 <- exp(-(int12+int13))}
  return(p11)}


#################################################################################################################### 
### For reproducing data.table p11alpha12_linear_gvhd_table. Takes approximately 5 minutes ###
#Table of p11alpha12 with linear effect of duration in state 1.
#p11alpha12_table_raw2 <- matrix(NA, nrow = length(time01), ncol = length(time12))
#colnames(p11alpha12_table_raw2) <- time12
#rownames(p11alpha12_table_raw2) <- time01
#for (i in 1:length(time01)){
#  for (j in 1:length(time12)){
#    if (time01[i] > time12[j]){p11alpha12_table_raw2[i,j] <- NA}
#    else {p11alpha12_table_raw2[i,j] <- p11_linear_gvhd(time01[i],time12[j])*haz12(s = time01[i], x = time12[j])}}}
#p11alpha12_linear_gvhd_table <- as.data.table(p11alpha12_table_raw2, keep.rownames = TRUE)
###########################################################################################################

p11alpha12_linear_gvhd_table <- read.csv("mellemregninger/p11alpha12-linear-gvhd.csv")
colnames(p11alpha12_linear_gvhd_table) <- c("rn",time12)
p11alpha12_linear_gvhd_table <- as.data.table(p11alpha12_linear_gvhd_table)

#P_11(u,x-)alpha12(x) with linear effect of duration in state 1.
p11alpha12_lingvhd <- function(start,end){
  return(unlist(p11alpha12_linear_gvhd_table[(match(start,time01)), (match(end, time12) + 1), with = FALSE]))}

#S2 = \Pi_{u <= t} (1 - dN2(u)/Y2(u)), Kaplan-Meier
S2_markov <- function(t){
  exits <- bmt_dt[rel ==1 & dead ==1 & intxsurv <= t]$intxsurv
  if (length(exits) == 0) {prob <- 1}
  else{
    exits_unique <- data.frame(table(exits))
    Y2 <- vapply(as.numeric(as.character(exits_unique$exits)), patients_in_2, numeric(1))
    frac <- (1 - exits_unique$Freq/Y2)
    prob <- prod(frac)}
  return(prob)}

#####################################################################################################################################
### For reproducing S2_markov_table. Takes approximately 10 minutes. ###
#Data.table with 1st column times, t, and 2nd column S2_markov(t)
#times_for_S2_markov <- seq(0,149.34, by = 0.01)
#S2_markov_of_time <- sapply(times_for_S2_markov, S2_markov)
#S2_markov_table <- as.data.table(cbind(times_for_S2_markov, S2_markov_of_time))
#####################################################################################################################################

s2_markov_table <- read.csv("mellemregninger/s2_table.csv")

#P22(u,t) = S2(t)/S2(u)
p22_markov <- function(start,end){
  s2_start <- s2_markov_table[start*100 + 1,2]
  s2_end <- s2_markov_table[end*100 + 1,2]
  return(as.numeric(s2_end/s2_start))}

######################################################### s = 3 #####################################################################
p02_direct_lingvhd_s3 <- function(t){
  time <- p00alpha02_s3[, 1][p00alpha02_s3[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha2 <- p00alpha02_s3[, 4][p00alpha02_s3[, 1] <= t]
    p22 <- sapply(time, p22_markov, end = t)
    prob02 <- sum(p00alpha2*p22)}
  return(prob02)}

#P_{0 -> 1 -> 2}(9,t) = \int_{9}^{t} (P_00(9,u-)alpha01(u) \int_{u}^{t} (P_11(u,x-)alpha12(x)P_22(x,t)))
p02_non_direct_lingvhd_s3 <- function(t){
  time <- p00alpha01_s3[, 1][p00alpha01_s3[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha1 <- p00alpha01_s3[, 4][p00alpha01_s3[, 1] <= t]
    p12 <- rep(0, length(time))
    possible_time12 <- time12[time12 < t]
    for (i in 1:length(time)){
      new_time12 <- possible_time12[possible_time12 > time[i]]
      if (length(new_time12) == 0){p12[i] <- 0}
      else{
        p11alpha2 <- sapply(new_time12, p11alpha12_lingvhd, start = time[i])
        p22 <- sapply(new_time12, p22_markov, end = t)
        p12[i] <- sum(p11alpha2*p22)}}
    prob02 <- sum(p00alpha1*p12)}
  return(prob02)}

#P_02(s,t) = P_{0->2}(s,t) + P_{0 -> 1 -> 2}(s,t)
p02_lingvhd_s3 <- function(t){
  return(p02_direct_lingvhd_s3(t) + p02_non_direct_lingvhd_s3(t))}

######################################################### s = 9 #####################################################################
p02_direct_lingvhd_s9 <- function(t){
  time <- p00alpha02_s9[, 1][p00alpha02_s9[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha2 <- p00alpha02_s9[, 4][p00alpha02_s9[, 1] <= t]
    p22 <- sapply(time, p22_markov, end = t)
    prob02 <- sum(p00alpha2*p22)}
  return(prob02)}

#P_{0 -> 1 -> 2}(9,t) = \int_{9}^{t} (P_00(9,u-)alpha01(u) \int_{u}^{t} (P_11(u,x-)alpha12(x)P_22(x,t)))
p02_non_direct_lingvhd_s9 <- function(t){
  time <- p00alpha01_s9[, 1][p00alpha01_s9[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha1 <- p00alpha01_s9[, 4][p00alpha01_s9[, 1] <= t]
    p12 <- rep(0, length(time))
    possible_time12 <- time12[time12 < t]
    for (i in 1:length(time)){
      new_time12 <- possible_time12[possible_time12 > time[i]]
      if (length(new_time12) == 0){p12[i] <- 0}
      else{
        p11alpha2 <- sapply(new_time12, p11alpha12_lingvhd, start = time[i])
        p22 <- sapply(new_time12, p22_markov, end = t)
        p12[i] <- sum(p11alpha2*p22)}}
    prob02 <- sum(p00alpha1*p12)}
  return(prob02)}

#P_02(s,t) = P_{0->2}(s,t) + P_{0 -> 1 -> 2}(s,t)
p02_lingvhd_s9 <- function(t){
  return(p02_direct_lingvhd_s9(t) + p02_non_direct_lingvhd_s9(t))}

###########################################################################################################################################
###########################################################################################################################################

###################################### PLUG-IN WITH LINEAR EFFECT OF DURATION OF RELAPSE ###############################################

rel_subset <- bmt_dt[rel == 1]

cox_23_linrel <- coxph(Surv(intxrel, intxsurv, dead) ~ tt(intxrel), tt = function(x,t,...) t - x, data = rel_subset, method = "breslow")
(beta23 <- cox_23_linrel$coefficients)

timej2 <- unique(sort(bmt_dt[rel ==1]$intxrel)) #Time of j -> 2 transitions, j = {0,1}
time23 <- unique(sort(bmt_dt[rel ==1 & dead ==1]$intxsurv)) #Time of 2 -> 3 transitions

#Baseline hazard for alpha23
basehaz23 <- function(t){
  dN23 <- nrow(bmt_dt[rel ==1 & dead == 1 & intxsurv == t])
  at_risk <- bmt_dt[rel == 1 & intxrel < t & intxsurv > (t - 0.000001)] 
  weigthed_at_risk <- sum(exp(beta23*(t - at_risk$intxrel)))
  haz23 <- dN23/weigthed_at_risk
  return(haz23)}
basehaz23_table <- as.data.table(cbind(time23, sapply(time23,basehaz23)))


haz23_linrel <- function(s,x){ #alpha23(x | T_j2 = s) = basehaz23(x)*exp(beta23(x-s))
  if (x %in% time23){
    haz <- basehaz23_table[basehaz23_table$time23 == x]$V2*exp(beta23*(x - s))}
  else {haz <- 0}
  return(haz)}

p22_linrel_slow <- function(start,end){
  if (start > end){p22 <- NA}
  else{
    times <- time23[time23 > start & time23 <= end]
    if (length(times) == 0){p22 <- 1}
    else {
      terms <- sapply(times, haz23_linrel, s = start)
      int <- sum(terms)
      p22 <- exp(-int)}}
  return(p22)}

###############################################################################################################################
### For reproducing data.table p22_linrel_table. Takes approximately 10 minutes ###
#Table of p22 with linear effect of relapse as covariate.
#p22_linrel_table_raw <- matrix(NA, nrow = length(timej2), ncol = length(time23))
#rownames(p22_linrel_table_raw) <- timej2
#colnames(p22_linrel_table_raw) <- time23
#for (i in 1:length(timej2)){
#  for (j in 1:length(time23)){
#    if (timej2[i] > time23[j]){p22_linrel_table_raw[i,j] <- NA}
#    else {p22_linrel_table_raw[i,j] <- p22_linrel_slow(timej2[i],time23[j])}}}
#p22_linrel_table <- as.data.table(p22_linrel_table_raw, keep.rownames = TRUE)
###############################################################################################################################

p22_linrel_table <- read.csv("mellemregninger/p22_linrel_table.csv")
colnames(p22_linrel_table) <- c("rn", time23)
p22_linrel_table <- as.data.table(p22_linrel_table)

#P_22(s,t)with linear effect of relapse as covariate.
p22_linrel <- function(start,end){
  end_point <- length(time23[time23 <=end]) + 1
  start_point <- length(timej2[timej2 <= start])
  p22 <- as.numeric(p22_linrel_table[start_point, end_point, with = FALSE])
  if (is.na(p22)) {p22 <- 1}
  return(p22)}


##################################### S = 3 ########################################################################################################### 

#P_{0->2}(3,t) = \int_{3}^{t} P_00(3,u-)alpha02(u)P_22(u,t)
p02_direct_linrel_s3 <- function(t){
  time <- p00alpha02_s3[, 1][p00alpha02_s3[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha2 <- p00alpha02_s3[, 4][p00alpha02_s3[, 1] <= t]
    p22 <- sapply(time, p22_linrel, end = t)
    prob02 <- sum(p00alpha2*p22)}
  return(prob02)}

#P_{0->1->2}(3,t) = \int_{3}^{t} (P_00(3,u-)alpha01(u) \int_{u}^{t} (P_11(u,x-)alpha12(x)P_22(x,t)))
p02_non_direct_linrel_s3 <- function(t){
  time <- p00alpha01_s3[, 1][p00alpha01_s3[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha1 <- p00alpha01_s3[, 4][p00alpha01_s3[, 1] <= t]
    p12 <- rep(0, length(time))
    possible_time12 <- time12[time12 < t]
    for (i in 1:length(time)){
      time12_sub <- possible_time12[possible_time12 > time[i]]
      if (length(time12_sub) == 0){p12[i] <- 0}
      else{
        p11alpha2 <- vapply(time12_sub, p11alpha12, start = time[i], numeric(1))
        p22 <- sapply(time12_sub, p22_linrel, end = t)
        p12[i] <- sum(p11alpha2*p22)}}
    prob02 <- sum(p00alpha1*p12)}
  return(prob02)}

#P_02(s,t) = P_{0->2}(s,t) + P_{0 -> 1 -> 2}(s,t)
p02_linrel_s3 <- function(t){
  return(p02_direct_linrel_s3(t) + p02_non_direct_linrel_s3(t))}



########################################### s = 9 ##############################################################################################



#P_{0->2}(9,t) = \int_{9}^{t} P_00(9,u-)alpha02(u)P_22(u,t)
p02_direct_linrel_s9 <- function(t){
  time <- p00alpha02_s9[, 1][p00alpha02_s9[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha2 <- p00alpha02_s9[, 4][p00alpha02_s9[, 1] <= t]
    p22 <- sapply(time, p22_linrel, end = t)
    prob02 <- sum(p00alpha2*p22)}
  return(prob02)}

#P_{0 -> 1 -> 2}(9,t) = \int_{9}^{t} (P_00(9,u-)alpha01(u) \int_{u}^{t} (P_11(u,x-)alpha12(x)P_22(x,t)))
p02_non_direct_linrel_s9 <- function(t){
  time <- p00alpha01_s9[, 1][p00alpha01_s9[, 1] <= t]
  if (length(time) == 0){prob02 <- 0}
  else{
    p00alpha1 <- p00alpha01_s9[, 4][p00alpha01_s9[, 1] <= t]
    p12 <- rep(0, length(time))
    possible_time12 <- time12[time12 < t]
    for (i in 1:length(time)){
      time12 <- possible_time12[possible_time12 > time[i]]
      if (length(time12) == 0){p12[i] <- 0}
      else{
        p11alpha2 <- sapply(time12, p11alpha12, start = time[i])
        p22 <- sapply(time12, p22_linrel, end = t)
        p12[i] <- sum(p11alpha2*p22)}}
    prob02 <- sum(p00alpha1*p12)}
  return(prob02)}

#P_02(s,t) = P_{0->2}(s,t) + P_{0 -> 1 -> 2}(s,t)
p02_linrel_s9 <- function(t){
  return(p02_direct_linrel_s9(t) + p02_non_direct_linrel_s9(t))}


###########################################################################################################################################
###########################################################################################################################################

###################################################### Figure 5.10 ##############################################################################

#For reproducing the data.frame semi_s3_data - takes approximately 1.5 hour
#Data.frame for plotting semi-Markov
#time_semi_s3 <- seq(3.01, 60, by = 0.03)
#semi_p02_s3 <- sapply(time_semi_s3, p02_semi_s3)
#semi_s3_data <- as.data.frame(cbind(time_semi_s3, semi_p02_s3))
#colnames(semi_s3_data) <- c("time", "pstate3")

semi_s3_data <- read.csv("mellemregninger/semi_s3_data.csv")

#LM Pepe data for plot with s = 3
time_pepe_s3 <- seq(3.01,60,0.01)
est_pepe_s3 <- sapply(time_pepe_s3, LMPepe_s3)
data_pepe_s3 <- as.data.frame(cbind(time_pepe_s3, est_pepe_s3))

#### For reproducing the data.frame lingvhd_s3_data. Takes approximately 1 hour
#Plug-in (linear effect of duration of gvhd) data for plot with s = 3
#time_lingvhd_s3 <- seq(3.01,60,by = 0.03)
#lingvhd_p02_s3 <- sapply(time_lingvhd_s3, p02_lingvhd_s3)
#lingvhd_s3_data <- as.data.frame(cbind(time_lingvhd_s3,lingvhd_p02_s3))
#colnames(lingvhd_s3_data) <- c("time", "pstate3")

lingvhd_s3_data <- read.csv("mellemregninger/data_lingvhd_s3.csv")


#For reproducing the data.frame linrel_s3_data. Takes approximately 1.5 hour
#time_linrel_s3 <- AJ3_max60$time
#time_linrel_s3 <- seq(3.01, 60, by = 0.03)
#linrel_p02_s3 <- sapply(time_linrel_s3, p02_linrel_s3)
#linrel_s3_data <- as.data.frame(cbind(time_linrel_s3, linrel_p02_s3))
#colnames(linrel_s3_data) <- c("time", "pstate3")

linrel_s3_data <- read.csv("mellemregninger/linrel_s3_data.csv")

#Plot with AaJ, LM AaJ and LM Pepe, s = 3
threemonths_ggplot_all <- ggplot(AJ3_max60,aes(time*(1/12),pstate3)) + 
  geom_step(data=semi_s3_data,aes(color="Semi")) +
  geom_step(data=lingvhd_s3_data,aes(color="LinGvHD")) +
  geom_step(data=linrel_s3_data,aes(color="LinRel")) +
  geom_step(aes(color="AaJ"))+
  geom_step(data=LMAJ3_max60,aes(color="LM AaJ")) +
  geom_step(data = data_pepe_s3, aes(x = time_pepe_s3*(1/12), y = est_pepe_s3, color = "LM Pepe")) +
  xlab("Time (years)") +
  ylab("Probability") +
  ylim(0,0.035)  +
  scale_color_manual(values = c("blue", "red", "green4", "black", "orange", "magenta")) + 
  theme_minimal() +
  labs(color="Estimator", title = "P(V(t) = 2 | V(0.25) = 0)") + border(size = 0.7) +
  theme(legend.position = c(0.92,0.85), legend.title = element_blank()) 
threemonths_ggplot_all

fig5.10 <- ggplot(AJ3_max60,aes(time,pstate3)) + 
  geom_step(data=semi_s3_data,aes(linetype="Semi")) +
  geom_step(data=lingvhd_s3_data,aes(linetype="LinGvHD")) +
  geom_step(data=linrel_s3_data,aes(linetype="LinRel")) +
  geom_step(aes(linetype="AaJ"))+
  geom_step(data=LMAJ3_max60,aes(linetype="LM AaJ")) +
  geom_step(data = data_pepe_s3, aes(x = time_pepe_s3, y = est_pepe_s3, linetype = "LM Pepe")) +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Probability") +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 60), breaks = seq(0, 60, by = 12)) + 
  scale_linetype_manual(values = c("solid", "twodash", "dotted", "longdash", "dotdash", "dashed")) +
  theme_general + 
  theme(legend.key.width = unit(1.5,"cm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 22))

fig5.10


ggsave("figures/h_BMT-all-s3.pdf", plot = fig5.10,
       width = 29.7, height = 21, units = "cm")


######################################### Figure 5.11 #################################################################################################

#For reproducing the data.frame semi_s9_data - takes approximately 15 minutes
#Data.frame for plotting semi-Markov
#time_semi_s9 <- seq(9.01, 60, by = 0.03)
#semi_p02_s9 <- sapply(time_semi_s9, p02_semi_s9)
#semi_s9_data <- as.data.frame(cbind(time_semi_s9, semi_p02_s9))
#colnames(semi_s9_data) <- c("time", "pstate3")

semi_s9_data <- read.csv("mellemregninger/semi_s9_data.csv")

#LM Pepe data for plot with s = 9
time_pepe_s9 <- seq(9.01,60,0.01)
est_pepe_s9 <- sapply(time_pepe_s9, LMPepe_s9)
data_pepe_s9 <- as.data.frame(cbind(time_pepe_s9,est_pepe_s9))

# For reproducing lingvhd_s9_data. Takes approximately 10 minutes.
#Plug-in (linear effect of duration of gvhd) data for plot with s = 9
#time_lingvhd_s9 <- seq(9.01,60,by = 0.03)
#lingvhd_p02_s9 <- sapply(time_lingvhd_s9, p02_lingvhd_s9)
#lingvhd_s9_data <- as.data.frame(cbind(time_lingvhd_s9,lingvhd_p02_s9))
#colnames(lingvhd_s9_data) <- c("time", "pstate3")

lingvhd_s9_data <- read.csv("mellemregninger/lingvhd_s9_data.csv")

#For reproducing the data.frame linrel_s9_data. Takes approximately 30 minutes
#time_linrel_s9 <- seq(9.01, 60, by = 0.03)
#linrel_p02_s9 <- sapply(time_linrel_s9, p02_linrel_s9)
#linrel_s9_data <- as.data.frame(cbind(time_linrel_s9, linrel_p02_s9))
#colnames(linrel_s9_data) <- c("time", "pstate3")

linrel_s9_data <- read.csv("mellemregninger/linrel_s9_data.csv")


#Plot with AaJ, LM AaJ, LM Pepe, s = 9
ninemonths_ggplot_all <- ggplot(AJ9_max60,aes(time*(1/12),pstate3)) + 
  geom_step(data=semi_s9_data,aes(color="Semi")) +
  geom_step(data=lingvhd_s9_data,aes(color="LinGvHD")) +
  geom_step(data=linrel_s9_data,aes(color="LinRel")) +
  geom_step(data = data_pepe_s9, aes(x = time_pepe_s9*(1/12), y = est_pepe_s9, color = "LM Pepe")) +
  geom_step(aes(color="AaJ"))+
  geom_step(data=LMAJ9_max60,aes(color="LM AaJ")) +
  xlab("Time (years)") +
  ylab("Probability") +
  ylim(0,0.035)  +
  scale_color_manual(values = c("blue", "red", "green4","black", "orange", "magenta")) + theme_minimal() +
  labs(color="Estimator", title = "P(V(t) = 2 | V(0.75) = 0)") + border(size = 0.7) +
  theme(legend.position = c(0.92,0.85), legend.title = element_blank(),legend.margin = margin(t = 0,r = 0,l = 0, b = 0)) 
ninemonths_ggplot_all


fig5.11 <- ggplot(AJ9_max60,aes(time,pstate3)) + 
  geom_step(data=semi_s9_data,aes(linetype="Semi")) +
  geom_step(data=lingvhd_s9_data,aes(linetype="LinGvHD")) +
  geom_step(data=linrel_s9_data,aes(linetype="LinRel")) +
  geom_step(data = data_pepe_s9, aes(x = time_pepe_s9, y = est_pepe_s9, linetype= "LM Pepe")) +
  geom_step(aes(linetype="AaJ"))+
  geom_step(data=LMAJ9_max60,aes(linetype="LM AaJ")) +
  xlab("Time since bone marrow transplantation (months)") +
  ylab("Probability") +
  ylim(0,0.035)  +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.05)), 
                     limits = c(0, 60), breaks = seq(0, 60, by = 12)) + 
  scale_linetype_manual(values = c("solid", "twodash", "dotted", "longdash", "dotdash", "dashed")) +
  theme_general + 
  theme(legend.key.width = unit(1.5,"cm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 22))

fig5.11

ggsave("figures/h_BMT-all-s9.pdf", plot = fig5.11,
       width = 29.7, height = 21, units = "cm")

####################################################### Table 5.5 #################################################################

cox12 <- coxph(formula = Surv(tgvhd, intxrel, rel) ~ tt(tgvhd), data = subset(bmt, gvhd ==1), tt = function(x, t, ...) t - x, method = "breslow")
summary(cox12)

cox13 <- coxph(formula = Surv(tgvhd, intxsurv, dead) ~ tt(tgvhd), data = subset(bmt, gvhd ==1), tt = function(x, t, ...) t - x, method = "breslow")
summary(cox13)

cox23 <- coxph(formula = Surv(intxrel, intxsurv, dead) ~ tt(intxrel), data = subset(bmt, rel == 1), tt = function(x, t, ...) t - x, method = "breslow")
summary(cox23)
