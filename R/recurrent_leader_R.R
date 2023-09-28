###################################################################################
###################################################################################
################## LEADER ARTICLE #################################################
###################################################################################
###################################################################################

# Load required packages 
require(survival)
require(mstate)
require(ggplot2)
require(survminer)
require(haven)


# Set working directory
setwd("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf")


# Read leader data 
recur <- read_sas("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/recurrent_leader.sas7bdat")

recur <- as.data.frame(recur); head(recur)
# View(recur)

table(recur$eventno, recur$TRT01P) # also censorings
colSums(table(recur$eventno, recur$TRT01P)[-1,]) # Number of events 


with(subset(recur, status == 1), table(eventno, TRT01P)) # Events - not taking censoring into account
colSums(with(subset(recur, status == 1), table(eventno, TRT01P)))

#### NON PARAMETRIC MODELS #### 

np_model1 <- survfit(Surv(starttime, stoptime, status2) ~ factor(TRT01P),
                     type = "fh2", 
                     data = recur)
summary(np_model1)
ggsurvplot(np_model1, risk.table = TRUE) # survival
ggsurvplot(np_model1, fun = function(x){-log(x)}, risk.table = TRUE) + ylab("Cumulative hazard") #cumhaz

# By previous no. events. 
recur$trt <- as.factor(recur$TRT01P)
recur$enum <- as.factor(recur$eventno)

np_model2 <- survfit(Surv(starttime, stoptime, status2) ~ trt + enum,
                     type = "fh2", 
                     data = recur)
summary(np_model2)

ggsurvplot(np_model2, risk.table = TRUE)$plot + facet_grid( ~ trt) # Probabilities

ggsurvplot(np_model2, fun = function(x){-log(x)}, risk.table = TRUE)$plot + facet_grid( ~ trt) + ylab("Cumulative hazard") # Cum haz


#### semi PARAMETRIC MODELS #### 
recur$trt2 <- relevel(recur$trt, ref = "Placebo")

# AG model - not stratified by no. events
sp_model1 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, # + cluster(SUBJID), 
                   robust = TRUE, method = "breslow", 
                   data = recur)
summary(sp_model1)

# Create survival curves and hazards by hand.. 
data <- data.frame(haz = c(basehaz(sp_model1)$hazard, basehaz(sp_model1)$hazard * exp(coef(sp_model1)[[1]])),
                   time = rep(basehaz(sp_model1)$time, 2), 
                   trt = c(rep("Lira", length(basehaz(sp_model1)$time)), rep("Placebo", length(basehaz(sp_model1)$time)))
                           ) 

ggplot(aes(x = time, y = haz), data = data) + geom_step(aes(color = trt))

# AG model - stratified by no. events
sp_model2 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + strata(eventno), 
                   subset = (eventno < 2),
                   robust = TRUE, method = "breslow", 
                   data = recur)
summary(sp_model2)

# predict(sp_model2, type = "lp")

# Not sure how to extract stratified hazards ... 

data <- data.frame (haz = coxph.detail(sp_model2)$hazard,
                    time = coxph.detail(sp_model2)$time,
                    strata = c(rep("noevents=1", coxph.detail(sp_model2)$strata[2]),
                               rep("noevents=2", coxph.detail(sp_model2)$strata[3]),
                               rep("noevents=3", coxph.detail(sp_model2)$strata[4]),
                               rep("noevents=4", coxph.detail(sp_model2)$strata[5]),
                               rep("noevents=5", coxph.detail(sp_model2)$strata[6]),
                               rep("noevents=6", coxph.detail(sp_model2)$strata[7]),
                               rep("noevents=7", coxph.detail(sp_model2)$strata[8]),
                               rep("noevents=8", coxph.detail(sp_model2)$strata[9])

                    )
)

ggplot(aes(x=time, y=haz), data = data) + geom_step(aes(color=strata))





# Marginal mean model
sp_model3 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + cluster(SUBJID), 
                   method = "breslow", 
                   data = recur)
summary(sp_model3)


# PWP total-time model
sp_model4 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + strata(enum) + cluster(SUBJID), 
                   method = "breslow",
                   data = recur)
summary(sp_model4)

# One model per transition 
sp_model4_1 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, 
                     method = "breslow", subset = (eventno == 1),
                     data = recur)
summary(sp_model4_1)

sp_model4_2 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, 
                     method = "breslow", subset = (eventno == 2),
                     data = recur)
summary(sp_model4_2)

sp_model4_3 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, 
                     method = "breslow", subset = (eventno == 3),
                     data = recur)
summary(sp_model4_3)

sp_model4_4 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, 
                     method = "breslow", subset = (eventno == 4),
                     data = recur)
summary(sp_model4_4)

sp_model4_5 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, 
                     method = "breslow", subset = (eventno == 5),
                     data = recur)
summary(sp_model4_5)




# Frailty model
sp_model6 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + frailty(SUBJID), 
                   method = "breslow", 
                   data = recur)
summary(sp_model6)


###########################################################################################
###########################################################################################

# In order to prepare for the WLW model, we need to restructure the data.. 
# Pull it from SAS

# Read leader data 
wlw_recur <- read_sas("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/wlw_recurrent_leader.sas7bdat")

wlw_recur <- as.data.frame(wlw_recur); head(wlw_recur)

wlw_recur$trt <- as.factor(wlw_recur$TRT01P)
wlw_recur$trt2 <- relevel(wlw_recur$trt, ref = "Placebo")


# Overall - just stratifying
cox_wlw <- coxph(Surv(stoptime, status2) ~ strata(wlw_eventno) + trt2*factor(wlw_eventno) + cluster(SUBJID), 
                 data = wlw_recur, method = "breslow")
summary(cox_wlw)
# You need to change the parametrization to get the same numbers as below.. 
# Similar as to multiplying the effects together.. 
# Consider: is it possible just to get one treatment effect estimate? (That would make it non marginal)


## YOU PROB DON'T NEED CLUSTER(ID) here
# T1: Time to entry in state 1 (1 event)
cox_wlw_event1 <- coxph(Surv(stoptime, status2) ~ trt2, 
                        data = subset(wlw_recur, wlw_eventno == 1), method = "breslow")
summary(cox_wlw_event1)

# T2: Time to entry in state 2 (2 events)
cox_wlw_event2 <- coxph(Surv(stoptime, status2) ~ trt2, 
                        data = subset(wlw_recur, wlw_eventno == 2), method = "breslow")
summary(cox_wlw_event2)

# T3: Time to entry in state 3 (3 events)
cox_wlw_event3 <- coxph(Surv(stoptime, status2) ~ trt2, 
                        data = subset(wlw_recur, wlw_eventno == 3), method = "breslow")
summary(cox_wlw_event3)

# T4: Time to entry in state 4 (4 events)
cox_wlw_event4 <- coxph(Surv(stoptime, status2) ~ trt2, 
                        data = subset(wlw_recur, wlw_eventno == 4), method = "breslow")
summary(cox_wlw_event4)

# T5: Time to entry in state 5 (5 events)
cox_wlw_event4 <- coxph(Surv(stoptime, status2) ~ trt2, 
                        data = subset(wlw_recur, wlw_eventno == 5), method = "breslow")
summary(cox_wlw_event4)



################################################################################################
################################################################################################
### MAKE A FUNCTION - outputs per endpoint #####################################################
################################################################################################


# Load required packages 
require(survival)
require(mstate)
require(ggplot2)
require(survminer)
require(haven)


# Set working directory
setwd("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf")


# Read leader data 
recur_all <- read_sas("P:/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/recurrent_leader_all.sas7bdat")

recur_all <- as.data.frame(recur_all); head(recur_all)


# By previous no. events. 
recur_all$trt <- as.factor(recur_all$TRT01P)
recur_all$enum <- as.factor(recur_all$eventno)
recur_all$trt2 <- relevel(recur_all$trt, ref = "Placebo")


################## SEMI-PARAMETRIC MODELS #################################

model_fitting <- function(endpointdat){


# AG model - not stratified by no. events
sp_model1 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, # + cluster(SUBJID), 
                   robust = TRUE, method = "breslow", 
                   data = endpointdat)
#summary(sp_model1)

# # AG model - stratified by no. events
# sp_model2 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + strata(eventno), 
#                    robust = TRUE, method = "breslow", 
#                    data = endpointdat)
# summary(sp_model2)


# Marginal mean model
sp_model3 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + cluster(SUBJID), 
                   method = "breslow", 
                   data = endpointdat)
#summary(sp_model3)


# PWP total-time model
sp_model4 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + strata(enum) + cluster(SUBJID), 
                   method = "breslow",
                   data = endpointdat)
#summary(sp_model4)

# One model per transition 
sp_model4_1 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, 
                     method = "breslow", subset = (eventno == 1),
                     data = endpointdat)
#summary(sp_model4_1)

sp_model4_2 <- coxph(Surv(starttime, stoptime, status2) ~ trt2, 
                     method = "breslow", subset = (eventno == 2),
                     data = endpointdat)
#summary(sp_model4_2)

sp_model4_3 <- coxph(Surv(starttime, stoptime, status2) ~ trt2,
                     method = "breslow", subset = (eventno == 3),
                     data = endpointdat)
# #summary(sp_model4_3)

sp_model4_4 <- coxph(Surv(starttime, stoptime, status2) ~ trt2,
                     method = "breslow", subset = (eventno == 4),
                     data = endpointdat)
# #summary(sp_model4_4)
#
sp_model4_5 <- coxph(Surv(starttime, stoptime, status2) ~ trt2,
                     method = "breslow", subset = (eventno == 5),
                     data = endpointdat)
# #summary(sp_model4_5)


# Frailty model
sp_model6 <- coxph(Surv(starttime, stoptime, status2) ~ trt2 + frailty(SUBJID), 
                   method = "breslow", 
                   data = endpointdat)
#summary(sp_model6)

mylist <- list()
mylist[[1]] <- summary(sp_model1) # AG
mylist[[2]] <- summary(sp_model3) # Marg mean
mylist[[3]] <- summary(sp_model4) # PWP
mylist[[4]] <- summary(sp_model4_1) # PWP trans 1
mylist[[5]] <- summary(sp_model4_2) # PWP trans 2 
mylist[[6]] <- summary(sp_model4_3) # PWP trans 3
mylist[[7]] <- summary(sp_model4_4) # PWP trans 4
mylist[[8]] <- summary(sp_model4_5) # PWP trans 5
mylist[[9]] <- summary(sp_model6) # Frailty

return(mylist)

}


####################
### NON-FATAL MI ###
####################

endpointdat_nfmi <- subset(recur_all, id == "recurrent_nfmi")

nf_mi <- model_fitting(endpointdat = endpointdat_nfmi)

nf_mi[[1]]
nf_mi[[2]]
nf_mi[[3]]
nf_mi[[4]]
nf_mi[[5]]
nf_mi[[6]]
nf_mi[[7]]
nf_mi[[8]]
nf_mi[[9]]



# Non-parametric fits

fit_1 <- survfit(Surv(starttime, stoptime, status2) ~ trt2 + eventno, data = endpointdat_nfmi,
                   #subset = ((from == 1) & (to == 2)), #& (trt2 == "Placebo")), 
                   type = "fh2")

ggsurvplot(fit_1, fun = "cumhaz")

fit_2 <- survfit(Surv(starttime, stoptime, status2) ~ trt2, data = endpointdat_nfmi,
                   #subset = ((from == 1) & (to == 2)), #& (trt2 == "Placebo")), 
                   type = "fh2")

ggsurvplot(fit_2, fun = "cumhaz")


########################
### NON-FATAL STROKE ###
########################

endpointdat_nfstroke <- subset(recur_all, id == "recurrent_nfstroke")
nf_stroke <- model_fitting(endpointdat = endpointdat_nfstroke)

nf_stroke[[1]]
nf_stroke[[2]]
nf_stroke[[3]]
nf_stroke[[4]]
nf_stroke[[5]]
# nf_stroke[[6]]
# nf_stroke[[7]]
nf_stroke[[9]]

fit_1 <- survfit(Surv(starttime, stoptime, status2) ~ trt2, data = endpointdat_nfstroke,
                 type = "fh2")

ggsurvplot(fit_1, fun = "cumhaz")

fit_2 <- survfit(Surv(starttime, stoptime, status2) ~ trt2 + eventno, data = endpointdat_nfstroke,
                 type = "fh2")

ggsurvplot(fit_2, fun = "cumhaz")



####################
### MI #############
####################

endpointdat_mi <- subset(recur_all, id == "recurrent_mi")

mi <- model_fitting(endpointdat = endpointdat_mi)

mi[[1]]
mi[[2]]
mi[[3]]
mi[[4]]
mi[[5]]
mi[[6]]
mi[[7]]
mi[[8]]
mi[[9]]

fit_1 <- survfit(Surv(starttime, stoptime, status2) ~ trt2, data = endpointdat_mi,
                 #subset = ((from == 1) & (to == 2)), #& (trt2 == "Placebo")), 
                 type = "fh2")

ggsurvplot(fit_1, fun = "cumhaz")



########################
### Stroke #############
########################

endpointdat_stroke <- subset(recur_all, id == "recurrent_stroke")

stroke <- model_fitting(endpointdat = endpointdat_stroke)

stroke[[1]]
stroke[[2]]
stroke[[3]]
stroke[[4]]
stroke[[5]]
stroke[[6]]
stroke[[9]]


#############################################
### Composite: NF stroke, NF MI #############
#############################################

endpointdat_comb_nfstr_nfmi <- subset(recur_all, id == "recurrent_comb_nfstr_nfmi")

comb_nfstr_nfmi <- model_fitting(endpointdat = endpointdat_comb_nfstr_nfmi)

comb_nfstr_nfmi[[1]]
comb_nfstr_nfmi[[2]]
comb_nfstr_nfmi[[3]]
comb_nfstr_nfmi[[4]]
comb_nfstr_nfmi[[5]]
comb_nfstr_nfmi[[6]]
comb_nfstr_nfmi[[7]]
comb_nfstr_nfmi[[8]]
comb_nfstr_nfmi[[9]]


#############################################
### Composite: stroke, MI ###################
#############################################

endpointdat_comb_str_mi <- subset(recur_all, id == "recurrent_comb_str_mi")

comb_str_mi <- model_fitting(endpointdat = endpointdat_comb_str_mi)

comb_str_mi[[1]]
comb_str_mi[[2]]
comb_str_mi[[3]]
comb_str_mi[[4]]
comb_str_mi[[5]]
comb_str_mi[[6]]
comb_str_mi[[7]]
comb_str_mi[[8]]
comb_str_mi[[9]]

fit_1 <- survfit(Surv(starttime, stoptime, status2) ~ trt2, data = endpointdat_comb_str_mi,
                 type = "fh2")

ggsurvplot(fit_1, fun = "cumhaz")

fit_2 <- survfit(Surv(starttime, stoptime, status2) ~ trt2 + eventno, data = endpointdat_comb_str_mi,
                 type = "fh2")

ggsurvplot(fit_2, fun = "cumhaz")


#############################################
### Composite: NF stroke, NF MI, CV death ###
#############################################

endpointdat_comb_nfstr_nfmi_cvdth <- subset(recur_all, id == "recurrent_comb_nfstr_nfmi_cvdth")

comb_nfstr_nfmi_cvdth <- model_fitting(endpointdat = endpointdat_comb_nfstr_nfmi_cvdth)

comb_nfstr_nfmi_cvdth[[1]]
comb_nfstr_nfmi_cvdth[[2]]
comb_nfstr_nfmi_cvdth[[3]]
comb_nfstr_nfmi_cvdth[[4]]
comb_nfstr_nfmi_cvdth[[5]]
comb_nfstr_nfmi_cvdth[[6]]
comb_nfstr_nfmi_cvdth[[7]]
comb_nfstr_nfmi_cvdth[[8]]

comb_nfstr_nfmi_cvdth[[9]]

#############################################
### Composite: Stroke, MI, CV death ###
#############################################

endpointdat_comb_str_mi_cvdth <- subset(recur_all, id == "recurrent_comb_str_mi_cvdth")

comb_str_mi_cvdth <- model_fitting(endpointdat = endpointdat_comb_str_mi_cvdth)

comb_str_mi_cvdth[[1]]
comb_str_mi_cvdth[[2]]
comb_str_mi_cvdth[[3]]
comb_str_mi_cvdth[[4]]
comb_str_mi_cvdth[[5]]
comb_str_mi_cvdth[[6]]
comb_str_mi_cvdth[[7]]
comb_str_mi_cvdth[[8]]

comb_str_mi_cvdth[[9]]






##############################################
## ADJUSTMENT FOR COMPETING RISK, NON CV DEATH 


fit_1 <- survfit(Surv(starttime, stoptime, status2) ~ trt2, data = endpointdat_comb_nfstr_nfmi_cvdth,
                 #subset = ((from == 1) & (to == 2)), #& (trt2 == "Placebo")), 
                 type = "fh2")

ggsurvplot(fit_1, fun = "cumhaz")

data <- data.frame(cumhaz = fit_1$cumhaz, 
                   time = fit_1$time, 
                   treat = c(rep("Placebo", fit_1$strata[[1]]), rep("Lira", fit_1$strata[[2]])),
                   type = rep("Nelson-Aalen", length(fit_1$cumhaz))
)


ggplot(aes(x = time / 365.25, y = cumhaz), data = data) + geom_step(aes(color = treat))

# With death -> cause specific 
require("cmprsk")
require("prodlim")
fit_rec <- survfit(Surv(starttime, stoptime, status == 1) ~ trt2, data = endpointdat_comb_nfstr_nfmi_cvdth,
                   #subset = ((from == 1) & (to == 2)), #& (trt2 == "Placebo")), 
                   type = "fh2")
fit_dea <- survfit(Surv(starttime, stoptime, status == 2) ~ trt2, data = endpointdat_comb_nfstr_nfmi_cvdth,
                   #subset = ((from == 1) & (to == 2)), #& (trt2 == "Placebo")), 
                   type = "fh2")

fit_adj <- prodlim(Surv(starttime, stoptime, status == 1) ~ trt2, data = endpointdat_comb_nfstr_nfmi_cvdth)
fit_adj_d <- prodlim(Surv(starttime, stoptime, status == 2) ~ trt2, data = endpointdat_comb_nfstr_nfmi_cvdth)
fit_adj$hazard

fit_adj <- cifreg(Event(starttime, stoptime, status == 1) ~ trt2, data = endpointdat_comb_nfstr_nfmi_cvdth, cause = 1)

plot(fit_adj, cause = 2) #; plot(fit_adj_d)

ggsurvplot(fit_rec, fun = "cumhaz")
ggsurvplot(fit_dea, fun = "cumhaz")

data_adj <- data.frame(cumhaz = fit_adj$cumhaz, 
                       time = fit_adj$time, 
                       treat = c(rep("Placebo", fit_adj$strata[[1]]), rep("Lira", fit_adj$strata[[2]])), 
                       type = rep("Aalen-Johansen", length(fit_adj$cumhaz))
)


# 
# data_adj <- data.frame(#cumhaz = fit_rec$cumhaz, 
#                        cumhaz = -log(fit_adj$surv),
#                        time = fit_adj$time, 
#                        treat = c(rep("Placebo", fit_adj$size.strata[[1]]), rep("Lira", fit_adj$size.strata[[2]])), 
#                        type = rep("Aalen-Johansen", length(fit_adj$surv))
# )
#             

combdata <- rbind(data, data_adj)

ggplot(aes(x = time, y = cumhaz), data = combdata) + geom_step(aes(color = treat, linetype = type))



table(endpointdat_comb_nfstr_nfmi_cvdth$status)
aj <- prodlim(Hist(stoptime, status) ~ trt2, data = endpointdat_comb_nfstr_nfmi_cvdth)
plot(aj, cause = 1)
plot(aj, cause = 2)





