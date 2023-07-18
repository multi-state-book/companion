## ---- Getdata
bissau <- data.frame(read.csv("data/bissau.csv"))

## ---- Table-3.1

table(bissau$bcg, bissau$dtp)
100*table(bissau$bcg, bissau$dtp) / rowSums(table(bissau$bcg, bissau$dtp))

## ---- Table-3.2

# Add extra variables
bissau$agein  <- with(bissau, age/(365.24/12))
bissau$ageout <- with(bissau, agein+fuptime/(365.24/12))
bissau$dtpany <- 1*with(bissau, dtp>0)

library(survival)
coxph(Surv(agein,ageout,dead!=0)~bcg,data=bissau,method="breslow",timefix=F)
coxph(Surv(agein,ageout,dead!=0)~dtpany,data=bissau,method="breslow",timefix=F)
coxph(Surv(agein,ageout,dead!=0)~bcg+dtpany,data=bissau,method="breslow",timefix=F)
coxph(Surv(agein,ageout,dead!=0)~bcg*dtpany,data=bissau,method="breslow",timefix=F)

