*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- CHAPTER 2 - EXERCISE SOLUTIONS ---------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We first load the data;

proc import out = chs_data
    datafile = 'C:/HRfinal/holter/cphholter.csv'
	dbms= csv replace;
	getnames=yes;
run;

* We will convert the time variables ('timeafib', 'timestroke', and 'timedeath') from days to years;

data chs_data;
	set chs_data;
	timeafib = timeafib/365.25;
	timestroke = timestroke/365.25;
	timedeath = timedeath/365.25;
run;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.1 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*------------------------------------------- 2.1.1 --------------------------------------------------------------------------------;

* We will estimate the cumulative hazards of death for subjects with or without ESVEA non-parametrically with the Nelson-Aalen 
  estimator. This can be done with the 'phreg' procedure. Since the event of interest is death the time variable used in the model 
  statement is 'timedeath' and the censoring variable is 'death'. Furthermore, 'esvea' is specified in the strata statement to 
  obtain different estimates of the hazards for subjects with or without ESVEA. The cumulative hazard are the stored in the
  'hazdata' file created with the baseline statement;

proc phreg data=chs_data noprint; 
	model timedeath*death(0)=;
	strata esvea;
	baseline out=hazdata cumhaz=naa;
run;

* The estimates are then plotted using the gplot procedure;

title "2.1.1: Nelson-Aalen estimates for subjects with or without ESVEA";
proc gplot data=hazdata;
	plot naa*timedeath=esvea /haxis=axis1 vaxis=axis2;
	axis1 order=0 to 15.19 by 1 label=('Time (Years)');
	axis2 order=0 to 1 by 0.2 label=(a=90 'Cumulative hazard');
	symbol1 i=stepjl c=red;
	symbol2 i=stepjl c=blue;
run;

*------------------------------------------- 2.1.2 --------------------------------------------------------------------------------;

* The hazards can be compared with a log rank test which is implemented in the lifetest procedure. The relevant time variable 
  'timedeath' and censoring variable 'death' are specified in the time statement and 'esvea' is used in the strata statement.;

title "2.1.2: Log rank test comparing the hazards for the subjects with or without ESVEA";
proc lifetest data=chs_data notable plots=none; 
	time timedeath*death(0);
	strata esvea;
run; 

* We get a Chi-square statistic of 17.66 on 1 degree of freedom with a corresponding p-value <0.0001;

*------------------------------------------- 2.1.3 --------------------------------------------------------------------------------;

* To estimate the cumulative hazard under the assumption that the hazard is constant within 5-year intervals we must first split 
  each record into subrecords with cut times at 5 and 10 years. We will create a new data set 'chs_pch213' where 'cens' indicates if the
  subject died during each 5-year interval, 'risktime' is the time at risk during each 5-year interval, logrisk is the logarithm of
  the risk time and 'interval' is a categorical variable with 6 levels corresponding to the different combinations of interval and 
  ESVEA status;
  
data chs_pch213; 
	set chs_data;
	cens=(timedeath<=5)*(death = 1);
	risktime=min(5,timedeath);
	logrisk = log(risktime);
	interval=1;
	if esvea = 1 then do; interval =4; end; output;
	if timedeath>5 then do;
	cens=(timedeath<=10)*(death = 1);
	risktime=min(5,timedeath-5);
	logrisk = log(risktime);
	interval=2;
	if esvea = 1 then do; interval = 5; end; output; end;
	if timedeath>10 then do;
	cens= (death = 1); 
	risktime=timedeath-10;
	logrisk = log(risktime);
	interval=3;
	if esvea = 1 then do; interval = 6; end; output; end;
run;


* We can then estimate the hazard within each of the intervals [0,5), [5,10), [10,15.19). 
  This can be done with the sql procedure;

title "2.1.3: Estimate of the hazards for subjects with or without ESVEA under the piece-wise constant hazards assumption";
proc sql;
    select interval,
    sum(cens) as sum_death, 
    sum(risktime) as sum_risktime,
    calculated sum_death / calculated sum_risktime as hazard
    from chs_pch213
    group by interval;
quit;

* Alternatively, all estimates can be obtained at once through a Poisson regression using the genmod procedure. Since 'interval' is 
  a categorical variable it must be included in the class statement. 'cens' is included in the model statement as the dependent
  variable whilst interval is included as the only explanatory variable. Furthermore the probability distribution is specified through 
  'dist = poi' and the logarithm of the offset is included as well as the argument 'noint' which is added to obtain a model without an 
  intercept. The estimates for the hazard for each of interval and status of ESVEA are then given using the estimate statement. These 
  estimates corresponds to the exponential of the estimates from the table of the maximum likelihood parameter estimation;

proc genmod data=chs_pch213;
	class interval;
	model cens=interval / dist=poi offset=logrisk noint;
	estimate '0-5 years, ESVEA = 0' interval 1 0 0 0 0 0;
	estimate '5-10 years, ESVEA = 0' interval 0 1 0 0 0 0;
	estimate '10+ years, ESVEA = 0' interval 0 0 1 0 0 0;
	estimate '0-5 years, ESVEA = 1' interval 0 0 0 1 0 0;
	estimate '5-10 years, ESVEA = 1' interval 0 0 0 0 1 0;
	estimate '10+ years, ESVEA = 1' interval 0 0 0 0 0 1;
run;

* To investigate if the hazards differs between the two groups we must fit a Poisson model where we do not stratify by ESVEA 
  and then compare the two Poisson models using a likelihood ratio test (LRT). Thus, we repeat the procedure and first create
  a new data set 'chs_pch_null213' where no distinction is made between subjects with or without ESVEA;

data chs_pch_null213;
  set chs_pch213;
  if interval in (4,5,6) then interval = interval - 3;
run;

* Then, the Poisson model is fitted;

title "2.1.3: Estimate of the hazard under the piece-wise constant hazards assumption with no distinction between subjects with or 
       without ESVEA";
proc genmod data=chs_pch_null213;
	class interval;
	model cens=interval/dist=poi offset=logrisk noint;
	estimate '0-5 years' interval 1 0 0;
	estimate '5-10 years' interval 0 1 0;
	estimate '10+ years' interval 0 0 1;
run;

* Lastly, the two models are compared with a likelihood ratio test. This is done by comparing the max log likelihoods for the two models 
  and then calculating the p-value using the Chi-square distribution with 3 degrees of freedom. The likelihood scores are found in the 'Criteria For
  Accessing Goodness Of Fit' table;
  
title "2.1.3: Likelihood ratio test comparing the two piece-wise constant hazards models";
data p;
	chi2=1444.6862-1428.4662;
	p=1-probchi(chi2,3);
proc print;
run;

* We get a Chi-square statistics of 16.22 on 3 degrees of freedom with a corresponding p-value of 0.001
  We conclude once more that the hazards seem different for the subjects with or without ESVEA also under this model.; 

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.2 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*------------------------------------------- 2.2.1 --------------------------------------------------------------------------------;

* To enable the analysis of the composite end-point stroke or death we will create a new time variable with the smallest value of
  'timestroke' and 'timedeath' called 'timestrokeordeath' and a censoring variable called 'strokeordeath' which is 1 if stroke and/
  or death occured and 0 otherwise for each subject.;

title "2.2.1";
data chs_data;
	set chs_data;
	timestrokeordeath = timedeath;
	if stroke = 1 then timestrokeordeath = timestroke;
	strokeordeath = death;
	if stroke = 1 then strokeordeath = 1;
run;

*------------------------------------------- 2.2.2 --------------------------------------------------------------------------------;

* We then repeat the procedure from exercise 2.1.1 to obtain the Nelson-Aalen estimate of the cumulative hazards.;

title "2.2.2";
proc phreg data=chs_data; 
	model timestrokeordeath*strokeordeath(0)=;
	strata esvea;
	baseline out=hazdat cumhaz=naa;
run;

* The estimates are then plotted;

title "2.2.2: Nelson-Aalen estimate of the cumulative hazards for stroke-free survival for subjects with or without ESVEA ";
proc gplot data=hazdat;
	plot naa*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 15.19 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;

*------------------------------------------- 2.2.3 --------------------------------------------------------------------------------;

* The hazards are compared with a log rank test which is implemented in the lifetest procedure;

title "2.2.3: Log rank test comparing the hazards for the subjects with or without ESVEA";
proc lifetest data=chs_data notable plots=none; 
	time timestrokeordeath*strokeordeath(0);
	strata esvea;
run; 

* With a Chi-square statistic of 18.6 and a corresponding p-value of <0.0001 we conclude that hazards are different
  for the groups with or without ESVEA under this model; 

*------------------------------------------- 2.2.4 --------------------------------------------------------------------------------;

* To estimate the cumulative hazard under the piece-wise constant hazards assumption we need to split each record into subrecords 
  with cut times at 5 and 10 years as we did in exercise 2.1.3.;

title "2.2.3";
data chs_pch; 
	set chs_data;
	cens=(timestrokeordeath<=5)*(strokeordeath=1);
	risktime=min(5,timestrokeordeath);
	logrisk=log(risktime); interval=1;
	if esvea=1 then do; interval=4; end; output;
	if timestrokeordeath>5 then do;
	cens=(timestrokeordeath<=10)*(strokeordeath=1);
	risktime=min(5,timestrokeordeath-5);
	logrisk=log(risktime); 
	interval=2;
	if esvea=1 then do; 
	interval=5; end; output; end;
	if timestrokeordeath>10 then do;
	cens=(strokeordeath=1); 
	risktime=timestrokeordeath-10;
	logrisk=log(risktime); interval=3;
	if esvea=1 then do; 
	interval = 6; end; output; end;
run;

* Then, the Poisson model is fitted using the genmod procedure.;

proc genmod data=chs_pch;
	class interval;
	model cens=interval/dist=poi offset=logrisk noint;
	estimate '0-5 years, ESVEA = 0'  interval 1 0 0 0 0 0;
	estimate '5-10 years, ESVEA = 0' interval 0 1 0 0 0 0;
	estimate '10+ years, ESVEA = 0' interval 0 0 1 0 0 0;
	estimate '0-5 years, ESVEA = 1' interval 0 0 0 1 0 0;
	estimate '5-10 years, ESVEA = 1' interval 0 0 0 0 1 0;
	estimate '10+ years, ESVEA = 1' interval 0 0 0 0 0 1;
run;

* To investigate if the hazards are different for the two groups we will fit a model where we do not stratify by ESVEA. 
  Once again we will first make a new version of the data set.;
  
data chs_pch_null;
  set chs_pch;
  if interval in (4,5,6) then interval = interval - 3;
run;

* Then, the model which does not distinguish between subject with or without ESVEA is fitted using the genmod procedure;

proc genmod data=chs_pch_null;
	class interval;
	model cens=interval/dist=poi offset=logrisk noint;
	estimate '0-5 years' interval 1 0 0;
	estimate '5-10 years' interval 0 1 0;
	estimate '10+ years' interval 0 0 1;
run;

* Lastly, the two models are compared with a likelihood ratio test. The deviance scores are found in the table 'Criteria for 
  assesing goodness of fit'.;

title "2.2.4: Likelihood ratio test comparing the two piece-wise constant hazards models";
data p;
	chi2=1558.4777-1534.7568;
	p=1-probchi(chi2,3);
proc print;
run;

* We get a Chi-square statistic of 23.72 on 3 degrees of freedom with a very small p-value. Thus, we conclude once
  more that the hazards are not the same for the two groups;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.3 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*------------------------------------------- 2.3.1 --------------------------------------------------------------------------------;

* A Cox model can be fitted with the phreg procedure. The model statement requires the time varible and a censoring variable to
  be specified on the left hand side of '=' and the explanatory variables on the right hand side. The censoring variable have to be followed 
  by a parenthesis containing the censoring value(s);

title "2.3.1";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea;
run;

* We get an estimated coefficient for ESVEA of 0.62849 which corresponds to a hazard ratio of exp(0.62849) = 1.874778. Thus, having 
  ESVEA almost doubles the rate of experiencing the composite end-point.

* We can conclude that the hazards for the subjects with or without ESVEA are different under the proportional hazards
  assumption as well, since we get a Chi-square statistic of 17.99 on 1 degree of freedom with a correpsonding p-value <0.0001.

*------------------------------------------- 2.3.2 --------------------------------------------------------------------------------;

* Likewise, we will investigate the influence of ESVEA under the Poisson model. We will thus add the variable ESVEA as
  an explanatory variable to our model statement;

title "2.3.2";
proc genmod data=chs_pch_null;
	class interval;
	model cens=interval esvea/dist=poi offset=logrisk;
run;

* We get a hazard ration of exp(0.6209) = 1.860602 for ESVEA under the assumption of piece-wise constant hazards model. This 
  difference in hazards between subjects with or without ESVEA is significant, since we get a Chi-square statistic of 17.56 on 1 
  degree of freedom with a corresponding p-value <0.0001.;

*------------------------------------------- 2.3.3 --------------------------------------------------------------------------------;

* The hazard ratio between subjects with or without ESVEA found using the Cox model is 1.87 while it was 1.86 assuming the 
  piece-wise constant hazards model. Thus, the estimates are nearly identical both predicting an increase of the rate for subjects 
  with ESVEA of a factor 1.86-1.87;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.4 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*------------------------------------------- 2.4.1 --------------------------------------------------------------------------------;

* We add sex, age, and systolic blood pressure as explanatory variables to our Cox model from before;

title "2.4.1";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp;
run;

* We now get a hazard ratio of exp(0.31830) = 1.374789 for ESVEA. The effect of ESVEA is still significant with Chi-square statistic
  of 4.3514 on 1 degree of freedom and a corresponding p-value of 0.0370;

*------------------------------------------- 2.4.2 --------------------------------------------------------------------------------;

* To obtain an estimate of the hazard ratio for ESVEA under the piece-wise constant hazards model where we also adjust for sex, age,
  and systolic blood pressure we must include these as explanatory variables in the model statement.;

title "2.4.2";
proc genmod data=chs_pch_null;
	class interval (ref = '1');
	model cens=interval esvea sex age sbp/dist=poi offset=logrisk;
run;

* We find a hazard ratio of ESVEA of exp(0.3111) = 1.364926. The hazards seems to be different between the subjects with or without ESVEA
  since we get a Chi-square statistic of 4.16 on 1 degree of freedom with a corresponding p-value of 0.0414;

*------------------------------------------- 2.4.3 --------------------------------------------------------------------------------;

* Again, the hazard ratios estimated with the Cox model and the piece-wise constant hazards model are almost identical. The 
  estimated hazard ratio for ESVEA is 1.36-1.37. This is however smaller than the hazard ratio estimated when we did not also adjust
  for sex, age, and systolic blood pressure.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.5 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*------------------------------------------- 2.5.1 --------------------------------------------------------------------------------;

* To check the assumption of proportional hazards for subjects with or without ESVEA under the Cox model we will fit seperate baseline
  hazards. This can be done by including a strata statement where ESVEA is specified. Furthermore, we specify the values of the 
  other covariates (sex, age, and systolic blood pressure) in the data frame 'covstr'. The hazard is then plotted against each 
  other. If the proportional hazards assumption holds the curve should be close to a straight line through (0,0) with a slope equal
  to the estimated hazard ratio of ESVEA.;

data covstr;
	sex=0; age=0; sbp=0;
run;

* The Cox model is fitted and a data set with the baseline hazards is saved as 'breslowstr';

title "2.5.1";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=sex age sbp;
	strata esvea;
	baseline out=breslowstr cumhaz=breslow covariates=covstr;
run;

* Now the baseline hazards for esvea = 0,1 are extracted and merged into one data set called 'breslow01';

data breslow0;
	set breslowstr; 
	if esvea=0; 
	a00=breslow; 
run;

data breslow1; 
	set breslowstr; 
	if esvea=1; 	
	a01=breslow; 
run;

data breslow01; 
	merge breslow0 breslow1; 
	by timestrokeordeath; 
run;

* Then, the empty cells of the cumulative hazards are replaced with the previous values;

data breslow01; 
set breslow01;
	by timestrokeordeath;
	retain last1 last2;
	if a00=. then cumhaz0=last1; if a00 ne . then cumhaz0=a00; 
	if a01=. then cumhaz1=last2; if a01 ne . then cumhaz1=a01;
	output;
	last1=cumhaz0; last2=cumhaz1; 
run;

* The straight line with a slope which equals the predicted hazard ratio for esvea from the model in 2.4.1, is added to the data set;

data breslow01; 
	set breslow01;
	line=1.3649*cumhaz0;
run;

* Finally, the cumulative baseline hazards for esvea = 0,1 are plotted against each other together with the straight line;

title "2.5.1: Checking the proportionality assumption of the hazards for ESVEA";
proc gplot data=breslow01;
	plot cumhaz1*cumhaz0 line*cumhaz0/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 0.0013 by 0.0001 minor=none label=('Cumulative baseline hazard: ESVEA = 0');
	axis2 order=0 to 0.0016 by 0.0001 minor=none label=(a=90 'Cumulative baseline hazard: ESVEA = 1');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=rl c=blue;
run;

* The cumulative hazard deviates a bit from the straight line. This suggests that the assumption of proportional hazards for
  subjects with or without ESVEA is questionable;

* To examine the assumption of proportional hazards for sex we will repeat the procedure. Thus, we will first specify the 
  other covariates (ESVEA, age, and systolic blood pressure) and the fit a Cox model including the strata statement.;

data covstr;
	esvea=0; age=0; sbp=0;
run;

title "2.5.1";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea age sbp;
	strata sex;
	baseline out=breslowstr cumhaz=breslow covariates=covstr;
run;

* Now the baseline hazards for sex = 0,1 are extracted and merged into one data set called 'breslow01';

data breslow0;
	set breslowstr; 
	if sex=0; 
	a00=breslow; 
run;

data breslow1; 
	set breslowstr; 
	if sex=1; 	
	a01=breslow; 
run;

data breslow01; 
	merge breslow0 breslow1; 
	by timestrokeordeath; 
run;

* Then, the empty cells of the cumulative hazards are replaced with the previous values;

data breslow01; 
set breslow01;
	by timestrokeordeath;
	retain last1 last2;
	if a00=. then cumhaz0=last1; if a00 ne . then cumhaz0=a00; 
	if a01=. then cumhaz1=last2; if a01 ne . then cumhaz1=a01;
	output;
	last1=cumhaz0; last2=cumhaz1; 
run;

* The straight line is then added to the data set;

data breslow01; 
	set breslow01;
	line=1.7840*cumhaz0;
run;

* Finally, the cumulative baseline hazards for sex = 0,1 are plotted against each other together with the line;

title "2.5.1: Examining the assumption of proportional hazards for sex";
proc gplot data=breslow01;
	plot cumhaz1*cumhaz0 line*cumhaz0/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 0.0014 by 0.0001 minor=none label=('Cumulative baseline hazard: Sex = 0');
	axis2 order=0 to 0.0022 by 0.0001 minor=none label=(a=90 'Cumulative baseline hazard: Sex = 1');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=rl c=blue;
run;

* The curve coincides nicely with the straight line. Thus, the proportional hazard assumption seems reasonable for sex.;

*------------------------------------------- 2.5.2 --------------------------------------------------------------------------------;

* We can check the linearity assumptions on the log(hazards)-scale for age and systolic blood pressure by including non-linear effects.;

* The linearity assumption of systolic blood pressure can be investigated by including a linear spline. We select 135 as a knot 
  since this value is typically the cut-point for medication for hypertension.;
  
title "2.5.2";
data chs_data;
   set chs_data;
   if sbp > 135 then
   hypertension = sbp - 135;  
   else hypertension = 0;
run;

* Then we include this new explanatory variable in our Cox model;

proc phreg data = chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp hypertension;
run;

* We compare the effect of adding a linear spline with the Cox model from 2.4.1 using a likelihood ratio test;

title "2.5.2: Checking the linearity on the log(hazard)-scale for systolic blood pressure";
data p;
	chi2=3457.955-3457.330;
	p=1-probchi(chi2,1);
proc print;
run;

* We get a Chi-square statistic of 0.625 on 1 degree of freedom with a correpsonding p-value of 0.4292. Thi is close
  to the Wald test. Thus, we conclude that
  including a linear spline in our model does not seem to describe data better than our model from exercise 2.4.1.;

* We will test the linearity assumption of age by including a quadratic term to the Cox model.;

title "2.5.2";
data chs_data;
   set chs_data;
   age2 = age*age;
run;

proc phreg data = chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age age2 sbp;
run;

*Then, we compare this model to the Cox model forom exercise 2.4.1. using a likelihood ratio test.;

title "2.5.2: Checking the linearity on the log(hazard)-scale for age";
data p;
	chi2=3457.955-3457.923;
	p=1-probchi(chi2,1);
proc print;
run;

* We get a Chi-square statistic of 0.032 on 1 degree of freedom and a corresponding p-value of 0.85803. Thus, we do not have
  evidence against the linearity assumption of age. Note, once more, the simlarity with the Wald test.;
  
*------------------------------------------- 2.5.3 --------------------------------------------------------------------------------;

* We can investigate the proportionality assumption of sex in the piece-wise constant hazards model by including an interaction term 
  of time and sex and compare this model with the model from exercise 2.4.2 using a likelihood ratio test (directly
  iobtainable using the type3 option).;

title "2.5.3";  
proc genmod data=chs_pch_null;
	class interval (ref = '1');
	model cens=interval esvea sex age sbp sex*interval /dist=poi offset=logrisk type3;
run;


title "2.5.3: Examining the assumption of proportional hazards for sex";
data p;
	chi2=1449.6828-1448.7745;
	p=1-probchi(chi2,2);
proc print;
run;

* We get a Chi-square statistic of 0.9083 on 2 degrees of freedom with a corresponding p-value of 0.63499. We conclude 
  that we do not have enough evidence against the proportional hazards assumption for sex.;

* To investigate the proportionality assumption of ESVEA we include an interaction term of time and ESVEA and compare this model with
  the model from exercise 2.4.2 using a likelihood ratio test.;
  
proc genmod data=chs_pch_null;
	class interval (ref = '1');
	model cens=interval esvea sex age sbp esvea*interval /dist=poi offset=logrisk type3;
run;

title "2.5.3: Examining the assumption of proportional hazards for ESVEA";
data p;
	chi2=1449.6828-1443.1417;
	p=1-probchi(chi2,2);
proc print;
run;

* We get a Chi-square statistic of 6.5411 on 2 degrees of freedom with a corresponding p-value of 0.037986. This suggest that 
  the proportional hazard assumption does not hold for ESVEA on a 0.05 significance level.;

* We can investigate the linearity assumption on the log(hazard)-scale for systolic blood pressure by including a variable for 
  hypertension to the piece-wise constant hazards model. Then, we can compare this model with the model from 2.4.2 using a likelihood 
  ratio test.;

title "2.5.3";
data chs_pch_null;
   set chs_pch_null;
   if sbp > 135 then
   hypertension = sbp - 135;  
   else hypertension = 0;
run;

proc genmod data=chs_pch_null;
	class interval;
	model cens=interval esvea sex age sbp hypertension /dist=poi offset=logrisk noint;
run;

title "2.5.3: Checking the linearity on the log(hazard)-scale for systolic blood pressure";
data p;
	chi2=1449.6828-1449.0699;
	p=1-probchi(chi2,1);
proc print;
run;

* We get a Chi-square statistic of 0.6129 on 1 degree of freedom with a corresponding p-value of 0.4337 when comparing the Poisson models. 
  Thus, for both the Cox model and the Poisson model we cannot reject the hypothesis that the effect of systolic blood pressure 
  is linear.;

* We can then investigate the linearity assumption for age by including a quadratic term to our model;

title "2.5.3";
data chs_pch_null;
   set chs_pch_null;
   age2 = age*age;
run;

proc genmod data=chs_pch_null;
	class interval;
	model cens=interval esvea sex age age2 sbp /dist=poi offset=logrisk noint;
run;

* Then this model is compared with the one from exercise 2.4.2.;

title "2.5.3: Checking the linearity on the log(hazard)-scale for age";
data p;
	chi2=1449.6828-1449.6446;
	p=1-probchi(chi2,1);
proc print;
run;

* We get a Chi-square statistic of 0.0382 on 1 degree of freedom and a corresponding p-value of 0.84504. We once more reach the 
  conclusion that we do not have evidence against the linearity assumption of age.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.6 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*------------------------------------------- 2.6.1 --------------------------------------------------------------------------------;

* To examine how ESVEA affects the mortality rate after stroke we will only consider the subjects who have a stroke during the 
  study. Thus, we will create a new data frame called 'stroke_data'.;
  
* Two observations (id 24 and 368) die on the same days as they have a stroke. We will add 0.5 day to their death time such that
  it is not identical to time of stroke. Otherwise, they will not be used in the regression.;

title "2.6.1";
data stroke_data;
   set chs_data;
   if stroke = 1;
   if timestroke = timedeath then timedeath = timedeath + 0.5/365.25;
run;

* If the time variable is time since recruitment we need to take delayed entry into account. 

* We will estimate the cumulative hazard non-parametrically with the Nelson-Aalen estimator which is implemented in the phreg
  procedure. We must take delayed entry into account. Thus, both the entry time 'timestroke' and exit time 'timedeath' is included
  in the model statement.;

proc phreg data=stroke_data; 
	model (timestroke,timedeath)*death(0)=;
	strata esvea;
	baseline out=hazdata cumhaz=naa;
run;

* Then, the estimates can be plotted;

title "2.6.1: Nelson-Aalen estimate of the cumulative hazards of death after stroke";
proc gplot data=hazdata;
	plot naa*timedeath=esvea / haxis=axis1 vaxis=axis2;
	axis1 order=0 to 15.19 by 1 minor=none label=('Years');
	axis2 order=0 to 3.3 by 0.1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;

*------------------------------------------- 2.6.2 --------------------------------------------------------------------------------;

*Similarly, we fit a Cox model where both entry and exit are included in the model statement.;

title "2.6.2";
proc phreg data = stroke_data;
	model (timestroke,timedeath)*death(0)=esvea;
run;

* We get a hazard ratio of exp(0.26242) = 1.3 for ESVEA with a Chi-square statistic of 0.7282 and a corresponding p-value of 0.3935. 
  This suggests that having ESVEA does not significantly increase the rate of death after stroke when the time variable is time 
  since recruitment.; 

*------------------------------------------- 2.6.3 --------------------------------------------------------------------------------;

* In the model where the time variable is time since stroke, a new variable 'timesincestroke' is added to the stroke_data data set ;

title "2.6.3";
data stroke_data;
   set stroke_data;
   timesincestroke = timedeath - timestroke;
run;

* The Nelson Aalen estimates of the cumulative hazards are then calculated and plotted as before;

proc phreg data=stroke_data; 
	model timesincestroke*death(0)=;
	strata esvea;
	baseline out=hazdata cumhaz=naa;
run;

title "2.6.3: Nelson-Aalen estimate of the cumulative hazards of death after stroke";
proc gplot data=hazdata;
	plot naa*timesincestroke = esvea / haxis=axis1 vaxis=axis2;
	axis1 order=0 to 15.19 by 1 minor=none label=('Years');
	axis2 order=0 to 3.5 by 0.1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;


* Likewise, the hazard ratio is estimated assuming proportional hazards with the Cox model; 

title "2.6.3";
proc phreg data = stroke_data;
	model timesincestroke*death(0)=esvea;
run;

* We get a hazard ratio of exp(0.15385) = 1.1663 for ESVEA with a  χ2 statistic of 0.2458 and a corresponding p-value of 0.62. Thus,
  ESVEA does not increase the risk of death after stroke significantly when the time variable is time since stroke as well.


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.7 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We fit a Cox model to the cause-specific hazard for the outcome stroke;

title '2.7';
proc phreg data = chs_data;
	model timestrokeordeath*stroke(0)= esvea sex age sbp;
run;

* We get a hazard ratio of 2.02 and a Chi-square statistic of 6.77 on 1 degree of freedom with a corresponding p-value of 0.009. 
  This suggest that having ESVEA does change the hazard significantly when the outcome is stroke.;

* To fit a Cox model for the outcome death wihtout stroke we must first add the censoring variable death without stroke to our data 
  set. Afterwards, the model can be fitted with the phreg procedure;
  
data chs_data;
	set chs_data;
	death_wo_stroke = death;
	if stroke = 1 then death_wo_stroke = 0;
run;

proc phreg data = chs_data;
	model timestrokeordeath*death_wo_stroke(0)= esvea sex age sbp;
run;

* We get a hazard ratio of 1.17 and a Chi-square statistic of 0.7344 on 1 degree of freedom with a corresponding p-value of 0.391. 
  This suggests that the status of ESVEA does not change the hazard significantly when the outcome is death without stroke.;

* 2.7.2: It is seen that (male) sex and age are both significantly associated with increased cause-specific hazards 
  whereas ESVEA and SysBP are only associated with the rate of stroke. The coefficients for the composite end-point 
  are between those for the cause-specific hazards. Please note that the Cox model for the composite end-point and 
  those for the separate cause-specific hazards are mathematically incompatible.
