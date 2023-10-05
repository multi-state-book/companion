*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- CHAPTER 3 - EXERCISE SOLUTIONS ---------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We first load the data;

proc import out = chs_data
    datafile = 'C:\HRfinal\holter\cphholter.csv'
	dbms= csv replace;
	getnames=yes;
run;

* We will convert the time variables, timeafib, timestroke, and timedeath, from days to years;

data chs_data;
	set chs_data;
	timeafib = timeafib/365.25;
	timestroke = timestroke/365.25;
	timedeath = timedeath/365.25;
run;

* Furthermore, we add variables for the composite end-point of stroke or death without stroke;

data chs_data;
	set chs_data;
	timestrokeordeath = timedeath;
	if stroke = 1 then timestrokeordeath = timestroke;
	strokeordeath = death;
	if stroke = 1 then strokeordeath = 1;
run;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 3.6 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We fit the Cox proportional hazards model for stroke-free survival including the covariates ESVEA, sex, age, and 
  systolic blood pressure like we did in exercise 2.4.1.;
  
title "3.6: Cox model from exercise 2.4.1";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp;
run;

* To investigate if the covariates may be described as time-constant hazard ratios, we will compare our model from exercise 2.4.1 
  with models where one of the covariates, Z, is assumed time-dependent, i.e. Z(t) = Z*f(t) for some function t of time.;

* We will investigate the assumption of a time-constant effect of ESVEA by fitting a Cox model with a time-dependent effect of ESVEA 
  described by Z(t) = x1*I(t > 5 years) + x2. Multiple other forms of the function could be chosen but we will only illustrate one 
  choice of f(t) for each covariate.;

title "3.6: Cox model with time-dependent effect of ESVEA";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp esveat0;
	esveat0 = (esvea=1)*(timestrokeordeath > 5);
run;

* We estimate the coefficient of the time-dependent effect of ESVEA to be Z(t) = 0.7788 - 0.0734*I(t > 5). The model thus predicts a decrease 
  in the effect of having ESVEA on the stroke-free survival hazard after 5 years.;

* We can then compare the model with the one from exercise 2.4.1 using a likelihood ratio test. This is done by comparing the likelihood 
  scores for the two models and then calculating the p-value using the Chi-square distribution with 1 degrees of freedom. The likelihood 
  scores are found in the 'Model Fit Statistics' table;

title "3.6: Likelihood ratio test investigating the assumption of a time-constant effect of ESVEA";
data p;
	chi2=3457.955-3452.967;
	p=1-probchi(chi2,1);
proc print;
run;

*  We get a Chi-square statistic of 4.988 on 1 degree of freedom and a correpsonding p-value of 0.0255. Thus, including ESVEA as a time-
   dependent covariate is significant on a 0.05 significance level. The same conclusion appears when studying the Wald test.; 

* Likewise, we can investigate the assumption of time-constant effect of sex by fitting a Cox model with Z(t) = x1*log(t) + x2.; 

title "3.6: Cox model with a time-dependent effect of sex";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp sexlogtime;
	sexlogtime = (sex=1)*log(timestrokeordeath);
run;

* We estimate the coefficient of the time-dependent effect of sex to be Z(t) = -0.0536*log(t) + 0.6783.;

*Then, we can compare this model with the one from exercise 2.4.1 with a likelihood ratio test.;

title "3.6: Likelihood ratio test investigating the assumption of a time-constant effect of sex";
data p;
	chi2=3457.955-3457.836;
	p=1-probchi(chi2,1);
proc print;
run;

* We get a Chi-square statistic of 0.119 on 1 degree of freedom and a correpsonding p-value of 0.7301. We conclude that we do not have 
  evidence against the assumption of a time-constant hazard ratio for sex. Same conclusion for the Wald test.;

* We will check the assumption of a time-constant effect of age by including an indicator of whether the current age is greater than 70, 
  i.e. Z(t) = x1 + x2*I((age + t) > 70), in our Cox model.;

title "3.6: Cox model with a time-dependent effect of age";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp aget0;
	aget0 = ((age + timestrokeordeath) > 70);
run;

* We estimate the coefficient of the time-dependent effect of age to be  Z(t) = 0.0585  +  0.3753*I(age + t > 70).;

* Then, we can compare this model with the one from exercise 2.4.1 using a likelihood ratio test.;

data p;
	chi2=3457.955-3455.135;
	p=1-probchi(chi2,1);
proc print;
run;

* We get a Chi-square statistic of 2.82 on 1 degree of freedom and a correpsonding p-value of 0.0931. Thus, our inclusion of age as a
  time-dependent covariate does not provide evidence against the assumption of a time-constant hazard ratio of age.
  Same conclusion for the Wald test.;

* Lastly, we will check the assumption of time-constant effect of systolic blood pressure. We will investigate a model where 
  Z(t) = x1* + x2.;

proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp sbptime;
	sbptime = sbp*timestrokeordeath;
run;

* We get an estimated time-dependent coefficient for the effect of systolic blood pressure of Z(t) = 5.32*10^{-3} - 2.11*10^{-5}*t.;

* Then, we can compare this model with the one from exercise 2.4.1 using a likelihood ratio test.;

data p;
	chi2=3457.955-3457.954;
	p=1-probchi(chi2,1);
proc print;
run;

* We get a Chi-square statistic of 0.001 on 1 degree of freedom and a correpsonding p-value of 0.9748.Thus, we do not have evidence 
  against the assumption of a time-constant effect of systolic blood pressure.  Same conclusion for the Wald test.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 3.7 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* Since the following exercises investigates the effect of AF on stroke-free survival, we will neglect diagnoses of AF 
  after a stroke. ;

data chs_data;
	set chs_data;
	if afib = 1 and timeafib > timestrokeordeath then afib = 0;
run;

* We will include a time-dependent covariate I(AF <= t) to our Cox model.;

title "3.7: Cox model including the time-dependent covariate I(AF < t)";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0) = esvea sex age sbp timeaf;
	timeaf = 0;
	if afib = 1 and timeafib <= timestrokeordeath then timeaf=1;
run;

* We get a hazard ratio of 1.254 for ESVEA with a Chi-square statistic of 2.1277 on 1 degree of freedom and a correpsonding p-value of 
  0.1447. Thus, having ESVEA is associated with shorter stroke-free survival. However, this effect is not significant!; 

* The effect of being diagnosed with AF before time t is highly significant and associated with a much smaller rate of
  stroke-free survival compared to those not diagnosed with AF at time t.;


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 3.8 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We will include AF as a time-fixed covariate to our Cox model which will lead to so-called immortal time bias.;

title "3.8: Cox model including AF as a time-fixed covariate";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)= esvea sex age sbp afib;
run;

* We get a hazard ratio of 0.858 for being diagnosed with AF. The effect of being diagnosed with AF during the study is associated 
  with a slightly shorter expected stroke-free survival compared to those not diagnosed. However we get a Chi-square statistic of 
  0.06751 on 1 degree and the p-value is 0.4113. Thus, the effect is not significant.;

* This contradicts the result from exercise 3.7 where AF was (correctly) included as an time-dependent covariate.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 3.9 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------- 3.9.1 ----------------------------------------------------------------------------;

* We will fit two separate Cox models for 0->2/3 and 1->2/3 transitions including ESVEA, sex, age, and systolic blood pressure using 
  the procedure described in section 3.8. Thus, we will first make the duplicated data set.

  We must first add two new time variables, 'timeoutof0' and 'timeoutof1' which will indicate time last seen in state 0 and state 1,
  respectively, and a new censoring variable, 'strokeordeath0' indicating whether a subject transitioned to state 2 or 3 directly 
  from state 0.;

data chs_data;
	set chs_data;
	timeoutof0 = timestrokeordeath;
	if afib = 1 then timeoutof0 = timeafib;
	strokeordeath0 = strokeordeath;
	if afib = 1 then strokeordeath0 = 0;
	timeoutof1 = timestrokeordeath;
	if afib = 1 and timeafib = timestrokeordeath then timeoutof1 = timestrokeordeath + 0.5/365.25;
run;


* Then, we can make the duplicated data set. All potential 0->2/3 transitions will be represented by a row where 'stratum' is 0
  whereas potential 1->2/3 transitions will have 'stratum' equal to 1. 'entrytime' and 'time' will denote the time of entry 
  and the last time seen in state 0 or 1 while 'status' indicates whether a stroke or death occurred at time 'time'. Furthermore, 
  the covariates associated with the 0->2/3 transition are denoted 'esvea0', 'sex0', age0', and 'sbp0' while the covariates 
  associated with the 1->2/3 transition are denoted 'esvea1', 'sex1', 'age1', and 'sbp1'.;

data double0;
	set chs_data;
	entrytime = 0;
	time = timeoutof0;
	status = strokeordeath0;
	esvea0 = esvea;
	sex0 = sex;
	age0 = age;
	sbp0 = sbp;
	esvea1 = 0;
	sex1 = 0;
	age1 = 0;
	sbp1 = 0;
	stratum = 0;
run;

data double1;
	set chs_data;
	if afib = 1;
	entrytime = timeafib;
	time = timeoutof1;
	status = strokeordeath;
	esvea0 = 0;
	sex0 = 0;
	age0 = 0;
	sbp0 = 0;
	esvea1 = esvea;
	sex1 = sex;
	age1 = age;
	sbp1 = sbp;
	stratum = 1;
run;

data double;
	set double0 double1;
run;

*Finally, we will use **phreg** to fit the two separate Cox models. 'entrytime', 'time' and 'status' are given to the left of '=' and
and 'esvea0', 'sex0', 'age0', 'sbp0', 'esvea1', 'sex1', 'age1', and 'sbp1' are included on the right side of '='. Furthermore, 'stratum' 
is given in the strata statement;

title "3.9.1: Separate Cox models for stroke-free survival for subjects
with or without AF";
proc phreg data = double;
	strata stratum;
	model (entrytime, time)*status(0)= esvea0 sex0 age0 sbp0 esvea1 sex1 age1 sbp1; 
run;


*----------------------------------------------- 3.9.2 ----------------------------------------------------------------------------;

* We will examine whether a model that assumes the same covariate effect for the two transitions describes our data just as good as 
  the model we fitted in exercise 3.9.1. ;

* Thus, we will fit a model where 'esvea0', 'sex0', 'age0', 'sbp0', 'esvea1', 'sex1', 'age1', and 'sbp1' are replaced with 'esvea', 
  'sex', 'age', and 'sbp'.;

title "3.9.2: Cox model assuming the same effect of the covariates for 0->2/3 and 1->2/3";
proc phreg data = double;
	strata stratum;
	model (entrytime, time)*status(0) = esvea sex age sbp;
run;

* Then we can compare this model with the one from exercise 3.9.1 with a likelihood ratio test using the scores from the tables 
  'Model Fit Statistics'.;

title "3.9.2: Examining the assumption of same covariate effect";
data p;
	chi2=3253.556-3248.745;
	p=1-probchi(chi2,4);
proc print;
run;

* We get a Chi-square statistics of 4.811 on 4 degrees of freedom and a corresponding p-value of 0.3072. Thus, it seems like the 
  new model with common regression coefficients describes our data as good as the one from exercise 3.9.1.;

* We will now examine the assumption of proportional hazards. We must extract the baseline hazards for subjects with or without AF 
  in the model where the hazards are assumed non-proportional. Thus, we fit a model from exercise 3.9.1 once more and specify that 
  all covariates should be 0.;  

data covstr;
	esvea0 = 0; sex0 = 0; age0 = 0; sbp0 = 0;
	esvea1 = 0; sex1 = 0; age1 = 0; sbp1 = 0;
run;

proc phreg data = double;
	strata stratum;
	model (entrytime, time)*status(0)= esvea0 sex0 age0 sbp0 esvea1 sex1 age1 sbp1;
	baseline out=hazdata cumhaz = breslow covariates = covstr;
run;

data breslow023;
	set hazdata; 
	if stratum=0; 
	a023=breslow; 
run;

data breslow123;
	set hazdata;
	if stratum=1;
	a123=breslow;
run;

data breslow; 
	merge breslow023 breslow123; 
	by time; 
run;

* Then, the empty cells of the cumulative hazards are replaced with the previous values;

data breslow; 
set breslow;
	by time;
	retain last1 last2;
	if a023=. then cumhaz0=last1; if a023 ne . then cumhaz0=a023; 
	if a123=. then cumhaz1=last2; if a123 ne . then cumhaz1=a123;
	output;
	last1=cumhaz0; last2=cumhaz1; 
run;

* Then, we must fit a Cox model where the hazards for subjects with or without AF are assumed proportional. This is done by 
  including AF as an time-dependent covariate like we did in exercise 3.7.;

proc phreg data = double;
	model (entrytime, time)*status(0) = esvea0 sex0 age0 sbp0 esvea1 sex1 age1 sbp1 timedepaf;
	timedepaf = 0;
	if afib = 1 and time > timeafib then timedepaf = 1;
run;

* Lastly, the predicted hazard ratio for the time-dependent covariate is added to the data as the slope of a straigth line 
  through (0,0);

data breslow; 
	set breslow;
	line=1.222*cumhaz0;
run;


* Finally, we can check the assumption of proportional hazards by plotting the hazards from the model assuming non-
  proportional hazards against each other together with the straigth line through (0,0) with a slope equal to the hazard ratio of 
  AF.;

title "3.9.2: Examining the assumption of proportional hazards";
proc gplot data=breslow;
	plot cumhaz1*cumhaz0 line*cumhaz0/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 0.0014 by 0.0001 minor=none label=('Cumulative baseline hazard for subjects without AF');
	axis2 order=0 to 0.0022 by 0.0001 minor=none label=(a=90 'Cumulative baseline hazard for subjects with AF');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=rl c=blue;
run;

*The model assuming proportional hazards does not seem to fit the data well, since the straigth blue line does not coincide nicely 
 with the red curve.; 

*The proportional hazards assumption may also be investigated using a time-dependent covariate: duration since af.;

proc phreg data = double;
	model (entrytime, time)*status(0) = esvea0 sex0 age0 sbp0 esvea1 sex1 age1 sbp1 timedepaf duration;
	timedepaf = 0;
	if afib = 1 and time > timeafib then timedepaf = 1;
	duration=0;
	if timedepaf = 1 then duration = time - timeafib;
run;

*This, once more, suggests borderline deviations from proportional hazards.;
