*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- CHAPTER 4 - EXERCISE SOLUTIONS ---------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We must first load the data;
proc import out = chs_data
    datafile = 'data/cphholter.csv'
	dbms= csv replace;
	getnames=yes;
run;

* We will convert the time variables (timeafib, timestroke, and timedeath) from days to years;
* Furthermore, we add variables for the composite end-point of stroke or death without stroke;
data chs_data;
	set chs_data;
	timeafib = timeafib/365.25;
	timestroke = timestroke/365.25;
	timedeath = timedeath/365.25;
	timestrokeordeath = timedeath;
	if stroke = 1 then timestrokeordeath = timestroke;
	strokeordeath = death;
	if stroke = 1 then strokeordeath = 1;
run;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.1 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We estimate the Kaplan-Meier survival function for subjects with or without ESVEA with the phreg procedure where 'esvea' is added
  in the strata statement. The result is saved as 'survdat'.;

title "4.1: Stroke-free survival probabilities estimated with the Kaplan-Meier estimator";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=;
	strata esvea;
	baseline out=survdat survival=km;
run;

* Then the estimates are plotted using the gplot procedure;

proc gplot data=survdat;
plot km*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 label=('Years');
	axis2 order=0 to 1 by 0.1 label=(a=90 'Survival probability');
	symbol1 i=stepjl c=red;
	symbol2 i=stepjl c=blue;
run;
quit;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.2 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------- 4.2.1 ----------------------------------------------------------------------------;

* To estimate the stroke-free survival functions for a 65-year old woman with a systolic blood pressure of 150mmHg with or without
  ESVEA we will first create a data frame 'cov' with the desired values of the covariate.;

data cov;
	esvea = 0; sex = 0; age = 65; sbp = 150; output;
	esvea = 1; sex = 0; age = 65; sbp = 150; output;
run;

* Then, a Cox model including ESVEA, sex, age, and systolic blood pressure is fitted with the phreg procedure and the stroke-free
  survival functions for subjects with values according to 'cov' are saved as 'survdata'.;

title "4.2: Stroke-free survival for a 65-year old woman with sbp = 150mmHg";
proc phreg data=chs_data;
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp;
	baseline out=survdata survival=surv covariates = cov;
run;

* Finally, the survival functions are plotted using the gplot procedure;

proc gplot data=survdata;
	plot surv*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 label=('Years');
	axis2 order=0 to 1 by 0.1 label=(a=90 'Stroke-free survival probability');
	symbol1  i=stepjl c=blue;
	symbol2  i=stepjl c=red;
run;
quit;

*----------------------------------------------- 4.2.2 ----------------------------------------------------------------------------;

* A Cox model including ESVEA, sex, age, and systolic blood pressure is fitted and 'diradj group = esvea' is added to obtain the
  predicted survival functions for patients with or without ESVEA using the g-formula. The data is saved as 'gsurv'.;

title "4.2: Cox model for the outcome stroke-free survival including ESVEA, sex, age, and systolic blood pressure";
proc phreg data=chs_data;
	class esvea (ref = '0');
	model timestrokeordeath*strokeordeath(0)=esvea sex age sbp;
	baseline out=gsurv survival=surv / diradj group=esvea;
run;

* The survival functions are then plotted using the gplot procedure;

title "4.2: Stroke-free survival probabilities estimated using the G-formula";
proc gplot data=gsurv;
	plot surv*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Estimated survival function (g-formula)');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;
quit;
*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.3 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;


* We will estimate the 3-year restricted mean time survival to the composite end-point strokke or death including ESVEA, sex, age, 
  and systolic blood pressure using the rmstreg procedure. We specify 'tau = 3' in the rmstreg statement to obtain a 3 year time 
  limit and 'link = linear' in the model statement to get a linear model. NB: requires SAS STAT 15.1;

title "4.3";
/* please note: method=ipcw(strata=esvea) to match rmst2() in R */
proc rmstreg data=chs_data tau=3;
   model timestrokeordeath*strokeordeath(0)=esvea sex age sbp /
			link=linear method=ipcw(strata=esvea);
run;

* Thus, we obtain the following model for the 3-year restricted mean time to the composite end-point stroke or death
  epsilon(3|Z) = 3.3552 - 0.0274*Z1 - 0.0540*Z2 - 0.0073*Z3 + 0.0005*Z4, where (Z1,Z2,Z3,Z4) are ESVEA, age, sex, and systolic blood 
  pressure;

* We will also present the non-parametric estimates. We restrict the data set at tau=3.  NB: requires SAS STAT 15.1;

proc lifetest data=chs_data rmst(tau=3);
time timestrokeordeath*strokeordeath(0);
strata esvea;
run;
  
*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.4 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We must first create a variable for the competing risks which we will call 'event'. 0 is censored, 1 is stroke, and 2 is death 
  without stroke;

data chs_data;
	set chs_data;
	death_wo_stroke = death;
	if stroke = 1 then death_wo_stroke = 0;
	event = 0;
	if stroke = 1 then event = 1;
	if death_wo_stroke = 1 then event = 2;
run;

* Then we will fit a Cox model returning the predicted cumulative incidence functions with the specified covariates. 
  This is done by adding the argument eventcode(cox) to the model statement and adding the 'cif' argument in the baseline 
  statement.  NB: requires SAS STAT 15.1;

proc phreg data = chs_data noprint; 
	model timestrokeordeath*event(0) =esvea sex age sbp / eventcode(cox) = 1;
	baseline covariates = cov out=cif44_stroke cif = cif;
run;

* Finally, the cumulative incidence functions are plotted using the gplot procedure;

title '4.4: CIF for the outcome stroke (based on Cox model)';
proc gplot data=cif44_stroke;
	plot cif*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 label=('Years');
	axis2 order=0 to 0.2 by 0.02 label=(a=90 'CIF for stroke');
	symbol1  i=stepjl c=blue;
	symbol2  i=stepjl c=red;
run;

* Then, we repeat the procedure for the outcome death without stroke;

proc phreg data = chs_data noprint; 
	model timestrokeordeath*event(0) =esvea sex age sbp / eventcode(cox) = 2;
	baseline covariates = cov out=cif44_death cif = cif;
run;

* We can now plot the cumulative incidence functions gpt death without stroke using the gplot procedure;

title '4.4: CIF for the outcome death without stroke (based on Cox model)';
proc gplot data=cif44_death;
	plot cif*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 label=('Years');
	axis2 order=0 to 0.3 by 0.03 label=(a=90 'CIF for death w/o stroke');
	symbol1  i=stepjl c=blue;
	symbol2  i=stepjl c=red;
run;
quit;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.5 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------- 4.5.1 ----------------------------------------------------------------------------;

* To obtain the CIF for the outcomes stroke and death without stroke using the Fine-Gray model, 'eventcode' is added instead of
  'eventcode(cox)' in the model statement of the phreg procedure;

* We will first estimate the CIF for the outcome stroke;

proc phreg data = chs_data; 
	model timestrokeordeath*event(0) =esvea sex age sbp / eventcode = 1;
	baseline covariates = cov out=cif451_stroke cif = cif;
run;

title '4.5.1: CIF for the outcome stroke (based on Fine-Gray model)';
proc gplot data=cif451_stroke;
	plot cif*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 label=('Years');
	axis2 order=0 to 0.125 by 0.0125 label=(a=90 'CIF for stroke');
	symbol1  i=stepjl c=blue;
	symbol2  i=stepjl c=red;
run;

* Then we will estimate the CIF for the outcome death without stroke;

proc phreg data = chs_data; 
	model timestrokeordeath*event(0) =esvea sex age sbp / eventcode = 2;
	baseline covariates = cov out=cif451_death cif = cif;
run;

title '4.5.1: CIF for the outcome stroke (based on Fine-Gray model)';
proc gplot data=cif451_death;
	plot cif*timestrokeordeath=esvea/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 label=('Years');
	axis2 order=0 to 0.25 by 0.025 label=(a=90 'CIF for death w/o stroke');
	symbol1  i=stepjl c=blue;
	symbol2  i=stepjl c=red;
run;
quit;

*----------------------------------------------- 4.5.2 ----------------------------------------------------------------------------;

* We will then make a macro function taking the cause of interest (stroke or death without stroke) and value of ESVEA as arguments;

%macro cif452(cause, esvea);
* Creating modified covariates. Observed values for sex, age, and sbp, while ESVEA is either 0 or 1 for all subjects;
data covar_temp;
	set chs_data;
	esvea = &esvea;
	keep id esvea sex age sbp;
run;
* Fine-Gray model returning the predicted cumulative incidence functions for the modified covariates for all patients;
proc phreg data = chs_data noprint; 
	model timestrokeordeath*event(0) =esvea sex age sbp / eventcode = &cause;
	baseline covariates = covar_temp out=ciftest cif = cif;
run;
* Splitting the data set ciftest and returning one data set per patient containing the predicted cumulative incidence function;
%do i = 1 %to 678;
	data cif&i;
		set ciftest;
		if id = &i;
		cif&i = cif;
		keep id timestrokeordeath cif&i;  
	run;	
%end;
* Calculating the average of the cumulative incidence functions;
data res&esvea;
	set cif1-cif678 ;
	merge cif1-cif678;
	by timestrokeordeath;
	ESVEA&esvea = mean(of cif1-cif678);
	keep timestrokeordeath ESVEA&esvea;
run;
%mend;

*We call the function for the outcome stroke and ESVEA = 0,1;

%cif452(1,1);
%cif452(1,0);

* Then we create a data set containing the cumulative incidence function for both ESVEA = 0 and ESVEA = 1 and plot the result.;

data stroke452;
	set res0 res1;
	merge res0 res1;
	by timestrokeordeath;
run;

proc sgplot data=stroke452;
title1'4.5.2 - Cumulative incidence function for stroke predicted with Fine-Gray models using the g-formula';
   step y=ESVEA0 x=timestrokeordeath;
   step y=ESVEA1 x=timestrokeordeath;
   xaxis label= "Time (Years since randomization)";
   yaxis label= "Cumulative incidence";
run;

* We repeat for the outcome death without stroke. First we call our macro function;

%cif452(2,1);
%cif452(2,0);

*Then the data is collected in one data set and afterwards we plot the result;

data death452;
	set res0 res1;
	merge res0 res1;
	by timestrokeordeath;
run;

proc sgplot data=death452;
title1'4.5.2 - Cumulative incidence function for death predicted with Fine-Gray models using the g-formula';
   step y=ESVEA0 x=timestrokeordeath;
   step y=ESVEA1 x=timestrokeordeath;
   xaxis label= "Time (Years since randomization)";
   yaxis label= "Cumulative incidence";
run;


/* For comparison, we also estimate the cumulative incidences non-parametrically. */

proc lifetest data=chs_data plots=(cif);
time timestrokeordeath*event(0)/eventcode=1;
strata esvea;
run;


proc lifetest data=chs_data plots=(cif);
time timestrokeordeath*event(0)/eventcode=2;
strata esvea;
run;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.6 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* There is currently not an implementation of a procedure for linear models of the number of years lost due to a specific cause in 
  SAS.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.7 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  4.7.1 ----------------------------------------------------------------------------;

* The prevalence of AF is given by Q_1 / (Q_0 + Q_1). Thus, we must first obtain estimates of Q_0 and Q_1.

* We will first fill in the empty spots of timeafib and replace paths where AF happens after stroke;
* We will also add an censoring variable for leaving state 0 called 'outof0'.;

title '4.7.1';
data chs_data;
	set chs_data;
	if afib = 1 and timeafib > timestrokeordeath then afib = 0;
	if afib = 0 then timeafib = timestrokeordeath;
	outof0 = 0;
	if afib = 1 or strokeordeath = 1 then outof0 = 1;
run;


* Then, we will fit a model for being in state 0 by using 'timeafib' as our time variable and 'outof0' as our censoring variable;

proc phreg data=chs_data; /* Q0(t) */
	model timeafib*outof0(0)=;
	baseline out=q0 survival=km0;
run;

* We will also fit a model for being in state 0 or state 1 by using 'timestrokeordeath' as our timevariable and 'strokeordeath' 
  as our censoring variable.;

proc phreg data=chs_data; /* Q0(t) + Q1(t) */
	model timestrokeordeath*strokeordeath(0)=;
	baseline out=q0andq1 survival = km01;
run;


* We will now estimate the probability of being in state 1 as the difference between the two models specified above. We must 
  first merge the data frames 'q0' and 'q0andq1' by the joint 'time' variable and fill the empty cells for the survival 
  probabilities with the last observed value. Then, 'q1' and the prevalence 'prev' is added to the dataframe 'allrev'.;

data q0; set q0; time=timeafib; run;
data q0andq1; set q0andq1; time=timestrokeordeath; run;
data all; merge q0 q0andq1; by time; run;

data allrev; 
set all;
	by time;
	retain last1 last2;
	if km0=. then q0=last1; if km0 ne . then q0=km0; *AF-free survival, Q0; 
	if km01=. then q01=last2; if km01 ne . then q01=km01; *Q0 + Q1;
	q1 = q01 - q0;
	prev = q1/q01;
	output;
	last1=q0; last2=q01; 
run;

*Finally, we plot the result;

proc gplot data=allrev;
	title'4.7.1: Prevalence of AF';
	plot prev*time / haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 5 minor=none label=('Years');
	axis2 order=0 to 0.12 by 0.02 minor=none label=(a=90 'Prevalence of AF');
	symbol  v=none i=stepjl c=blue;
run;

*----------------------------------------------  4.7.2 ----------------------------------------------------------------------------;


* To estimate the expected length of stay in state 0 and state 1 up till 3 years we will integrate the functions for being in state 
  0 or 1 from 0 t0 3 years. 

* We must first add the time point for 3 years to the data set 'allrev'.;

data end_point; time = 3; run;
data allrev; merge allrev end_point; by time; run;

* Then we will calculate the product of the length of each time period and the value of q0 or q1 and then sum these products to 
  obtain estimates of the expected length of state in state 0 or 1.;

data allrev;
	set allrev;
	retain elos0 elos1;
	dq0 = dif(time)*lag(q0);
	dq1 = dif(time)*lag(q1);
	elos0 + dq0;
	elos1 + dq1;
run;

* Finally, we will print the result;

title'4.7.2';
proc print data = allrev;
	var time elos0 elos1;
	where time = 3;
run;

*----------------------------------------------  4.7.3 ----------------------------------------------------------------------------;

* We will first make our bootstrap data frames. We will make 200 bootstrap samples and each bootstrap sample contains 678 rows which
  are sampled with replacement from our originqal data.;

title'4.7.3';
data boot_chs;
	do sampnum = 1 to 200; /* nboot=200*/
	do i = 1 to 678; /*nobs=678*/
	x=round(ranuni(0)*678); /*nobs=678*/
	set chs_data
	point=x;
	output;
	end;
	end;
	stop;
run;

* Then we will estimate the expected length of stay in state 0 before 3 years based on each of the 200 bootstrap samples;

* We will first estimate the probability of being in state 0 and the probability of being in either state 0 or state 1.;

proc phreg data=boot_chs noprint; /* Q0(t) */
	by sampnum;
	model timeafib*outof0(0)=;
	baseline out=q0_boot survival=km0;
run;

proc phreg data=boot_chs noprint; /* Q0(t) + Q1(t) */
	by sampnum;
	model timestrokeordeath*strokeordeath(0)=;
	baseline out=q0andq1_boot survival = km01;
run;

* Then, we will merge the two datasets by time and add the time point for 3 years to all bootstrap samples.;

data q0_boot; set q0_boot; time=timeafib; run;
data q0andq1_boot; set q0andq1_boot; time=timestrokeordeath; run;

data end_point; do sampnum = 1 to 200; time = 3; output; end; run;
data all_boot; merge q0_boot q0andq1_boot end_point; by sampnum time; run;


* Then, we will estimate q1 as the difference between the probability of being in either state 0 or state 1 and the 
  probability of being in state 0.;

data allrev_boot; 
set all_boot;
	by sampnum time;
	retain last1 last2;
	if km0=. then q0=last1; if km0 ne . then q0=km0; *AF-free survival, Q0; 
	if km01=. then q01=last2; if km01 ne . then q01=km01; *Q0 + Q1;
	q1 = q01 - q0;
	output;
	last1=q0; last2=q01; 
run;

* The expected length of stay in state 0 or 1 is then estimated as the sum of the products of the length of the time intervals 
  and the value of q0 or q1;

data allrev_boot;
	set allrev_boot;
	if time > 3 then delete;
	retain elos0 elos1;
	dq0 = dif(time)*lag(q0);
	dq1 = dif(time)*lag(q1);
	if time = 0 then do; elos0 = 0; elos1 = 0; dq0 = .; dq1 = .; end;
	elos0 + dq0;
	elos1 + dq1;
run;

* We will only keep the result for the value of the expected lengths of stay up at 3 years in the data frame 'elos_est';

data elos_est;
	set allrev_boot;
	keep time elos0 elos1;
	where time = 3;
run;

* Finally, the mean and standard deviation is calculated;

proc means data=elos_est stddev mean;
	var elos0 elos1;
run;


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.10 ---------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  4.10.1 ---------------------------------------------------------------------------;

*We will estimate the censoring distribution G(t) with the Kaplan-Meier estimator where 'censoring' is the event and 'death' acts 
 as censoring.;


title "4.10.1";
proc phreg data=chs_data;
	model timedeath*death(1)=;
	baseline out=survdat survival=km;
run;

* Then the estimates are plotted using the gplot procedure;

proc gplot data=survdat;
plot km*timedeath/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 16 by 2 label=('Years');
	axis2 order=0 to 1 by 0.1 label=(a=90 'Censoring probability');
run;
quit;

*----------------------------------------------  4.10.2 ---------------------------------------------------------------------------;

*To examine to what extent the censoring distribution depends on the variables ESVEA, sex, age, and systolic blood pressure we will
 fit a Cox model including these covariates and where 'censoring' is the event and 'death' acts as censoring.;

title "4.10.2";
proc phreg data=chs_data;
	model timedeath*death(1)=esvea sex age sbp;
run;

*The p-values indicate that censoring depends strongly on ESVEA, sex, and age but not on systolic blood pressure.;


*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;
