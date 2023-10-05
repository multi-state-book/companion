*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- CHAPTER 5 - EXERCISE SOLUTIONS ---------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We first load the data;

proc import out = chs_data
    datafile = 'C:\\Users\\htk287\\Desktop\\CHS\\holt.csv'
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
*---------------------------------------- EXERCISE 5.3 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* For the process of an illness-death model to be Markovian alpha13(.) must only depend on time since start of the study, 
  t. We will test the Markov assumption by adding the following time-dependent covariate in a Cox model for alpha13(.),
  
   Z_i1(t) = I(d_i(t) > 30 days),
   
   where d_i(t) = t - T_1i  and T_1i is the time of diagnosis of AF for patient i.;

* The model is fitted using the 'phreg' procedure. The only covariate in the model is 'wait30' which is designed as described above. 
  To account for delayed entry we must add the argument '/ entry = timeafib' in the model statement.;

title "5.3";
proc phreg data = chs_data;
	model timestrokeordeath*strokeordeath(0)=wait30 / entry = timeafib;
	wait30 = 0;
	if ((timestrokeordeath - timeafib) > 30/365.25) then wait30 = 1;
run;
	
* The time-dependent covariate is highly significant (p-value of 0.0016). Thus, the Markov assumption does seem to be violated, 
  i.e. there seems to be an effect of the duration of AF.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 5.4 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  5.4.1 ----------------------------------------------------------------------------;

* We will now fit separate land-mark models at times 3, 6, and 9 years for the mortality rate, including AF, stroke, ESVEA, sex, 
  age, and systolic blood pressure. AF and stroke will enter as the time-dependent covariates I(AF < s) and I(stroke < s) for s = 3, 6,
  9 years;

* We will first make a new data set called 'landmark' where we record the status of AF and stroke for all subjects at risk at the 
  three landmark times;

data landmark; set chs_data;
	if timedeath >=3 then do;
		landmark=3; entry=3;
		af = 0; if afib = 1 and timeafib <= 3 then af = 1;
		str =0; if stroke = 1 and timestroke <=3 then str=1;
	output; end;
	if timedeath >=6 then do;
		landmark=6; entry=6;
		af = 0; if afib = 1 and timeafib <= 6 then af = 1;
		str =0; if stroke = 1 and timestroke <=6 then str=1;
	output; end;
	if timedeath>=9 then do;
		landmark=9; entry=9;
		af = 0; if afib = 1 and timeafib <= 9 then af = 1;
		str =0; if stroke = 1 and timestroke <=9 then str=1;
	output; end;
run;

* We will fit all three models within one 'phreg' procedure. Therefore, we must use three different names for the 
  covariates at each land-mark time. The names are specified in the data frame 'covar541';

data cov541;
	afib3 = 0; afib6 = 0; afib9 = 0;
	stroke3 = 0; stroke6 = 0; stroke9 = 0;
	esvea3 = 0; esvea6 = 0; esvea9 = 0;
	sex3 = 0; sex6 = 0; sex9 = 0;
	age3 = 0; age6 = 0; age9 = 0;
	sbp3 = 0; sbp6 = 0; sbp9 = 0;
run;

* We must include the argument '/entry = entry' in the model statement and 'landmark' in the strata statement.;

title "5.4.1";
proc phreg data=landmark;
	model timedeath*death(0)= afib3 stroke3 esvea3 sex3 age3 sbp3
						  	  afib6 stroke6 esvea6 sex6 age6 sbp6
						  	  afib9 stroke9 esvea9 sex9 age9 sbp9 /entry=entry;
	afib3 = af*(landmark = 3); stroke3 = str*(landmark = 3); esvea3 = esvea*(landmark = 3); 
	sex3 = sex*(landmark = 3); age3 = age*(landmark = 3); sbp3 = sbp*(landmark = 3);
	afib6 = af*(landmark = 6); stroke6 = str*(landmark = 6);esvea6 = esvea*(landmark = 6); 
	sex6 = sex*(landmark = 6); age6 = age*(landmark = 6); sbp6 = sbp*(landmark = 6);
	afib9 = af*(landmark = 9); stroke9 = str*(landmark = 9); esvea9 = esvea*(landmark = 9); 
	sex9 = sex*(landmark = 9); age9 = age*(landmark = 9); sbp9 = sbp*(landmark = 9);
	strata landmark;
run;

*----------------------------------------------  5.4.2 ----------------------------------------------------------------------------;

* We will now fit land-mark ‘super-models’ where the coefficients vary smoothly among land-marks but with separate baseline hazards 
  at each land-mark. We will do as suggested in the book, m = 2 and f1(s) = (s - s1)/(sL - s1) = (s - 3)/6 and f2(s) = f1(s)^2.;

* Once again, we specify the names of the covariates in a data step;

data cov542;
	afib_ = 0; afibtime = 0; afibtime2 = 0;
	stroke_ = 0; stroketime = 0; stroketime2 = 0;
	esvea_ = 0; esveatime = 0; esveatime2 = 0;
	sex_ = 0; sextime = 0; sextime2 = 0;
	age_ = 0; agetime = 0; agetime2 = 0;
	sbp_ = 0; sbptime = 0; sbptime2 = 0;
run;

* Then, we fit the models using 'phreg' just as before.;

title "5.4.2";
proc phreg data=landmark;
	model timedeath*death(0)= afib_ afibtime afibtime2
							  stroke_ stroketime stroketime2
							  esvea_ esveatime esveatime2
							  sex_ sextime sextime2
							  age_ agetime agetime2
							  sbp_ sbptime sbptime2 / entry=entry;
	afib_ = af; afibtime = af*(landmark - 3)/6; afibtime2 = af*((landmark-3)/6)**2;
	stroke_ = str; stroketime = str*(landmark - 3)/6; stroketime2 = str*((landmark-3)/6)**2;
	esvea_ = esvea; esveatime = esvea*(landmark - 3)/6; esveatime2 = esvea*((landmark-3)/6)**2;
	sex_ = sex; sextime = sex*(landmark - 3)/6; sextime2 = sex*((landmark-3)/6)**2;
	age_ = age; agetime = age*(landmark - 3)/6; agetime2 = age*((landmark-3)/6)**2;
	sbp_ = sbp; sbptime = sbp*(landmark - 3)/6; sbptime2 = sbp*((landmark-3)/6)**2;
	strata landmark;
run;


*----------------------------------------------  5.4.3 ----------------------------------------------------------------------------;

* Finally, we will fit a land-mark ‘super-model’ where both the coefficients and the baseline hazards vary smoothly among 
  land-marks. We choose m = 2 and g_l = f_l as suggested in the book.;

* We thus add g_1 ('strtime') and  g_2 ('strtime2') to our data frame;

data cov543;
	afib_ = 0; afibtime = 0; afibtime2 = 0;
	stroke_ = 0; stroketime = 0; stroketime2 = 0;
	esvea_ = 0; esveatime = 0; esveatime2 = 0;
	sex_ = 0; sextime = 0; sextime2 = 0;
	age_ = 0; agetime = 0; agetime2 = 0;
	sbp_ = 0; sbptime = 0; sbptime2 = 0;
	strtime=0; strtime2=0;
run;

*  We once again fit the model using the 'phreg' procedure (but this time we must exclude the strata statement);

title "5.4.3";
proc phreg data=landmark;
	model timedeath*death(0)= afib_ afibtime afibtime2
							  stroke_ stroketime stroketime2
							  esvea_ esveatime esveatime2
							  sex_ sextime sextime2
							  age_ agetime agetime2
							  sbp_ sbptime sbptime2 
							  strtime strtime2 / entry=entry;
	afib_ = af; afibtime = af*(landmark - 3)/6; afibtime2 = af*((landmark-3)/6)**2;
	stroke_ = str; stroketime = str*(landmark - 3)/6; stroketime2 = str*((landmark-3)/6)**2;
	esvea_ = esvea; esveatime = esvea*(landmark - 3)/6; esveatime2 = esvea*((landmark-3)/6)**2;
	sex_ =sex; sextime = sex*(landmark - 3)/6; sextime2 = sex*((landmark-3)/6)**2;
	age_ = age; agetime = age*(landmark - 3)/6; agetime2 = age*((landmark-3)/6)**2;
	sbp_ = sbp; sbptime = sbp*(landmark - 3)/6; sbptime2 = sbp*((landmark-3)/6)**2;
	strtime=(landmark-3)/6; strtime2=strtime**2;
run;


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 5.5 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* There is currently no implementation of direct binomial regression in SAS.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 5.6 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  5.6.1 ----------------------------------------------------------------------------;

* We will use cumulative Schoenfeld residuals to test whether the assumption of time-constant hazard ratios of ESVEA, sex, age, and
  systolic blood pressure is reasonable. This can be done by adding the 'assess' statement followed by 'ph'. 
  Furthermore, '\ resample' is added to generate 1000 paths from the approximate asymptotic distribution which is used to calculate
  a goodness-of-fit test.;
  
* For the assumption of time-constant hazard ratio to hold we expect the observed Schoenfeld residuals to vary randomely around 0 
  and that the p-value is non-significant.;

title "5.6.1";
proc phreg data=chs_data;
   model timestrokeordeath*strokeordeath(0)=esvea sex age sbp;
   assess ph / resample;
run;

* The assumption of time-constant hazard ratios of ESVEA and age may not be reasonable since the p-values are borderline significant
  and the observed Schoenfeld residuals do not seem to vary randomly around 0.

* In contrast, the effect of sex and systolic blood pressure does seem to fulfill the assumption of time-constant hazard ratios, 
  since the curves of the observed Schoenfeld residuals vary randomely around 0 and the p-values are far from significant.;

*----------------------------------------------  5.6.2 ----------------------------------------------------------------------------;

* To test whether the assumption of linearity (on the log(hazard) scale) is reasonable for the effect of age and systolic blood 
  pressure we will consider the cumulative martingale residuals.;

* This can be done by adding the 'assess' statement to the 'phreg' procedure followed by the argument 'var' which specifies the
  covariates of interest. '\ resample' is added to generate 1000 paths from the approximate asymptotic distribution which is used to 
  calculate a goodness-of-fit test.;  

* Once again, we expect non-significant p-values and curves that vary randomly around 0 for the assumption of linearity to hold.;

title "5.6.2";
proc phreg data=chs_data;
   model timestrokeordeath*strokeordeath(0)=esvea sex age sbp;
   assess var=(age sbp) /  resample;
run;

* The assumption of lineraty seems to hold for both age and systolic blood pressure based on the result of both the goodness-of-fit 
  test as well as the graphical assessment.;
