proc genmod data=pbc3mult;
  class tment (ref='0') interv;
  model fail= tment interv alb  bili / dist=poi offset=logrisk ;
run;
*------------------------------------------------------------------;
*------- Chapter 2, SAS code, PBC3 data ---------------------------;
*------------------------------------------------------------------;

* Load pbc data; 
proc import out=pbc3
	datafile="data/pbc3.csv"
	dbms=csv replace;
run;

*---------------------------------------------------------------;
*--------------------- Figure 2.1 ------------------------------;
*---------------------------------------------------------------;
 
proc phreg data=pbc3; 
	model days*status(0)=;
	strata tment;
	baseline out=hazdat cumhaz=naa stdcumhaz=sdnaa;
run;
data hazdat; 
	set hazdat; 
	daysyears = days/365.25; 
run;
proc gplot data=hazdat;
	plot naa*daysyears=tment/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 0.7 by 0.1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;

* In-text: Nelson-Aalen estimates at 2 years;
proc print data=hazdat;
	where 1.5<daysyears<=2; 
run;

* Logrank for tment; 
proc lifetest data=pbc3 notable; 
	time days*status(0);
	strata tment;
run; 

*---------------------------------------------------------------;
*--------------------- Table 2.1 -------------------------------;
*---------------------------------------------------------------; 

/* prepare for poisson model */

data pbc3mult; 
	set pbc3;
	fail=(days<=2 * 365.25)*(status ne 0);
	risktime=min(2 * 365.25, days);
	logrisk=log(risktime); interv=1; output;  
	if days>2 * 365.25 then do;
	fail=(days<=4 * 365.25)*(status ne 0);
	risktime=min(2 * 365.25 ,days-2 * 365.25);
	logrisk=log(risktime); interv=2; output; end;
	if days>4 * 365.25  then do;
	fail=status ne 0; 
	risktime=days-4 * 365.25;
	logrisk=log(risktime); interv=3; output; end;
run;


proc means sum  data=pbc3mult;
	class tment interv;
	var fail risktime;
run;


proc genmod data=pbc3mult;
	where tment=1;
	class interv;
	model fail=interv/dist=poi offset=logrisk;
	estimate '0-2 years' intercept 1  interv 1 0 0/exp;
	estimate '2-4 years' intercept 1  interv 0 1 0/exp;
	estimate '4-6 years' intercept 1  interv 0 0 1/exp;
run;

proc genmod data=pbc3mult;
	where tment=0;
	class interv;
	model fail=interv/dist=poi offset=logrisk;
	estimate '0-2 years' intercept 1  interv 1 0 0/exp;
	estimate '2-4 years' intercept 1  interv 0 1 0/exp;
	estimate '4-6 years' intercept 1  interv 0 0 1/exp;
run;

ods output modelfit=full;
proc genmod data=pbc3mult;
	class interv tment;
	model fail=interv|tment /dist=poi offset=logrisk;
run;
data full;
  set full;
	 where Criterion="Deviance";
	 full_like=value;
	 full_df=df;
	 keep full_like full_df;
run;
proc print;run;
ods output modelfit=reduced;
proc genmod data=pbc3mult;
	class interv tment;
	model fail=interv /dist=poi offset=logrisk;
run;
data reduced;
  set reduced;
	 where Criterion="Deviance";
	 reduced_like=value;
	 reduced_df=df;
	 keep reduced_like reduced_df;
run;
data lrt_pval;
 merge full reduced;
  lrt = reduced_like - full_like;
	df  = reduced_df - full_df;
  p_value = 1 - probchi(lrt,df);
run;
proc print data=lrt_pval;
  title "LR test statistic and p-value";
run;
*---------------------------------------------------------------;
*--------------------- Figure 2.3 ------------------------------;
*---------------------------------------------------------------; 
 
data hazdat; 
	set hazdat;
	if days<=2 * 365.25 then 
		pwch=days*(27.0000000/104856*(1-tment)+24.0000000/107931.5*tment);

	if 2 * 365.25 <days<=4 * 365.25 then
		pwch=2* 365.25*(27.0000000/104856*(1-tment)+24.0000000/107931.5*tment)
		+(days-2* 365.25)*(17.0000000/49673*(1-tment)+18.0000000/50284*tment);

	if 4 * 365.25 <days then
		pwch=2* 365.25 *(27.0000000/104856*(1-tment)+24.0000000/107931.5*tment)
		+(2* 365.25)*(17.0000000/49673*(1-tment)+18.0000000/50284*tment)
		+(days-4* 365.25)*(2.0000000/8642*(1-tment)+2.0000000/7599*tment);
	daysyears = days/365.25; 
run; 


proc gplot data=hazdat; where tment=1;
	plot (naa pwch)*daysyears/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 0.7 by 0.1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl r=1 c=red;
	symbol2 v=none i=join r=1 c=blue;
run;


*---------------------------------------------------------------;
*--------------------- In-text Cox model -----------------------;
*---------------------------------------------------------------; 
 
proc phreg data=pbc3;
	model days*status(0)=tment/rl;
run;


*---------------------------------------------------------------;
*--------------------- Figure 2.5 ------------------------------;
*---------------------------------------------------------------; 
 
data cov;
	tment=0;
run;
proc phreg data=pbc3;
	model days*status(0)=tment/rl;
	baseline out=breslow cumhaz=breslow covariates=cov;
run;
data breslow; 
	set breslow; 
	daysyears = days/365.25; 
run; 
proc gplot data=breslow;
	plot breslow*daysyears=tment/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 0.7 by 0.1 minor=none label=(a=90 'Cumulative baseline hazard');
	symbol1  v=none i=stepjl c=blue;
run;



*---------------------------------------------------------------;
*--------------------- Table 2.3 -------------------------------;
*---------------------------------------------------------------; 
proc means data=pbc3 mean;
  class tment;
	var alb bili;
run; 

*---------------------------------------------------------------;
*--------------------- Table 2.4 -------------------------------;
*---------------------------------------------------------------; 

proc phreg data=pbc3;
	model days*status(0)=tment alb bili / rl;
run;

 
*---------------------------------------------------------------;
* --In-text: Poisson model with treatment only (and time)-------;
*---------------------------------------------------------------;

proc genmod data=pbc3mult;
	class interv tment;
	model fail = interv tment / dist=poi offset=logrisk;
	estimate 'CyA  vs placebo' tment -1 1 / exp;
run;

*---------------------------------------------------------------;
*--------------------- Table 2.5 -------------------------------;
*---------------------------------------------------------------;
options ps=200;
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= tment interv alb  bili / dist=poi offset=logrisk ;
run;


*---------------------------------------------------------------;
*--------------------- Table 2.6 -------------------------------;
*---------------------------------------------------------------; 
 
* Transform albumin and bilirubin values; 
data pbc3; 
	set pbc3;
	albnorm=(alb-35)*(alb>35);
	alb10=alb/10; 
	alb2=alb*alb;
	bilihigh=(bili-17.1)*(bili>17.1);
	bilitoohigh=(bili-34.2)*(bili>34.2);
	bilimuchtoohigh=(bili-51.3)*(bili>51.3);
	bili100=bili/100;
	bili2=bili*bili;
	log2bili=log2(bili);
	logbilihigh=(log2bili-log2(17.1))*(bili>17.1);
	logbilitoohigh=(log2bili-log2(34.2))*(bili>34.2);
	logbilimuchtoohigh=(log2bili-log2(51.3))*(bili>51.3);
	log2bili2=log2bili*log2bili;
run;

** Cox models **;

* Base models for LR tests; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb bili/rl ties=breslow CONVERGELIKE=1E-8 type3(lr);
run;
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili/rl ties=breslow CONVERGELIKE=1E-8 type3(lr);
run;

* Splines Cox 1; 
* LRT can here be read of type3 test result; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb albnorm bili/rl ties=breslow CONVERGELIKE=1E-8 type3(lr);
run;

* Quadratic Cox 1; 
* LRT can here be read of type3 test result; 
proc phreg data=pbc3;
  alb2=alb*alb;
	class tment (ref='0');
	model days*status(0)=tment alb10 alb2 bili100/rl type3(LR);
run;

* Splines Cox 2;
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb bili bilihigh bilitoohigh 
	                         bilimuchtoohigh/rl type3(lr) ties=breslow CONVERGELIKE=1E-8;
run;
* LRT; 
data p;
	chi2=826.830-802.434;
	p=1-probchi(chi2,3);
proc print;
run;

* Quadratic Cox 2; 
* LRT can here be read of type3 test result; 
proc phreg data=pbc3;
  bili2=bili*bili;
	class tment (ref='0');
	model days*status(0)=tment alb10 bili bili2 / rl type3(LR);
run;

* Splines Cox 3;
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili logbilihigh 
	                         logbilitoohigh logbilimuchtoohigh/rl ties=breslow CONVERGELIKE=1E-8;
*linbili: test  logbilihigh=logbilitoohigh=logbilimuchtoohigh=0;
run;
* LRT; 
data p;
	chi2=805.881-804.267;
	p=1-probchi(chi2,3);
proc print;
run;

* Quadratic Cox 3; 
* LRT can here be read of type3 test result; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb10 log2bili log2bili2/rl type3(LR);
run;

* Split data into timegroups for Poisson models; 
data pbc3mult; 
	set pbc3;
	fail=(days<=2 * 365.25)*(status ne 0);
	risktime=min(2 * 365.25,days);
	logrisk=log(risktime); interv=1; output;  
	if days>2* 365.25 then do;
	fail=(days<=4* 365.25)*(status ne 0);
	risktime=min(2* 365.25,days-2* 365.25);
	logrisk=log(risktime); interv=2; output; end;
	if days>4* 365.25 then do;
	fail=status ne 0; 
	risktime=days-4* 365.25;
	logrisk=log(risktime); interv=3; output; end;
run;


** Possion models **;

* Models to compare with for the LR tests; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb bili/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb log2bili/ dist=poi offset=logrisk type3;
run;

* Splines Poisson 1; 
* LRT can here be read of type3 test result; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb albnorm bili/ dist=poi offset=logrisk type3;
run;

* Quadratic Poisson 1;
* LRT can here be read of type3 test result; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb10 alb2 bili100/ dist=poi offset=logrisk type3;
run;

* Splines Poisson 2; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb bili bilihigh bilitoohigh  
	             bilimuchtoohigh/ dist=poi offset=logrisk type3;
run;
* LRT; 
data p;
	chi2=350.7253-326.1805;
	p=1-probchi(chi2,3);
proc print;
run;

* Quadratic Poisson 2;
* LRT can here be read of type3 test result; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb10 bili100 bili2/ dist=poi offset=logrisk type3;
run;

* Splines Poisson 3;
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb log2bili logbilihigh logbilitoohigh 
	             logbilimuchtoohigh / dist=poi offset=logrisk type3;
run;
* LRT; 
data p;
	chi2=329.2509-327.5375;
	p=1-probchi(chi2,3);
proc print;
run;

* Quadratic Poisson 3;
* LRT can here be read of type3 test result; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb10 log2bili log2bili2/ dist=poi offset=logrisk type3;
run;



*---------------------------------------------------------------;
*--------------------- Figure 2.6 ------------------------------;
*---------------------------------------------------------------; 


* Below is estimates from the following models collected; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb bili/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb albnorm bili/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb bili bilihigh bilitoohigh 
             logbilimuchtoohigh / dist=poi offset=logrisk type3;
run;



* Collect linear predictor information; 

data lin2; set pbc3;
	lp1=1.6076-0.1123*alb; * bilirubin = 0; 
	lp2=0.7828-0.0864*alb-0.0474*albnorm;
	lp3=1.6076+0.0085*bili-38.7*0.1123; *also an effect of albumin at 38.7; 
	lp4=-0.5534+0.0617*bili-0.0168*bilihigh+0.0027*bilitoohigh
	    -0.0428*bilimuchtoohigh-0.0865*38.7;
run;

proc sort data=lin2; by alb; run;

proc gplot data=lin2;
	plot (lp1 lp2)*alb/overlay haxis=axis1 vaxis=axis2;
	axis1 order=20 to 60 by 10 minor=none label=('Se-albumin');
	axis2 order=-6 to 0 by 1 minor=none label=(a=90 'Linear predictor');
	symbol1  v=circle i=join r=1 c=red;
	symbol2 v=none i=join r=1 c=blue;
run; 


*---------------------------------------------------------------;
*--------------------- Figure 2.7 ------------------------------;
*---------------------------------------------------------------; 

proc sort data=lin2; 
	by bili; 
run;


proc gplot data=lin2;
	plot (lp3 lp4)*bili/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 100 minor=none label=('Se-bilirubin');
	axis2 order=-5 to 2 by 1 minor=none label=(a=90 'Linear predictor');
	symbol1  v=circle i=join r=1 c=red;
	symbol2 v=none i=join r=1 c=blue;
run;



*---------------------------------------------------------------;
*--------------------- Figure 2.8 ------------------------------;
*---------------------------------------------------------------; 

* Linear predictors from the following models; 

proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb log2bili / dist=poi offset=logrisk type3;
run;


proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb log2bili logbilihigh logbilitoohigh 
	             logbilimuchtoohigh / dist=poi offset=logrisk type3;
run;

* Create data with linear predictors as a function of bilirubin values; 
data log2; 
	set pbc3;
	lp3=-2.0162+0.6469*log2bili-38.7*0.087;
	lp4=-0.6194+0.198*log2bili+0.8815*logbilihigh-0.2336*logbilitoohigh
	    -0.3139*logbilimuchtoohigh-0.0844*38.7;
run; 


proc sort data=log2; 
	by bili; 
run;


proc gplot data=log2;
	plot (lp3 lp4)*log2bili/overlay haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(Se-bilirubin)');
	axis2 order=-5 to 1 by 1 minor=none label=(a=90 'Linear predictor');
	symbol1 v=circle i=join r=1 c=red;
	symbol2 v=none i=join r=1 c=blue;
run;



*---------------------------------------------------------------;
*--------------------- Table 2.7 -------------------------------;
*---------------------------------------------------------------;

* Cox model; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili / rl;
run;


* Poisson model;
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb log2bili / dist=poi offset=logrisk;
run;



*---------------------------------------------------------------;
*--------------------- Table 2.8 -------------------------------;
*---------------------------------------------------------------;

* Cox models; 
* Cox model 1; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb tment*alb log2bili/rl type3 (LR);
	estimate 'alb, Z=0' alb 1;
	estimate 'alb, Z=1' alb 1 tment*alb 1;
run;

* Cox model 2; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili tment*log2bili/rl  type3 (LR);
	estimate 'log2bili, Z=0' log2bili 1;
	estimate 'log2bili, Z=1' log2bili 1 tment*log2bili 1;
run;


* Poisson models;
* Poisson model 1; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb log2bili alb*tment/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='1') interv;
	model fail= interv tment alb log2bili alb*tment/ dist=poi offset=logrisk type3;
run;

* Poisson model 2; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv;
	model fail= interv tment alb log2bili log2bili*tment/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='1') interv;
	model fail= interv tment alb log2bili log2bili*tment/ dist=poi offset=logrisk type3;
run;


*---------------------------------------------------------------;
*--------------------- Table 2.9 -------------------------------;
*---------------------------------------------------------------;


* Column 1; 
proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='1');
	model fail= interv tment interv*tment alb log2bili/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='2');
	model fail= interv tment interv*tment alb log2bili/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='3');
	model fail= interv tment interv*tment alb log2bili/ dist=poi offset=logrisk type3;
run;


* Column 2;
proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='1');
	model fail= interv tment alb log2bili alb*interv/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='2');
	model fail= interv tment alb log2bili alb*interv/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='3');
	model fail= interv tment alb log2bili alb*interv/ dist=poi offset=logrisk type3;
run;


* Column 3;
proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='1');
	model fail= interv tment alb log2bili log2bili*interv/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='2');
	model fail= interv tment alb log2bili log2bili*interv/ dist=poi offset=logrisk type3;
run;

proc genmod data=pbc3mult;
	class tment (ref='0') interv(ref='3');
	model fail= interv tment alb log2bili log2bili*interv/ dist=poi offset=logrisk type3;
run;


*---------------------------------------------------------------;
*--------------------- Figure 2.9 ------------------------------;
*---------------------------------------------------------------;

data covstr;
alb=0; log2bili=0;
run;

proc phreg data=pbc3;
	model days*status(0)= alb log2bili/rl;
	strata tment;
	baseline out=breslowstr cumhaz=breslow covariates=covstr;
run;

data breslow0;
	set breslowstr; 
	if tment=0; 
	a00=breslow; 
run;

data breslow1; 
	set breslowstr; 
	if tment=1; 	
	a01=breslow; 
run;

data breslow01; 
	merge breslow0 breslow1; 
	by days; 
run;

data breslowrev; 
set breslow01;
	by days;
	retain last1 last2;
	if a00=. then cumhaz0=last1; if a00 ne . then cumhaz0=a00; 
	if a01=. then cumhaz1=last2; if a01 ne . then cumhaz1=a01;
	output;
	last1=cumhaz0; last2=cumhaz1; 
run;

data breslowrev; 
set breslowrev;
	line=exp(-0.574)*cumhaz0;
run;

proc gplot data=breslowrev;
	plot cumhaz1*cumhaz0 line*cumhaz0/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 1.5 by 0.5 minor=none label=('Cumulative baseline hazard: placebo');
	axis2 order=0 to 1 by 0.5 minor=none label=(a=90 'Cumulative baseline hazard: CyA');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=rl c=blue;
run;




*---------------------------------------------------------------;
*--------------------- Figure 2.10 -----------------------------;
*---------------------------------------------------------------;


* Made in R with timereg package ; 


*---------------------------------------------------------------;
*--------------------- Figure 2.11 -----------------------------;
*---------------------------------------------------------------;


* Made in R with timereg package ; 


*---------------------------------------------------------------;
*--------------------- Table 2.11 ------------------------------;
*---------------------------------------------------------------;


data pbc3add; 
set pbc3mult;
	time1=(interv=1)*risktime; 
	time2=(interv=2)*risktime; 
	time3=(interv=3)*risktime;
	tment0=(tment=0)*risktime; 
	tment1=(tment=1)*risktime;
	albny=(alb-35)/100*risktime;
	biliny=(bili-50)/1000*risktime;
run;


proc genmod data=pbc3add;
model fail=time1 time2 time3 tment1/dist=poi link=id noint;
run;

proc genmod data=pbc3add;
model fail=time1 time2 time3 tment1 albny biliny/dist=poi link=id noint;
run;


