*------------------------------------------------------------------;
*------- Chapter 2, SAS code, Bissau data -------------------------;
*------------------------------------------------------------------;

* Load data; 
options ps=200;

proc import out=bissau
	datafile="data/bissau.csv"
	dbms=csv replace;
run;
data bissau; 
	set bissau; 
	agem = int(age/30.44);
run;

*---------------------------------------------------------------;
*--------------------- Table 2.12 ------------------------------;
*---------------------------------------------------------------;

* Age in months; 

* Covariates; 
data cov1;
	input bcg agem ;
	datalines;
	1 3
;
run;

* Cox model fit - column 1; 
proc phreg data=bissau;
	model fuptime*dead(0)=bcg agem/rl;
	baseline out=hazfuptime covariates=cov1 cumhaz=basefup;
run;


* Make age the time variable instead;
data bissau; 
	set bissau;
	age12=age/(365.24/12);
	ageout=age12+fuptime/(365.24/12);
run;

data cov2;
input bcg;
datalines;
1
;
run;

* Cox model fit - column 2; 
proc phreg data=bissau;
	model ageout*dead(0)=bcg/entry=age12 rl;
	baseline out=hazage covariates=cov2 cumhaz=baseage;
run;


*---------------------------------------------------------------;
*--------------------- Figure 2.12 -----------------------------;
*---------------------------------------------------------------;

data hazfuptime_months; 
	set hazfuptime; 
	fuptime_m = fuptime / (365.24/12); 
run;

proc gplot data=hazfuptime_months;
	plot basefup*fuptime_m/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 7 by 1 minor=none label=('Follow-up time (months)');
	axis2 order=0 to 0.06 by 0.01 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=blue;
run;



*---------------------------------------------------------------;
*--------------------- Figure 2.13 -----------------------------;
*---------------------------------------------------------------;

proc gplot data=hazage;
	plot baseage*ageout/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 15 by 5 minor=none label=('Age (months)');
	axis2 order=0 to 0.15 by 0.05 minor=none label=(a=90 'Cumulative baseline hazard');
	symbol1  v=none i=stepjl c=blue;
run;
