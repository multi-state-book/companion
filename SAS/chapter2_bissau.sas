*------------------------------------------------------------------;
*------- Chapter 2, SAS code, Bissau data -------------------------;
*------------------------------------------------------------------;
options ps=200;

* Load data; 
proc import out=bissau
	datafile="data/bissau.csv"
	dbms=csv replace;
run;
data bissau; 
	set bissau; 
	agem = int(age/30.4); * Age in months; 
run;

*---------------------------------------------------------------;
*--------------------- Table 2.12 ------------------------------;
*---------------------------------------------------------------;

proc phreg data=bissau;
	model fuptime*dead(0)=bcg / rl;
run;

proc phreg data=bissau;
	model fuptime*dead(0)=bcg agem / rl;
run;

* Make age the time variable instead;
data bissau; 
	set bissau;
	agein=age/(365.24/12);
	ageout=agein+fuptime/(365.24/12);
run;

* Cox model fit - column 2; 
proc phreg data=bissau;
	model ageout*dead(0)=bcg / entry=agein rl;
run;


*---------------------------------------------------------------;
*--------------------- Figure 2.12 -----------------------------;
*---------------------------------------------------------------;

* Covariates; 
data cov1;
	input bcg agem ;
	datalines;
	0 3
;
run;

* Cox model fit and get cumulative baseline hazard ; 
proc phreg data=bissau;
	model fuptime*dead(0)=bcg agem/rl;
	baseline out=hazfuptime covariates=cov1 cumhaz=basefup;
run;
data hazfuptime_months; 
	set hazfuptime; 
	fuptime_m = fuptime / (365.24/12); 
run;
proc gplot data=hazfuptime_months;
	plot basefup*fuptime_m/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 7 by 1 minor=none
		label=('Follow-up time (months)');
	axis2 order=0 to 0.06 by 0.01 minor=none
		label=(a=90 'Cumulative hazard');
	symbol1 v=none i=stepjl c=blue;
run;
quit;
* or using proc sgplot;
ods graphics on;
proc sgplot data=hazfuptime_months;
  step x=fuptime_m y=basefup / justify=left;
	xaxis grid values=(0 to 7 by 1);
	yaxis grid values=(0 to 0.06 by 0.01);
	label fuptime_m="Follow-up time (months)";
	label basefup="Cumulative hazard";
run;


*---------------------------------------------------------------;
*--------------------- Figure 2.13 -----------------------------;
*---------------------------------------------------------------;

data cov2;
input bcg;
datalines;
0
;
run;
* Cox model fit and get cumulative baseline hazard ; 
proc phreg data=bissau;
	model ageout*dead(0)=bcg / entry=agein rl;
	baseline out=hazage covariates=cov2 cumhaz=baseage;
run;
proc gplot data=hazage;
	plot baseage*ageout/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 15 by 5 minor=none
		label=('Age (months)');
	axis2 order=0 to 0.15 by 0.05 minor=none 
		label=(a=90 'Cumulative baseline hazard');
	symbol1  v=none i=stepjl c=blue;
run;
quit;
