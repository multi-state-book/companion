*------------------------------------------------------------------;
*------- Chapter 1, SAS code, Holter data  ------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import 
datafile="data/cphholter.csv"
	out=holter
	dbms = csv
	replace;
run;

* Summarise data set; 
proc contents 
	data=holter; 
run;

*---------------------------------------------------------------;
*--------------------- Table 1.5 -------------------------------;
*---------------------------------------------------------------;

proc freq data=holter; 
  title 'Total';
	tables esvea / nocol norow nopercent; 
run;

proc freq data=holter; 
	title '';
	tables afib*stroke*death*esvea / nocol norow nopercent; 
run;

proc freq data=holter; 
	title '0 -> AF -> Stroke -> dead (yes/no)';
	where .<timeafib <= timestroke; 
	tables afib*stroke*death*esvea / nocol norow nopercent; 
run;

proc freq data=holter; 
	title '0 -> Stroke -> AF -> dead (yes/no)';
	where .<timestroke < timeafib; 
	tables afib*stroke*death*esvea / nocol norow nopercent; 
run;
