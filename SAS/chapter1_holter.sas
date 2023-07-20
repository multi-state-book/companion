*------------------------------------------------------------------;
*------- Chapter 1, SAS code, Holter data  ------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import 
datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/cphholter.csv"
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
	tables esvea; 
run;

proc freq data=holter; 
	tables afib*stroke*death*esvea; 
run;


proc freq data=holter; 
where timeafib < timestroke; 
	tables afib*stroke*death*esvea; 
run;


proc freq data=holter; 
where timestroke < timeafib; 
	tables afib*stroke*death*esvea; 
run;
