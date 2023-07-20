*------------------------------------------------------------------;
*------- Chapter 1, SAS code, Bissau data -------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import 
datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/bissau.csv"
	out=bissau
	dbms = csv
	replace;
run;

* Summarise data set; 
proc contents 
	data=bissau; 
run;

*---------------------------------------------------------------;
*--------------------- Table 1.1 -------------------------------;
*---------------------------------------------------------------;

proc freq data=bissau; 
	table bcg*dead;
run;