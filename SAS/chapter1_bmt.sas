*------------------------------------------------------------------;
*------- Chapter 1, bmt data, SAS code ----------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=bmt
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/bmt.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=bmt; 
run;

*---------------------------------------------------------------;
*--------------------- Table 1.4 -------------------------------;
*---------------------------------------------------------------;
proc freq data=bmt; 
	table rel death rel*death gvhd gvhd*rel gvhd*death; 
run; 