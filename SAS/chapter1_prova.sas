*------------------------------------------------------------------;
*------- Chapter 1, SAS code, PROVA data -------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=provany
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/prova.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=provany; 
run;


*---------------------------------------------------------------;
*--------------------- Table 1.2 -------------------------------;
*---------------------------------------------------------------;

data provany; 
	set provany; 
	beh = scle*2 + beta; 
run; 

proc freq data=provany; 
	table beh beh*bleed beh*death beh*death*bleed; 
run;