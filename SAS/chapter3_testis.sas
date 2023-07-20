*------------------------------------------------------------------;
*------- Chapter 3, testis data, SAS code -------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=testis
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/testis.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=testis; 
run;

* Add extra variables; 
data testis; 
	set testis; 
	lpyrs = log(pyrs); 
	par2 = (parity>=2); 
run;

*---------------------------------------------------------------;
*--------------------- Table 3.5 -------------------------------;
*---------------------------------------------------------------;

* Column 1; 
proc genmod;
	class age (ref='20') par2 (ref='0');
	model cases=par2 age/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;

* Column 2; 
proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='0');
	model cases=par2 age motherage cohort/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;


* In-text tests; 
proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='0');
	model cases=par2 age motherage cohort par2*age/dist=poi offset=lpyrs type3;
run;

proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='0');
	model cases=par2 age motherage cohort par2*cohort/dist=poi offset=lpyrs type3;
run;

proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='0');
	model semi=par2 age motherage cohort/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;

proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='0');
	model nonsemi=par2 age motherage cohort/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;
