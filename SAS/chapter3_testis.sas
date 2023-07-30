*------------------------------------------------------------------;
*------- Chapter 3, testis data, SAS code -------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=testis
	datafile="data/testis.csv"
	dbms=csv replace;
run;
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
	class age (ref='20') par2 (ref='1');
	model cases=par2 age/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;

* Column 2; 
proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='1');
	model cases=par2 age motherage cohort/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;


* In-text interaction test; 
proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='1');
	model cases=par2 age motherage cohort par2*age/dist=poi offset=lpyrs type3;
run;

* seminomas;
proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='1');
	model semi=par2 age motherage cohort/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;

* non-seminomas;
proc genmod;
	class age (ref='20') motherage(ref='30') cohort(ref='1973') par2 (ref='1');
	model nonsemi=par2 age motherage cohort/dist=poi offset=lpyrs type3;
	estimate 'RR' par2 1 -1/exp;
run;
