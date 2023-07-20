*------------------------------------------------------------------;
*------- Chapter 3, SAS code, Bissau data -------------------------;
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
*--------------------- Table 3.1 -------------------------------;
*---------------------------------------------------------------;

data bissau; 
	set bissau;
	dtpany=(dtp>0);
	drop var1; 
run; 

proc freq;
	tables bcg*dtp bcg*dtpany/ nocol nopercent;
run;


*---------------------------------------------------------------;
*--------------------- Table 3.2 -------------------------------;
*---------------------------------------------------------------;

data bissau; 
	set bissau;
	age12=age/(365.24/12);
	ageout=age12+fuptime/(365.24/12);
run;

proc phreg data=bissau;
	class bcg(ref="0");
	model ageout*dead(0)=bcg/entry=age12 rl;
run;

proc phreg data=bissau;
	class bcg(ref="0");
	model ageout*dead(0)=dtpany/entry=age12 rl;
run;

proc phreg data=bissau;
	class bcg(ref="0");
	model ageout*dead(0)=bcg dtpany/entry=age12 rl;
run;

proc phreg data=bissau;
	class bcg(ref="0");
	model ageout*dead(0)=bcg dtpany bcg*dtpany/entry=age12 rl;
run;
