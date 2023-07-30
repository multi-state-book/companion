*------------------------------------------------------------------;
*------- Chapter 3, SAS code, Bissau data -------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import 
	datafile="data/bissau.csv" out=bissau
	dbms = csv replace;
run;

*---------------------------------------------------------------;
*--------------------- Table 3.1 -------------------------------;
*---------------------------------------------------------------;

data bissau; 
	set bissau;
	dtpany=(dtp>0);
	drop var1; 
run; 
proc freq data=bissau;
	tables bcg*dtp bcg*dtpany/ nocol nopercent;
run;


*---------------------------------------------------------------;
*--------------------- Table 3.2 -------------------------------;
*---------------------------------------------------------------;

data bissau; 
	set bissau;
	agein=age/(365.24/12);
	ageout=agein+fuptime/(365.24/12);
run;

proc phreg data=bissau;
	model ageout*dead(0) = bcg / entry=agein rl;
run;
proc phreg data=bissau;
	model ageout*dead(0) = dtpany / entry=agein rl;
run;
proc phreg data=bissau;
	model ageout*dead(0) = bcg dtpany / entry=agein rl;
run;
proc phreg data=bissau;
	model ageout*dead(0) = bcg dtpany bcg*dtpany / entry=agein rl;
run;
