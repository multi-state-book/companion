*------------------------------------------------------------------;
*------- Chapter 4, SAS code, PROVA data --------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=prova
	datafile="data/prova.csv"
	dbms=csv replace;
run;
data cens;
  set prova;
	if timebleed=. then time=timedeath;
	else time= timebleed + timedeath;
	beh = scle + beta*2; 
	log2bili = log2(bili);
run;

*---------------------------------------------------------------;
*--------------------- Figure 4.22 -----------------------------;
*---------------------------------------------------------------;

* distribution of cens.;
proc phreg data=cens atrisk noprint;
	model time*death(1)=;
	baseline out=survcens survival=kmc / method=pl;
run;
data survcens; 
	set survcens; 
	timey = time/365.25; 
run;
proc gplot data=survcens;
	plot kmc*timey/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Probability of no censoring');
	symbol1  v=none i=stepjl c=black;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Table 4.12 ------------------------------;
*---------------------------------------------------------------;

proc phreg data=cens;
	class beh (ref='0');
	model time*death(1)=beh;
run;
proc phreg data=cens;
	class varsize (ref='1');
	model time*death(1)=varsize;
run;
proc phreg data=cens;
	class sex (ref='1');
	model time*death(1)=sex;
run;
proc phreg data=cens;
	model time*death(1)=coag;
run;
proc phreg data=cens;
	model time*death(1)=log2bili;
run;
proc phreg data=cens;
	model time*death(1)=age;
run;
