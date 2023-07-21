*------------------------------------------------------------------;
*------- Chapter 4, SAS code, PROVA data --------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=prova
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/prova.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=prova; 
run;

data prova; 
	set prova; 
	beh = scle + beta*2; 
	log2bili = log2(bili);
	if bleed = 1 then wait = timedeath - timebleed;
run; 

data provany; 
	set prova;
	if bleed=1 then do; btime=timebleed; d0time=timebleed; dead0=0; outof0=1; 
	bdtime=timedeath; deadb=death; wait=bdtime-timebleed; 
	end;
	if bleed=0 then do; btime=timedeath; d0time=timedeath; dead0=death; outof0=death; 
	bdtime=.; deadb=.; wait=.; end;
	log2bili=log2(bili);
run;


*---------------------------------------------------------------;
*--------------------- Figure 4.24 -----------------------------;
*---------------------------------------------------------------;

* distribution of cens.;
data provany; 
	set provany;
	if wait>0 then time=btime+bdtime;
	if wait = . then time=btime;
run;

proc phreg data=provany atrisk;
	model time*death(1)=;
	baseline out=survcens survival=kmc / method=pl;
run;

data survcens; 
	set survcens; 
	timey = time / 365.25; 
run;

proc gplot data=survcens;
	plot kmc*timey/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Probability of no censoring');
	symbol1  v=none i=stepjl c=black;
run;

*---------------------------------------------------------------;
*--------------------- Table 4.13 ------------------------------;
*---------------------------------------------------------------;

proc phreg data=provany;
	class beh (ref='0');
	model time*death(1)=beh/rl;
run;

proc phreg data=provany;
	class varsize (ref='1');
	model time*death(1)=varsize/rl;
run;

proc phreg data=provany;
	class sex (ref='1');
	model time*death(1)=sex/rl;
run;


proc phreg data=provany;
	model time*death(1)=coag/rl;
run;


proc phreg data=provany;
	model time*death(1)=log2bili/rl;
run;

proc phreg data=provany;
	model time*death(1)=age/rl;
run;
