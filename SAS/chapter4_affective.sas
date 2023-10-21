*------------------------------------------------------------------;
*------- Chapter 4, SAS code, affective data  ---------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=affective
	datafile="data/affective.csv"
	dbms=csv replace;
run;
data affective; 
	set affective; 
	wait = stop - start; 
run; 
data angstprev; 
	set affective;
	by id;
	retain prev;
	if first.id then prev=0; 
	output; 
	if state=1 then prev=start; if state=0 then prev=stop;
run;


*---------------------------------------------------------------;
*--------------------- Table 4.7 -------------------------------;
*---------------------------------------------------------------;

* Mortality treated as censoring; 
proc phreg covs(aggregate) data=angstprev;
	where state=0 or status=2 or status=3;
	class bip (ref='0');
	model stop*status(2 3)=bip/entry=prev rl;
	id id;
run;

* Mortality treated as competing risk; 
proc phreg data=angstprev;
	where state=0 or status=2 or status=3;
	class bip (ref='0');
	model stop*status(3)=bip/entry=prev eventcode=1 rl;
run;

*---------------------------------------------------------------;
*In-text, p. 145: Cox for mortality;
*---------------------------------------------------------------;

proc phreg data=affective;
  model (start,stop)*status(0,1,3) = bip / rl;
run;

*---------------------------------------------------------------;
*--------------------- Figure 4.18 -----------------------------;
*---------------------------------------------------------------;

ods graphics on; 
proc phreg plots(overlay=row)=mcf covs(aggregate) data=angstprev;
	where state=0 or status=2 or status=3;
	class bip;
	model stop*status(2 3)=/entry=prev;
	id id;
	strata bip;
	baseline out=mcfdata cmf=naa;
run;
data mcfdata; set mcfdata;
	years=stop/12;
run;
proc gplot data=mcfdata;
	plot naa*years=bip/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 30 by 5 minor=none label=('Years');
	axis2 order=0 to 12 by 2 minor=none label=(a=90 'Expected number of episodes');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Figure 4.19 -----------------------------;
*---------------------------------------------------------------;

/* Using "fine-gray model" in PHREG gives an alternative solution to 
  the estimator for CMF using the Breslow type estimator for 
  the baseline mean function (see p. 199 in book). The estimator is not
	exactly the same as Cook-Lawless because of a different procedures 
	for ties of terminating events and censorings. If no ties 
	(or no censorings) it equals Cook & Lawless */

proc phreg data=angstprev;
	where state=0 or status=2 or status=3;
	model stop*status(3)=/entry=prev eventcode=1;
	strata bip;
	baseline out=mcfdata1 cif=naa1;
run;
data mcfdata1;
	set mcfdata1;
	cmf=-log(1-naa1);
	years=stop/12;
run;
proc gplot data=mcfdata1;
	plot cmf*years=bip/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 30 by 5 minor=none label=('Years');
	axis2 order=0 to 12 by 2 minor=none label=(a=90 'Expected number of episodes');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;
quit;


/*** Calc Cook & Lawless or (Ghosh & Lin (GL)) estimator for CMF 'by hand' ***/
/* First create KM data for death */
proc phreg data=angstprev noprint;
	where state=0 or status=2 or status=3;
  model stop*status(1 3)= / entry=prev; /* status=2=death */
  strata bip;
  baseline out=kmdata survival=km / method=pl ;
run;
/* Second create NAa data */
proc phreg data=angstprev noprint;
	where state=0 or status=2 or status=3;
  model stop*status(2 3)= / entry=prev;/* status=1=event */
  strata bip;
  baseline out=nadata cumhaz=na;
run;
/* Use NA data to calculate dA(u), i.e., increments in NAa */
data na;
  set nadata;
  dAu=na-lag(na);
  if stop=0 then dAu=0;
  keep bip stop dAu na;
run;
/* merge NAa and KM data */
data merged;
  merge na kmdata;
  by bip stop;
run;
/* multiply S(u-) and dA(u) */
data fill;
   set merged;
   retain _km;
   if not missing(km) then _km=km;
   else km=_km;
   /* S(u-) */
   S_uminus=lag(km);
   if stop=0 then S_uminus=1;

   if dAu=. then dAu=0;
   GLfactor=S_uminus*dAu;
   keep bip stop na dAu S_uminus GLfactor;
run;
data GLdata;
  set fill;
  by bip;
  if first.bip then GL=0;
  else GL+GLfactor;
run;
proc sgplot data=GLdata;
  step x=stop y=GL / group=bip;
  step x=stop y=na / group=bip;
run;

*---------------------------------------------------------------;
*--------------------- Table 4.9 ------------------------------;
*---------------------------------------------------------------;

data angstwlw; 
	set affective;
	if episode<5 and (state=0 or status=2 or status=3);
run;
proc sort data=angstwlw; 
	by id; 
run;
proc freq data=angstwlw;
	tables episode*status; 
run;

*---------------------------------------------------------------;
*--------------------- Table 4.10 ------------------------------;
*---------------------------------------------------------------;

proc print data=affective;
  var id state start stop status episode;
run;
data angstwlw; 
	set affective;
	where episode<5 and (state=0 or status=2 or status=3);
run;
proc freq;
  table id*episode/nocol norow nopercent;
	run;
proc sort data=angstwlw; 
	by id; 
run;
proc print;
  var id state start stop status episode;
run;
data angstwlw4; 
	set angstwlw;
	by id;
	time=stop; dc=status; stratum=episode;
	output; 
	/* if last episode is not #4 then later episodes are either
		 censored (1 or 3) or the 'end in death' (2) */
	if last.id then do;
		if episode=3 then do;
			time=stop;
  		if status=1 or status=3 then dc=0; 
			if status=2 then dc=2;
			stratum=4;
			output;
  	end;
		if episode=2 then do;
			time=stop; 
			if status=1 or status=3 then dc=0; 
			if status=2  then dc=2;
			stratum=3;
			output; 
			time=stop;
			if status=1 or status=3 then dc=0; 
			if status=2 then dc=2;
			stratum=4;
			output; 
		end;
		if episode=1 then do; 
			time=stop;
			if status=1 or status=3 then dc=0; 
			if status=2 then dc=2;
			stratum=2;
			output; 
			time=stop; if status=1 or status=3 then dc=0; 
			if status=2  then dc=2;
			stratum=3;
			output; 
			time=stop;
			if status=1 or status=3 then dc=0; 
			if status=2 then dc=2;
			stratum=4;
			output; 
		end;
	end;
run;
proc print;
run;
/*
proc export data=angstwlw4
	outfile="data/affectivewlw.csv"
	dbms=csv replace;
run;
*/

data angstwlw4; set angstwlw4;
	bip1=bip*(stratum=1); bip2=bip*(stratum=2);
	bip3=bip*(stratum=3); bip4=bip*(stratum=4);
run;

/* composite end point */ 
proc phreg data=angstwlw4 covs(aggregate);
	model time*dc(0 3)=bip1 bip2 bip3 bip4;
	strata stratum;
	id id;
	bip: test bip1=bip2=bip3=bip4;
run;

proc phreg data=angstwlw4 covs(aggregate);
	model time*dc(0 3)=bip;
	strata stratum;
	id id;
run;

/* Same analyses of cause-spec. hazards for 1.,2.,3.,4. event */
proc phreg data=angstwlw4 covs(aggregate);
	model time*dc(0 2 3)=bip1 bip2 bip3 bip4;
	strata stratum;
	id id;
	bip: test bip1=bip2=bip3=bip4;
run;
proc phreg data=angstwlw4 covs(aggregate);
	model time*dc(0 2 3)=bip;
	strata stratum;
	id id;
run;


*---------------------------------------------------------------;
*--------------------- Figure 4.23 -----------------------------;
*---------------------------------------------------------------;

data cens; 
	set affective;
	by id;
	if last.id;
run;
proc phreg data=cens atrisk noprint;
	model stop*status(2)=;
	baseline out=angstcens survival=kmc / method=pl;
run;
data angstcens; 
	set angstcens; 
	years=stop/12; 
run;
proc gplot data=angstcens;
	plot kmc*years/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 30 by 5 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none
	label=(a=90 'Probability of no censoring');
	symbol1  v=none i=stepjl c=black;
run;
quit;

* in-text p-value;
proc phreg data=cens;
	model stop*status(2)=bip/rl;
run;

proc freq data=cens;table status;run;
