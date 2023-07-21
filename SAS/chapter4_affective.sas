*------------------------------------------------------------------;
*------- Chapter 4, SAS code, affective data  ---------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=affective
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/affective.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=affective; 
run;

data affective; 
	set affective; 
	wait = stop - start; 
run; 


*---------------------------------------------------------------;
*--------------------- Figure 4.21 -----------------------------;
*---------------------------------------------------------------;

data angstprev; 
	set affective;
	by id;
	retain prev;
	if first.id then prev=0; 
	output; 
	if state=1 then prev=start; if state=0 then prev=stop;
run;


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

*---------------------------------------------------------------;
*--------------------- Table 4.8 -------------------------------;
*---------------------------------------------------------------;

* Mortality not taken into account; 
proc phreg covs(aggregate) data=angstprev;
	where state=0 or status=2 or status=3;
	class bip (ref='0');
	model stop*status(2 3)=bip/entry=prev rl;
	id id;
run;

* Mortality taken into account; 
proc phreg data=angstprev;
	where state=0 or status=2 or status=3;
	class bip (ref='0');
	model stop*status(3)=bip/entry=prev eventcode=1 rl;
run;

*---------------------------------------------------------------;
*--------------------- Table 4.10 ------------------------------;
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
*--------------------- Table 4.11 ------------------------------;
*---------------------------------------------------------------;

data angstwlw4; 
set angstwlw;
	by id;
	time=stop; dc=status; stratum=episode; output; 

	/* if last episode is not #4 then later episodes are either censored (1 or 3) or they
	   'end in death' (2) */

	if last.id then do;
	if episode=3 then do;
	time=stop; if status=1 or status=3 then dc=0; 
	if status=2 then dc=2; stratum=4;  output; end;


	if episode=2 then do;
	time=stop; if status=1 or status=3 then dc=0; 
	if status=2  then dc=2; stratum=3; output; 
	time=stop; if  status=1 or status=3 then dc=0; 
	if status=2  then dc=2; stratum=4; output; end;

	if episode=1 then do; 
	time=stop; if status=1 or status=3 then dc=0; 
	if status=2 then dc=2; stratum=2; output; 
	time=stop; if status=1 or status=3 then dc=0; 
	if status=2  then dc=2; stratum=3; output; 
	time=stop; if status=1 or status=3 then dc=0; 
	if status=2 then dc=2; stratum=4; output; 
	end;

	end;
run;
/*
proc export data=angstwlw4
	outfile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/affectivewlw.csv"
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
*--------------------- Figure 4.22 -----------------------------;
*---------------------------------------------------------------;


proc phreg data=angstprev;
where state=0 or status=2 or status=3;
model stop*status(3)=/entry=prev eventcode=1;
strata bip;
baseline out=mcfdata1 cif=naa1;
run;

data mcfdata1; set mcfdata1;
cmf=-log(1-naa1);
years=stop/12;
run;

* fig 4.22; 
proc gplot data=mcfdata1;
plot cmf*years=bip/haxis=axis1 vaxis=axis2;
axis1 order=0 to 30 by 5 minor=none label=('Years');
axis2 order=0 to 12 by 2 minor=none label=(a=90 'Expected number of episodes');
symbol1  v=none i=stepjl c=red;
symbol2  v=none i=stepjl c=blue;
run;

*---------------------------------------------------------------;
*--------------------- Figure 4.26 -----------------------------;
*---------------------------------------------------------------;

data cens; 
	set affective;
	by id;
	if last.id;
run;

proc phreg data=cens atrisk;
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
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Probability of no censoring');
	symbol1  v=none i=stepjl c=black;
run;

proc phreg data=cens;
	model stop*status(2)=bip/rl;
run;
