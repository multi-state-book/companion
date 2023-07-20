*------------------------------------------------------------------;
*------- Chapter 4, bmt data, SAS code ----------------------------;
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

* Add extra variables;
data bmt; 
set bmt;
	intxsurv=timedeath; dead=death;
	if rel=1 then intxrel=timerel; if rel=0 then intxrel=timedeath;
	trm=0; if rel=0 and death=1 then trm=1;
	state0=rel+2*trm;
	if gvhd=1 then tgvhd=timegvhd; if gvhd=0 then tgvhd=intxrel;
run;

*---------------------------------------------------------------;
*--------------------- Figure 4.17 -----------------------------;
*---------------------------------------------------------------;

proc phreg data=bmt; /* Relapse-free surv */
	model intxrel*state0(0)=;
	baseline out=surv survival=km;
run;

proc phreg data=bmt; /* Relapse */
	model intxrel*state0(0)=/eventcode=1;
	baseline out=cif1 cif=cif1;
run;

proc phreg data=bmt; /* Death in remission */
	model intxrel*state0(0)=/eventcode=2;
	baseline out=cif2 cif=cif2;
run;

proc phreg data=bmt; /* Overall surv. */
	model intxsurv*dead(0)=/eventcode=1;
	baseline out=dead cif=cif23;
run;

/* We need the same time variable for all prob's */
data dead; set dead; time=intxsurv; run;
data surv; set surv; time=intxrel; run;
data cif1; set cif1; time=intxrel; run;
data cif2; set cif2; time=intxrel; run;
data all; merge surv cif1 cif2 dead; by time; run;

data allrev; 
set all;
	by time;
	retain last1 last2 last3 last4;
	if km=. then rfs=last1; if km ne . then rfs=km; 
	if cif1=. then c1=last2; if cif1 ne . then c1=cif1;
	if cif2=. then c2=last3; if cif2 ne . then c2=cif2;
	if cif23=. then c23=last4; if cif23 ne . then c23=cif23;
	output;
	last1=rfs; last2=c1; last3=c2; last4=c23;
run;

data allrev; 
set allrev;
	q0=rfs; q2=c2; q3=c23-c2; q1=c1-q3; sum=q0+q1+q2+q3; prev=q1/(q0+q1); tment=0;
run;

proc gplot data=allrev;
	plot prev*time q1*time/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 120 by 10 minor=none label=('Months');
	axis2 order=0 to 0.05 by 0.01 minor=none label=(a=90 'Relapse prev. and prob.');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;

*---------------------------------------------------------------;
*--------------------- Table 4.9 -------------------------------;
*---------------------------------------------------------------;

/* Relapse, relapse-free and overall survival
   without and with adjustment for center */

proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0");
	model intxrel*rel(0)=bmonly all age;
run;

proc phreg data=bmt covs(aggregate);
	class bmonly(ref="0") all(ref="0") team;
	model intxrel*rel(0)=bmonly all age;
	id team;
run;


proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0");
	model intxrel*state0(0)=bmonly all age;
run;

proc phreg data=bmt covs(aggregate);
	class bmonly(ref="0") all(ref="0");
	model intxrel*state0(0)=bmonly all age;
	id team;
run;


proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0");
	model intxsurv*dead(0)=bmonly all age;
run;

proc phreg data=bmt covs(aggregate);
	class bmonly(ref="0") all(ref="0");
	model intxsurv*dead(0)=bmonly all age;
	id team;
run;


*---------------------------------------------------------------;
*--------------------- Table 4.11 -------------------------------;
*---------------------------------------------------------------;

data double02; set bmt; 
	/* joint analysis of relapse-free and overall survival */
	version=1; dc=state0>0; time=intxrel; gsource0=bmonly; gsource2=0;
	disease0=all; disease2=0; age0=age; age2=0; output;
	version=2; dc=dead; time=intxsurv; gsource2=bmonly; gsource0=0;
	disease2=all; disease0=0; age2=age; age0=0; output;
run;

proc phreg data=double02 covs(aggregate) covout outest=params2;
	/* NB bmonly and all are now binary quantitative */
	model time*dc(0)=gsource0 gsource2 disease0 disease2 age0 age2;
	strata version;
	id id;
	gs: test gsource0=gsource2;
	dis: test disease0=disease2;
	a: test age0=age2;
run;

data corr;
	gsource02=0.0059020581/sqrt(0.0058582143*0.0061954571);
	disease02=0.0058782867/sqrt(0.0059598232*0.0062648812);
	age02=6.6520543/sqrt(6.7253187*7.0003067);
proc print;
run;

proc phreg data=double02 covs(aggregate);
	model time*dc(0)=bmonly disease0 disease2 age;
	strata version;
	id id;
	dis: test disease0=disease2;
run;

data bmt; set bmt;
	gvhdny=gvhd; 
	if gvhdny=1 then nytgvhd=tgvhd;
	if gvhdny=0 then do; nytgvhd=intxsurv; if dead=1 then gvhdny=1; end; 
	/* NB her tæller både GvHD og død uden GvHD med */
run;

data triple02G; set bmt; 
	/* joint analysis of relapse-free, GvHD-free and overall survival */
	version=1; dc=state0>0; time=intxrel; gsource0=bmonly; 
	gsource2=0; gsourceG=0;
	disease0=all; disease2=0; diseaseG=0; age0=age; age2=0; ageG=0; output;
	version=2; dc=dead; time=intxsurv; gsource2=bmonly; 
	gsource0=0; gsourceG=0;
	disease2=all; disease0=0; diseaseG=0; age2=age; age0=0; ageG=0; output;
	version=3; dc=gvhdny; time=nytgvhd; gsourceG=bmonly; gsource0=0; gsource2=0;
	diseaseG=all; disease0=0; disease2=0; ageG=age; age0=0; age2=0; output;
run;

proc phreg data=triple02G covs(aggregate) covout outest=params3;
	/* bmonly and all binary quatitative*/
	model time*dc(0)=gsource0 gsource2 gsourceG disease0 disease2 diseaseG
	age0 age2 ageG;
	strata version;
	id id;
run;
