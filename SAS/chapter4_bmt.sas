*------------------------	------------------------------------------;
*------- Chapter 4, bmt data, SAS code ----------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=bmt
	datafile="data/bmt.csv"
	dbms=csv replace;
run;
data bmt; 
set bmt;
	intxsurv=timedeath; dead=death;
	if rel=1 then intxrel=timerel; if rel=0 then intxrel=timedeath;
	trm=0; if rel=0 and death=1 then trm=1;
	state0=rel+2*trm;
	if gvhd=1 then tgvhd=timegvhd; if gvhd=0 then tgvhd=intxrel;
run;

*---------------------------------------------------------------;
*--------------------- Figure 4.15 -----------------------------;
*---------------------------------------------------------------;

proc phreg data=bmt noprint; /* Relapse-free surv */
	model intxrel*state0(0)=;
	baseline out=surv survival=km;
run;

proc phreg data=bmt noprint; /* Relapse */
	model intxrel*state0(0)=/eventcode=1;
	baseline out=cif1 cif=cif1;
run;

proc phreg data=bmt noprint; /* Death in remission */
	model intxrel*state0(0)=/eventcode=2;
	baseline out=cif2 cif=cif2;
run;

proc phreg data=bmt noprint; /* Overall surv. */
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
	axis1 order=0 to 150 by 10 minor=none label=('Months');
	axis2 order=0 to 0.05 by 0.01 minor=none label=(a=90 'Relapse prev. and prob.');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;
quit;


*---------------------------------------------------------------;
*--------------------- In-text Q_1(t) etc ----------------------;
*---------------------------------------------------------------;

* Macro for area under CIF; 
%macro areastep(data,trt,grp,time,y,tau);
	data select; set &data;	where &trt=&grp;
	run;
	data select; set select;by &trt;
		retain mu;
		if first.&trt then mu=0;
		if &time>&tau then do;
		&time=&tau; &y=lag(&y); end;
		mu+lag(&y)*(&time-lag(&time));
		if last.&trt then do;
		if &time<&tau then mu+(&tau-&time)*&y; end;
	run;
	data last;
	  set select;
		by &trt;
		if last.&trt;
	run;
	proc print; var mu; run;
%mend areastep;

%areastep(allrev,tment,0,time,q1,120);
%areastep(allrev,tment,0,time,q0,120);
%areastep(allrev,tment,0,time,c23,120);
%areastep(allrev,tment,0,time,q2,120);
%areastep(allrev,tment,0,time,q3,120);

/* Bootstrap */
data bootbmt;
	do sampnum = 1 to 1000; /* nboot=1000*/
	do i = 1 to 2009; /*nobs=2009*/
	x=round(ranuni(0)*2009); /*nobs=2009*/
	set bmt
	point=x;
	output;
	end;
	end;
	stop;
run;

%macro areastepby(data,byvar,beh,grp,tid,y,tau);
	data select;
		set &data;
		where &beh=&grp;
	run;
	data select;
		set select;
		by &byvar;
		retain mu oldt oldy;
		if first.&byvar then do oldt=0; oldy=1; mu=0;  end;
		if &tid>&tau then do;
		&tid=&tau; &y=oldy; end;
		if not first.&byvar then mu+oldy*(&tid-oldt);
		if last.&byvar then do;
		if &tid<&tau then mu+(&tau-&tid)*&y; end;
		oldy=&y; oldt=&tid;
	run;
	data last;
		set select;
		by  &byvar;
		if last.&byvar;
	run;
%mend areastepby;

proc phreg data=bootbmt noprint; /* Relapse-free surv */
by sampnum;
model intxrel*state0(0)=;
baseline out=surv survival=km;
run;

proc phreg data=bootbmt noprint; /* Relapse */
by sampnum;
model intxrel*state0(0)=/eventcode=1;
baseline out=cif1 cif=cif1;
run;

proc phreg data=bootbmt noprint; /* Death in remission */
by sampnum;
model intxrel*state0(0)=/eventcode=2;
baseline out=cif2 cif=cif2;
run;

proc phreg data=bootbmt noprint; /* Overall surv. */
by sampnum;
model intxsurv*dead(0)=/eventcode=1;
baseline out=dead cif=cif23;
run;

data dead; set dead; time=intxsurv; drop intxsurv;  run;  
data surv; set surv; time=intxrel; drop intxrel; run;
data cif1; set cif1; time=intxrel; drop intxrel; run;
data cif2; set cif2; time=intxrel; drop intxrel; run;
data all; merge surv cif1 cif2 dead ; by sampnum time; run;

data allrev; set all;
by sampnum time; 
retain last1 last2 last3 last4;
if km=. then rfs=last1; if km ne . then rfs=km; 
if cif1=. then c1=last2; if cif1 ne . then c1=cif1;
if cif2=. then c2=last3; if cif2 ne . then c2=cif2;
if cif23=. then c23=last4; if cif23 ne . then c23=cif23;
output;
last1=rfs; last2=c1; last3=c2; last4=c23;
run;

data allrev; set allrev;
q0=rfs; q2=c2; q3=c23-c2; q1=c1-q3; sum=q0+q1+q2+q3; prev=q1/(q0+q1); tment=0;
run;

%areastepby(allrev,sampnum,tment,0,time,q0,120);
proc means data=last mean stddev;
var mu;
run;

/* macro need to be changed for cuminc (start in 0) */

%macro areastepby0(data,byvar,beh,grp,tid,y,tau);
	data select;
		set &data;
		where &beh=&grp;
	run;
	data select;
		set select;
		by &byvar;
		retain mu oldt oldy;
		if first.&byvar then do oldt=0; oldy=0; mu=0;  end;
		if &tid>&tau then do;
		&tid=&tau; &y=oldy; end;
		if not first.&byvar then mu+oldy*(&tid-oldt);
		if last.&byvar then do;
		if &tid<&tau then mu+(&tau-&tid)*&y; end;
		oldy=&y; oldt=&tid;
	run;
	data last;
		set select;
		by  &byvar;
		if last.&byvar;
	run;
%mend areastepby;

%areastepby0(allrev,sampnum,tment,0,time,q1,120);
proc means data=last mean stddev;
var mu;
run;
%areastepby0(allrev,sampnum,tment,0,time,c23,120);
proc means data=last mean stddev;
var mu;
run;
%areastepby0(allrev,sampnum,tment,0,time,q2,120);
proc means data=last mean stddev;
var mu;
run;
%areastepby0(allrev,sampnum,tment,0,time,q3,120);
proc means data=last mean stddev;
var mu;
run;


*---------------------------------------------------------------;
*--------------------- Table 4.8 -------------------------------;
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
proc phreg data=double02 covs(aggregate);
	/* NB bmonly and all are now binary quantitative */
	model time*dc(0)=gsource0 gsource2 disease0 disease2 age0 age2 / corrb;
	strata version;
	id id;
	gs:  test gsource0=gsource2;
	dis: test disease0=disease2;
	a:   test age0=age2;
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
