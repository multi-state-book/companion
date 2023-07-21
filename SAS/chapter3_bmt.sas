*------------------------------------------------------------------;
*------- Chapter 3, bmt data, SAS code ----------------------------;
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

* Extra variables; 
data bmt; 
set bmt;
	intxsurv=timedeath; dead=death;
	if rel=1 then intxrel=timerel; if rel=0 then intxrel=timedeath;
	trm=0; if rel=0 and death=1 then trm=1;
	state0=rel+2*trm;
	if gvhd=1 then tgvhd=timegvhd; if gvhd=0 then tgvhd=intxrel;
run;

*---------------------------------------------------------------;
*--------------------- Figure 3.5 ------------------------------;
*---------------------------------------------------------------;

* Model fit; 
proc phreg data=bmt;
	model tgvhd*gvhd(0)=;
	baseline out=alfagvh cumhaz=naagvh;
run;

* Make plot; 
proc gplot data=alfagvh;
	plot naagvh*tgvhd/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 156 by 12 minor=none label=('Months');
	axis2 order=0 to 0.8 by 0.1 minor=none label=(a=90 'Cumulative GvHD hazard');
	symbol1  v=none i=stepjl c=blue;
run;


*---------------------------------------------------------------;
*--------------------- Figure 3.6 ------------------------------;
*---------------------------------------------------------------;

* Reformatting; 
data bmt; 
set bmt; /* Censor at GvHD */
	if tgvhd<intxrel then do;
	nyreltime=tgvhd; nyrel=0; nytrm=0; end;
	if tgvhd=intxrel then do;
	nyreltime=intxrel; nyrel=rel; nytrm=trm; end;
run;

* Cumulative relapse rate without gvhd; 
proc phreg data=bmt;
	model nyreltime*nyrel(0)=;
	baseline out=alfa0rel cumhaz=naa0rel;
run;

* Cumulative relapse rate after gvhd; 
proc phreg data=bmt;
	where gvhd=1;
	model intxrel*rel(0)=/entry=tgvhd;
	baseline out=alfagvhrel cumhaz=naagvhrel;
run;

data alfa0rel; 
	set alfa0rel; 
	reltime=nyreltime; 
run;

data alfagvhrel; 
	set alfagvhrel; 
	reltime=intxrel; 
run;

data rel; 
	merge alfa0rel alfagvhrel; 
	by reltime; 
run;

data relrev; 
	set rel;
	by reltime;
	retain last1 last2;
	if naa0rel=. then a02=last1; if naa0rel ne . then a02=naa0rel; 
	if naagvhrel=. then a13=last2; if naagvhrel ne . then a13=naagvhrel;
	output;
	last1=a02; last2=a13; 
run;

proc gplot data=relrev;
	plot a02*reltime a13*reltime/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 156 by 12 minor=none label=('Months');
	axis2 order=0 to 0.2 by 0.1 minor=none label=(a=90 'Cumulative relapse hazard');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;

*---------------------------------------------------------------;
*--------------------- Figure 3.7 ------------------------------;
*---------------------------------------------------------------;

data relrev; 
	set relrev;
	line=0.858*a02;
run;


proc gplot data=relrev;
	plot a13*a02 line*a02/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 0.2 by 0.05 minor=none label=('Cumulative relapse hazard: no GvHD');
	axis2 order=0 to 0.2 by 0.05 minor=none label=(a=90 'Cumulative relapse hazard: GvHD');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=rl c=blue;
run;

*---------------------------------------------------------------;
*--------------------- Figure 3.8 ------------------------------;
*---------------------------------------------------------------;

* Cumulative death rate without gvhd; 
proc phreg data=bmt; 
	model nyreltime*nytrm(0)=;
	baseline out=alfa0dead cumhaz=naa0dead;
run;

* Cumulative death rate with gvhd; 
proc phreg data=bmt;
	where gvhd=1;
	model intxrel*trm(0)=/entry=tgvhd;
	baseline out=alfagvhdead cumhaz=naagvhdead;
run;

* Book-keeping for plot;
data alfa0dead; 
	set alfa0dead; 
	deadtime=nyreltime; 
run;

data alfagvhdead; 
	set alfagvhdead; 
	deadtime=intxrel; 
run;

data dead; 
	merge alfa0dead alfagvhdead; 
	by deadtime; 
run;

data deadrev; 
	set dead;
	by deadtime;
	retain last1 last2;
	if naa0dead=. then a02=last1; if naa0dead ne . then a02=naa0dead; 
	if naagvhdead=. then a13=last2; if naagvhdead ne . then a13=naagvhdead;
	output;
	last1=a02; last2=a13; 
run;

* Add coefficient; 
data deadrev; 
	set deadrev;
	line=3.113*a02;
run;

* Make plot; 
proc gplot data=deadrev;
	plot a13*a02 line*a02/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 0.2 by 0.05 minor=none label=('Cumulative death hazard: no GvHD');
	axis2 order=0 to 1 by 0.2 minor=none label=(a=90 'Cumulative death hazard: GvHD');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=rl c=blue;;
run;


*---------------------------------------------------------------;
*--------------------- Figure 3.4 ------------------------------;
*---------------------------------------------------------------;

* Cumulative death hazards wrt relapse ;
proc phreg data=bmt; 
	model intxrel*trm(0)=;
	baseline out=alfa02 cumhaz=naa02;
run;

data bmt; 
	set bmt; 
	if intxrel eq intxsurv and rel eq 1 then intxsurv = intxsurv + 0.01; 
run; 

proc phreg data=bmt;
	where rel=1;
	model intxsurv*dead(0)=/entry=intxrel;
	baseline out=alfa13 cumhaz=naa13;
run;

* Book-keeping; 
data alfa02; 
	set alfa02; 
	deadtime=intxrel; 
run;

data alfa13; 
	set alfa13; 
	deadtime=intxsurv; 
run;

data alfadead; 
	merge alfa02 alfa13; 
	by deadtime; 
run;

data alfarev; set alfadead;
	by deadtime;
	retain last1 last2;
	if naa02=. then a0rel=last1; if naa02 ne . then a0rel=naa02; 
	if naa13=. then arel=last2; if naa13 ne . then arel=naa13;
	output;
	last1=a0rel; last2=arel; 
run;

* 3 mortality rates simultaneously;
data alfa3dead; 
	merge deadrev alfarev; 
	by deadtime; 
run;

* Make plot;
proc gplot data=alfa3dead;
plot arel*deadtime a02*deadtime a13*deadtime
                 / overlay haxis=axis1 vaxis=axis2;
axis1 order=0 to 120 by 10 minor=none label=('Months');
axis2 order=0 to 9 by 1 minor=none label=(a=90 'Cumulative hazard of death');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
symbol3  v=none i=stepjl c=black;
run;


*---------------------------------------------------------------;
*--------------------- Table 3.11 ------------------------------;
*---------------------------------------------------------------;

** Relapse **; 
* Column 1;
proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0");
	model intxrel*rel(0)=tdcgvhd age bmonly all/rl;
	tdcgvhd=0;
	if gvhd=1 and intxrel>tgvhd then tdcgvhd=1;
run;

* Column 2;
proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0");
	model intxrel*rel(0)=tdcanc500 age bmonly all tdcgvhd/rl;
	tdcanc500=0;
	if anc500=1 and intxrel>timeanc500 then tdcanc500=1;
	tdcgvhd=0;
	if gvhd=1 and intxrel>tgvhd then tdcgvhd=1;
run;


** Death in remission **;
* Column 1; 
proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0");
	model intxrel*trm(0)=tdcgvhd age bmonly all/rl;
	tdcgvhd=0;
	if gvhd=1 and intxrel>tgvhd then tdcgvhd=1;
run;

* Column 2;
proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0");
	model intxrel*trm(0)=tdcanc500 age bmonly all tdcgvhd/rl;
	tdcanc500=0;
	if anc500=1 and intxrel>timeanc500 then tdcanc500=1;
	tdcgvhd=0;
	if gvhd=1 and intxrel>tgvhd then tdcgvhd=1;
run;



*---------------------------------------------------------------;
*--------------------- Table 3.13 ------------------------------;
*---------------------------------------------------------------;

* Make "double" data set; 
data doubledod; 
	set bmt; 
	nygsource=bmonly; 
	nydisease=all; 
	age10=age/10;
	time=intxrel; 
	status=trm; 
	entry=0;

if gvhd=1 then do; 
	time=tgvhd; 
	status=0; 
end;
	age1=age10; 
	age2=0; 
	gsource1=nygsource; 
	gsource2=0;
	disease1=nydisease; 
	disease2=0; 
	stratum=1; 
output;

if gvhd=1 then do;
	time=intxrel; 
	status=trm; 
	entry=tgvhd;
	age1=0; 
	age2=age10; 
	gsource1=0; 
	gsource2=nygsource;
	disease1=0; 
	disease2=nydisease; 
	stratum=2; 
	output; end;
run;

* Row 1; 
data covar; 
	input disease1 disease2 age1 age2; 
	datalines; 
	0 0 0 0
; 
run; 
proc phreg data=doubledod; 
	model time*status(0)=disease1 disease2 
	age1 age2/entry=entry;
	strata stratum;
	baseline out=trm cumhaz=breslow covariates=covar;
run;

* Row 2; 
proc phreg data=doubledod; 
	model time*status(0)=disease1 disease2 
	age1 age2 type/entry=entry;
	type=stratum-1;
run;

* Row 3; 
proc phreg data=doubledod; 
	class stratum (ref='1') nydisease (ref='0');
	model time*status(0)=nydisease age10 type/entry=entry;
	type=stratum-1;
run;


*---------------------------------------------------------------;
*--------------------- Figure 3.10 -----------------------------;
*---------------------------------------------------------------;

proc gplot data=trm;
	plot breslow*time=stratum/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 140 by 10 minor=none label=('Months');
	axis2 order=0 to 0.7 by 0.1 minor=none label=(a=90 'Cumulative hazard of death in remission');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;;
run;

*---------------------------------------------------------------;
*--------------------- Table 3.14 ------------------------------;
*---------------------------------------------------------------;

* Log-normal; 
proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0") team;
	model intxrel*state0(0)=bmonly all age/rl;
	random team;
run;

* Gamma; 
proc phreg data=bmt;
	class bmonly(ref="0") all(ref="0") team;
	model intxrel*state0(0)=bmonly all age/rl;
	random team/dist=gamma;
run;