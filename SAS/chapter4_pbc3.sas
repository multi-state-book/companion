*------------------------------------------------------------------;
*------- Chapter 4, SAS code, PBC3 data ---------------------------;
*------------------------------------------------------------------;

* Load pbc data; 
proc import out=pbc3
datafile="data/pbc3.csv"
dbms=csv replace;
run;

data pbc3; 
set pbc3;
albnorm=(alb-35)*(alb>35);
alb10=alb/10; alb2=alb10*alb10;
bilihigh=(bili-17.1)*(bili>17.1);
bilitoohigh=(bili-34.2)*(bili>34.2);
bilimuchtoohigh=(bili-51.3)*(bili>51.3);
bili100=bili/100; bili2=bili100*bili100;
log2bili=log2(bili);
logbilihigh=(log2bili-log2(17.1))*(bili>17.1);
logbilitoohigh=(log2bili-log2(34.2))*(bili>34.2);
logbilimuchtoohigh=(log2bili-log2(51.3))*(bili>51.3);
log2bili2=log2bili*log2bili;
followup=days/365.25;
run;


*---------------------------------------------------------------;
*--------------------- Figure 4.2 ------------------------------;
*---------------------------------------------------------------;

* Kaplan-Meier plot per treatment; 

* Using proc lifetest;
proc lifetest data=pbc3 notable plots=(survival(nocensor));
time followup*status(0);
strata tment;
run;

* Using proc phreg;
proc phreg data=pbc3;
model days*status(0)=;
strata tment;
baseline out=survdat survival=km stderr=sd / method=pl;
run;

data survdat;
set survdat; 
daysyears = days/365.25; 
run; 

proc gplot data=survdat;
plot km*daysyears=tment/haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Survival probability');
symbol1  v=none i=stepjl c=red;
symbol2  v=none i=stepjl c=blue;
run;

* In-text: Kaplan-Meier estimates at 2 years;
proc lifetest data=pbc3 timelist=2;
time followup*status(0);
strata tment;
run;

proc print data=survdat;
where days<=2*365.25;
run;

*---------------------------------------------------------------;
*--------------------- Figure 4.4 ------------------------------;
*---------------------------------------------------------------;

data survdat; 
set survdat;
if days<=2 * 365.25 then
pwch=exp(-(days*(27.0000000/104856*(1-tment)+24.0000000/107931.5*tment)));
if 2 * 365.25 <days<=4 * 365.25 then
pwch=exp(-(2* 365.25*(27.0000000/104856*(1-tment)+24.0000000/107931.5*tment)
           +(days-2* 365.25)*(17.0000000/49673*(1-tment)+18.0000000/50284*tment)));
if 4 * 365.25 <days then
pwch=exp(-(2* 365.25 *(27.0000000/104856*(1-tment)+24.0000000/107931.5*tment)
           +(2* 365.25)*(17.0000000/49673*(1-tment)+18.0000000/50284*tment)
           +(days-4* 365.25)*(2.0000000/8642*(1-tment)+2.0000000/7599*tment)));
run;

data survdat;
set survdat; 
daysyears = days/365.25; 
run; 

proc gplot data=survdat; 
where tment=1;
plot (km pwch)*daysyears/overlay haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Survival probability');
symbol1  v=none i=stepjl r=1 c=red;
symbol2 v=none i=join r=1 c=blue;
run;
quit;



*---------------------------------------------------------------;
*--------------------- Figure 4.5 ------------------------------;
*---------------------------------------------------------------;

data cov;
tment=0; alb=38; log2bili=log2(45); output;
tment=1; alb=38; log2bili=log2(45); output;
run;

proc phreg data=pbc3;
model days*status(0)=tment alb log2bili/rl;
baseline out=predsurv survival=surv covariates=cov/ method=breslow;
run;

data predsurv;
set predsurv; 
daysyears = days/365.25; 
run; 

proc gplot data=predsurv;
plot surv*daysyears=tment/haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Estimated survival function');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
run;
quit;


*---------------------------------------------------------------;
*--------------------- Figure 4.6 ------------------------------;
*---------------------------------------------------------------;

data predsurv2; 
set predsurv; 
logminlogsurv = log(-log(surv)); 
run; 

data predsurv2;
set predsurv2; 
daysyears = days/365.25; 
run; 

proc gplot data=predsurv2;
plot logminlogsurv*daysyears=tment/haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=-7 to 0 by 1 minor=none label=(a=90 'log(-log(survival function))');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
run;
quit;




*---------------------------------------------------------------;
*--------------------- Figure 4.7 ------------------------------;
*---------------------------------------------------------------;

proc phreg data=pbc3;
class tment;
model days*status(0)=tment alb log2bili/rl;
baseline out=gsurv survival=surv stderr=sd / method=breslow diradj group=tment;
run;

data gsurv;
set gsurv; 
daysyears = days/365.25; 
run; 

proc gplot data=gsurv;
plot surv*daysyears=tment/haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Estimated survival function (g-formula)');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
run;
quit;

proc print data=gsurv;
where days=714;
run;

*---------------------------------------------------------------;
*--------------------- Figure 4.8 ------------------------------;
*---------------------------------------------------------------;

proc phreg data=pbc3;
*class tment;
model days*status(0)=tment alb log2bili/rl;
baseline out=predsurvpl survival=surv covariates=cov / method=pl;
run;

data predsurvpl;
set predsurvpl; 
daysyears = days/365.25; 
run; 

proc gplot data=predsurvpl;
plot surv*daysyears=tment/haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Estimated survival function');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
run;
quit;



*---------------------------------------------------------------;
*--------------------- Table 4.1 -------------------------------;
*---------------------------------------------------------------;

proc sort data=pbc3 out=pbc3sorted;
by tment;
run;
proc rmstreg data=pbc3sorted tau=3;
by tment;
model followup*status(0)= / 
  link=linear method=ipcw(strata=tment);
run;
* Bootstrapping using 'point='; 
data bootpbc;
do sampnum = 1 to 1000; /* nboot=1000*/
  do i = 1 to 349;      /*nobs=349*/
  x=round(ranuni(0)*349); 
set pbc3
point=x;
output;
end;
end;
stop;
run;
proc sort data=bootpbc out=boot;
by sampnum tment;
run;
proc rmstreg data=boot tau=3;
by sampnum tment;
model followup*status(0)= / 
  link=linear method=ipcw(strata=tment);
ods output parameterestimates=pe;
run;
proc means data=pe mean stddev;
class tment;
var estimate;
run;

* AUC under stepcurves; 
%macro areastepby(data,byvar,trt,grp,time,y,tau);
data select;
set &data;
where &trt=&grp;
run;

data select;
set select;
by &byvar;
retain mu oldt oldy;
if first.&byvar then do;
oldt=0; oldy=1; mu=0;
end;
if &time>&tau then do;
&time=&tau; &y=oldy;
end;
if not first.&byvar then mu+oldy*(&time-oldt);
if last.&byvar then do;
if &time<&tau then mu+(&tau-&time)*&y; end;
oldy=&y; oldt=&time;
run;

data last; set select;
by  &byvar;
if last.&byvar;
run;
%mend areastepby;

proc phreg data=bootpbc noprint;
by sampnum;
model followup*status(0)=;
strata tment;
baseline out=survdat survival=km / method=pl;
run;

* Restricted mean;
%areastepby(survdat,sampnum,tment,0,followup,km,3);
title "Placebo";
proc means data=last mean stddev;
var mu;
run;
%areastepby(survdat,sampnum,tment,1,followup,km,3);
title "CyA";
proc means data=last mean stddev;
var mu;
run;

/* pred surv given covariates */
  data cov;
tment=0; alb=38; log2bili=log2(45); output;
tment=1; alb=38; log2bili=log2(45); output;
run;
proc phreg data=bootpbc noprint;
by sampnum;
model followup*status(0)=tment alb log2bili/rl;
baseline out=predsurv survival=surv covariates=cov/ method=breslow;
run;

%areastepby(predsurv,sampnum,tment,0,followup,surv,3);
title "Placebo";
proc means data=last mean stddev;
var mu;
run;
%areastepby(predsurv,sampnum,tment,1,followup,surv,3);
title "CyA";
proc means data=last mean stddev;
var mu;
run;

data cov;
tment=0; alb=20; log2bili=log2(90); output;
tment=1; alb=20; log2bili=log2(90); output;
run;
proc phreg data=bootpbc noprint;
by sampnum;
model followup*status(0)=tment alb log2bili/rl;
baseline out=predsurv2 survival=surv covariates=cov/ method=breslow;
run;

%areastepby(predsurv2,sampnum,tment,0,followup,surv,3);
title "Placebo";
proc means data=last mean stddev;
var mu;
run;
%areastepby(predsurv2,sampnum,tment,1,followup,surv,3);
title "CyA";
proc means data=last mean stddev;
var mu;
run;

/* g-formel */
  
  proc phreg data=bootpbc noprint;
by sampnum;
class tment (ref='0');
model followup*status(0)=tment alb log2bili/rl;
baseline out=gsurv survival=surv stderr=se/ method=breslow diradj group=tment;
run;

%areastepby(gsurv,sampnum,tment,0,followup,surv,3);
title "Placebo";
proc means data=last mean stddev;
var mu;
run;
%areastepby(gsurv,sampnum,tment,1,followup,surv,3);
title "CyA";
proc means data=last mean stddev;
var mu;
run;


*---------------------------------------------------------------;
*--------------------- Figure 4.10 -----------------------------;
*---------------------------------------------------------------;

proc phreg data=pbc3 noprint;
model days*status(0)=;
baseline out=overallsurv survival=surv / method=;
run;
proc sort;by tment days;run;
proc phreg data=pbc3 noprint;
model days*status(0)=/eventcode=1;
strata tment;
baseline out=cuminc1 cif=cif1;
run;
proc sort;by tment days;run;
proc phreg data=pbc3 noprint;
model days*status(0)=/eventcode=2;
strata tment;
baseline out=cuminc2 cif=cif2;
run;
proc sort;by tment days;run;
data plot0; 
merge overallsurv cuminc1 cuminc2; 
by tment days;
if tment=0;
fail=1-surv;
run;

data plot0ok; 
set plot0;
by days;
retain last0 last1 last2;
if fail=. then c0=last0; if fail ne . then c0=fail;
if cif1=. then c1=last1; if cif1 ne . then c1=cif1;
if cif2=. then c2=last2; if cif2 ne . then c2=cif2;
output;
last0=c0; last1=c1; last2=c2;
run;
data plot0ok; 
set plot0ok;
cum1=c1; 
cum2=c1+c2; 
cum3=c1+c2+(1-c0);
run;
proc print;run;

data plot0ok; 
set plot0ok; 
daysyears = days/365.25; 
run; 

proc gplot 
data=plot0ok;
plot cum1*daysyears cum2*daysyears cum3*daysyears /overlay haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Stacked cumulative incidence and survival');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
symbol3  v=none i=stepjl c=black;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Figure 4.11 -----------------------------;
*---------------------------------------------------------------;

proc phreg data=pbc3 noprint;
model days*status(0 2)=;
strata tment;
baseline out=cuminc1wrng survival=s1wrng;
run;

proc phreg data=pbc3 noprint;
model days*status(0 1)=;
baseline out=cuminc2wrng survival=s2wrng;
strata tment;
run;

data plot0wrng; 
merge overallsurv cuminc1wrng cuminc2wrng; 
where tment=0; by days;
run;

data plot0wrng; 
set plot0wrng;
by days;
retain last1 last2;
if s1wrng=. then c1=last1; if s1wrng ne . then c1=1-s1wrng;
if s2wrng=. then c2=last2; if s2wrng ne . then c2=1-s2wrng;
output;
last1=c1; last2=c2;
run;

data plot0wrng; 
set plot0wrng;
cum1=c1; cum2=c1+c2; cum3=c1+c2+surv; one=1;
run;

data plot0wrng; 
set plot0wrng; 
daysyears = days/365.25; 
run; 

proc gplot 
data=plot0wrng;
plot cum1*daysyears cum2*daysyears cum3*daysyears one*daysyears
/overlay haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=0 to 1.1 by 0.1 minor=none label=(a=90 'Stacked cumulative incidence and survival');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
symbol3  v=none i=stepjl c=black;
symbol4  v=none i=stepjl l=2 c=black;
run;
quit;


*---------------------------------------------------------------;
*--------------------- Figure 4.12 -----------------------------;
*---------------------------------------------------------------;
filename cumincpr url 'https://biostat.ku.dk/pka/CUMINC/cuminc.sas';
%inc cumincpr;

* Macro and double data set;
%macro CumInc(Data,Strata,Time,Surv);

/* Number of Stratas */
  
  proc sort data=&Data;
by &Strata;

data nstrat;
set &Data end=last;
by &Strata;
firstS=first.&Strata;
retain nStrata 0;
nStrata+firstS;
if last then call symput('nStrata',nStrata);
drop firstS;
run;


/* Observations in ciData are deleted */
  
  
  data newData; 
set &Data;
start=(&Time eq 0);   
retain komb 0;
komb+start;
med=0;
%do k=0 %to (&nStrata-1) %by 1;
med=med+(komb eq (1+&k*(&nStrata+1)));
%end;
if (med eq 1);  
drop start med komb;
run;

/* Time vector */
  
  proc sort data=newData;
by &Time;

data Time (keep=&Time); 
set newData;
%do i=1 %to &nStrata %by 1;
if &Time=lag(&Time) then delete;
%end;


/* Strata vector */
  
  proc sort data=newData;
by &Strata &Time;

data temp1;
set newData; by &Strata;
retain stratum 0;
stratum+first.&Strata;

/*  A=sum(A_1,...,A_nStrata) */
  
  %do i=1 %to &nStrata %by 1;
data data&i;
set temp1;
if stratum=&i;

data data&i (keep = &Time A&i);
merge data&i Time; by &Time;  
retain Surv&i;
if not (&Surv=.) then Surv&i=&Surv;
A&i=-log(Surv&i);

%end;


data A_and_S;
A=0;
%do i=1 %to &nStrata %by 1;
merge data&i; by &Time;
A+A&i;
%end;

dA=A-lag(A);
if dA eq . then dA=0;

retain S 1;
S+S*(1-dA)-S;
lagS=lag(S);
run;



/* Output */
  
  data temp2 (keep = &Strata &Time);
set temp1; 
proc sort data=temp2; by &Time;

data temp3;
merge temp2 A_and_S; by &Time;
proc sort data=temp3; by &Strata &Time;


data data0 (keep = &Time p stratum);
set temp3; by &Strata;
*   lagS=lag(S);

if first.&Strata then do;
lagS=1; end;

retain stratum 0;
stratum+first.&Strata;

%do i=1 %to &nStrata %by 1;
lagA&i=lag(A&i);
if ((stratum eq &i) and (first.&Strata)) then lagA&i=0;
if (stratum ne &i) then lagA&i=0;

l&i=(stratum=&i)*(lagS*(A&i-lagA&i));
retain p&i 0;
p&i+l&i; 
if (stratum ne &i) then p&i=0;
%end;

p=0;
%do i=1 %to &nStrata %by 1;
p=p+p&i; %end;
run;





/* put data into the right form */
  
  %do i=1 %to &nStrata;
data data&i;
set data0;
if (stratum eq &i);
p0&i=p; drop p stratum;

proc sort data=data&i;
by &Time;
%end;

data data;
set data1;
%do i=1 %to &nStrata;
data data;
merge data data&i;
by &Time;
%end;


/* ... and complete the p's */

proc sort data=data;
 by &Time;

data data; 
 set data end=last;
  by &Time;
   n=_N_;
  if last then call symput('nobs',n);
  drop n;


data data;
 set data;
  %do i=1 %to &nStrata;
   %do j=1 %to &nobs;
    dummy=lag(p0&i);
    if (p0&i eq .) then p0&i=dummy;
   %end;
  %end;
  drop dummy;
  p00=1;
  %do i=1 %to &nStrata;
   p00=p00-p0&i;
  %end;
run;
%mend CumInc;



data pbc32; 
	set pbc3 pbc3;
	h=1+(_N_ gt 349);
	time=days;
	d=(status=1)*(h=1)+(status=2)*(h=2);
	tment1=tment*(h=1); tment2=tment*(h=2);
	sex1=sex*(h=1); sex2=sex*(h=2);
	age1=age*(h=1); age2=age*(h=2);
	alb1=alb*(h=1); alb2=alb*(h=2);
	log2bili1=log2bili*(h=1); log2bili2=log2bili*(h=2);
run;


*---------------------------------------------------------------;
*--------------------- Figure 4.12a ----------------------------;
*---------------------------------------------------------------;

data cov;
input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
datalines; 
0 0 0 0 40 0 38 0 45 0
0 0 0 0 0 40 0 38 0 45
;
run;

data cov; 
set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;
proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	                log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data data; set data; tment=0; run;
data plot1; 
	set data; 
	cum1=p01; cum2=p01+p02; cum3=1;
	p = p01+p02+p00;
run;
data plot1; 
	set plot1; 
	daysyears = time/365.25; 
run; 
proc gplot data=plot1;
	plot cum1*daysyears cum2*daysyears cum3*daysyears/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Stacked cumulative incidence and survival');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
run;
quit;


*---------------------------------------------------------------;
*--------------------- Figure 4.12b ----------------------------;
*---------------------------------------------------------------;


data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	0 0 0 0 40 0 20 0 90 0
	0 0 0 0 0 40 0 20 0 90
	;
run;
data cov; 
set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;
proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	                log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data data; set data; tment=0; run;
data plot2plac; set data; 
	cum1=p01; cum2=p01+p02; cum3=1;
run;
data plot2plac; 
	set plot2plac; 
	daysyears = time/365.25; 
run; 
proc gplot data=plot2plac;
	plot cum1*daysyears cum2*daysyears cum3*daysyears/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Stacked cumulative incidence and survival');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Figure 4.12c ----------------------------;
*---------------------------------------------------------------;

data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	0 0 0 0 60 0 38 0 45 0
	0 0 0 0 0 60 0 38 0 45
	;
run;
data cov; set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;
proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	                log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data data; set data; tment=0; run;
data plot3plac; set data; 
	cum1=p01; cum2=p01+p02; cum3=1;
run;
data plot3plac; 
	set plot3plac; 
	daysyears = time/365.25; 
run; 
proc gplot data=plot3plac;
	plot cum1*daysyears cum2*daysyears cum3*daysyears/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Stacked cumulative incidence and survival');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
run;
quit;


/* USED? */
proc phreg data=pbc3;
model followup*status(0)=;
strata tment;
baseline out=overallsurv survival=surv / method=pl;
run;

*---------------------------------------------------------------;
*--------------------- Table 4.2 -------------------------------;
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
	proc print; run;
%mend areastep;

* No adjustment: Transplantation; 
proc phreg data=pbc3 noprint;
	model followup*status(0)=/eventcode=1;
	strata tment;
	baseline out=cuminc1 cif=cif1 stdcif=std1;
run;
title "Placebo";
%areastep(cuminc1,tment,0,followup,cif1,3);
title "CyA";
%areastep(cuminc1,tment,1,followup,cif1,3);

* No adjustment: Death wo transplant;
proc phreg data=pbc3 noprint;
	model followup*status(0)=/eventcode=2;
	strata tment;
	baseline out=cuminc2 cif=cif2 stdcif=std2;
run;
title "Placebo";
%areastep(cuminc2,tment,0,followup,cif2,3);
title "CyA";
%areastep(cuminc2,tment,1,followup,cif2,3);

* Cuminc macro requires duplication; 
data pbc32; 
	set pbc3 pbc3;
	h=1+(_N_ gt 349);
	time=followup;
	d=(status=1)*(h=1)+(status=2)*(h=2);
	tment1=tment*(h=1); tment2=tment*(h=2);
	sex1=sex*(h=1); sex2=sex*(h=2);
	age1=age*(h=1); age2=age*(h=2);
	alb1=alb*(h=1); alb2=alb*(h=2);
	log2bili1=log2bili*(h=1); log2bili2=log2bili*(h=2);
run;


*** age=40 alb=38 bili=45 ***; 
* CyA; 
data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	1 0 0 0 40 0 38 0 45 0
	0 1 0 0 0 40 0 38 0 45
	;
run;

data cov; set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;

proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data dataC; 
	set data; 
	tment=1; 
run;
%areastep(dataC,tment,1,time,p01,3);
%areastep(dataC,tment,1,time,p02,3);

* Placebo; 
data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	0 0 0 0 40 0 38 0 45 0
	0 0 0 0 0 40 0 38 0 45
	;
run;
data cov; set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;
proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
									log2bili1 log2bili2;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data dataP; 
set data; 
tment=0; 
run;
%areastep(dataP,tment,0,time,p01,3);
%areastep(dataP,tment,0,time,p02,3);

*** age=40 alb=20 bili=90 ***; 
* CyA; 
data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	1 0 0 0 40 0 20 0 90 0
	0 1 0 0 0 40 0 20 0 90
	;
run;

data cov; set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;

proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;

%cuminc(cidata,h,time,surv);

data dataC; 
	set data; 
	tment=1; 
run;

%areastep(dataC,tment,1,time,p01,3);
%areastep(dataC,tment,1,time,p02,3);

* Placebo; 
data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	0 0 0 0 40 0 20 0 90 0
	0 0 0 0 0 40 0 20 0 90
	;
run;

data cov; set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;
proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data dataP; 
	set data; 
	tment=0; 
run;
%areastep(dataP,tment,0,time,p01,3);
%areastep(dataP,tment,0,time,p02,3);

*** age=60 alb=38 bili=45 ***; 
* CyA; 
data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	1 0 0 0 60 0 38 0 45 0
	0 1 0 0 0 60 0 38 0 45
	;
run;
data cov; set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;
proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data dataC; 
	set data; 
	tment=1; 
run;
%areastep(dataC,tment,1,time,p01,3);
%areastep(dataC,tment,1,time,p02,3);

* Placebo; 
data cov;
	input tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 log2bili1 log2bili2;
	datalines; 
	0 0 0 0 60 0 38 0 45 0
	0 0 0 0 0 60 0 38 0 45
	;
run;
data cov; set cov;
	if log2bili1>0 then log2bili1=log2(log2bili1);
	if log2bili2>0 then log2bili2=log2(log2bili2);
run;
proc phreg data=pbc32 noprint;
	model time*d(0)=tment1 tment2 sex1 sex2 age1 age2 alb1 alb2 
	log2bili1 log2bili2/rl;
	strata h;
	baseline out=cidata covariates=cov survival=surv/nomean method=ch;
run;
%cuminc(cidata,h,time,surv);
data dataP; 
	set data; 
	tment=0; 
run;
%areastep(dataP,tment,0,time,p01,3);
%areastep(dataP,tment,0,time,p02,3);



*---------------------------------------------------------------;
*--------------------- Table 4.4 -------------------------------;
*---------------------------------------------------------------;

proc rmstreg data=pbc3 tau=3;
  model followup*status(0)=tment  / 
        link=linear method=ipcw(strata=tment);
run;

proc rmstreg data=pbc3 tau=3;
  model followup*status(0)=tment log2bili / 
        link=linear method=ipcw(strata=tment);
run;

proc rmstreg data=pbc3 tau=3;
  model followup*status(0)=tment alb log2bili / 
        link=linear method=ipcw(strata=tment);
run;



*---------------------------------------------------------------;
*--------------------- Table 4.5 -------------------------------;
*---------------------------------------------------------------;

* Death without transplantation; 
proc phreg data=pbc3;
	class sex (ref='1') tment (ref='0');
	model days*status(0)=sex tment age log2bili alb / rl eventcode=2;
run;


* Transplantation; 
proc phreg data=pbc3;
	class sex (ref='1') tment (ref='0');
	model days*status(0)=sex tment age log2bili alb / rl eventcode=1;
run;


*---------------------------------------------------------------;
*--------------------- Figure 4.17a -----------------------------;
*---------------------------------------------------------------;

* Covariates; 
data cov0;
	input sex tment age alb log2bili;
	log2bili=log2(log2bili);
	datalines;
	0 0 40 38 45
	;
run;
data cov1;
	input sex tment age alb log2bili;
	log2bili=log2(log2bili);
	datalines;
	0 1 40 38 45
	;
run;

* For death without transplantation; 
proc phreg data=pbc3;
	model days*status(0)=sex tment age log2bili alb/eventcode=2;
	baseline out=cuminc20 covariates=cov0 cif=cif20;
run;
proc phreg data=pbc3;
	model days*status(0)=sex tment age log2bili alb/eventcode=2;
	baseline out=cuminc21 covariates=cov1 cif=cif21;
run;
data cuminc2; 
	merge cuminc20 cuminc21; 
	by days; 
run;
data cuminc2; 
	set cuminc2; 
	daysyears = days / 365.25; 
run; 
proc gplot data=cuminc2;
	plot cif20*daysyears cif21*daysyears/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Cumulative incidences for death');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;
quit;



*---------------------------------------------------------------;
*--------------------- Figure 4.17b -----------------------------;
*---------------------------------------------------------------;

* For transplantation; 
proc phreg data=pbc3;
	*class sex tment (ref='0');
	model days*status(0)=sex tment age log2bili alb/eventcode=1;
	baseline out=cuminc10 covariates=cov0 cif=cif10;
run;
proc phreg data=pbc3;
	*class sex tment (ref='0');
	model days*status(0)=sex tment age log2bili alb/eventcode=1;
	baseline out=cuminc11 covariates=cov1 cif=cif11;
run;
data cuminc1; 
	merge cuminc10 cuminc11; 
	by days; 
run;
data cuminc1; 
	set cuminc1; 
	daysyears = days / 365.25; 
run; 
proc gplot data=cuminc1;
	plot cif10*daysyears cif11*daysyears/overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none label=(a=90 'Cumulative incidence for transpl.');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;
quit;

*---------------------------------------------------------------;
*---------------- Figure 4.21 and cox models for censoring------;
*---------------------------------------------------------------;

proc phreg data=pbc3 atrisk noprint;
	model followup*status(1 2)=;
	baseline out=survcens survival=kmc / method=pl;
run;
proc phreg data=pbc3 atrisk noprint;
	model followup*status(0)=;
	baseline out=survdat survival=kms / method=pl;
run;
data plotsurv; 
	merge survcens survdat; 
	by followup;
run;
data plotfin; set plotsurv;
	by followup;
	retain lastc lasts;
	if kmc=. then do c=lastc; s=kms; end; 
	if kms=. then do s=lasts; c=kmc; end;
	if kmc ne . and kms ne . then do; c=kmc; s=kms; end;
	output;
	lastc=c; lasts=s;
run;
proc gplot data=plotfin;
	plot c*followup s*followup/haxis=axis1 vaxis=axis2 overlay;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 1 by 0.1 minor=none
	label=(a=90 '"Survival" probability');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;
quit;

proc phreg data=pbc3;
	class tment (ref='0');
	model followup*status(1 2)=tment/rl;
run;
proc phreg data=pbc3;
	model followup*status(1 2)=alb/rl;
run;
proc phreg data=pbc3;
	model followup*status(1 2)=bili/rl;
run;
