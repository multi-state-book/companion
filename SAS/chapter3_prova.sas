*------------------------------------------------------------------;
*------- Chapter 3, SAS code, PROVA data -------------------------;
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
*--------------------- Table 3.3 -------------------------------;
*---------------------------------------------------------------;

* Table 3.3 column 1;
* Variceal bleeding; 
proc phreg data=provany;
	class beh (ref='0');
	model btime*bleed(0)=beh;
run;

* Death without bleeding; 
proc phreg data=provany;
	class beh (ref='0');
	model d0time*dead0(0)=beh;
run;

* Table 3.3 column 2;
* Variceal bleeding;
proc phreg data=provany;
	class beh (ref='0') varsize (ref='1');
	model btime*bleed(0)=beh sex coag log2bili varsize;
run;

* Death without bleeding; 
proc phreg data=provany;
	class beh (ref='0') varsize (ref='1');
	model d0time*dead0(0)=beh sex coag log2bili varsize;
run;

* In text; 
proc phreg data=provany;
	class beh (ref='0');
	model d0time*dead0(0)=beta scle;
run;

proc phreg data=provany;
	class beh (ref='0');
	model d0time*dead0(0)=scle;
run;

* LR test; 
data p;
LRT1=468.379-466.752;
p1=1-probchi(LRT1,1); 
LRT2=468.727-468.379;
p2=1-probchi(LRT2,1);
proc print;
run;


*---------------------------------------------------------------;
*--------------------- Table 3.4 -------------------------------;
*---------------------------------------------------------------;

proc phreg data=provany;
	class beh (ref='0') varsize (ref='1');
	model btime*outof0(0)=beh sex coag log2bili varsize;
run;

*---------------------------------------------------------------;
*--------------------- Table 3.7 -------------------------------;
*---------------------------------------------------------------;

** Column 1; 

data covar; 
	input beh sex log2bili; 
	datalines; 
	0 0 0
; 
run; 

* Time since randomisation;
proc phreg data=provany atrisk;
	class beh (ref='0');
	model bdtime*deadb(0)=beh sex log2bili/entry=btime rl;
	baseline out=cumhaztime cumhaz=breslowtime covariates=covar;
run;

* Duration; 
proc phreg data=provany;
	class beh (ref='0');
	model wait*deadb(0)=beh sex log2bili;
	baseline out=cumhazwait cumhaz=breslowwait covariates=covar;
run;


** Column 2; 

* Time since randomisation;
proc phreg  data=provany;
	class beh (ref='0');
	model bdtime*deadb(0)=beh sex wait1 wait2 log2bili/entry=btime rl;
	wait1=0; wait2=0;
	if (bdtime-btime<5) then wait1=1;
	if (5<=bdtime-btime<10) then wait2=1;
	duration: test wait1=0, wait2=0;
run;

* Duration;
proc phreg data=provany;
	class beh (ref='0');
	model wait*deadb(0)=beh sex time1 time2 log2bili;
	time1=0; time2=0;
	if (btime+wait<365.25) then time1=1;
	if (365.25<=btime+wait<2*365.25) then time2=1;
	timeeff: test time1=0, time2=0;
run;


*---------------------------------------------------------------;
*--------------------- Figure 3.2 ------------------------------;
*---------------------------------------------------------------;

data cumhaztime;
	set cumhaztime; 
	bdtimeyears = bdtime / 365.25; 
run;

proc gplot data=cumhaztime;
	plot breslowtime*bdtimeyears/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 4 by 1 minor=none 
	      label=('Time since randomization (Years)');
	axis2 order=0 to 4 by 1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=blue;
run;


*---------------------------------------------------------------;
*--------------------- Figure 3.3 ------------------------------;
*---------------------------------------------------------------;

data cumhazwait;
	set cumhazwait; 
	waityears = wait / 365.25; 
run;


proc gplot data=cumhazwait;
	plot breslowwait*waityears/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 3 by 1 minor=none label=('Duration (Years)');
	axis2 order=0 to 1.5 by 0.5 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=blue;
run;



*---------------------------------------------------------------;
*--------------------- Table 3.8 -------------------------------;
*---------------------------------------------------------------;

data provasplit1;
	set provany;
	where bleed=1;
	fail=(wait<5)*(deadb ne 0);
	risktime=min(5,wait);
	logrisk=log(risktime); wint=1; 
	start=btime; slut=btime+min(5,wait); output;  
	if wait>=5 then do;
	fail=(wait<10)*(deadb ne 0);
	risktime=min(5,wait-5);
	logrisk=log(risktime); wint=2; 
	start=btime+5; slut=btime+min(10,wait); output; end;
	if wait>10 then do;
	fail=deadb ne 0; 
	risktime=wait-10;
	logrisk=log(risktime); wint=3; 
	start=btime+10; slut=btime+wait; output; end;
run;

data provasplit2; 
	set provasplit1;
	if start<365.25 then do; risktime2=min(slut,365.25)-start;
	fail2=fail*(slut<365.25); logrisk2=log(risktime2); tint=1; output;
	if slut>365.25 then do; risktime2=min(slut,2*365.25)-365.25; logrisk2=log(risktime2);
	fail2=fail*(slut<2*365.25); tint=2; output; end;
	if slut>2*365.25 then do; risktime2=slut-2*365.25; logrisk2=log(risktime2);
	fail2=fail; tint=3; output; end;
	end;
	if 365.25<=start<2*365.25 then do; risktime2=min(slut,2*365.25)-start;
	fail2=fail*(slut<2*365.25); logrisk2=log(risktime2); tint=2; output;
	if slut>2*365.25 then do; risktime2=slut-2*365.25; logrisk2=log(risktime2);
	fail2=fail; tint=3; output; end;
	end;
	if start>=2*365.25 then do; risktime2=slut-start; logrisk2=log(risktime2);
	fail2=fail; tint=3; output; 
end;
run;

* Table 3.8 output; 
data provasplit2; 
	set provasplit2;
	risktime2ys=risktime2/365.25;
run;

proc means data=provasplit2 sum; 
	class wint tint;
	var fail2 risktime2ys;
run;


*---------------------------------------------------------------;
*--------------------- Table 3.9 -------------------------------;
*---------------------------------------------------------------;

* Column 1; 
proc genmod data=provasplit1;
	class beh (ref='0') wint;
	model fail=beh wint sex log2bili/dist=poi offset=logrisk type3;
run;

* Column 2; 
proc genmod data=provasplit2;
	class beh (ref='0') tint;
	model fail2=beh tint sex log2bili/dist=poi offset=logrisk2 type3;
run;

* Column 3; 
proc genmod data=provasplit2;
	class beh (ref='0') wint tint;
	model fail2=beh wint tint sex log2bili/dist=poi offset=logrisk2 type3;
run;

* Interaction model, in-text;
proc genmod data=provasplit2;
	class beh (ref='0') wint tint;
	model fail2=beh wint tint wint*tint sex log2bili/dist=poi offset=logrisk2 type3;
run;


*---------------------------------------------------------------;
*--------------------- Table 3.12 ------------------------------;
*---------------------------------------------------------------;

* Prepare data as described; 
data double; 
set provany; 
	time=d0time; 
	status=dead0; 
	entrytime=0; 
	sex1=sex; 
	sex2=0;
	age1=age; 
	age2=0; 
	bili1=log2bili; 
	bili2=log2bili*0; 
	stratum=1; 
		output;
	if bleed=1 then do;
	time=bdtime; 
	status=deadb; 
	entrytime=btime; 
	sex1=0; 
	sex2=sex;
	age1=0; 
	age2=age; 
	bili1=log2bili*0; 
	bili2=log2bili; 
	stratum=2; 
		output; end;
run;

* Row 1; 
data covar; 
	input sex1 sex2 bili1 bili2; 
	datalines; 
	0 0 0 0
; 
run; 
	
proc phreg data=double;
	model time*status(0)=sex1 sex2 bili1 bili2/entry=entrytime;
	strata stratum;
	baseline out=mort cumhaz=breslow covariates=covar;
run;

* Row 2; 
proc phreg data=double; 
	model time*status(0)=sex bili1 bili2/entry=entrytime;
	strata stratum;
run;

* Row 3; 
proc phreg data=double; 
	model time*status(0)=sex bili1/entry=entrytime;
	strata stratum;
run;


*---------------------------------------------------------------;
*--------------------- Figure 3.9 ------------------------------;
*---------------------------------------------------------------;

data mort; 
	set mort; 
	timeyears = time /365.25; 
run;


proc gplot data=mort; 
	plot breslow*timeyears=stratum/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 4 by 1 minor=none 
	      label=('Time since randomization (Years)');
	axis2 order=0 to 4 by 1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;



