*------------------------------------------------------------------;
*------- Chapter 3, SAS code, PROVA data -------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=prova
	datafile="data/prova.csv"
	dbms=csv replace;
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
* logrank test is the score test;
proc phreg data=provany;
	class beh (ref='0');
	model btime*bleed(0)=beh  / type3(lr);
run;

* Death without bleeding; 
* logrank test is the score test;
proc phreg data=provany;
	class beh (ref='0');
	model d0time*dead0(0)=beh / type3(lr);
run;
* Death without bleeding - in text LRT for additive model; 
proc phreg data=provany;
	model d0time*dead0(0)=scle|beta / type3(lr);
	estimate 'both' scle 1 beta 1 scle*beta 1;
run;
proc phreg data=provany;
	model d0time*dead0(0)=scle beta / type3(lr);
run;
* Death without bleeding - remove propranolol;
proc phreg data=provany;
	model d0time*dead0(0)=scle / type3(lr);
run;

* Table 3.3 column 2;
* Variceal bleeding;
proc phreg data=provany;
	class beh (ref='0') varsize (ref='1');
	model btime*bleed(0)=beh sex coag log2bili varsize / type3(lr);
run;

* Death without bleeding; 
proc phreg data=provany;
	class beh (ref='0') varsize (ref='1');
	model d0time*dead0(0)=beh sex coag log2bili varsize / type3(lr);
run;

*---------------------------------------------------------------;
*--------------------- Table 3.4 -------------------------------;
*---------------------------------------------------------------;

proc phreg data=provany;
	class beh (ref='0') varsize (ref='1');
	model btime*outof0(0)=beh sex coag log2bili varsize  / type3(lr);
run;

*---------------------------------------------------------------;
*--------------------- Table 3.8 -------------------------------;
*---------------------------------------------------------------;

* Time since randomisation;
* Column 1; 
proc phreg data=provany atrisk;
	class beh (ref='0');
	model bdtime*deadb(0)=beh sex log2bili / entry=btime rl type3(lr);
run;
* LRT interaction scle*beta;
proc phreg data=provany atrisk;
	class beh (ref='0');
	model bdtime*deadb(0)=scle|beta sex log2bili / entry=btime rl type3(lr);
run;
* Column 2; 
proc phreg  data=provany;
	class beh (ref='0');
	model bdtime*deadb(0)=beh sex log2bili wait1 wait2 / entry=btime rl type3(lr);
	wait1=0; wait2=0;
	if (bdtime-btime<5) then wait1=1;
	if (5<=bdtime-btime<10) then wait2=1;
	duration: test wait1=0, wait2=0;
run;

* In text: Linear effect of time-dependent covariate;
proc phreg  data=provany;
	class beh (ref='0');
	model bdtime*deadb(0)=beh sex log2bili lin / entry=btime rl type3(lr);
	lin=bdtime-btime;
run;


* Duration; 
* Column 1; 
proc phreg data=provany;
	class beh (ref='0');
	model wait*deadb(0)=beh sex log2bili / type3(lr);
	baseline out=cumhazwait cumhaz=breslowwait covariates=covar;
run;

* Column 2;
proc phreg data=provany;
	class beh (ref='0');
	model wait*deadb(0)=beh sex log2bili time1 time2 / type3(lr);
	time1=0; time2=0;
	if (btime+wait<365.25) then time1=1;
	if (365.25<=btime+wait<2*365.25) then time2=1;
	timeeff: test time1=0, time2=0;
run;


*---------------------------------------------------------------;
*--------------------- Figure 3.2 ------------------------------;
*---------------------------------------------------------------;

data covar; 
	input beh sex log2bili; 
	datalines; 
	0 0 0
; 
run; 
proc phreg data=provany atrisk;
	class beh (ref='0');
	model bdtime*deadb(0)=beh sex log2bili/entry=btime rl type3(lr);
	baseline out=cumhaztime cumhaz=breslowtime covariates=covar;
run;
data cumhaztime;
	set cumhaztime; 
	bdtimeyears = bdtime / 365.25; 
run;
proc gplot data=cumhaztime;
	plot breslowtime*bdtimeyears/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 4 by 1 minor=none 
	      label=('Time since randomization (Years)');
	axis2 order=0 to 1.1 by 0.1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=blue;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Figure 3.3 ------------------------------;
*---------------------------------------------------------------;

* Duration; 
data covar; 
	input beh sex log2bili; 
	datalines; 
	0 0 0
; 
run; 
proc phreg data=provany;
	class beh (ref='0');
	model wait*deadb(0)=beh sex log2bili / type3(lr);
	baseline out=cumhazwait cumhaz=breslowwait covariates=covar;
run;
data cumhazwait;
	set cumhazwait; 
	waityears = wait / 365.25; 
run;
proc gplot data=cumhazwait;
	plot breslowwait*waityears/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 4 by 1 minor=none label=('Duration (Years)');
	axis2 order=0 to 0.2 by 0.05 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=blue;
run;
quit;


*---------------------------------------------------------------;
*--------------------- Table 3.9 -------------------------------;
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
data provasplit2; 
	set provasplit2;
	risktime2ys=risktime2/365.25;
run;
proc means data=provasplit2 sum; 
	class wint tint;
	var fail2 risktime2ys;
run;


*---------------------------------------------------------------;
*--------------------- Table 3.10 -------------------------------;
*---------------------------------------------------------------;

* part (a);
proc genmod data=provasplit1;
	class beh (ref='0') wint;
	model fail=beh wint sex log2bili/dist=poi offset=logrisk type3;
run;

* part (b);
proc genmod data=provasplit2;
	class beh (ref='0') tint;
	model fail2=beh tint sex log2bili/dist=poi offset=logrisk2 type3;
run;

* part (c);
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
*--------------------- Table 3.13 ------------------------------;
*---------------------------------------------------------------;

* Prepare data set for analysis - double; 
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
	bleeding=1; 
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
		bleeding=2; 
		output; 
	end;
run;

* part (a);
proc phreg data=double;
	model time*status(0)=sex1 sex2 bili1 bili2 /entry=entrytime type3(lr);
	strata bleeding;
	test sex1=sex2; /* wald tests instead of LRT */
	test bili1=bili2;
run;

* part (b);
proc phreg data=double; 
	model time*status(0)=sex bili1 bili2 /entry=entrytime type3(lr);
	strata bleeding;
run;

* part (c);
proc phreg data=double; 
	model time*status(0)=sex bili1 /entry=entrytime type3(lr);
	strata bleeding;
run;

* In-text: LRT proportionality;
proc phreg data=double; 
  bleedinglogt=bleeding*log(time);
	model time*status(0)=sex bili1 bleeding bleedinglogt /entry=entrytime type3(lr);
run;

*---------------------------------------------------------------;
*--------------------- Figure 3.9 ------------------------------;
*---------------------------------------------------------------;

data covar; 
	input sex1 sex2 bili1 bili2; 
	datalines; 
	0 0 0 0
; 
run; 
* part (a);
proc phreg data=double;
	model time*status(0)=sex1 sex2 bili1 bili2 /entry=entrytime type3(lr);
	strata bleeding;
	baseline out=mort cumhaz=breslow covariates=covar;
run;
data mort; 
	set mort; 
	timeyears = time /365.25; 
run;
proc gplot data=mort; 
	plot breslow*timeyears=bleeding/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 4 by 1 minor=none 
	      label=('Time since randomization (Years)');
	axis2 order=0 to 4 by 1 minor=none label=(a=90 'Cumulative hazard');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
run;
quit;
