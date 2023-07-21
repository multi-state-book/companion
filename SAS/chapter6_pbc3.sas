*------------------------------------------------------------------;
*------- Chapter 6, SAS code, PBC3 data ---------------------------;
*------------------------------------------------------------------;

* Load pbc data; 
proc import out=pbc3
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/pbc3.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=pbc3; 
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
run;


*---------------------------------------------------------------;
*--------------------- Pseudo-macro ----------------------------;
*---------------------------------------------------------------;

* From: http://192.38.117.59/~linearpredictors/datafiles/pseudosurv.sas; 


%macro pseudosurv(indata,time,dead,howmany,datatau,outdata);

/* The subsequent SAS code is adapted from that described by J.P.Klein,
    M.Gerster, P.K.Andersen, S.Tarima, M.Pohar Perme (2008): "SAS and R
    functions to compute pseudo-values for censored data regression"
    Comp. Meth. Progr. Biomed., vol. 89,pp. 289-300. */ 

/*    MACRO COMPUTES PSEUDOVALUES BASED ON THE KAPLAN-MEIER ESTIMATOR AT EACH
      TIME
      INPUTS:
      INDATA---INPUT DATA SET
      TIME--TIME VARIABLE
      DEAD---EVENT INDICATOR (1-EVENT, 0-CENSORED)
      HOWMANY---SAMPLE SIZE
      DATATAU---SUBSET OF INPUT DATA SET AT WHICH PSEUDO VALUES ARE COMPUTED DATA
                SET HAS SINGLE VARIABLE TIME
     
      OUTDATA---OUTPUT DATA SET WHICH CONTAINS PSUK,K=1,...,HOWMANY THE PSEUDO
      VALUES AT EACH TIME POINT (Note output data set includes orginal data sorted
                           by time
      
*/
 
proc sort data=&indata;
by &time;

data keep;  set &datatau;
find=1;
proc sort data=keep;  by time;

data point;  set &indata;
time=&time;
keep=1;
 

data point;  merge point keep;  by time;
keep time find keep;
 
data useme;  set point;
retain temp -1;
if keep = 1 then temp=time;
tuse=temp; 
if find ne 1 then delete;
&time=tuse;
 run;
/* WORKING DATA SET THAT INCLUDE A SET OF N INDICATORS WHERE FOR THE KTH INDICATOR
 THE EVENT IS MISSING */
data newdat;  set &indata;
id+1;
array iobs(&howmany) dead1-dead&howmany;
do j=1 to &howmany;
iobs(j)=&dead;
if j=id then iobs(j)=.;
end;

data out;  set newdat;

/* COMPUTE KME FOR FULL SAMPLE */

proc lifetest data=newdat noprint;
time &TIME*&dead(0);
survival out=sall;
data sall;  set sall;
sall=survival;
 keep &time sall;
 
%do ip=1 %to &howmany;

/* COMPUTE KME FOR REDUCED SAMPLE */
proc lifetest data=newdat noprint;
time &time*dead&ip(0);
survival out=stmp;
data stmp;  set stmp;
s&ip=survival;
keep &time s&ip;
 
/*merge KMEs AND COMPUTE PSEUDOVALUES FOR OBSERVATION IP*/
data pstmp; merge sall stmp;  by &time;
retain stemp 1;
if s&ip=. then s&ip=stemp;
stemp=s&ip;
psu&ip=&howmany*sall- (&howmany-1)*s&ip;
 
data out; merge out pstmp useme ;  by &time;
if find ne 1 then delete;
keep &time psu1-psu&ip;
%end;

data out;  set out;
retain dup -1;
if &time=dup then delete;
 else dup=&time;
 jd+1;
 
data io;  set out;
array ps psu1-psu&howmany;
do id=1 to &howmany;
pseudo=ps(id);
time=&time;
output;
end;
keep time id pseudo jd;
proc sort data=io ;
  by id;
  
data a;  set &indata;
id+1;
data io;  merge a io;  by id;
 proc sort data=io;  by jd;

proc sort data=&datatau;
by time;

data taus;  set &datatau ;
 jd+1;
 tpseudo=time;
 keep jd tpseudo;
 
 
  data &outdata;  merge io taus;  by jd;
  drop jd id &time &dead;
%mend;



*---------------------------------------------------------------;
*--------------------- Figure 6.1 ------------------------------;
*---------------------------------------------------------------;

proc sort data=pbc3; 
	by days; 
run;

data pbc3; 
	set pbc3; 
	fail=(status>0); 
run; 

proc print data=pbc3; 
	var days fail; 
run;

data timepoints;
	input time;
	datalines;
	366
	743
	1105
	;
run;

%pseudosurv(pbc3, days, fail, 349, timepoints, outdata);

* Pseudo-values for S(t) at ~1 2 3 years are computed and stored in permanent data set;
proc sort data=outdata; 
	by tpseudo; 
run;

data outdata; 
	set outdata; 
	fail=(status>0); 
run;

libname h "/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/";

data h.pbcpseu123; 
	set outdata; 
run;

data outdata;
  set h.pbcpseu123;
run;


data minus1; 
	set pbc3; /* drop failure with t~1 yr */
	if ptno=415 then delete;
run;

proc lifetest 
	data=pbc3;
	time days*fail(0);
	survival out=sall;
run;

data sall; 
	set sall; 
	survall=survival; 
run;

proc lifetest data=minus1;
	time days*fail(0);
	survival out=sminus1;
run;

data ekstra;
	input days survival;
	datalines;
	366 0.9225440684
	;
run;

data sminus1ny; 
	set sminus1 ekstra; 
	survny=survival;
run; 

proc sort; 
	by days;
run;

data final1; 
	merge sall sminus1ny; 
	by days; 
	followup=days/365;
	pseudo1=349*survall-348*survny;
run;

proc gplot data=final1;
	plot pseudo1*followup;
run;

data minus2; 
	set pbc3; /* drop cens with t~1 yr */
	if ptno=458 then delete;
run;

data ekstra;
	input days survival;
	datalines;
	365 0.9225440684
	;
run;

proc lifetest data=minus2;
	time days*fail(0);
	survival out=sminus2;
run;


data sminus2ny; 
	set sminus2 ekstra; 
	survny=survival;
run; 

proc sort; 
	by days;
run;

data final2; 
	merge sall sminus2ny; 
	by days; 
	followup=days/365;
	pseudo2=349*survall-348*survny;
run;

proc gplot data=final2;
	plot pseudo2*followup;
run;

data final; 
	merge final1 final2; 
	by days; 
	en=1; nul=0;
run;

proc gplot data=final;
	plot pseudo1*followup pseudo2*followup en*followup nul*followup
	      /overlay haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=-0.1 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=none i=join c=blue;
	symbol2  v=none i=join c=red;
	symbol3 v=none i=join c=black;
	symbol4 v=none i=join c=black;
run;

*---------------------------------------------------------------;
*--------------------- Figure 6.2 ------------------------------;
*---------------------------------------------------------------;

proc gplot data=outdata;
	where time=366;
	plot pseudo*followup=fail/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=-0.9 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;


proc gplot data=outdata;
	where time=743;
	plot pseudo*followup=fail/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=-0.9 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;

proc gplot data=outdata;
	where time=1105;
	plot pseudo*followup=fail/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=-0.9 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;


*---------------------------------------------------------------;
*--------------------- Figure 6.3 ------------------------------;
*---------------------------------------------------------------;

data outdata; 
	set outdata;  
	log2bili=log2(bili); 
run;

proc sort data=outdata; 
	by bili;
run;

* Top plot; 
proc gplot data=outdata;
	where time=743;
	plot pseudo*bili/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 50 minor=none label=('Bilirubin');
	axis2 order=-0.4 to 1.2 by 0.2 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=x i=sm70;
run;

* Bottom plot; 
proc loess data=outdata;
	where time=743;
	model pseudo=bili/smooth=0.7;
	output out=smbili p=smooth;
run;

data smbili; 
	set smbili;
	line=log(-log(smooth));
run;

proc sort data=smbili; 
	by bili; 
run;

proc gplot data=smbili;
	plot line*bili/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 50 minor=none label=('Bilirubin');
	axis2 order=-5 to 2 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=none i=join;
run;


*---------------------------------------------------------------;
*--------------------- Figure 6.4 ------------------------------;
*---------------------------------------------------------------;


proc gplot data=outdata;
	where time=743;
	plot pseudo*log2bili/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(bilirubin)');
	axis2 order=-0.4 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=x i=sm70;
run;

proc loess data=outdata;
	where time=743;
	model pseudo=log2bili/smooth=0.7;
	output out=smlogbili p=smooth;
run;

data smlogbili; 
	set smlogbili;
	line=log(-log(smooth));
run;

proc sort data=smlogbili; 
	by bili; 
run; 

proc gplot data=smlogbili;
	plot line*log2bili/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(Bilirubin)');
	axis2 order=-5 to 2 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=none i=join;
run;


*check; 
proc gplot data=smlogbili;
	plot smooth*log2bili/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(Bilirubin)');
	axis2 order=-5 to 2 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=none i=join;
run;


*---------------------------------------------------------------;
*--------------------- Table 6.1 -------------------------------;
*---------------------------------------------------------------;

/*
proc genmod data=outdata;
	where time=743;
	class ptno;
	fwdlink link=log(-log(_mean_));
	invlink ilink=exp(-exp(_xbeta_));
	model pseudo=log2bili/dist=normal noscale;
	repeated subject=ptno/corr=ind;
run;
*/

proc genmod data=outdata;
	where time=743;
	class ptno;
	fwdlink link=log(-log(_mean_));
	invlink ilink=exp(-exp(_xbeta_));
	model pseudo=log2bili alb tment/dist=normal noscale;
	repeated subject=ptno/corr=ind;
run;

proc genmod data=outdata;
	class ptno time;
	fwdlink link=log(-log(_mean_));
	invlink ilink=exp(-exp(_xbeta_));
	model pseudo=time log2bili alb tment/dist=normal noscale noint;
	repeated subject=ptno/corr=ind;
run;

*---------------------------------------------------------------;
*--------------------- Figure 6.5 ------------------------------;
*---------------------------------------------------------------;

data outdata; set outdata;
	if time=366 then do;
	linpred=-2.4746+0.6841*log2bili-0.0939*alb-0.5985*tment;
	pred=exp(-exp(linpred)); res=pseudo-pred; end;
	if time=743 then do;
	linpred=-1.5540+0.6841*log2bili-0.0939*alb-0.5985*tment;
	pred=exp(-exp(linpred)); res=pseudo-pred; end;
	if time=1105 then do;
	linpred=-1.1232+0.6841*log2bili-0.0939*alb-0.5985*tment;
	pred=exp(-exp(linpred)); res=pseudo-pred; end;
run;

proc gplot data=outdata;
	plot res*log2bili=time/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(Bilirubin)');
	axis2 order=-2 to 1 by 1 minor=none label=(a=90 'Pseudo-residuals');;
	symbol1 v=x i=sm50;
	symbol2 v=o i=sm50;
	symbol3 v=+ i=sm50;
run;

*---------------------------------------------------------------;
*--------------------- Table 6.2 -------------------------------;
*---------------------------------------------------------------;

proc genmod data=outdata;
	class ptno time;
	*fwdlink link=log(_mean_);
	*invlink ilink=exp(_xbeta_);
	fwdlink link=-log(_mean_);
	invlink ilink=exp(-_xbeta_);
	model pseudo=time bili alb tment/dist=normal noscale noint;
	repeated subject=ptno/corr=ind;
run;


*---------------------------------------------------------------;
*--------------------- Figure 6.6 ------------------------------;
*---------------------------------------------------------------;

proc loess data=outdata;
	where time=743;
	model pseudo=bili/smooth=0.7;
	output out=smbili p=smooth;
run;

data smbili; 
	set smbili;
	line=-log(smooth);
run;

proc sort data=smbili; 
	by bili; 
run;

proc gplot data=smbili;
	plot line*bili/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 50 minor=none label=('Bilirubin');
	axis2 order=-1 to 3 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=none i=join;
run;

proc sort data=outdata; 
	by bili; 
run;

proc gplot data=outdata;
	where time=743;
	plot pseudo*bili/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 50 minor=none label=('Bilirubin');
	axis2 order=-0.4 to 1.2 by 0.2 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=x i=;
run;


*---------------------------------------------------------------;
*--------------------- Figure 6.7 ------------------------------;
*---------------------------------------------------------------;

data outdata; 
set outdata;
	if time=366 then do;
	linpred2=0.3403+0.0042*bili-0.0097*alb-0.0484*tment;
	pred2=exp(-linpred2); res2=pseudo-pred2; end;
	if time=743 then do;
	linpred2=0.4120+0.0042*bili-0.0097*alb-0.0484*tment;
	pred2=exp(-linpred2); res2=pseudo-pred2; end;
	if time=1105 then do;
	linpred2=0.5075+0.0042*bili-0.0097*alb-0.0484*tment;
	pred2=exp(-linpred2); res2=pseudo-pred2; end;
run;

proc gplot data=outdata;
	plot res2*bili=time/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 100 minor=none label=('Bilirubin');
	axis2 order=-2 to 1 by 1 minor=none label=(a=90 'Pseudo-residuals');;
	symbol1 v=x i=sm50;
	symbol2 v=o i=sm50;
	symbol3 v=+ i=sm50;
run;


*---------------------------------------------------------------;
*--------------------- Table 6.3 -------------------------------;
*---------------------------------------------------------------;


* Macro for computing pseudo-obs of the restricted mean; 
%macro pseudomean(indata,time,dead,howmany,tmax,outdata);
 /* MACRO ARGUMENTS
      INDATA--NAME OF INPUT DATA SET
      TIME--NAME OF TIME VARIABLE
      DEAD--NAME OF EVENT INDICATOR VARIABLE--(1-DEAD,0-CENSORED)
      HOWMANY--SAMPLE SIZE
      TMAX--UPPER LIMIT OF INTEGRATION FOR RESTRICTED MEAN
      OUTDATA--NAME OF OUTPUT DATA SET WITH PSEUDOVALUES FPR RESTRICTED
               MEAN IN VARIABLE "PSUMEAN"  */
 
 /* CREATE A DATA SET WHERE EVERYTHING ABOVE TMAX IS CENSORED */
 DATA work; SET &indata;
  restime = MIN(&tmax, &time);
  resdead = &dead;
  IF restime EQ &tmax THEN resdead = 0;

/* CREATE DATA SET WITH SET OF INDICATORS DEADK THAT HAS MISSING VALUE
     FOR KTH OBSERVATION, K=1,...,HOWMANY*/
 DATA work;  SET work;
   id+1;
   ARRAY iobs(&howmany) dead1-dead&howmany;
   DO j = 1 TO &howmany;
    iobs(j) = resdead;
    IF j = id THEN iobs(j) = .;
   END;
 
 /* COMPUTE RESTRICTED MEAN FOR COMPLETE SAMPLE USING PROC LIFETEST */
 
 PROC LIFETEST DATA = work OUTSURV = km;
   TIME restime*resdead(0);
   ODS SELECT MEANS;
   ODS OUTPUT MEANS = mall;
 RUN;

  DATA km; SET km;
    IF _CENSOR_ EQ 0;
  PROC SORT DATA=km;
    BY restime;
  RUN; 
  DATA km; SET km END=LAST;
    IF NOT(LAST) THEN DELETE;
    area = (&tmax - restime)*survival;
    KEEP area;

  DATA psu; MERGE km mall;
    meanall = mean + area;
    KEEP meanall;
    
 %DO ip = 1 %TO &howmany;
 
   /* COMPUTE RESTRICTED MEAN FOR SAMPLE WITH IPTH OBSERVATION DELETED
      USING PROC LIFETEST */
 
   PROC LIFETEST DATA = work OUTSURV = km1;
     TIME restime*dead&ip(0);
     ODS SELECT means;
     ODS OUTPUT MEANS = m1;
   RUN;

   DATA km1; SET km1;
     IF _CENSOR_ EQ 0;
   PROC SORT DATA = km1;
     BY restime;
   RUN;
   DATA km1; SET km1 END=LAST;
     IF NOT(LAST) THEN DELETE;
     area = (&tmax - restime)*survival;
     KEEP area;
   DATA km1; MERGE km1 m1;
     mean = mean + area;
     KEEP mean;

   /* COMPUTE PSEUDOVALUE FOR IPTH OBSERVATION*/
   DATA psu; MERGE psu km1;
     psu&ip=&howmany*meanall-(&howmany-1)*mean;
 %END;
 
 /* TRANSPOSE DATASET AND MERGE WITH RAW DATA*/
 DATA out;  SET psu;
   ARRAY y(&howmany) psu1-psu&howmany;
   DO j = 1 TO &howmany;
     psumean=y(j);
     OUTPUT;
   END;
 KEEP psumean;
 
 DATA &outdata; MERGE &indata out;
 run; 
%MEND;

%pseudomean(pbc3,followup,fail,349,3,outmean3);


data outmean3;
set outmean3;
	log2bili=log2(bili);
run;


*GEE model fit; 
proc genmod data=outmean3;
	class ptno;
	model psumean=log2bili alb tment/dist=normal;
	repeated subject=ptno/corr=ind;
run;


*---------------------------------------------------------------;
*--------------------- Figure 6.8 ------------------------------;
*---------------------------------------------------------------;

proc sort data=outmean3; 
	by followup; 
run;

proc gplot data=outmean3;
	plot psumean*followup=fail/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 4 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;


*---------------------------------------------------------------;
*--------------------- Figure 6.9 ------------------------------;
*---------------------------------------------------------------;

proc sort data=outmean3; 
	by bili; 
run;

proc gplot data=outmean3;
	plot psumean*log2bili/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(bilirubin)');
	axis2 order=0 to 4 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=x i=sm70;
run;


*---------------------------------------------------------------;
*--------------------- Table 6.4 -------------------------------;
*---------------------------------------------------------------;