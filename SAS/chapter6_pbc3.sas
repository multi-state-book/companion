*------------------------------------------------------------------;
*------- Chapter 6, SAS code, PBC3 data ---------------------------;
*------------------------------------------------------------------;

* Load pbc data; 
proc import out=pbc3
	datafile="data/pbc3.csv"
	dbms=csv replace;
run;
* NB: Inside the macros below a variable called `id` is created which will interfere 
  with your identification variable if it is also called `id` in your data set,
  in which case you will need to rename your `id` variable. We have to do it for pbc3;
data pbc3; 
	set pbc3;
	fail=(status>0); 
	log2bili=log2(bili);
	years=days/365.25;
	rename id=ptno;
run;

*---------------------------------------------------------------;
*--------------------- Pseudo-macro ----------------------------;
*---------------------------------------------------------------;

* From: http://192.38.117.59/~linearpredictors/datafiles/pseudosurv.sas; 

* 'noprint plots = none' are included twice in the proc lifetest statement of the MACRO to avoid none important outputs being 
  printed;

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

proc lifetest data=newdat noprint plots = none;
time &TIME*&dead(0);
survival out=sall;
data sall;  set sall;
sall=survival;
 keep &time sall;
 
%do ip=1 %to &howmany;

/* COMPUTE KME FOR REDUCED SAMPLE */
proc lifetest data=newdat noprint plots = none;
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
	run;
%mend;



*---------------------------------------------------------------;
*--------------------- Figure 6.1 ------------------------------;
*---------------------------------------------------------------;

proc sort data=pbc3; 
	by days; 
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
quit;

*---------------------------------------------------------------;
*--------------------- Figure 6.2 ------------------------------;
*---------------------------------------------------------------;

* Load pbc data; 
proc import out=pbc3
	datafile="data/pbc3.csv"
	dbms=csv replace;
run;
data pbc3; 
	set pbc3;
	fail=(status>0); 
	log2bili=log2(bili);
	years=days/365.25;
	rename id=ptno;
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

data outdata; 
	set outdata;
	fail=(status>0); 
run;

proc gplot data=outdata;
	where time=366;
	plot pseudo*years = fail / haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=-0.9 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;
quit;

proc gplot data=outdata;
	where time=743;
	plot pseudo*years = fail / haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=-0.9 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;
quit;

proc gplot data=outdata;
	where time=1105;
	plot pseudo*years = fail / haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=-0.9 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Figure 6.3 ------------------------------;
*---------------------------------------------------------------;

proc sort data=outdata; 
	by bili;
run;

* Left plot; 
proc gplot data=outdata;
	where time=743;
	plot pseudo*bili/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 50 minor=none label=('Bilirubin');
	axis2 order=-0.4 to 1.2 by 0.2 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=x i=sm70; * i=sm70 specifies that a smooth line is fit to data;
run;
quit;

* Right plot; 
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
quit;


*---------------------------------------------------------------;
*--------------------- Figure 6.4 ------------------------------;
*---------------------------------------------------------------;


proc gplot data=outdata;
	where time=743;
	plot pseudo*log2bili/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(bilirubin)');
	axis2 order=-0.4 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=x i=sm70;* i=sm70 specifies that a smooth line is fit to data;
run;
quit;
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
quit;


*check; 
proc gplot data=smlogbili;
	plot smooth*log2bili/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(Bilirubin)');
	axis2 order=-5 to 2 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=none i=join;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Table 6.1 -------------------------------;
*---------------------------------------------------------------;

* At days 743 ----------------------------------------------------------;
proc genmod data=outdata;
	where time=743;
	class ptno;
	fwdlink link=log(-log(_mean_));
	invlink ilink=exp(-exp(_xbeta_));
	model pseudo=tment alb log2bili / dist=normal noscale;
	repeated subject=ptno/corr=ind;
run;

* At days 366 , 743, 1105 ---------------------------------------------;
proc genmod data=outdata;
	class ptno time;
	fwdlink link=log(-log(_mean_));
	invlink ilink=exp(-exp(_xbeta_));
	model pseudo=tment alb log2bili time / dist=normal noscale noint;
	repeated subject=ptno/corr=ind;
run;

*---------------------------------------------------------------;
*--------------------- Figure 6.5 ------------------------------;
*---------------------------------------------------------------;

data fig6_5; set outdata;
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

proc gplot data=fig6_5;
	plot res*log2bili=time/haxis=axis1 vaxis=axis2;
	axis1 order=1 to 9 by 1 minor=none label=('log2(Bilirubin)');
	axis2 order=-2 to 1 by 1 minor=none label=(a=90 'Pseudo-residuals');;
	symbol1 v=x i=sm50;
	symbol2 v=o i=sm50;
	symbol3 v=+ i=sm50;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Table 6.2 -------------------------------;
*---------------------------------------------------------------;

*  At days 366 , 743, 1105 using link function "-log"; 
proc genmod data=outdata;
	class ptno time;
	fwdlink link=-log(_mean_);
	invlink ilink=exp(-_xbeta_);
	model pseudo=tment alb bili time  /dist=normal noscale noint;
	repeated subject=ptno / corr=ind;
run;


*---------------------------------------------------------------;
*--------------------- Figure 6.6 ------------------------------;
*---------------------------------------------------------------;

* Right plot; 
proc sort data=outdata; 
	by bili; 
run;
proc gplot data=outdata;
	where time=743;
	plot pseudo*bili/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 50 minor=none label=('Bilirubin');
	axis2 order=-0.4 to 1.2 by 0.2 minor=none label=(a=90 'Pseudo-values');
	symbol1 v=x i=sm70;
run;
quit;

* Left plot; 
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
quit;

*---------------------------------------------------------------;
*--------------------- Figure 6.7 ------------------------------;
*---------------------------------------------------------------;

data fig6_7; 
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

proc gplot data=fig6_7;
	plot res2*bili=time/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 500 by 100 minor=none label=('Bilirubin');
	axis2 order=-2 to 1 by 1 minor=none label=(a=90 'Pseudo-residuals');;
	symbol1 v=x i=sm50;
	symbol2 v=o i=sm50;
	symbol3 v=+ i=sm50;
run;
quit;

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
 
 PROC LIFETEST DATA = work OUTSURV = km ;
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
 
   PROC LIFETEST DATA = work OUTSURV = km1 ;
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

%pseudomean(pbc3,years,fail,349,3,outmean3);

proc genmod data=outmean3;
	class ptno;
	model psumean = tment alb log2bili /dist=normal;
	repeated subject=ptno / corr=ind;
run;

proc rmstreg data=pbc3 tau=3;
   model years*status(0)=tment alb log2bili / link=linear;
run;

*---------------------------------------------------------------;
*--------------------- Figure 6.8 ------------------------------;
*---------------------------------------------------------------;

proc sort data=outmean3; 
	by years; 
run;
proc gplot data=outmean3;
	plot psumean*years=fail/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 4 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;
quit;

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
quit;

*---------------------------------------------------------------;
*--------------------- Table 6.4 -------------------------------;
*---------------------------------------------------------------;

* The pseudo observations for the cumulative incidences of stroke and death without stroke at three years are computed using the 
  pseudoci MACRO. To use this MACRO we must remember to include the cuminc MACRO.;

%macro cuminc(datain,x,re,de,dataout,cir,cid);
/*  THIS MACRO COMPUTES THE CUMULATIVE INCIDENCE FUNCTIONS FOR
    BOTH COMPETING RISKS USING PROC PHREG OUTPUT
    INPUTS TO MACRO
    DATAIN--NAME OF INPUT DATA SET CONTAINING 
    	X--TIME TO EVENT
    	RE--INDICATOR OF FIRST COMPETING RISK (1-YES, 0-NO)
    	DE--INDICATOR OF SECOND COMPETING RISK
    DATAOUT--NAME OF OUTPUT DATA SET CONTAINING
    	CIR--CUMULATIVE INCIDENCE FUNCTION FOR 1ST COMPETING RISK
    	CID--CUMULATIVE INCIDENCE FUNCTION FOR 2ST COMPETING RISK

*/

data work;  set &datain;
t=&x;
r=&re;
d=&de;
zero=0;

/* COMPUTE CRUDE CUMUALTIVE HAZARD FOR FIRST COMPETING RISK */
proc phreg data=work noprint; 
model t*r(0)=zero;
output out=rel  logsurv=chr  /method=emp;
 
 /* COMPUTE CRUDE CUMUALTIVE HAZARD FOR SECOND COMPETING RISK */
proc phreg data=work noprint; 
model t*d(0)=zero;
output out=dead  logsurv=chd  /method=emp;
 
 
 /* COMPUTE cumualtive incidence */
data both;  merge rel dead;  by t;
retain s 1
retain cr 0;
retain cd 0;
retain cumincr 0;
retain cumincd 0;
hr=-(cr+chr);
hd=-(cd+chd);

/* NOTE HR AND HD ARE THE JUMPS IN THE CUMUALTIVE CRUDE HAZARDS AT THIS TIME */

cr=-chr;
cd=-chd;
cir=cumincr+hr*s;
cumincr=cir;
cid=cumincd+hd*s;
cumincd=cid;
s=s*(1-hr-hd);
/* NOTE S IS KAPLAN-MEIER ESTIMATE IGNORING CAUSE OF FAILURE */
data &dataout;  set both;
&x=t;
&cir=cir;  &cid=cid;
keep &x &cir &cid;
run;
%mend;


%macro pseudoci(datain,x,r,d,howmany,datatau,dataout);

/*    MACRO COMPUTES PSEUDOVALUES BASED ON THE CUMUALTIVE INCIDENCE FUNCTION
      FOR BOTH OF TWO COMPETING RISKS  
      TIME
      INPUTS:
      DATAIN---INPUT DATA SET
      X--TIME VARIABLE
      R--INDICATOR OF FIRST COMPETING RISK (1-YES, 0-NO)
      D--INDICATOR OF SECOND COMPETING RISK
      HOWMANY---SAMPLE SIZE
     
      DATATAU---SUBSET OF INPUT DATA SET AT WHICH PSEUDO VALUES ARE COMPUTED 
                DATA SET HAS SINGLE VARIABLE "TIME"
                
      DATAOUT---OUTPUT DATA SET WHICH CONATINS PSUK,K=1,...,HOWMANY THE PSEUDO
                VALUES AT EACH TIME POINT (Note output data set
                 includes orginal data sorted by time)
      
*/

proc sort data=&datain;  by &x;

data keep;  set &datatau;
find=1;

proc sort data=keep;  by time;

data point;  set &datain;
time=&x;
keep=1;
data point;  merge point keep;  by time;
keep time find keep;
 
data useme;  set point;
retain temp -1;
if keep = 1 then temp=time;
tuse=temp; 
if find ne 1 then delete;
&x=tuse;
proc print;

/* PREPARE DATA SET WITH MISSING VALUES FOR DEADK AND RELAPSEK TO BE USED IN COMPUTING
   ESTIMATED CUMULATIVE INCIDENCE WITH KTH OBSERVATION DELETED*/
proc sort data=&datain;
by &x;
data newdat;  set &datain ;
id+1;
array iobsd(&howmany) dead1-dead&howmany;
array iobsr(&howmany) relapse1-relapse&howmany;
do j=1 to &howmany;
iobsd(j)=&d;
iobsr(j)=&r;
if j=id then do; iobsr(j)=.; iobsd(j)=.; end;
end;

data out;  set newdat;
drop dead1-dead&howmany relapse1-relapse&howmany;
/* COMPUTE CI FOR 1ST (CIRALL) AND 2ND (CIDALL) FOR FULL SAMPLE, STORE IN SALL*/
%cuminc(newdat,&x,&r,&d,sall,cirall,cidall);

%do ip=1 %to &howmany;

/* COMPUTE CI FOR 1ST (CIRALL) AND 2ND (CIDALL) FOR REDUCED SAMPLE, STORE IN SIP*/
%cuminc(newdat,&x,relapse&ip,dead&ip,stemp,cir1,cid1);

/* COMPUTE PSEUDOVALUES FOR BOTH RISK AT EVERY DATA POINT AND ADD TO FILE */ 
data ps; merge sall stemp;  by &x;
retain cirtemp 0;
retain cidtemp 0;
if cir1=. then cir1=cirtemp;
cirtemp=cir1;
rpsu&ip=&howmany*cirall- (&howmany-1)*cir1;
 if cid1=. then cid1=cidtemp;
cidtemp=cid1;
dpsu&ip=&howmany*cidall- (&howmany-1)*cid1;

data out; merge out ps useme; by &x;
if find ne 1 then delete;
keep time rpsu1-rpsu&ip dpsu1-dpsu&ip &x;
run;
%end;
 
 data &dataout;  set newdat; 
 drop dead1-dead&howmany relapse1-relapse&howmany;
 
 data all;  set out;
  
 array yr(&howmany) rpsu1-rpsu&howmany;
array yd(&howmany) dpsu1-dpsu&howmany;
do j=1 to &howmany;
rpseudo=yr(j);
dpseudo=yd(j);
id=j;
output;
end;
keep id time rpseudo dpseudo;
 proc sort data=all;  by id;
 data &dataout; merge &dataout all;
 by id;
 retain otime -1;
 retain oid -1;
 if id eq oid and otime=time then delete;
 else do; oid=id; otime=time; end;
 run;
 %mend;

* NB: Inside the macros below a variable called `id` is created which will interfere 
  with your identification variable is called `id` in your data set,
  in which case you will need to rename your `id` variable; 

* create indicator variables for each competing risk as required by the macro;
proc import out=pbc3
	datafile="data/pbc3.csv"
	dbms=csv replace;
run;
data pbc3; 
	set pbc3;
	log2bili=log2(bili);
	years=days/365.25;
	trans=status=1;
	death=status=2;
	rename id=ptno;
run;

data timepoint;
	input time;
	datalines;
	2
	;
run;

%pseudoci(pbc3,years,trans,death,349,timepoint,cumincpv1);

* rpseudo: transplant; 
* dpseudo: death;

* logit link function;
proc genmod data=cumincpv1;
	class ptno;
	model dpseudo = tment / dist=normal noscale link=logit; 
	repeated subject=ptno / corr=ind;
run;
proc genmod data=cumincpv1;
	class ptno;
	model dpseudo = tment alb log2bili / dist=normal noscale link=logit; 
	repeated subject=ptno / corr=ind;
run;

* cloglog link function;
proc genmod data=cumincpv1;
	class ptno;
	fwdlink link = log(-log(1-_mean_));
	invlink ilink = 1 - exp(-exp(_xbeta_));
	model dpseudo = tment / dist=normal noscale ; 
	repeated subject=ptno / corr=ind;
run;
proc genmod data=cumincpv1;
	class ptno;
	fwdlink link = log(-log(1-_mean_));
	invlink ilink = 1 - exp(-exp(_xbeta_));
	model dpseudo = tment alb log2bili / dist=normal noscale ; 
	repeated subject=ptno / corr=ind;
run;


* In-text, p. 233: cloglog link function - year 1,2,3;
data timepoints;
	input time;
	datalines;
	1
	2
	3
	;
run;
%pseudoci(pbc3,years,trans,death,349,timepoints,cumincpv3);

* cloglog link function - year 1,2,3;
proc genmod data=cumincpv3;
	class ptno time;
	fwdlink link = log(-log(1-_mean_));
	invlink ilink = 1 - exp(-exp(_xbeta_));
	model dpseudo = tment alb log2bili time / dist=normal noscale ; 
	repeated subject=ptno / corr=ind;
run;

*---------------------------------------------------------------;
*--------------------- Table 6.5 -------------------------------;
*---------------------------------------------------------------;

* cloglog link function - 3 time points;
proc genmod data=cumincpv3;
	class ptno time;
	fwdlink link = log(-log(1-_mean_));
	invlink ilink = 1 - exp(-exp(_xbeta_));
	model dpseudo = tment alb log2bili sex age time / dist=normal noscale ; 
	repeated subject=ptno / corr=ind;
run;

* cloglog link function - 10 time points;
data time10points;
	input time;
	datalines;
	.5
	1
	1.5
	2
	2.5
	3
	3.5
	4
	4.5
	5
	;
run;
%pseudoci(pbc3,years,trans,death,349,time10points,cumincpv10);

proc genmod data=cumincpv10;
	class ptno time;
	fwdlink link = log(-log(1-_mean_));
	invlink ilink = 1 - exp(-exp(_xbeta_));
	model dpseudo = tment alb log2bili sex age time / dist=normal noscale ; 
	repeated subject=ptno / corr=ind;
run;

libname h ".";
data h.cumincpv1;set cumincpv1;run;
data h.cumincpv3;set cumincpv3;run;
data h.cumincpv10;set cumincpv10;run;

*---------------------------------------------------------------;
*--------------------- Table 6.6 -------------------------------;
*---------------------------------------------------------------;

* No macro available;
