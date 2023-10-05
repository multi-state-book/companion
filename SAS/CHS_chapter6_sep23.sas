*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- CHAPTER 6 - EXERCISE SOLUTIONS ---------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We first load the data;

proc import out = chs_data
    datafile = 'C:\HRfinal\holter\cphholter.csv'
	dbms= csv replace;
	getnames=yes;
run;

* We will convert the time variables, timeafib, timestroke, and timedeath, from days to years;

data chs_data;
	set chs_data;
	timeafib = timeafib/365.25;
	timestroke = timestroke/365.25;
	timedeath = timedeath/365.25;
run;

* Furthermore, we add variables for the composite end-point of stroke or death without stroke;

data chs_data;
	set chs_data;
	timestrokeordeath = timedeath;
	if stroke = 1 then timestrokeordeath = timestroke;
	strokeordeath = death;
	if stroke = 1 then strokeordeath = 1;
run;


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 6.3 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  6.3.1 ----------------------------------------------------------------------------;

* We will compute the pseudo observations for the incompletely observed survival indicators I(T_i > 3 years), i = 1,...,678 using 
  the pseudosurv MACRO.;

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
%mend;

* Inside the MACRO a variable called 'id' is created which interferes with the 'id' variable in our data set, chs_data. Therefore,
  we will first rename the 'id' variable in the chs_data;

data chs_data; 
	set chs_data; 
	ptno = id;
	drop id;
run; 

* Then, we make a data set named 'timepoint' which contains one column 'time' with the timepoint where we wish the pseudo observations to
  be computed.;

data timepoint;
	input time;
	datalines;
	3
	;
run;

* The pseudosurv MACRO is then called with 'chs_data', 'timestrokeordeath', 'strokeordeath', '678'. 'timepoint' and 'outdatadata631' 
  as arguments;

* The MACRO returns a data set 'outdata631' where the pseudo observations are stored in the column 'pseudo'.;

%pseudosurv(chs_data, timestrokeordeath, strokeordeath, 678, timepoint, outdata631);

* The output of the MACRO 'outdata631' is then used as data when fitting the cloglog model using the genmod procedure.;

* We must specify the link function and the inverse link function in the fwdlink and invlink statements;

title '6.3.1';
proc genmod data=outdata631;
	class ptno;
	fwdlink link=log(-log(_mean_));
	invlink ilink=exp(-exp(_xbeta_));
	model pseudo= esvea sex age sbp/dist=normal noscale;
	repeated subject=ptno/corr=ind;
run;

* Thus, we obtain the following model for stroke-free survival at three years including the covariates ESVEA, sex, age, and systolic 
  blood pressure
  S(3|Z) = exp(-exp(-8.4177 + 0.2358*Z1 + 0.5525*Z2 +  0.0810*Z3 - 0.0019*Z4)),
  where Z1, Z2, Z3, and Z4 are ESVEA, sex, age, and systolic blood pressure, respectively.;


*----------------------------------------------  6.3.2 ----------------------------------------------------------------------------;

* In Exercise 2.4 we obtained the following estimates of the coefficients for ESVEA, sex, age, and systolic blood pressure assuming 
  the Cox model and almost identical estimates was found assuming a Poisson model;
   
* beta_esvea = 0.31830, beta_sex = 0.57759, beta_age = 0.07666, beta_sbp = 0.00515;

* Thus, the models agree in the direction and in some degree the magnitude of the effect of ESVEA, sex, and age. The effect of 
  systolic blood pressure is on the other hand estimated to have the opposite effect in the model using pseudo observations 
  compared with the Cox and Poisson models. However, the effect of systolic blood pressure (and also ESVEA and sex) is not 
  significant at a 0.05 level in the model based on pseudo values.;


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 6.4 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  6.4.1 ----------------------------------------------------------------------------;

* The pseudo observations for the 3 year restricted mean time to the composite end-point are computed using the pseudomean MACRO;

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

* The pseudomean MACRO is now called with 'chs_data', 'timestrokeordeath', 'strokeordeath', '678'. '3' and 'outdatadata641' 
  as arguments;

* The MACRO returns a data set 'outdata641' where the pseudo observations are stored in the column 'psumean'.;

title '6.4.1';
%pseudomean(chs_data, timestrokeordeath, strokeordeath, 678, 3, outdata641);

* The model is then fitted with the genmod procedure; 

proc genmod data=outdata641;
	class ptno;
	model psumean=esvea sex age sbp /dist=normal;
	repeated subject=ptno/corr=ind;
run;

* We get the following model for the 3-year restricted mean time to the composite end-point stroke-free survival
  epsilon(3|Z) = 3.3552 - 0.0274*Z1 - 0.0540*Z2 - 0.0073*Z3 + 0.0005*Z4,
  where Z1, Z2, Z3, and Z4 are ESVEA, sex, age, and systolic blood pressure, respectively.;


*----------------------------------------------  6.4.2 ----------------------------------------------------------------------------;

* The model obtained in Exercise 4.3 is identical to the model fitted using pseudo observations. This is because none of the subjects 
  are censored before 3 years and the pseudo observation for subject i is thus min(3,T_i) where T_i is the time of stroke or death 
  without stroke for subject i.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 6.5 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  6.5.1 ----------------------------------------------------------------------------;

* The pseudo observations for the cumulative incidences of stroke and death without stroke at three years are computed using the 
  pseudoci MACRO. To use this MACRO we must remember to include the cuminc MACRO.;

%macro cuminc(datain,x,re,de,dataout,cir,cid);
/*  THIS MACRO COMPUTES THE CUMULATIVE INCIDENCE FUNCTIONS FOR
    BOTH COMPETING RISKS USING PROC PHREG OUTPUT
    INPUTS TO MACRO
    DATAIN--NAME OF INPUT DATA SET CONTAINING 
    	X--TIME TO EVENT
    	RE--INDICCATOR OF FIRST COMPETING RISK (1-YES, 0-NO)
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
%mend;


%macro pseudoci(datain,x,r,d,howmany,datatau,dataout);

/*    MACRO COMPUTES PSEUDOVALUES BASED ON THE CUMUALTIVE INCIDENCE FUNCTION
      FOR BOTH OF TWO COMPETING RISKS  
      TIME
      INPUTS:
      DATAIN---INPUT DATA SET
      X--TIME VARIABLE
      R--INDICCATOR OF FIRST COMPETING RISK (1-YES, 0-NO)
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
 %mend;


* We must create a variable for the end-point death without stroke;

data chs_data;
	set chs_data;
	death_wo_stroke = death;
	if stroke = 1 then death_wo_stroke = 0;
run;

* The pseudoci MACRO is called with 'chs_data', 'timestrokeordeath', 'stroke', 'death_wo_stroke', '678'. 'timepoint' and 'outdatadata661'
  as arguments;

* The MACRO returns a data set 'outdata661' where the pseudo observations for the end-point stroke are stored in the column 'rpseudo'
  and the pseudo observations for death without stroke are stored in the column 'dpseudo'.;

title '6.6.1';
 %pseudoci(chs_data, timestrokeordeath, stroke, death_wo_stroke, 678, timepoint, outdata661);

* We fit a cloglog model for the cumulative incidence of stroke with the genmod procedure where we specify the 'fwdlink' and 'invlink' 
  statements;

title '6.5.1 - stroke';
proc genmod data=outdata661;
	class ptno;
	fwdlink link=log(-log(1-_mean_));
	invlink ilink= 1 - exp(-exp(_xbeta_));
	model rpseudo= esvea sex age sbp /dist=normal  noscale ; 
	repeated subject=ptno/corr=ind;
run;

* Thus, we obtain the following model for the cumulative incidence of stroke at 3 years
  F_stroke(3|Z) = 1 - exp(-exp(-13.2294 + 0.0312*Z1 + 0.9711*Z2 + 0.1199*Z3  + 0.0040*Z4)),  
  where Z1, Z2, Z3, and Z4 are ESVEA, sex, age, and systolic blood pressure, respectively.;

* Likewise, we will fit a model for the cumulative incidence of death without stroke at 3 years;

title '6.5.1 - death without stroke';
proc genmod data=outdata661;
	class ptno;
	fwdlink link=log(-log(1-_mean_));
	invlink ilink= 1 - exp(-exp(_xbeta_));
	model dpseudo= esvea sex age sbp /dist=normal noscale ; 
	repeated subject=ptno/corr=ind;
run;

* We obtain the following model for the cumulative incidence of death without stroke at 3 years
  F_death(3|Z) = 1 - exp(-exp(-3.3924 + 0.2862*Z1 + 1.2618*Z2 + 0.0771*Z3 - 0.0442*Z4)), 
  where Z1, Z2, Z3, and Z4 are ESVEA, sex, age, and systolic blood pressure, respectively.;

*----------------------------------------------  6.6.2 ----------------------------------------------------------------------------;

* In Exercise 4.5 we used Fine-Gray models to estimate the cumulative incidence functions for stroke and death without stroke. 
  For the end-point stroke we obtained the following coefficients
  beta_esvea = 0.59398, beta_sex = 0.37919, beta_age = 0.06335, beta_sbp = 0.01063;

* and for death without stroke we obtained the following coefficients
  beta_esvea = -0.00596, beta_sex = 0.53042, beta_age = 0.06651, beta_sbp = 0.00160; 

* The estimates of the coefficients are quite different for the Fine-Gray model and the model based on pseudo 
observations for 
  both end-points. This is likely due to the fact that the Fine-Gray model assumes that the linear 
predictor has the same association with the cumulative incidence at all time points. We also notice that only one of the covariates in each model using pseudo values (age in the case of stroke 
  and systolic blood pressure in the case of death without stroke) are associated with p-values smaller than 0.05 which might explain 
  the poor agreement.;

* Comparison with Exercise 5.8 is not meaningful - sorry!;
