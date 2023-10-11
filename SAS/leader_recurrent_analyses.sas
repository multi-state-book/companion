* Load data; 
proc import out=leader
	datafile="c:/Users/hnrv/OneDrive - Novo Nordisk/Book/leader/data/leader_mi_3p.csv"
	dbms=csv replace;
run;
data leader_mi; 
	set leader; 
	where type = "recurrent_mi"; 
run; 


* Frailty model - slow!; 
proc phreg data=leader_mi covs(aggregate);
  class id;
  model stop*status(0 2) = treat / entry=start; 
  random id / dist=gamma;
  title1 'Frailty model for recurrent event data';
run;

/* 
------------------------------------------------------------------------------
------------------------------------------------------------------------------
 LEADER ARTICLE: Recurrent events
------------------------------------------------------------------------------
------------------------------------------------------------------------------
This program does all the analyses I can think of :-) 
Data sets created in leader_play_recurrent.sas

*/

*** Set-up access; 

%access_v01;

* LEADER CTR adam data - check with SNOL & SRRM which freeze to use; 
libname adam4738 "/projstat/ex2211/ex2211-3748/ctr_20160811_er/stats/data/ADaM";

libname myout "/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf";




* MI; 
data recurrent_mi; 
	set myout.recurrent_leader_all(where=(id eq "recurrent_mi")); 
run; 

data recurrent_wlw_mi; 
	set myout.recurrent_leader_wlw_all(where=(id eq "recurrent_mi")); 
run; 


** All models; 
* Need to use the convergelike=1E-9 to obtain same optimization as in R using N-R algorithm; 

* Time-to-first event Cox model (adjusted for treat);
proc phreg data=recurrent_mi(where=(eventno=1));
  class trt01p;
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow CONVERGELIKE=1E-9; 
  title1 'Cox model for time to first event';
run;

* Andersen-Gill model (AG); 
proc phreg data=recurrent_mi;
  class trt01p;
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow  CONVERGELIKE=1E-9; 
  title1 'A-G model for recurrent event data (AG 1982)';
run;

* PWP 1 model (one overall treatment estimate (a common treatment effect)); 
proc phreg data=recurrent_mi covs(aggregate)/*(where=(eventno<=5))*/;
  class trt01p;
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow CONVERGELIKE=1E-9; 
  strata eventno; 
  id subjid; 
  title1 'PWP 1 model for recurrent event data (PWP 1989)';
run;

* PWP 2 model (one per transition);
proc sort data=recurrent_mi; 
	by eventno; 
run; 
 
proc phreg data=recurrent_mi(where=(eventno <= 5));
  class trt01p;
  by eventno; 
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow CONVERGELIKE=1E-9; 
  title1 'PWP 2 model for recurrent event data (PWP 1989)';
run;

* WLW model; 
proc sort data=recurrent_wlw_mi; 
	by wlw_eventno; 
run; 

proc phreg data=recurrent_wlw_mi covs(aggregate);
  class trt01p;
  strata wlw_eventno;
  id subjid;  
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow CONVERGELIKE=1E-9; 
  title1 'WLW 1 model for recurrent event data';
run;

proc phreg data=recurrent_wlw_mi;
  class trt01p;
  by wlw_eventno; 
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow CONVERGELIKE=1E-9; 
  title1 'WLW 2 model for recurrent event data';
run;


* Marginal means and rates model (LWYY); 
proc phreg data=recurrent_mi covs(aggregate);
  class trt01p;
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow CONVERGELIKE=1E-9; 
  id subjid;
  title1 'Marginal means model for recurrent event data (LWYY 2000)';
  *baseline out=lwyy_mi cumhaz=cumhaz; 
run;


* Same as - treat death as a censoring; 
/*
proc phreg data=recurrent_mi covs(aggregate);
  class trt01p;
  model (starttime, stoptime)*status(0 2) = trt01p / ties=breslow; 
  id subjid;
  title1 'Marginal means model for recurrent time-to-event data';
  baseline out=lwyy_mi cumhaz=cumhaz; 
run;
*/

* Ghosh & Lin (2002); 
proc phreg data=recurrent_mi covs(aggregate);
  class trt01p;
  model (starttime, stoptime)*status(0) = trt01p / ties=breslow eventcode=1 CONVERGELIKE=1E-9; 
  id subjid;
  title1 'Marginal means model (adjusted for death) for recurrent event data (Ghosh-Lin 2002)';
run;






* Non-parametric two-sample estimates (Lawless & Nadeau); 
proc phreg data=recurrent_mi covs(aggregate);
  class trt01p;
  model (starttime, stoptime)*status2(0) = / ties=breslow CONVERGELIKE=1E-9; 
  strata trt01p; 
  id subjid;
  title1 'Marginal means model (non-parametric) for recurrent event data (Lawless & Nadeau)';
  baseline out=law_mi cumhaz=cumhaz; 
run;

* And the test - score test result; 
proc phreg data=recurrent_mi covs(aggregate);
  class trt01p;
  model (starttime, stoptime)*status2(0) = trt01p / ties=breslow CONVERGELIKE=1E-9; 
  id subjid;
  title1 'Marginal means model (non-parametric) for recurrent event data (Lawless & Nadeau)';
run;
/*
* Non-parametric two-sample estimates, adjusted for death (Cook & Lawless 1997, Ghosh & Lin 2000); 
proc phreg data=recurrent_mi covs(aggregate);
  class trt01p;
  model (starttime, stoptime)*status(0) = / ties=breslow eventcode=1 CONVERGELIKE=1E-9; 
  strata trt01p; 
  id subjid;
  title1 'Marginal means model (non-parametric), adjusted for death, for recurrent event data (Ghosh & Lin)';
  baseline out=gl_np_mi cumhaz=cumhaz; 
run;

* Two-sample test; 
proc phreg data=recurrent_mi covs(aggregate);
  class trt01p;
  model (starttime, stoptime)*status(0) = trt01p / ties=breslow eventcode=1 CONVERGELIKE=1E-9; 
  id subjid;
  title1 'Marginal means model (non-parametric), adjusted for death, for recurrent event data (Ghosh & Lin)';
  baseline out=gl_np_mi cumhaz=cumhaz; 
run;*/


* LWA; 
proc phreg data=recurrent_mi;
  class trt01p;
  strata eventno;
  model stoptime*status2(0) = trt01p / ties=breslow  CONVERGELIKE=1E-9; 
  title1 'LWA model for recurrent event data';
run;


**************** CHECK LEADER ARTICLE PLOT ***************************; 
* MCF plot for 3-p MACE; 
data recurrent_3pmace; 
	set myout.recurrent_leader_all(where=(id eq "recurrent_comb_str_mi_cvdth")); 
run; 

* NELSON AALEN; 
proc phreg data=recurrent_3pmace;
  strata trt01p;
  model (starttime, stoptime)*status(0 2) = / ties=breslow CONVERGELIKE=1E-9; 
  title1 'Nelson-Aalen estimates (censoring for death)';
  baseline out=na_data cumhaz=naa;
run;

data na_est;
	set na_data; 
	type = "Nelson-Aalen";
	cumevent = naa; 
	treat_type = trim(trt01p) || ", "  || type; 
run; 

* COOK & LAWLESS (GHOSH & LIN);
proc phreg data=recurrent_3pmace;
  strata trt01p;
  model (starttime, stoptime)*status(0) = /  eventcode=1 ties=breslow CONVERGELIKE=1E-9; 
  title1 'Ghosh Lin (2000) estimates (adjusting for non-CV death)';
  baseline out=gl_data cif=cuminc;
run;

data gl_est;
	set gl_data; 
	type = "Ghosh & Lin";
	cumevent = -log(1-cuminc); 
	treat_type = trim(trt01p) || ", " || type; 
run; 

data comb; 
	set na_est gl_est; 
	time = stoptime;
	drop naa cuminc;
run;
	
* Plot;
proc sgplot data=comb;
	step x=time y=cumevent / group=treat_type; 
run; 

/*
* OK; * What do I think they have done;
proc phreg data=recurrent_3pmace;
  strata trt01p;
  model (starttime, stoptime)*status(0 2) = /  ties=exact CONVERGELIKE=1E-9; 
  title1 'Puscv009 approach 1';
  baseline out=pudata1 cumhaz=cumhaz;
run;

data pudata11; 
	set pudata1; 
	cumevent = cumhaz; 
	treat_type = trim(trt01p) || ", " || "Puscv009 1";
run; 

proc phreg data=recurrent_3pmace;
  strata trt01p;
  model (starttime, stoptime)*status(0) = /  eventcode=1 CONVERGELIKE=1E-9; 
  title1 'Puscv009 approach 2';
  baseline out=pudata2 cif=cuminc;
run;

data pudata22; 
	set pudata2; 
	cumevent = -log(1-cuminc); 
	treat_type = trim(trt01p) || ", " || "Puscv009 2";
run; 

* CHECK; 
data comb; 
	set pudata11 pudata22; 
	time = stoptime;
	drop cumhaz cuminc;
run;
	

* Plot;
proc sgplot data=comb;
	step x=time y=cumevent / group=treat_type; 
run; 
*/


* Check data set; 

* Times to non-cv death, from ADTTE in LEADER; 
data noncvdeath; 
	set adam4738.adtte; 
	where paramcd='NONCVTM'; 
	keep ady cnsr subjid;
run;

data check; 
    merge recurrent_3pmace noncvdeath;
    by subjid;
    if status2=1 then event=1;
    else if status2=0 then do;
        if cnsr=1 then event=0;
        else if cnsr=0 and ady<=stoptime+1 then event=2;
		else if cnsr=0 and ady=stoptime then event=2;
        else if cnsr=0 and ady>stoptime then event=0;
    end;
run; 


data check2; 
	set check; 
	where status ne event; 
run; 

data in2;
    trt01p='Lira   '; trt=2; output;
    trt01p='Placebo'; trt=1; output;
run;

* Perform Andersen-Gill recurrent event analysis with graph output;
* covs only on SEs; 
proc phreg data=recurrent_3pmace covs(aggregate) covm;
    class trt01p;
      id subjid;
    model (starttime,stoptime)*status(0 2)=trt01p / ties=exact;
    baseline covariates=in2 out=out cmf=_all_ / nomean;
run;

*With competing risk;
proc phreg data=recurrent_3pmace;
    class trt01p;
    model stoptime*status(0)=trt01p / entry=starttime eventcode=1;
    strata trt01p; 
    baseline out=mcfdata cif=cuminc;
run;

data out2;
    length trt01p $40.;
    set out(in=out) mcfdata(in=mcfdata drop=trt01p);
    if out then trt01p=strip(trt01p);
    else if mcfdata then do;
        trt01p=strip(trt01p2) || ' (with competing risk)';
        cmf=-log(1-cuminc);
    end;
run;


proc sgplot data=out2;
	step x=stoptime y=cmf/ group=trt01p; 
run; 









/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 
/*************************************************************************************/ 
****** BLADDER CANCER DATA CHECK ; 

libname myout "/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf";

*data bladder1; 
*	set myout.bladder1; 
*run; 

proc import datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/bladder1.csv" out=bladder1
dbms=csv replace;
run;	


* LWYY model with treatment as covariate; 
proc phreg data=bladder1 covs(aggregate) covm;
    class Z(ref="0");
    id id;
    model (startt, stopt)*status(0 2 3)=Z / ties=breslow CONVERGELIKE=1E-9;
	*output out=resOut resmart=resmart;
run;

*proc print data=resOut; *where resmart is missing; *run;

* LWYY model with treatment, size and number as covariate; 
proc phreg data=bladder1 covs(aggregate) covm;
    class Z(ref="0");
    id id;
    model (startt, stopt)*status(0 2 3)=Z size number/ ties=breslow CONVERGELIKE=1E-9;
run;


*GL model with treatment as covariate; 
proc phreg data=bladder1 covs(aggregate) covm;
    class Z(ref="0");
    id id;
    model (startt, stopt)*status3(0)=Z / ties=breslow eventcode=1 CONVERGELIKE=1E-9;
run;


* GL model with treatment, size and number as covariates; 
proc phreg data=bladder1 covs(aggregate) covm;
    class Z(ref="0");
    id id;
    model (startt, stopt)*status3(0)=Z size number / ties=breslow eventcode=1 CONVERGELIKE=1E-9;
run;






*  Cox model on death with treat; 
proc phreg data=bladder1 covs(aggregate) covm;
    class Z(ref="0");
    id id;
    model (startt, stopt)*death(0)=Z / ties=breslow CONVERGELIKE=1E-9;
run;
