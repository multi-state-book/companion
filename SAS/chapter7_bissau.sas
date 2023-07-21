*------------------------------------------------------------------;
*------- Chapter 8, SAS code, Bissau data -------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import 
datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/bissau.csv"
	out=bissau
	dbms = csv
	replace;
run;

* Summarise data set; 
proc contents 
	data=bissau; 
run;

*---------------------------------------------------------------;
*--------------------- Table 8.1 -------------------------------;
*---------------------------------------------------------------;

/* 1.: Fit Cox model */
data bissau; 
	set bissau; 
	agem = int(age/30.44);
run; 

proc phreg data=bissau;
  class bcg agem;
  model fuptime*dead(0)=bcg agem / rl;
run;

/* 2.: Create a 12.5 pct sub-cohort */

data bissau; 
	set bissau;
	seed=260452;
	s=ranbin(seed,1,0.125);
run;

proc freq; 
	tables s*dead; 
run;

/* 3.: Fit Cox model to case-cohort data */

data casecoho; set bissau;
	epsilon=0.001;
	if dead=1 and s=0 then do;
	d=1; start=fuptime-epsilon; stop=fuptime; w=1; output; end;
	if dead=0 and s=1 then do;
	d=0; start=0; stop=fuptime; w=1/0.125; output; end;
	if dead=1 and s=1 then do;
	d=0; start=0; stop=fuptime-epsilon; w=1/0.125; output;
	d=1; start=fuptime-epsilon; stop=fuptime; w=1; output; end;
run;

proc phreg data=casecoho covsandwich(aggregate);
	class bcg agem;
	model stop*d(0)=bcg agem/rl entry=start;
	weight w; id id;
run;

/* NCC MACRO */

data source; set bissau;
	study_id=id;
	age_entry=0;
	age_dlo=fuptime;
	censor=dead;
run;

%macro caseset;
	%let sampling = 1;
	%let ratio = 3;
	/* Enumerate Cases */
	data cases;
	set source ;
	if censor = 1;
	run;
	data cases;
	set cases end = eof;
	if eof then call symput ('ncases', put(_n_,6.));
	run;
	/* Create Risk Set */
	%do iter = 1 %to &ncases;
	data temp_case;
	set cases;

	if _n_ = &iter ;
	call symput ('rs', put(_n_,6.));
	call symput ('age_rs', put(age_dlo,8.)); call symput ('case_id',
	put(study_id,8.));
	run;
	data temp_control;
	set source;
	if age_entry  <= &age_rs  <= age_dlo;
	/* Exclude Index Case */
	if study_id = &case_id then delete;
	number = ranuni(0);
	age_rs = &age_rs;
	censor = 0;
	run;
	/* Sample Controls */
	%if &sampling = 1 %then %do;
	proc sort data = temp_control;
	by number;
	data temp_control;
	set temp_control;
	by age_rs;
	retain m;
	if first.age_rs then m = 0;
	m=m+1;
	if m <= &ratio then output temp_control;
	run;
	%end; 
	/* End If Sampling = 1 */
	/* Combine Case with Controls */
	data rs_&iter;
	set temp_case
	temp_control;
	rs = &rs;
	age_rs = &age_rs;
	run;
	/* DM Output 'Clear'; Log 'Clear'; */

	%end; 
	/* End Loop Creating Risk Set */
	/* Append Risk Sets */
	%do j = 2 %to &ncases;
	proc append base = rs_1 data = rs_&j;
	run;
	%end;
	data final; set rs_1; run;
%mend ; 
/* End Macro */

/* Invoke Macro */
%caseset;

proc phreg data=final;
class bcg agem;
  model fuptime*censor(0)=bcg agem / rl;
  strata rs;
run;

