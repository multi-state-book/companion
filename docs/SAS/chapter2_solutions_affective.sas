*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- EXERCISE SOLUTIONS - AFFECTIVE DISORDERS -----------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;
* We first load the data;

proc import out=affective_data
	datafile='data/affective.csv' 
	dbms=csv replace;
	getnames=yes;
run;


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 2.8 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  2.8.1 ----------------------------------------------------------------------------;

* We will non-parametrically estimate the cumulative event intensities for unipolar and bipolar patients using the phreg procedure 
  where we include 'status(0,2,3)' as censoring variables in the model statement since the event of interest is 1 (admission 
  to hospital) and 'bip' in the strata. Furthermore, we include 'start' in the model statement because of delayed entry. Attention 
  is restricted to records where state=0 (out of hospital).;

title '2.8.1';
proc phreg data=affective_data; 
    where state=0;
	model (start,stop)*status(0,2,3)=;
	strata bip;
	baseline out=data281 cumhaz=naa;
run;

* We plot the result using the gplot procedure;

proc gplot data=data281;
	plot naa*stop=bip /haxis=axis1 vaxis=axis2;
	axis1 order=0 to 340 by 24 label=('Time since first admission');
	axis2 order=0 to 10 by 0.5 label=(a=90 'Cumulative event intensity');
	symbol1 i=stepjl c=red;
	symbol2 i=stepjl c=blue;
run;

* The cumulative event intensity is at all times higher for bipolar patients compared to unipolar patients.;

*----------------------------------------------  2.8.2 ----------------------------------------------------------------------------;

*We fit a simple AG-type model (like the one described in section 2.4.1) using the phreg procedure on the subset of data where 'state = 0', 
i.e. periods where patients are out of hospital.; 

title '2.8.2';
proc phreg data=affective_data; 
	model (start,stop)*status(0,2,3)= bip year;
	where state = 0;
run;

* Thus, we estimate hazard ratio of exp(0.35) = 1.4291 of admission to hospital for being diagnosed with bipolar disorder 
  compared to those diagnosed with unipolar disorder.;

*----------------------------------------------  2.8.3 ----------------------------------------------------------------------------;

* We fit a PWP model using the phreg procedure. The set up is the same as for the AG type model but we also include 'episode' in 
  the strata statement;

title '2.8.3';
proc phreg data=affective_data; 
	model (start,stop)*status(0,2,3)= bip year;
	where state = 0;
	strata episode;
run;

* We get an estimate of the hazard ratio of readmission for patients diagnosed with bipolar disorder of exp(0.2421) = 1.274 compared
  to patients diagnosed with unipolar disorder.;

*----------------------------------------------  2.8.4 ----------------------------------------------------------------------------;

* We get a smaller coefficient of 'bip' with the PWP model than with the AG model. This is because, by taking the number of 
  previous episodes into account (via stratification), some of the discrepancy between bipolar and unipolar patients disappears 
  since the occurrence of repeated episodes is itself affected by the initial diagnosis as discussed in the book.;
