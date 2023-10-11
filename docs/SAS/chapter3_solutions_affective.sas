*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 3.11 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We first load the data;

proc import out=affective_data
	datafile='data/affective.csv' 
	dbms=csv replace;
	getnames=yes;
run;

*----------------------------------------------  3.11.1 ----------------------------------------------------------------------------;

* We can fit a gamma frailty model using the phreg procedure where we include a random statement where we specify 
'dist = gamma'.; 

* We will first fit a frailty model for the subset of patients with unipolar disorder.;

title '3.11.1 - unipolar disorder';
proc phreg data=affective_data;
	class id;
	model (start,stop)*status(0,2,3)=episode;
	random id/dist=gamma;
	where state = 0 and bip = 0;
run;

* Thus, we obtain a hazard ratio of exp(0.028) = 1.028 for the number of previous episodes for patients with unipolar disorder.;

* Then, we will fit a frailty model for the patients with bipolar disorder;

title '3.11.1 - bipolar disorder';
proc phreg data=affective_data;
	class id;
	model (start,stop)*status(0,2,3)=episode;
	random id/dist=gamma;
	where state = 0 and  bip = 1;
run;

* We get an estimate of the hazard ratio for the number of episodes of exp(0.077) = 1.08 for patients diagnosed with bipolar 
  disorder.;

*----------------------------------------------  3.11.2 ----------------------------------------------------------------------------;

* The recurrence rates tend to increase with the number of previous episodes both for patients diagnosed with unipolar and bipolar 
  disorder.;
