*------------------------------------------------------------------;
*------- Chapter 2, SAS code, affective data  ---------------------;
*------------------------------------------------------------------;

* Load data; 

*---------------------------------------------------------------;
*--------------------- Table 2.13 ------------------------------;
*---------------------------------------------------------------;

proc sort data=affective out=state0;
  where state=0;
  by state episode;
run;

* Cox model for 1., 2., 3., 4. episode 'Markov': Column 1; 
proc phreg data=state0;
  where episode<=4;
	model stop*status (2 3)= bip / entry=start;
	by episode;
run;

* AG model, no past; 
proc phreg data=state0;
	model stop*status (2 3)= bip / entry=start;
run;

* PWP model; 
proc phreg data=state0;	
	model stop*status (2 3)= bip / entry=start;
	strata episode;
run; 

* Cox model for 1., 2., 3., 4. episode 'Gap time': Column 2; 
proc phreg data=state0;
  where episode<=4;
	model wait*status (2 3)= bip;
	by episode;
run;

* AG gap time model; 
proc phreg data=state0;
	model wait*status (2 3)= bip;
run;

* PWP gap time model; 
proc phreg data=state0;
	model wait*status (2 3)= bip;
	strata episode;
run;

