*------------------------------------------------------------------;
*------- Chapter 2, SAS code, affective data  ---------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=affective
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/affective.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=affective; 
run;

data affective; 
	set affective; 
	wait = stop - start; 
run; 

*---------------------------------------------------------------;
*--------------------- Table 2.13 ------------------------------;
*---------------------------------------------------------------;

* Cox model for 1., 2., 3., 4. episode 'Markov': Column 1; 
proc phreg data=affective;
	where  state=0 and episode=1;
	model stop*status (2 3)= bip/entry=start;
run;

proc phreg data=affective;
	where state=0 and episode=2;
	model stop*status (2 3)= bip/entry=start;
run;

proc phreg data=affective;
	where state=0 and episode=3;
	model stop*status (2 3)= bip/entry=start;
run;

proc phreg data=affective;
	where state=0 and episode=4;
	model stop*status (2 3)= bip/entry=start;
run;

* Cox model for 1., 2., 3., 4. episode 'Gap time': Column 2; 
proc phreg data=affective;
	where state=0 and episode=1;
	model wait*status (2 3)= bip;
run;

proc phreg data=affective;
	where state=0 and episode=2;
	model wait*status (2 3)= bip;
run;

proc phreg data=affective;
	where state=0 and episode=3;
	model wait*status (2 3)= bip;
run;

proc phreg data=affective;
	where state=0 and episode=4;
	model wait*status (2 3)= bip;
run;

* AG model, no past; 
proc phreg data=affective;
	where state=0; 
	model stop*status (2 3)= bip/entry=start;
run;

* PWP model; 
proc phreg data=affective;	
	where state=0; 
	model stop*status (2 3)= bip/entry=start;
	strata episode;
run; 

* AG gap time model; 
proc phreg data=affective;
	where state=0; 
	model wait*status (2 3)= bip;
run;

* PWP gap time model; 
proc phreg data=affective;
	where state=0; 
	model wait*status (2 3)= bip;
	strata episode;
run;



