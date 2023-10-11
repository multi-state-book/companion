*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 5.10 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;
* We first load the data;

proc import out=affective_data
	datafile='data/affective.csv' 
	dbms=csv replace;
	getnames=yes;
run;

* We make a data set correpsonding to the setting with cycles depicted in figure 1.5, i.e. the interval in hospital is 
  included in the time between events.;

title '5.10';
data data510; 
	set affective_data;
	by id;
	retain prev;
	if first.id then prev=0; 
	output; 
	if state=1 then prev=start; if state=0 then prev=stop;
run;
data data510;
	set data510;
	if state = 0 or status = 2 or status = 3;
* Thus the entry and exit time are now 'prev' and 'stop';
run;

* Death are defined as part of composite and add an extra record with (prev,stop) length 0.5 for
  death events as to be able to trick phreq to estimate Mao-Lin model:
  Mi and deaths count as events and death also as competing risk;
data angstML;
	set data510;
	if status ne 2 then output;
	if status = 2 then do;
		status=1; output; 
		prev=stop; stop=stop+0.5; status=2; output;
	end; 
run;
* Mao-Lin model;
proc phreg data=angstML covs(aggregate);
	model stop*status(0)= bip / entry=prev eventcode=1 rl;
	id id;
run;

