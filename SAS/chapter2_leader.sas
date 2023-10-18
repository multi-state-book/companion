*------------------------------------------------------------------;
*------- Chapter 2, SAS code, LEADER data  ------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=leader_mi
	datafile="c:/Users/hnrv/OneDrive - Novo Nordisk/Book/leader/data/leader_mi.csv"
	dbms=csv replace;
run;
*---------------------------------------------------------------;
*--------------------- Figure 2.15 -----------------------------;
*---------------------------------------------------------------;

* Nelson-Aalen estimates; 
proc phreg data=leader_mi noprint;
	model stop*status(0 2)=/entry=start;
	id id;
	strata treat;
  baseline out=na_data cumhaz=naa;
run;
data na_est;
	set na_data; 
	time = stop/(365.25/12);
run; 
proc sgplot data=na_est;
	step x=time y=naa / group=treat justify=left;
	xaxis grid values=(0 to 60 by 12);
	yaxis grid values=(0 to 0.12 by 0.02);
	label time="Time since randomisation (months)";
	label cumevent="Cumulative hazard"; 
run; 


*---------------------------------------------------------------;
*--------------------- Table 2.15 -----------------------------;
*---------------------------------------------------------------;

proc sort data=leader_mi;
  by eventno;
run;

* AG Cox type;
proc phreg data=leader_mi;
  model stop*status(0 2)= treat / entry=start; 
run;

* Cox model first event and PWP models 2nd, ..., 5th event;
proc phreg data=leader_mi;
  where eventno<=5;
  model stop*status(0 2)= treat / entry=start; 
	by eventno;
run;

* PWP model all events;
proc phreg data=leader_mi;
  model stop*status(0 2)= treat / entry=start; 
	strata eventno;
run;
