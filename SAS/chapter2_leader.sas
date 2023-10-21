*------------------------------------------------------------------;
*-	------ Chapter 2, SAS code, LEADER data  ------------------------;
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


* AG Cox type;
proc phreg data=leader_mi;
  model stop*status(0 2)= treat / entry=start; 
run;

* AG model, piece-wise constant hazards with 5 equally-sized intervals;
proc icphreg data=leader_mi;
	model stop*status(0 2) = treat / entry=start basehaz=pch(nintervals=5);
run;

* Alternatively, one can split data at time points and use proc genmod: 
* cuts: split by 5 equally-sized intervals;
* For the piece-wise constant model we use the macro `lexis` by Bendix Carstensen,
  The macro splits observation time into time interval. See also the book p. 85.
  For more information on the macro visit https://bendixcarstensen.com and look for Software; 

filename lexis url 'https://bendixcarstensen.com/Lexis/Lexis.sas';
%inc lexis;

/* We name the failure 'mi' and the time intervals 'timegroup' */
data leader_mi; 
	set leader_mi;
	mi=status=1;
run;
%Lexis(data   = leader_mi,
			 out    = pois,
			 entry  = start,
			 exit   = stop,
			 fail   = mi,
			 left   = timegroup,
       breaks = %str(0,379.2,758.4,1137.6,1516.8,1896.0));

proc genmod data=pois;
  class treat(ref="0") timegroup;
	model mi = treat timegroup / dist=poisson  offset=lrisk;
run;

* Cox model first event and PWP models 2nd, ..., 5th event;
proc sort data=leader_mi;
  by eventno;
run;
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

