*------------------------------------------------------------------;
*------- Chapter 5, SAS code, LEADER data  ------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=leader_mi_3p
	datafile="/projstat/ex2211/ex2211-exploratory/otsot006/stats/program/Draft/jukf/MSB/data/leader_mi_3p.csv"
	dbms=csv replace;
run;

* Summarise data set; 
proc contents 
	data=leader_mi_3p; 
run;

* Make a data set per endpoint; 
data leader_mi; 
	set leader_mi_3p; 
	where type = "recurrent_mi"; 
run; 

data leader_3p; 
	set leader_mi_3p; 
	where type = "recurrent_comb"; 
run; 

*---------------------------------------------------------------;
*----------- In-text, p. 213: Mao-Lin models -------------------;
*---------------------------------------------------------------;

data leader_mi; 
	set leader_mi; 

	if status = 2 then status = 1; 
run;

data leader_mi_ml; 
	set leader_mi;

	if death=1 then do; 
		id = id; 
	    start = stop; 
		stop = start +1; 
		status = 2; 
	output;
	end; 
run;

data leader_mi_ml_extra;
	set leader_mi leader_mi_ml; 
run;

* Mao-Lin model; 
proc phreg data=leader_mi_ml_extra covs(aggregate);
  class treat(ref="0");
  model (start, stop)*status(0) = treat/ 
			ties=breslow eventcode=1 CONVERGELIKE=1E-9; 
  id id;
run;