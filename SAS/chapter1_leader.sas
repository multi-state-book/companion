*------------------------------------------------------------------;
*------- Chapter 1, SAS code, LEADER data  ------------------------;
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
*--------------------- Table 1.3 -------------------------------;
*---------------------------------------------------------------;

proc freq data=leader_mi; 
	tables status*treat*eventno; 
run; 

proc freq data=leader_3p; 
	tables status*treat*eventno; 
run; 

