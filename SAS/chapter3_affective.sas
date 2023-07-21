*------------------------------------------------------------------;
*------- Chapter 3, SAS code, affective data  ---------------------;
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
*--------------------- Table 3.6 -------------------------------;
*---------------------------------------------------------------;

proc phreg data=affective;
	where state=0; 
	model stop*status(2 3)= bip/entry=start rl;
run;

proc phreg data=affective;
	where state=0; 
	model stop*status(2 3)= bip episode/entry=start rl;
run;

proc phreg data=affective;
	where state=0; 
	model stop*status(2 3)= bip episode episode*episode/entry=start rl;
run;


proc phreg data=affective;
	where state=0; 
	model stop*status(2 3)= bip period1 period2 period3 period4/entry=start rl;
	period=year+0.5+stop/12;
	period1=0; period2=0; period3=0; period4=0;
	if 71>period>=66 then period1=1;
	if 76>period>=71 then period2=1;
	if 81>period>=76 then period3=1;
	if period>=81 then period4=1;
run;
