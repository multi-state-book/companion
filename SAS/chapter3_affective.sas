*------------------------------------------------------------------;
*------- Chapter 3, SAS code, affective data  ---------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=affective
	datafile="data/affective.csv"
	dbms=csv replace;
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
	model stop*status(2 3)= bip / entry=start rl type3(lr);
run;

proc phreg data=affective;
	where state=0; 
	model stop*status(2 3)= bip episode / entry=start rl type3(lr);
run;

* Episode as categorical;
data affective2;
  set affective;
	if episode>10 the episode=10;
proc phreg data=affective2;
	where state=0; 
	class episode(ref="1");
	model stop*status(2 3)= bip episode / entry=start rl type3(lr);
run;

proc phreg data=affective;
	where state=0; 
	model stop*status(2 3)= bip episode episode*episode / entry=start rl type3(lr);
run;


*---------------------------------------------------------------;
*--------------------- Table 3.7 -------------------------------;
*---------------------------------------------------------------;

proc phreg data=affective;
	where state=0; 
	model stop*status(2 3)= bip period1 period2 period3 period4 / entry=start rl type3(lr);
	period=year+0.5+stop/12;
	period1=0; period2=0; period3=0; period4=0;
	if 71>period>=66 then period1=1;
	if 76>period>=71 then period2=1;
	if 81>period>=76 then period3=1;
	if period>=81 then period4=1;
run;
