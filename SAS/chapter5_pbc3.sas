proc import out=pbc3
	datafile="data/pbc3.csv"
	dbms=csv replace;
run;
data pbc3; 
	set pbc3;
	log2bili=log2(bili);
	years=days/365.25;
run;

*---------------------------------------------------------------;
*--------------------- Figure 5.12    --------------------------;
*---------------------------------------------------------------;

ods graphis on; 
proc phreg data=pbc3;
	class tment (ref='0');
	model years*status(0)=tment alb bili / rl;
	assess var=(alb) / resample=1000 npaths=50;
run;


*---------------------------------------------------------------;
*--------------------- Figure 5.13    --------------------------;
*---------------------------------------------------------------;

ods graphis on; 
proc phreg data=pbc3;
	class tment (ref='0');
	model years*status(0)=tment alb bili / rl;
	assess var=(bili) / resample=1000 npaths=50;
run;

*---------------------------------------------------------------;
*--------------------- Figure 5.14 -----------------------------;
*---------------------------------------------------------------;

ods graphis on; 
proc phreg data=pbc3;
	class tment (ref='0');
	model years*status(0)=tment alb log2bili /rl;
	assess var=(log2bili) / resample=1000 npaths=50;
run;

*---------------------------------------------------------------;
*--------------------- Figure 5.15-17 --------------------------;
*---------------------------------------------------------------;

ods graphis on; 
proc phreg data=pbc3;
	class tment (ref='0');
	model years*status(0)=tment alb log2bili /rl;
	assess ph / resample=1000 npaths=50;
run;
