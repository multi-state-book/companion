*------------------------------------------------------------------;
*------- Chapter 3, SAS code, PBC3 data ---------------------------;
*------------------------------------------------------------------;

* Load pbc data; 
proc import out=pbc3
	datafile="data/pbc3.csv"
	dbms=csv replace;
run;
data pbc3; 
	set pbc3;
	albnorm=(alb-35)*(alb>35);
	alb10=alb/10;
	alb2=alb*alb;
	bilihigh=(bili-17.1)*(bili>17.1);
	bilitoohigh=(bili-34.2)*(bili>34.2);
	bilimuchtoohigh=(bili-51.3)*(bili>51.3);
	bili100=bili/100;
	bili2=bili*bili;
	log2bili=log2(bili);
	logbilihigh=(log2bili-log2(17.1))*(bili>17.1);
	logbilitoohigh=(log2bili-log2(34.2))*(bili>34.2);
	logbilimuchtoohigh=(log2bili-log2(51.3))*(bili>51.3);
	log2bili2=log2bili*log2bili;
run;

*---------------------------------------------------------------;
*--------------------- Table 3.11 ------------------------------;
*---------------------------------------------------------------;

* Treatment; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili tmenttime/rl;
	tmenttime=(tment=1)*days;
run;

proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili tmentlogtime/rl;
	tmentlogtime=(tment=1)*log(days);
run;

proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili tmentt0/rl;
	tmentt0=(tment=1)*(days>2*365.25);
run;

* Log bilirubin;
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili bilitime/rl;
	bilitime=log2bili*days;
run;

proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili bililogtime/rl;
	bililogtime=log2bili*log(days);
run;

proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili bilit0/rl;
	bilit0=log2bili*(days>2*365.25);
run;

* Albumin; 
proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili albtime/rl;
	albtime=alb*days;
run;

proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili alblogtime/rl;
	alblogtime=alb*log(days);
run;

proc phreg data=pbc3;
	class tment (ref='0');
	model days*status(0)=tment alb log2bili albt0/rl;
	albt0=alb*(days>2*365.25);
run;


