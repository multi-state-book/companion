*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.9 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;
* We first load the data;

proc import out=affective_data
	datafile='data/affective.csv' 
	dbms=csv replace;
	getnames=yes;
run;

*----------------------------------------------  4.9.1 ----------------------------------------------------------------------------;

* We can estimate the mean number of episodes in [0,t] for unipolar and bipolar patients, taking the mortality into account 
  by equation (4.13).

* We make a data set correpsonding to the setting with cycles depicted in figure 1.5, i.e. the interval in hospital is 
  included in the time between events.;

title '4.9.1';
data data491; 
	set affective_data;
	by id;
	retain prev;
	if first.id then prev=0; 
	output; 
	if state=1 then prev=start; if state=0 then prev=stop;
run;

data data491;
	set data491;
	if state = 0 or status = 2 or status = 3;
run;

* Thus the entry and exit time are now 'prev' and 'stop';

* We can estimate S(X-) using the phreg procedure. Since the event of interest is death (status = 2) we include 'status(0 1 3)' as 
  censoring variables in the model statement. The result is saved as 'kmdata491';

proc phreg data=data491;
	class bip;
	model (prev,stop)*status(0 1 3)=;
	id id;
	strata bip;
	baseline out=kmdata491 survival=km;
run;

* Likewise, we can estimate dN(X)/Y(X) using the phreg procedure with censoring variables 'status(0 2 3)' in the model statement. The 
  result is saved as 'naadata491';

proc phreg data=data491;
	class bip;
	model (prev,stop)*status(0 2 3)=;
	id id;
	strata bip;
	baseline out=naadata491 cumhaz=naa;
run;

* We then create a data set for the unipolar patients and one for the bipolar patients containing the estimates of S(X-) and
  dN(X)/Y(X);

data naa_uni;
	set naadata491;
	if bip = 0;
run;

data km_uni; 
	set kmdata491;
	if bip = 0;
run;

data uni;
	merge naa_uni km_uni;
	by stop;
run;

data naa_bip;
	set naadata491;
	if bip = 1;
run;

data km_bip;
	set kmdata491;
	if bip = 1;
run;

data bip;
	merge naa_bip km_bip;
	by stop;
run;

* We then fill the empty cells in the data set with the previous value of S(X-) and dN(X)/Y(X);


data uni;
	set uni;
	retain _km _naa;
	if km ne . then _km = km;
	if naa ne . then _naa = naa;
	years = stop/12;
run;

data bip;
	set bip;
	retain _km _naa;
	if km ne . then _km = km;
	if naa ne . then _naa = naa;
	years = stop/12;
run;

* Finally, we estimate m(t) for unipolar and bipolar patients respectively;

data uni; 
	set uni;
	dA = dif(_naa); 
	if years = 0 then dA = 0;
	dmu = _km*dA;
	retain mu;
	mu + dmu;
	bip = 0;
	keep years _naa mu bip;
run;


data bip; 
	set bip;
	dA = dif(_naa); 
	if years = 0 then dA = 0;
	dmu = _km*dA;
	retain mu;
	mu + dmu;
	bip = 1;
	keep years _naa mu bip;
run;

* We merge the data sets for unipolar and bipolar patients and plot the data using the gplot procedure;

data plotdata491;
	set uni bip;
run;

title '4.9.1';
proc gplot data=plotdata491;
	plot mu*years=bip/ haxis=axis1 vaxis=axis2;
	axis1 order=0 to 30 by 5 minor=none label=('Years');
	axis2 order=0 to 10 by 2 minor=none label=(a=90 'Expected number of episodes');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;

* The expected number of episodes is larger for bipolar patients at all times compared to unipolar patients.;

* Less transparently, but in a way a lot easier, we can 'cheat' proc phreg to do the computations by
  fitting an empty Fine-Gray model and transform the cumulative sub-distribution hazard!;

proc phreg data=data491;
	model stop*status(0 3)=/entry=prev eventcode=1;
	strata bip;
	baseline out=mcfdata1 cif=naa1;
run;

data mcfdata1; set mcfdata1;
	cmf=-log(1-naa1);
	years=stop/12;
run;

proc gplot data=mcfdata1;
plot cmf*years=bip/haxis=axis1 vaxis=axis2;
axis1 order=0 to 30 by 5 minor=none label=('Years');
axis2 order=0 to 8 by 0.5 minor=none label=(a=90 'Expected number of episodes');
symbol1  v=none i=stepjl c=red;
symbol2  v=none i=stepjl c=blue;
run;
quit;

*----------------------------------------------  4.9.2 ----------------------------------------------------------------------------;

* We can estimate the expected number of episodes neglecting mortality by equation (4.11). This is 
  in fact the Nelson-Aalen estimate we saved in the data set 'naadata491';

* We use the gplot procedure to plot the two estimates of the expected number of episodes;

legend1 label = ('Accounting for mortality') value = ('Yes' 'No');

title '4.9.2 - unipolar';
proc gplot data=uni;
	plot mu*years _naa*years/ overlay haxis=axis1 vaxis=axis2 legend =legend1;
	axis1 order=0 to 30 by 5 minor=none label=('Years');
	axis2 order=0 to 6 by 2 minor=none label=(a=90 'Expected number of episodes');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;


title '4.9.2 - bipolar';
proc gplot data=bip;
	plot mu*years _naa*years/ overlay haxis=axis1 vaxis=axis2 legend =legend1;
	axis1 order=0 to 30 by 5 minor=none label=('Years');
	axis2 order=0 to 10 by 2 minor=none label=(a=90 'Expected number of episodes');
	symbol1  v=none i=stepjl c=red;
	symbol2  v=none i=stepjl c=blue;
run;


* For both unipolar and bipolar patients we get larger estimates of the mean number of episodes when we do not account for mortality.;
