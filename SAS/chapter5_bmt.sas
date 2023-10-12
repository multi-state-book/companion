*------------------------------------------------------------------;
*------- Chapter 5, bmt data, SAS code ----------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=bmt
	datafile="data/bmt.csv" 
	dbms=csv replace;
run;
* Add extra variables;
data bmt; 
  set bmt;
	intxsurv=timedeath;
	dead=death;
	if rel=1 then intxrel=timerel;
	if rel=0 then intxrel=timedeath;
	trm=0;
	if rel=0 and death=1 then trm=1;
	state0=rel+2*trm;
	if gvhd=1 then tgvhd=timegvhd;
	if gvhd=0 then tgvhd=intxrel;
	dytxanc5=timeanc500*30; 
run;

*---------------------------------------------------------------;
*--------------------- Table 5.1 -------------------------------;
*---------------------------------------------------------------;

* landmarks at 0.5 ...  2.5 mo. horizon 6 mo. ahead;
data landmark; set bmt;
	if intxrel>=0.5 then do;
		time=min(intxrel,6.5); if time<6.5 then status=state0; 
		if time>=6.5 then status=0;
		landmark=0.5; entry=0.5;
		anc=0; if anc500=1 and timeanc500<=0.5 then anc=1;
		gvh=0; if gvhd=1 and tgvhd<=0.5 then gvh=1;
	output; end;
	if intxrel>=1 then do;
		time=min(intxrel,7); if time<7 then status=state0; 
		if time>=7 then status=0;
		landmark=1; entry=1;
		anc=0; if anc500=1 and timeanc500<=1 then anc=1;
		gvh=0; if gvhd=1 and tgvhd<=1 then gvh=1;
	output; end;
		if intxrel>=1.5 then do;
		time=min(intxrel,7.5); if time<7.5 then status=state0; 
		if time>=7.5 then status=0;
		landmark=1.5; entry=1.5;
		anc=0; if anc500=1 and timeanc500<=1.5 then anc=1;
		gvh=0; if gvhd=1 and tgvhd<=1.5 then gvh=1;
	output; end;
	if intxrel>=2 then do;
		time=min(intxrel,8); if time<8 then status=state0; 
		if time>=8 then status=0;
		landmark=2; entry=2;
		anc=0; if anc500=1 and timeanc500<=2 then anc=1;
		gvh=0; if gvhd=1 and tgvhd<=2 then gvh=1;
	output; end;
		if intxrel>=2.5 then do;
		time=min(intxrel,8.5); if time<8.5 then status=state0; 
		if time>=8.5 then status=0;
		landmark=2.5; entry=2.5;
		anc=0; if anc500=1 and timeanc500<=2.5 then anc=1;
		gvh=0; if gvhd=1 and tgvhd<=2.5 then gvh=1;
	output; end;
run;
proc freq data=landmark; 
	tables anc*landmark gvh*landmark/ nocol norow nopercent; 
run;


*---------------------------------------------------------------;
*--------------------- Table 5.2 -------------------------------;
*---------------------------------------------------------------;

data cov1;
	anc05=0; anc1=0; anc15=0; anc2=0; anc25=0; 
	gvh05=0; gvh1=0; gvh15=0; gvh2=0; gvh25=0;
run;
proc phreg data=landmark covs(aggregate);
	class landmark;
	model time*status(0)=anc05 anc1 anc15 anc2 anc25 
	                     gvh05 gvh1 gvh15 gvh2 gvh25/entry=entry;
	anc05=anc*(landmark=0.5); gvh05=gvh*(landmark=0.5);
	anc1=anc*(landmark=1); gvh1=gvh*(landmark=1);
	anc15=anc*(landmark=1.5); gvh15=gvh*(landmark=1.5);
	anc2=anc*(landmark=2); gvh2=gvh*(landmark=2);
	anc25=anc*(landmark=2.5); gvh25=gvh*(landmark=2.5);
	strata landmark;
	id id;
	baseline out=base1 covariates=cov1 survival=km0;
run;


*---------------------------------------------------------------;
*-------------------- Figure 5.1 -------------------------------;
*---------------------------------------------------------------;

title 'Figure 5.1 (a)';
proc gplot data=base1;
	plot km0*time=landmark/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 9 by 1 minor=none label=('Months');
	axis2 order=0 to 1 by 0.1 minor=none 
	label=(a=90 'Conditional survival probability');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
	symbol4 v=none i=stepjl c=green;
	symbol5 v=none i=stepjl c=orange;
run;
quit;

data base1; set base1;
	if landmark=0.5 then km1=km0**exp(-0.33460+0.70277);
	if landmark=1 then km1=km0**exp(-0.60851+0.67945);
	if landmark=1.5 then km1=km0**exp(-1.68550+0.86339);
	if landmark=2 then km1=km0**exp(-2.96384+0.80243);
	if landmark=2.5 then km1=km0**exp(-3.35448+0.83065);
run;

title 'Figure 5.1 (b)';
proc gplot data=base1;
	plot km1*time=landmark/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 9 by 1 minor=none label=('Months');
	axis2 order=0 to 1 by 0.1 minor=none 
	label=(a=90 'Conditional survival probability');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
	symbol4 v=none i=stepjl c=green;
	symbol5 v=none i=stepjl c=orange;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Table 5.3 -------------------------------;
*---------------------------------------------------------------;
* Column 1; 
data cov2;
	anc=0; anctime=0; anctime2=0;  
	gvh=0; gvhtime=0; gvhtime2=0; 
run;
proc phreg data=landmark covs(aggregate);
	model time*status(0)=anc anctime anctime2 
	                     gvh gvhtime gvhtime2/entry=entry;
	anctime=anc*(landmark-0.5)/2; anctime2=anc*((landmark-0.5)/2)**2;
	gvhtime=gvh*(landmark-0.5)/2; gvhtime2=gvh*((landmark-0.5)/2)**2;
	strata landmark;
	id id;
	baseline out=base2 covariates=cov2 survival=km0;
run;
data regn;
	do i=1 to 5;
	b1=-0.33460; b2=-0.60851; b3=-1.68550; b4=-2.96384; b5=-3.35448; 
	c1=0.70277; c2=0.67945; c3=0.86339; c4=0.80243; c5=0.83065;
	t=0.5*(i-1)/2;
	d=-0.32202+(-1.19065)*t+(-2.25731)*t*t; 
	e=0.66257+0.39131*t+(-0.22640)*t*t;
	output;
	end;
proc print;
run; 

* Column 2; 
data cov3;
	anc=0; anctime=0; anctime2=0; 
	gvh=0; gvhtime=0; gvhtime2=0;
	strtime=0; strtime2=0;
run;

proc phreg data=landmark covs(aggregate);
	model time*status(0)=anc anctime anctime2 
	                     gvh gvhtime gvhtime2
	                     strtime strtime2/entry=entry;
	anctime=anc*(landmark-0.5)/2; anctime2=anc*((landmark-0.5)/2)**2;
	gvhtime=gvh*(landmark-0.5)/2; gvhtime2=gvh*((landmark-0.5)/2)**2;
	strtime=(landmark-0.5)/2; strtime2=strtime**2;
	id id;
	baseline out=base3 covariates=cov3 survival=km;
run;

*---------------------------------------------------------------;
*--------------------- Figure 5.2 ------------------------------;
*---------------------------------------------------------------;

proc gplot data=base2;
	plot km0*time=landmark/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 9 by 1 minor=none label=('Months');
	axis2 order=0 to 1 by 0.1 minor=none 
	label=(a=90 'Conditional survival probability');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
	symbol4 v=none i=stepjl c=green;
	symbol5 v=none i=stepjl c=orange;
run;
quit;

data base2; set base2;
	lp=-0.32202+(-1.19065)*(landmark-0.5)/2
	+(-2.25731)*((landmark-0.5)/2)**2
	+0.66257+0.39131*(landmark-0.5)/2+(-0.22640)*((landmark-0.5)/2)**2;
	km1=km0**exp(lp);
run;


proc gplot data=base2;
	plot km1*time=landmark/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 9 by 1 minor=none label=('Months');
	axis2 order=0 to 1 by 0.1 minor=none 
	label=(a=90 'Conditional survival probability');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
	symbol4 v=none i=stepjl c=green;
	symbol5 v=none i=stepjl c=orange;
run;
quit;

*---------------------------------------------------------------;
*--------------------- Figure 5.3 ------------------------------;
*---------------------------------------------------------------;
proc freq; 
	tables gvhtime * gvhtime2; 
run;

data base305; set base3;
	landmark=0.5;
	km0=km; 
	lpz=-0.29798+0.67374;
	km1=km0**exp(lpz);
	if time<6.5 then output;
run;


data base31; set base3;
	landmark=1;
	if time>1 then do;
	lpt=1.41386*0.5/2+1.94038*0.5*0.5/4;
	km0=(km/0.9784097802)**exp(lpt);
	lpz=-0.29798+0.67374+(-1.39310)*0.5/2+(-2.03750)/4*0.5*0.5
	   +0.33291*0.5/2+(-0.17548)/4*0.5*0.5; 
	km1=km0**exp(lpz);
	if time<7 then output;
	end;
run;

data base315; set base3;
	landmark=1.5;
	if time>1.5 then do;
	lpt=1.41386*1/2+1.94038*1*1/4;
	km0=(km/0.9554490221)**exp(lpt);
	lpz=-0.29798+0.67374+(-1.39310)*1/2+(-2.03750)/4*1*1
	   +0.33291*1/2+(-0.17548)/4*1*1;
	km1=km0**exp(lpz);
	if time<7.5 then output;
	end;
run;

data base32; set base3;
	landmark=2;
	if time>2 then do;
	lpt=1.41386*1.5/2+1.94038*1.5*1.5/4;
	km0=(km/0.9395794634)**exp(lpt);
	lpz=-0.29798+0.67374+(-1.39310)*1.5/2+(-2.03750)/4*1.5*1.5
	   +0.33291*1.5/2+(-0.17548)/4*1.5*1.5;
	km1=km0**exp(lpz);
	if time<8 then output;
	end;
run;


data base325; set base3;
	landmark=2.5;
	if time>=2.5 then do;
	lpt=1.41386*2/2+1.94038*2*2/4;
	km0=(km/0.9142153063)**exp(lpt);
	lpz=-0.29798+0.67374+(-1.39310)*2/2+(-2.03750)/4*2*2
	   +0.33291*2/2+(-0.17548)/4*2*2;
	km1=km0**exp(lpz);
	if time<8.5 then output;
	end;
run;

data base3slut; 
	set base305 base31 base315 base32 base325;
run;

proc gplot data=base3slut;
	plot km0*time=landmark/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 9 by 1 minor=none label=('Months');
	axis2 order=0 to 1 by 0.1 minor=none 
	label=(a=90 'Conditional survival probability');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
	symbol4 v=none i=stepjl c=green;
	symbol5 v=none i=stepjl c=orange;
run; 

proc gplot data=base3slut;
	plot km1*time=landmark/haxis=axis1 vaxis=axis2;
	axis1 order=0 to 9 by 1 minor=none label=('Months');
	axis2 order=0 to 1 by 0.1 minor=none 
	label=(a=90 'Conditional survival probability');
	symbol1  v=none i=stepjl c=blue;
	symbol2  v=none i=stepjl c=red;
	symbol3  v=none i=stepjl c=black;
	symbol4 v=none i=stepjl c=green;
	symbol5 v=none i=stepjl c=orange;
run; 
 
 
 
