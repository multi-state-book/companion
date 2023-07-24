/* Kaplan-Meier */
proc lifetest data=pbc3 plots=survival(nocensor) notable;
  time followup*fail(0);
run;
proc freq data=pbc3;
where status>0;
table days;
run;
/* Tian et al */
proc rmstreg data=pbc3 tau=3;
  model followup*status(0)=tment alb log2bili / link=linear method=ipcw;
run;
proc rmstreg data=pbc3 tau=3;
  model followup*status(0)=tment alb log2bili / 
        link=linear method=ipcw(strata=tment);
run;

/* Pseudo values */
proc rmstreg data=pbc3 tau=3 outpv=pv3;
  model followup*fail(0)=tment alb log2bili / link=linear;
run;
proc gplot data=pv3;
	plot pseudovalue*followup=fail /haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 4 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;
quit;


/* Tian et al */
proc rmstreg data=pbc3 tau=5;
  model followup*fail(0)=tment alb log2bili / link=linear method=ipcw  converge=0.000000001;
run;

/* Pseudo values */
proc rmstreg data=pbc3 tau=5 outpv=pv5 ;
  model followup*fail(0)=tment alb log2bili / link=linear converge=0.000001;
run;
proc gplot data=pv5;
	plot pseudovalue*followup=fail /haxis=axis1 vaxis=axis2;
	axis1 order=0 to 6 by 1 minor=none label=('Years');
	axis2 order=0 to 6 by 1 minor=none label=(a=90 'Pseudo-values');
	symbol1  v=x i=none c=black;
	symbol2  v=o i=none c=black;
run;
quit;


proc rmstreg data=pbc3 tau=2;
  class tment;
  model followup*fail(0)=tment / link=linear method=ipcw ;
run;
proc rmstreg data=pbc3 tau=2;
  model followup*fail(0)=tment / link=linear ;
run;

proc rmstreg data=pbc3 tau=3;
  model followup*fail(0)=tment / link=linear method=ipcw;
run;
proc rmstreg data=pbc3 tau=3;
  model followup*fail(0)=tment / link=linear ;
run;

proc rmstreg data=pbc3 tau=5;
  model followup*fail(0)=tment / link=linear;
run;

proc rmstreg data=pbc3 tau=3;
  model followup*fail(0)=tment  alb log2bili / link=linear method=ipcw(strata=tment);
run;
proc rmstreg data=pbc3 tau=3;
  model followup*fail(0)=tment  alb log2bili / link=linear method=ipcw;
run;

proc rmstreg data=pbc3 tau=3;
  model followup*fail(0)=tment  alb log2bili / link=linear method=pv;
run;

proc rmstreg data=pbc3 tau=5;
  model followup*fail(0)= / link=linear;
run;
