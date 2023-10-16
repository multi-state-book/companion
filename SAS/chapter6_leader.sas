*------------------------------------------------------------------;
*------- Chapter 6, SAS code, LEADER data  ------------------------;
*------------------------------------------------------------------;


proc import out=leader_mi
	datafile="c:/Users/hnrv/OneDrive - Novo Nordisk/Book/leader/data/leader_mi.csv"
	dbms=csv replace;
run;

*---------------------------------------------------------------;
*-------- Figure 6.11  ------------------------------------------;
*---------------------------------------------------------------;
/* Using "fine-gray model" in PHREG gives an alternative solution to 
  the estimator for CMF using the Breslow type estimator for 
  the baseline mean function (see p. 199 in book). The estimator is not
	exactly the same as Cook-Lawless because of a different procedures 
	for ties of terminating events and censorings. If no ties 
	(or no censorings) it equals Cook & Lawless */

* Left plot;
proc phreg data=leader_mi noprint;
  model (start, stop)*status(0)=/eventcode=1; 
  strata treat;
  baseline out=cmf cif=cuminc;
run;
data cmf;
	set cmf; 
	cumevent = -log(1-cuminc); 
	time = stop/(365.25/12);
run;
ods graphics on;
proc sgplot data=cmf;
	step x=time y=cumevent/group=treat justify=left;
	xaxis grid values=(0 to 60 by 12);
	yaxis grid values=(0 to 0.12 by 0.02);
	label time="Time since randomisation (months)";
	label cumevent="Expected number events per subject"; 
run; 


/*** Calc Cook & Lawless or (Ghosh & Lin (GL)) estimator for CMF 'by hand' ***/
/* First create KM data for death */
proc phreg data=leader_mi noprint;
  model stop*status(0 1)= / entry=start; /* status=2=death */
  strata treat;
  baseline out=kmdata survival=km / method=pl ;
run;
/* Second create NAa data */
proc phreg data=leader_mi noprint;
  model stop*status(0 2)= / entry=start; /* status=1=event */
  strata treat;
  baseline out=nadata cumhaz=na;
run;
/* Use NA data to calculate dA(u), i.e., increments in NAa */
data na;
  set nadata;
  dAu=na-lag(na);
  if stop=0 then dAu=0;
  keep treat stop dAu na;
run;
/* merge NAa and KM data */
data merged;
  merge na kmdata;
  by treat stop;
run;
/* multiply S(u-) and dA(u) */
data fill;
   set merged;
   retain _km;
   if not missing(km) then _km=km;
   else km=_km;
   /* S(u-) */
   S_uminus=lag(km);
   if stop=0 then S_uminus=1;

   if dAu=. then dAu=0;
   GLfactor=S_uminus*dAu;
   keep treat stop na dAu S_uminus GLfactor;
run;
data GLdata;
  set fill;
  by treat;
  if first.treat then GL=0;
  else GL+GLfactor;
	time = stop/(365.25/12);
run;
proc sgplot data=GLdata;
  step x=time y=GL / group=treat;
	xaxis grid values=(0 to 60 by 12);
	yaxis grid values=(0 to 0.12 by 0.02);
	label time="Time since randomisation (months)";
	label na="Expected number events per subject"; 
run;

* Right plot: Kaplan-Meier;
data kmplot;
  set kmdata;
	time = stop/(365.25/12);
run;
proc sgplot data=kmplot;
  step x=time y=km / group=treat;
	xaxis grid values=(0 to 60 by 12);
	yaxis grid values=(0 to 1 by 0.01);
	label time="Time since randomisation (months)";
	label na="Kaplan-Meier estimates"; 
run;
