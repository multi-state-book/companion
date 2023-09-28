*------------------------------------------------------------------;
*------- Chapter 5, SAS code, LEADER data  ------------------------;
*------------------------------------------------------------------;

* Load data; 
proc import out=leader_mi_3p
	datafile="c:/Users/hnrv/OneDrive - Novo Nordisk/Book/leader/data/leader_mi_3p.csv"
	dbms=csv replace;
run;
data leader_mi; 
	set leader_mi_3p; 
	where type = "recurrent_mi"; 
run; 

*---------------------------------------------------------------;
*--------------------- Figure 4.20 -----------------------------;
*---------------------------------------------------------------;

/* Using "fine-gray model" in PHREG gives an alternative solution to 
  the estimator for CMF using the Breslow type estimator for 
  the baseline mean function (see p. 199 in book). The estimator is not
	exactly the same as Cook-Lawless because of a different procedures 
	for ties of terminating events and censorings. If no ties 
	(or no censorings) it equals Cook & Lawless */

* NELSON-AALEN; 
proc phreg data=leader_mi noprint;
	model stop*status(0 2)=/entry=start;
	id id;
	strata treat;
  baseline out=na_data cumhaz=naa;
run;
data na_est;
	set na_data; 
	type = "Nelson-Aalen";
	cumevent = naa; 
	treat_type = trim(treat) || ", "  || type; 
run; 

* COOK & LAWLESS (GHOSH & LIN);
proc phreg data=leader_mi noprint;
  model (start, stop)*status(0)=/eventcode=1; 
  strata treat;
  baseline out=gl_data cif=cuminc;
run;
data gl_est;
	set gl_data; 
	type = "Ghosh & Lin";
	cumevent = -log(1-cuminc); 
	treat_type = trim(treat) || ", " || type; 
run; 
data comb; 
	set na_est gl_est; 
	time = stop/(365.25/12);
	drop naa cuminc;
run;
proc sgplot data=comb;
	step x=time y=cumevent/group=treat_type justify=left;
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
  step x=time y=na / group=treat;
  step x=time y=GL / group=treat;
	xaxis grid values=(0 to 60 by 12);
	yaxis grid values=(0 to 0.12 by 0.02);
	label time="Time since randomisation (months)";
	label na="Expected number events per subject"; 
run;

*---------------------------------------------------------------;
*----------- In-text: LWYY and Ghosh-Lin models ----------------;
*---------------------------------------------------------------;

title 'LWYY model';
proc phreg data=leader_mi covs(aggregate);
  class treat(ref="0");
  model (start, stop)*status(0 2) = treat / rl; 
  id id;
run;
title 'Ghosh-Lin model';
proc phreg data=leader_mi covs(aggregate);
  class treat(ref="0");
  model (start, stop)*status(0) = treat / rl
		eventcode=1 convergelike=1E-9; 
  id id;
run;
