*------------------------------------------------------------------;
*------- Chapter 5, SAS code, LEADER data  ------------------------;
*------------------------------------------------------------------;


proc import out=leader_mi
	datafile="c:/Users/hnrv/OneDrive - Novo Nordisk/Book/leader/data/leader_mi.csv"
	dbms=csv replace;
run;

proc import out=leader_3p
	datafile="c:/Users/hnrv/OneDrive - Novo Nordisk/Book/leader/data/leader_3p.csv"
	dbms=csv replace;
run;

*---------------------------------------------------------------;
*-------- Section 5.8.3, p. 213: Mao-Lin models ----------------;
*---------------------------------------------------------------;

* recurrent MI+death;
* Death are defined as part of composite and add an extra record with (start,stop)
  of length 0.5 for death events as to be able to trick phreq to estimate Mao-Lin model:
  Mi and deaths count as events and death also as competing risk;
data mi; 
	set leader_mi;
	event=status in (1 2);
  output;
	if status=2 then do;
		event=2;
    start=stop;
    stop=stop+0.5;
		output;
	end;
run;
* Mao-Lin model for recurrent MI+death ; 
proc phreg data=mi covs(aggregate);
  model (start, stop)*event(0) = treat / eventcode=1 rl  convergelike=1E-9; 
  id id;
run;

* recurrent 3p-MACE non-CV deaths incorrectly treated as censoring;
data macecens; 
	set leader_3p;
	event=status in (1 2 3);
  if status=4 then event=0;
  output;
	if status=3 then do;
		event=2;
    start=stop;
    stop=stop+0.5;
		output;
	end;
run;
* Mao-Lin model and non-CV deaths incorrectly treated as censoring; 
proc phreg data=macecens covs(aggregate);
  model (start, stop)*event(0) = treat / eventcode=1 rl  convergelike=1E-9; 
  id id;
run;

* recurrent 3p-MACE;
data mace; 
	set leader_3p;
	event=status in (1 2 3);
  if status=4 then event=2;
  output;
	if status=3 then do;
		event=2;
    start=stop;
    stop=stop+0.5;
		output;
	end;
run;
* Mao-Lin model for 3p-MACE and non-CV deaths treated as competing risks;
proc phreg data=mace covs(aggregate);
  model (start, stop)*event(0) = treat / eventcode=1 rl convergelike=1E-9;
  id id;
run;
