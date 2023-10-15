*------------------------------------------------------------------;
*------- Chapter 3, SAS code, LEADER data  ------------------------;
*------------------------------------------------------------------;

 
*---------------------------------------------------------------;
*--------------------- Table 3.16 ------------------------------;
*---------------------------------------------------------------;
proc import out=leader_mi
	datafile="c:/Users/hnrv/OneDrive - Novo Nordisk/Book/leader/data/leader_mi.csv"
	dbms=csv replace;
run; 

proc phreg data=leader_mi covs(aggregate);
  class id;
  model stop*status(0 2) = treat / entry=start; 
  random id / dist=gamma;
  title 'Table 3.16';
run;
