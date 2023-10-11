*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.8 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* We first load the data;

proc import out=bissau_data datafile='data/bissau.csv' 
		dbms=csv replace;
	getnames=yes;
run;

* To make our results comparable to table 2.12 we must first convert the age variable from days to months.;

data bissau_data;
	set bissau_data;
	agem = age/30.4;
run;

*----------------------------------------------  4.8.1 ----------------------------------------------------------------------------;

* We fit a marginal hazard model for the mortality rate adjusting for cluster and the variables BCG, DTP and age using the 'phreg' 
  procedure where we include 'covs(aggregate)' in the phreg statement to obtain robust SD estimates and 'cluster' in the id 
  statement;

title '4.8.1';
proc phreg data=bissau_data covs(aggregate);
	class cluster;
	model fuptime*dead(0)= bcg dtp agem;
	id cluster;
run;

* We obtain a coefficient of -0.41 for BCG, 0.073 for DTP and 0.04 for age.;

*----------------------------------------------  4.8.2 ----------------------------------------------------------------------------;

* The estimates of the coefficients are almost identical with the ones obtained with the gamma frailty model. The standard 
  deviations of the marginal hazards model are, however, a little bit smaller than the standard errors from the gamma frailty 
  model. 


