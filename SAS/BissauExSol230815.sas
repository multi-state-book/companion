*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- EXERCISE SOLUTIONS - GUINEA-BISSAU -----------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;
* We first load the data;

proc import out=bissau_data datafile='data/bissau.csv' 
		dbms=csv replace;
	getnames=yes;
run;


*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 3.11 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

*----------------------------------------------  3.11.1 ----------------------------------------------------------------------------;

* To make our results comparable to table 2.12 we must first convert the age variable from days to months.;

data bissau_data;
	set bissau_data;
	agem = age/(365.25/12);
run;

* Then, we will fit the gamma frailty model using the 'phreg' procedure where we include a random statement with 'dist = gamma' 
  included.;

title '3.11.1';
proc phreg data=bissau_data;
	class cluster;
	model fuptime*dead(0)=bcg dtp agem;
	random cluster/dist=gamma;
run;

* The coefficient for age (month) is estimated to be 0.04 and the coefficient for BCG vaccination is -0.41. In table 2.12 the 
  coefficient for age (month) was 0.055 and the coefficient for BCG was -0.353. Thus, an increasing age increase the hazard and 
  being vaccinated with BCG decrease the hazard in both models in both models.;

* We furthermore note that our model also adjust for the effect of DTP vaccination unlike the model from table 2.12. 
  We could also compare with Table 3.2 where dtp is included. Both models suggest an effect of BCG but not DTP.;


*----------------------------------------------  3.11.2 ----------------------------------------------------------------------------;

* We fit a Cox model stratified on cluster by adding a 'strata' statement in the 'phreg' procedure;

title '3.11.2';
proc phreg data=bissau_data;
	model fuptime*dead(0)=bcg dtp agem;
	strata cluster;
run;

* The Cox model stratified on cluster and the gamma frailty model agree in the direction BCG, DTP and age change the hazard. 
  However, the gamma frailty model estimates a bigger decrease by being vaccinated with BCG, a smaller increase by being vaccinated 
  with DTP and a larger increase of a positive change in age of the hazard compared to the stratified Cox model.;

*----------------------------------------------------------------------------------------------------------------------------------;
*---------------------------------------- EXERCISE 4.8 ----------------------------------------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

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


